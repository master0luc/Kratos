//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//                   Philipp Bucher
//

#if !defined(KRATOS_COMPOSITE_PROPERTY_ASSIGNMENT_HPP_INCLUDED )
#define  KRATOS_COMPOSITE_PROPERTY_ASSIGNMENT_HPP_INCLUDED

// System includes
#include <iterator>
#include <math.h>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

// Application includes
#include "structural_mechanics_application_variables.h"

namespace Kratos {
	///@name Kratos Globals
	///@{
	///@}
	///@name Type Definitions
	///@{
	typedef ModelPart::ElementsContainerType::iterator ElementsIteratorType;

	///@}
	///@name  Enum's
	///@{
	///@}
	///@name  Functions
	///@{
	///@}
	///@name Kratos Classes
	///@{
	/// TODO write a proper description
	/**
	* Example Python Code:
	* # Eigenvectors are first computed and stored in the nodal variable EIGENVECTOR_MATRIX.
	* for step in range(NumEigenvalues):
	*   main_model_part.ProcessInfo[TIME] = float(step+1)
	*   EigenvectorToSolutionStepVariableTransferUtility().Transfer(main_model_part,step,0)
	*   gid_output.PrintOutput()
	*/
	class CompositePropertyAssignment {
	public:
		///@name Type Definitions
		///@{
		KRATOS_CLASS_POINTER_DEFINITION(CompositePropertyAssignment);

		///@}
		///@name Life Cycle
		///@{
		CompositePropertyAssignment(ModelPart& rSubModelPart) : mrModelPart(rSubModelPart) {}

		virtual ~CompositePropertyAssignment() {}

		///@}
		///@name Operators
		///@{
		///@}
		///@name Operations
		///@{

        void Execute(Parameters MethodParameters)
        {
            // TODO do more checks on valid input
            // TODO if needed use the echo_level. If not we can remove it
            MethodParameters.ValidateAndAssignDefaults(mDefaultParameters);
            const std::string method_name = MethodParameters["method"].GetString();
            if (method_name == "")
            {
                KRATOS_ERROR << "Please specify a method name" << std::endl;
            }

            auto specific_parameters = MethodParameters["method_specific_settings"];

            if (method_name == "simple")
            {
                if (specific_parameters.Has("cs_axis_1"))
                {
                    specific_parameters.ValidateAndAssignDefaults(mSimpleMethodCSDefaultParameters);

                    Vector3 cs_axis_1;
                    Vector3 cs_axis_2;
    
                    CheckAndReadVectors(specific_parameters, "cs_axis_1", cs_axis_1);
                    CheckAndReadVectors(specific_parameters, "cs_axis_2", cs_axis_2);
                    // TODO check axis for orthogonality(?)

                    const double rotation_angle = specific_parameters["cs_rotation_angle"].GetDouble();
                    if (rotation_angle > 360.0) // TODO check for all directions
                    {
                        KRATOS_ERROR << "Rotation angle is too large" << std::endl;
                    }

                    if (specific_parameters["cs_normal_axis"].GetInt() == 0)
                    {
                        KRATOS_ERROR << "The normal axis has to be specified" << std::endl;
                    }
                    // TODO check if the input is valid (i.e. 1,2 or 3)

                    std::cout << "Using the simple method with user-defined cs" << std::endl;

                    const int normal_axis_number = specific_parameters["cs_normal_axis"].GetInt();

                    ExecuteCustomCS(cs_axis_1, cs_axis_2, 
                                    normal_axis_number, rotation_angle);
                }
                else 
                {
                    specific_parameters.ValidateAndAssignDefaults(mSimpleMethodDefaultParameters);
                    Vector3 global_fiber_direction;
                    Vector3 normal_vector;
    
                    CheckAndReadVectors(specific_parameters, "global_fiber_direction", global_fiber_direction);
                    CheckAndReadVectors(specific_parameters, "normal_vector", normal_vector);

                    std::cout << "Using the simple method" << std::endl;   
                    ExecuteOLD(global_fiber_direction, normal_vector, 4);                 
                }
            }
            else if (method_name == "advanced")
            {
                specific_parameters.ValidateAndAssignDefaults(mAdvancedMethodDefaultParameters);   
                
                Vector3 global_fiber_direction;
                Vector3 normal_vector;

                CheckAndReadVectors(specific_parameters, "global_fiber_direction", global_fiber_direction);
                CheckAndReadVectors(specific_parameters, "normal_vector", normal_vector);

                const int level = specific_parameters["smoothness_level"].GetInt();

                std::cout << "Using the advanced method" << std::endl;

                ExecuteOLD(global_fiber_direction, normal_vector, level);
            }
            else
            {
                KRATOS_ERROR << "Method \"" << method_name << "\" does not exist" << std::endl;
            }
        }


        void WriteFiberAngles(std::string FileName = "CompositeFiberAngles.txt")
        {
            std::ofstream element_data_file (FileName);

            if (element_data_file.is_open())
            {
                element_data_file << "Begin ElementalData FIBER_ANGLE" << std::endl;
                for (auto& element : mrModelPart.Elements())
                {
                    element_data_file << "\t" << element.Id() 
                                      << " " << element.GetValue(FIBER_ANGLE) 
                                      << std::endl;
                }
                element_data_file << "End ElementalData" << std::endl;
                
                element_data_file.close();
            }
            else std::cout << "Unable to open file";
        }


		///@}
		///@name Access
		///@{
		///@}
		///@name Inquiry
		///@{
		///@}
		///@name Friends
		///@{
		///@}

	protected:
		///@name Protected static Member Variables
		///@{
		///@}
		///@name Protected member Variables
		///@{
		///@}
		///@name Protected Operators
		///@{
		///@}
		///@name Protected Operations
		///@{
		///@}
		///@name Protected  Access
		///@{
		///@}
		///@name Protected Inquiry
		///@{
		///@}
		///@name Protected LifeCycle
		///@{
		///@}

	private:
		///@name Static Member Variables
		///@{
		///@}
		///@name Member Variables
		///@{

        ModelPart& mrModelPart;
            
        Parameters mDefaultParameters = Parameters( R"(
        {
            "method"                   : "",
            "method_specific_settings" : {},
            "echo_level"               : 0
        }  )" );
        
        Parameters mSimpleMethodDefaultParameters = Parameters( R"(
        {
            "global_fiber_direction" : [0,0,0],
            "normal_vector"   : [0,0,0]
        }  )" );
        
        Parameters mSimpleMethodCSDefaultParameters = Parameters( R"(
        {
            "cs_axis_1" : [0,0,0],
            "cs_axis_2" : [0,0,0],
            "cs_rotation_angle" : 0,
            "cs_normal_axis" : 0
        }  )" );
        
        Parameters mAdvancedMethodDefaultParameters = Parameters( R"(
        {
            "global_fiber_direction" : [0,0,0],
            "normal_vector"   : [0,0,0],
            "smoothness_level" : 2
        }  )" );
		///@}
		///@name Private Operators
		///@{
		///@}
		///@name Private Operations
		///@{

        void ExecuteOLD(Vector3 GlobalFiberDirection, Vector3 normalVector, const int Level)
        {
            // Check to see if the composite orientation assignment has already
            // been performed on the current modelPart
            // Just look at the first element to save time
            const ElementsIteratorType& firstElement = mrModelPart.ElementsBegin();
            Properties elementProperties = (*firstElement).GetProperties();

            if ((*firstElement).Has(FIBER_ANGLE))
            {
                // the composite orientation assignment has already been done
                std::cout << "Composite Assignment is skipped, because Fiber Angles are already Present in the first element" << std::endl;
            }
            else
            {
                // perform the composite orientation assignment
                compositeOrientationAssignment(GlobalFiberDirection, normalVector, Level);
            }
        }

        void ExecuteCustomCS(const Vector3 lc1, const Vector3 lc2, 
                             const int normalAxisNumber, const double normalRotationDegrees)
        {
            auto current_process_info = mrModelPart.GetProcessInfo();
            // Check to see if the composite orientation assignment has already
            // been performed on the current modelPart
            // Just look at the first element to save time
            const ElementsIteratorType& firstElement = mrModelPart.ElementsBegin();
            Properties elementProperties = (*firstElement).GetProperties();

            if ((*firstElement).Has(FIBER_ANGLE))
            {
                // the composite orientation assignment has already been done
                std::cout << "Composite Assignment is skipped, because Fiber Angles are already Present in the first element" << std::endl;             
            }
            else
            {
                // Check inputs
                if (abs(inner_prod(lc1,lc2)) > 1E-6)
                {
                    KRATOS_ERROR <<
                        "Local axis 1 and 2 of the user specified CS must be orthogonal" << std::endl;
                }

                // Create a copy of the User Coordinate System (UCS)
                Vector ucs1 = Vector(lc1);
                ucs1 /= std::sqrt(inner_prod(ucs1, ucs1));
                Vector ucs2 = Vector(lc2);
                ucs2 /= std::sqrt(inner_prod(ucs2, ucs2));
                Vector ucs3 = Vector(MathUtils<double>::CrossProduct(ucs1, ucs2));

                // Declare working variables
                Matrix LCSOrientation, R;
                Vector localGlobalFiberDirection, rotation_axis;
                Vector shellLocalAxis1 = ZeroVector(3);
                Vector shellLocalAxis2 = ZeroVector(3);
                Vector shellLocalAxis3 = ZeroVector(3);
                Properties::Pointer pElementProps;


                // Consider what UCS axis the user said was most normal to the shell
                Vector& ucsNormal = ucs3;
                Vector& ucsToBeProjectedOntoShell = ucs1;
                if (normalAxisNumber == 1)
                {
                    // ucs1 is normal to shell
                    // thus, ucs2 will be projected on to the shell
                    ucsNormal = ucs1;
                    ucsToBeProjectedOntoShell = ucs2;
                }
                else if (normalAxisNumber == 2)
                {
                    // ucs2 is normal to shell
                    // thus, ucs3 will be projected on to the shell
                    ucsNormal = ucs2;
                    ucsToBeProjectedOntoShell = ucs3;
                }
                else if (normalAxisNumber == 3)
                {
                    // already done
                }
                else
                {
                    KRATOS_ERROR <<
                        "Normal axis number must be 1, 2 or 3." << std::endl;
                }


                // Rotate the UCS by the given angle about the user specified 
                // normal-most axis
                double rotationAngleRad = normalRotationDegrees / 180.0 * KRATOS_M_PI;
                R = setUpRotationMatrix(rotationAngleRad, ucsNormal);
                ucsToBeProjectedOntoShell = prod(R, ucsToBeProjectedOntoShell);


                // Loop over all elements in part
                for (auto& element : mrModelPart.Elements())
                {
                    // get current element properties
                    pElementProps = element.pGetProperties();

                    // get local orientation of GlobalFiberDirection
                    element.Calculate(LOCAL_ELEMENT_ORIENTATION, LCSOrientation, current_process_info);

                    // get element local axis vectors (global cartesian)
                    for (size_t i = 0; i < 3; i++)
                    {
                        shellLocalAxis1[i] = LCSOrientation(0, i);
                        shellLocalAxis2[i] = LCSOrientation(1, i);
                        shellLocalAxis3[i] = LCSOrientation(2, i);
                    }


                    // normalise local axis vectors (global cartesian)
                    shellLocalAxis1 /= std::sqrt(inner_prod(shellLocalAxis1, shellLocalAxis1));
                    shellLocalAxis2 /= std::sqrt(inner_prod(shellLocalAxis2, shellLocalAxis2));
                    shellLocalAxis3 /= std::sqrt(inner_prod(shellLocalAxis3, shellLocalAxis3));



                    // TODO make this a warning
                    // Check that this user specified normal axis isn't actually 
                    // orthogonal to the shell normal
                    if (abs(inner_prod(ucsNormal, shellLocalAxis3)) < 1E-6)
                    {
                        KRATOS_ERROR << "The user axis (axis" << normalAxisNumber << "=" << ucsNormal << ") you said was normal to the shell is actually orthogonal to shell element "  << element.Id() << "(shell normal = " << shellLocalAxis3 << ")" << std::endl;
                    }
                    

                    // Project the vector onto the shell surface
                    // http://www.euclideanspace.com/maths/geometry/elements/plane/lineOnPlane/index.htm
                    // Projected vector = A
                    // Surface normal = B
                    Vector& A = ucsToBeProjectedOntoShell;
                    const Vector& B = shellLocalAxis3;
                    double B_length = std::sqrt(inner_prod(B, B));
                    Vector ACrossB = Vector(MathUtils<double>::CrossProduct(A, B));
                    ACrossB /= B_length;
                    Vector projectedResult = Vector(MathUtils<double>::CrossProduct(B, ACrossB));
                    projectedResult /= B_length;


                    // Find the angle between our projected result and the 
                    // current shell localAxis1
                    double cosTheta = inner_prod(shellLocalAxis1, projectedResult);
                    double theta = std::acos(cosTheta);
                    // make sure the angle is positively defined according to right
                    // hand rule
                    double dotCheck = inner_prod(shellLocalAxis2, projectedResult);
                    if (dotCheck < 0.0)
                    {
                        // theta is currently negative, flip to positive definition
                        theta *= -1.0;
                    }

                    // set required rotation in element
                    element.SetValue(FIBER_ANGLE, theta);

                }// sub-modelpart element loop



            }
        }



        void compositeOrientationAssignment(Vector3 GlobalFiberDirection, 
                                            Vector3 normalVector, const int Level)
		{
            auto current_process_info = mrModelPart.GetProcessInfo();
			
			// case 1
			// (strictly fulfills alignment via hard orthogonality)
			//
			// Iterative approach which rotates the vector in the element's plane to make sure 'element_fiber_dir <dot> <orthogonal_vector> = 0, where orthogonal_vector = 'global_fiber <cross prod> projection_dir'.


			// case 2 
			// (generally aligns with global fiber)
			// set perform_normal_alignment to FALSE for less accuracy but more smoothness
			bool perform_normal_alignment = true;

			// OPTIONAL (if perform_normal_alignment == true):
			// Computes the 3d rotation between element normal and projection_dir. 
			// 'local_fiber = 3d_rotation x global_fiber'
			
			// ALWAYS PERFORMED:
			// Determine angle between lc1 and local fiber


			// case 3
			// (Abaqus default projection)
			// http://130.149.89.49:2080/v6.8/books/gsa/default.htm?startat=ch05s03.html
			// Shell local axis 1 is the projection of Global X vector onto the shell surface.
			// If the Global X vector is normal to the shell surface, 
			// the shell local 1-direction is the projection of the 
			// Global Z vector onto the shell surface


			// case 4
			// (Abaqus custom projection)
			// http://130.149.89.49:2080/v6.8/books/gsa/default.htm?startat=ch05s03.html
			// Shell local axis 1 is the projection of Global X vector onto the shell surface.
			// If the Global X vector is normal to the shell surface, 
			// the shell local 1-direction is the projection of the 
            // Global Z vector onto the shell surface
            
            // Select approach -------------------------------------------------
            int caseId = 0;
            if (Level == 1)
            {
                caseId = 1; 
            }
            else if (Level == 2)
            {
                caseId = 2; 
            }
            else if (Level == 3)
            {
                caseId = 2; 
                perform_normal_alignment = false;
            }
            else if (Level == 4)
            {
                caseId = 3; 
            }
            else
            {
                KRATOS_ERROR << "Wrong Level!" << std::endl;
            }

			std::cout << "Using composite alignment method ";
			switch (caseId)
			{
			case 1:
				std::cout << caseId << " (strictly fulfills alignment via hard orthogonality)" << std::endl;
				break;
			case 2:
				std::cout << caseId << " (generally aligns with global fiber); normal_alignment=" << perform_normal_alignment << std::endl;
				break;
			case 3:
				std::cout << caseId << " (Abaqus default projection to Global Cartesian axis)" << std::endl;
				break;
			default:
				break;
			}
			
			// -----------------------------------------------------------------

			// Declare working variables
			Matrix LCSOrientation, R;
			Vector localGlobalFiberDirection, rotation_axis;
			Vector localAxis1 = ZeroVector(3);
			Vector localAxis2 = ZeroVector(3);
			Vector localAxis3 = ZeroVector(3);
			double cosTheta, theta, rotation_angle, dotCheck;
			Properties::Pointer pElementProps;
			Vector orthogonal_vector;
			bool debug_printout = true;

			// Optional printout of details for current part
			bool printDetails = true;
			if (printDetails)
			{
				std::cout << "Composite orientation assignment details for model part '" << mrModelPart.Name() << "': " << std::endl;
				std::cout << "\tGlobal fiber direction = \t" << GlobalFiberDirection << std::endl;
				std::cout << "\tProjection direction = \t\t" << normalVector << std::endl;
			}

			// Check incoming GlobalFiberDirection and normalVector are valid
			if (inner_prod(GlobalFiberDirection, GlobalFiberDirection) == 0.0)
			{
				KRATOS_ERROR <<
					"Defined global fiber direction for subModelPart " <<
					mrModelPart.Name() << " has zero length" << std::endl;
			}
			else if (inner_prod(normalVector, normalVector) == 0.0)
			{
				KRATOS_ERROR <<
					"Defined normal vector for subModelPart " <<
					mrModelPart.Name() << " has zero length" << std::endl;
			}

			// Normalize
			GlobalFiberDirection /= std::sqrt(inner_prod(GlobalFiberDirection, GlobalFiberDirection));
			normalVector /= std::sqrt(inner_prod(normalVector, normalVector));

			// Loop over all elements in part
			for (auto& element : mrModelPart.Elements())
			{
				// get current element properties
				pElementProps = element.pGetProperties();

				// get local orientation of GlobalFiberDirection
				element.Calculate(LOCAL_ELEMENT_ORIENTATION, LCSOrientation, current_process_info);

				// get element local axis vectors (global cartesian)
				for (size_t i = 0; i < 3; i++)
				{
					localAxis1[i] = LCSOrientation(0, i);
					localAxis2[i] = LCSOrientation(1, i);
					localAxis3[i] = LCSOrientation(2, i);
				}

				// normalise local axis vectors (global cartesian)
				localAxis1 /= std::sqrt(inner_prod(localAxis1, localAxis1));
				localAxis2 /= std::sqrt(inner_prod(localAxis2, localAxis2));
				localAxis3 /= std::sqrt(inner_prod(localAxis3, localAxis3));

				// Make deep copy of local fiber direction (global cartesian)
				localGlobalFiberDirection = Vector3(GlobalFiberDirection);

				// Flip projection vector such that is it in same dir as LC3
				dotCheck = inner_prod(normalVector, localAxis3);
				Vector correctedNormalVector = Vector(normalVector);
				if (dotCheck < 0.0)
				{
					correctedNormalVector *= -1.0;
				}


				// Perform the assignment method selected ----------------------
				switch (caseId)
				{
				case 1:
					// use hard iterative approach

					// create vector which we must be orthogonal to
					orthogonal_vector = Vector(MathUtils<double>::CrossProduct(GlobalFiberDirection, correctedNormalVector));
					//orthogonal_vector = Vector(MathUtils<double>::CrossProduct(GlobalFiberDirection, localAxis3));

					theta = iterativelyDetermineBestAngle(localAxis1, localAxis3, orthogonal_vector,GlobalFiberDirection);
					break;

				case 2:

					// OPTIONAL (if perform_normal_alignment == true):
					// Computes the 3d rotation between element normal and projection_dir. 
					// 'local_fiber = 3d_rotation x global_fiber'
								
					if (perform_normal_alignment)
					{
						// get rotation matrix to align element normal with projection vec (global cartesian)
						// and apply it to global fiber direction to get 'draped' global fiber direction
						// in the element's plane.
						//
						// Using this option increases alignment accuracy but can also give 'harder' transitions

						rotation_axis = MathUtils<double>::CrossProduct(correctedNormalVector, localAxis3);
						rotation_angle = inner_prod(correctedNormalVector, localAxis3);
						if (abs(rotation_angle) < (1.0 - 1E-6)) // skip if already co-linear
						{
							rotation_angle = std::acos(rotation_angle);
							R = setUpRotationMatrix(rotation_angle, rotation_axis);
							localGlobalFiberDirection = prod(R, localGlobalFiberDirection);
						}
					}


					// ALWAYS PERFORMED:
					// Determine angle between lc1 and local fiber

					// Put everything in local space (local cartesian)
					localGlobalFiberDirection = prod(LCSOrientation, localGlobalFiberDirection);
					localAxis1 = prod(LCSOrientation, localAxis1);
					localAxis2 = prod(LCSOrientation, localAxis2);

					// compute angle 'theta' between local axis 1 and localGlobalFiberDirection (local cartesian)
					cosTheta = inner_prod(localAxis1, localGlobalFiberDirection);
					theta = std::acos(cosTheta);

					// dot between lc2 and localFiberDir (local cartesian)
					dotCheck = inner_prod(localAxis2, localGlobalFiberDirection);
					if (dotCheck < 0.0)
					{
						// theta is currently negative, flip to positive definition
						theta *= -1.0;
					}

					break;

				case 3:
					// (Abaqus default projection)
					// http://130.149.89.49:2080/v6.8/books/gsa/default.htm?startat=ch05s03.html
					// Shell local axis 1 is the projection of Global X vector onto the shell surface.
					// If the Global X vector is normal to the shell surface, 
					// the shell local 1-direction is the projection of the 
					// Global Z vector onto the shell surface

					theta = defaultGlobalProjection(localAxis1, localAxis2, localAxis3);

				default:
					break;
				}

				// set required rotation in element
				element.SetValue(FIBER_ANGLE, theta);

				// add option to write out angles so they don't have to be computed next time
					// or maybe this should be a separate python call
			}// sub-modelpart element loop
		}

		double iterativelyDetermineBestAngle(Vector localAxis1, Vector localAxis3, Vector orthogonal_vector, Vector GlobalFiberDirection)
		{
			double tolerance = 1E-9;
			double steps = 16.0;
			double step_size = 2.0*KRATOS_M_PI / steps; // initially 45 degrees
			double central_angle = 0.0; // from current alignment
			double min_dot_prod = 10.0;
			double best_angle = 0.0;
			bool converged = false;
			int iteration_limit = 20;
			int iteration = 0;
			Vector tempFiber, bestFiber;

			while (converged == false)
			{
				for (size_t angle_step = 0; angle_step < steps; angle_step++)
				{
					double current_angle = best_angle + (angle_step - steps / 2.0)*step_size;
					Matrix R = setUpRotationMatrix(current_angle, localAxis3);
					tempFiber = prod(R, localAxis1);
					double current_dot_prod = inner_prod(tempFiber, orthogonal_vector);
					
					if (abs(current_dot_prod) < abs(min_dot_prod))
					{
						min_dot_prod = current_dot_prod;
						best_angle = current_angle;
						bestFiber = Vector(tempFiber);
					}
				}
				step_size /= (steps / 2);
				iteration++;
				if (abs(min_dot_prod) < tolerance)
				{
					converged = true;
				}
				if (iteration > iteration_limit)
				{
					std::cout << "iteration lim" << std::endl;
					converged = true;
				}
			}

			// Make sure we are pointing in the right direction.
			if (inner_prod(bestFiber, GlobalFiberDirection) < 0.0)
			{
				best_angle += KRATOS_M_PI;
			}

			return best_angle;
		}

		double defaultGlobalProjection(const Vector localAxis1, const Vector localAxis2, const Vector localAxis3)
		{
			// (Abaqus default projection)
			// http://130.149.89.49:2080/v6.8/books/gsa/default.htm?startat=ch05s03.html
			// Shell local axis 1 is the projection of Global X vector onto the shell surface.
			// If the Global X vector is normal to the shell surface, 
			// the shell local 1-direction is the projection of the 
			// Global Z vector onto the shell surface

			// Initially set up Global X vector as global vector
			Vector globalVector = ZeroVector(3);
			globalVector[0] = 1.0;


			// First, check if Global X vector is normal to the shell surface
			if (abs(inner_prod(globalVector, localAxis1)) < 1E-6)
			{
				// Now we use Global Z as the global vector
				globalVector[0] = 0.0;
				globalVector[2] = 1.0;
			}


			// Second, project the global vector onto the shell surface
			// http://www.euclideanspace.com/maths/geometry/elements/plane/lineOnPlane/index.htm
			// Projected vector = A
			// Surface normal = B
			Vector& A = globalVector;
			const Vector& B = localAxis3;
			double B_length = std::sqrt(inner_prod(B, B));
			Vector ACrossB = Vector(MathUtils<double>::CrossProduct(A, B));
			ACrossB /= B_length;
			Vector projectedResult = Vector(MathUtils<double>::CrossProduct(B, ACrossB));
			projectedResult /= B_length;


			// Third, find the angle between our projected result and the 
			// current shell localAxis1
			double cosTheta = inner_prod(localAxis1, projectedResult);
			double theta = std::acos(cosTheta);
			// make sure the angle is positively defined according to right
			// hand rule
			double dotCheck = inner_prod(localAxis2, projectedResult);
			if (dotCheck < 0.0)
			{
				// theta is currently negative, flip to positive definition
				theta *= -1.0;
			}

			return theta;
		}

		Matrix setUpRotationMatrix(double angle, Vector& rotation_axis)
		{
			Matrix rotationMatrix(3, 3, 0.0);

			double u = rotation_axis[0];
			double v = rotation_axis[1];
			double w = rotation_axis[2];

			double L = (u*u + v * v + w * w);
			double u2 = u * u;
			double v2 = v * v;
			double w2 = w * w;

			rotationMatrix(0, 0) = (u2 + (v2 + w2) * cos(angle)) / L;
			rotationMatrix(0, 1) = (u * v * (1 - cos(angle)) - w * sqrt(L) * sin(angle)) / L;
			rotationMatrix(0, 2) = (u * w * (1 - cos(angle)) + v * sqrt(L) * sin(angle)) / L;
			//rotationMatrix(0,3) = 0.0;

			rotationMatrix(1, 0) = (u * v * (1 - cos(angle)) + w * sqrt(L) * sin(angle)) / L;
			rotationMatrix(1, 1) = (v2 + (u2 + w2) * cos(angle)) / L;
			rotationMatrix(1, 2) = (v * w * (1 - cos(angle)) - u * sqrt(L) * sin(angle)) / L;
			//rotationMatrix(1,3) = 0.0;

			rotationMatrix(2, 0) = (u * w * (1 - cos(angle)) - v * sqrt(L) * sin(angle)) / L;
			rotationMatrix(2, 1) = (v * w * (1 - cos(angle)) + u * sqrt(L) * sin(angle)) / L;
			rotationMatrix(2, 2) = (w2 + (u2 + v2) * cos(angle)) / L;
			//rotationMatrix(2,3) = 0.0;

			//rotationMatrix(3,0) = 0.0;
			//rotationMatrix(3,1) = 0.0;
			//rotationMatrix(3,2) = 0.0;
			//rotationMatrix(3,3) = 1.0;

			return rotationMatrix;
        }
        
        void CheckAndReadVectors(Parameters ThisParameters, const std::string KeyName, Vector3& rVector)
        {
            // Note: The parameters have already been validated

            if (ThisParameters[KeyName].size() != 3)
            {
                KRATOS_ERROR << "\" "<< KeyName << "\" is not of size 3" << std::endl;
            }
            
            // This is a workaround until we have the "GetVector" Method
            rVector[0] = ThisParameters[KeyName][0].GetDouble();
            rVector[1] = ThisParameters[KeyName][1].GetDouble();
            rVector[2] = ThisParameters[KeyName][2].GetDouble();

            if (inner_prod(rVector, rVector) < 1E-6)
			{
				KRATOS_ERROR <<	"Vector \" "<< KeyName << "\" has zero length" << std::endl;
			}
        }

		///@}
		///@name Private  Access
		///@{
		///@}
		///@name Private Inquiry
		///@{
		///@}
		///@name Un accessible methods
		///@{
		///@}
	}; // class CompositePropertyAssignment

	   ///@}

	   ///@name Type Definitions
	   ///@{
	   ///@}
}
// namespace Kratos
#endif  // KRATOS_COMPOSITE_PROPERTY_ASSIGNMENT_HPP_INCLUDED defined