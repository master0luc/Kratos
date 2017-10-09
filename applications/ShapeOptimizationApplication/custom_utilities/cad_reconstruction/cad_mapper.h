// ==============================================================================
/*
 KratosShapeOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//==============================================================================
//
//   Project Name:        KratosShape                            $
//   Created by:          $Author:    daniel.baumgaertner@tum.de $
//   Last modified by:    $Co-Author: daniel.baumgaertner@tum.de $
//   Date:                $Date:                      Decem 2016 $
//   Revision:            $Revision:                         0.0 $
//
// ==============================================================================

#ifndef CAD_MAPPER_H
#define CAD_MAPPER_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <math.h>
#include <map>
#include <list>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/python.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "../../kratos/includes/define.h"
#include "../../kratos/processes/process.h"
#include "../../kratos/includes/node.h"
#include "../../kratos/includes/element.h"
#include "../../kratos/includes/model_part.h"
#include "../../kratos/includes/kratos_flags.h"
#include "../../kratos/spatial_containers/spatial_containers.h"
#include "../../kratos/utilities/binbased_fast_point_locator.h"
#include "linear_solvers/linear_solver.h"
#include "patch.h"
#include "brep_element.h"
#include "cad_model_reader.h"
#include "shape_optimization_application.h"
#include "data_point.h"
// ==============================================================================

namespace Kratos
{
class CADMapper
{
  public:
    ///@name Type Definitions
    ///@{

    // Fort vector / matrix operations
	typedef UblasSpace<double, CompressedMatrix, Vector> CompressedSpaceType;
	typedef UblasSpace<double, SparseMatrix, Vector> SparseSpaceType;
	typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;
    typedef typename CompressedSpaceType::MatrixType CompressedMatrixType;
	typedef typename SparseSpaceType::MatrixType SparseMatrixType;
    typedef std::vector<double> DoubleVector;
    typedef std::vector<int> IntVector;
    typedef std::vector<ControlPoint> ControlPointVector;
	typedef std::vector<Patch::Pointer> PatchVector;
	typedef std::vector<BREPElement> BREPElementVector;
	typedef std::vector<BREPGaussPoint> BREPGaussPointVector;
	typedef Node<3> NodeType;
    typedef std::vector<NodeType::Pointer> NodeVector;
    typedef LinearSolver<CompressedSpaceType, DenseSpaceType > CompressedLinearSolverType;
	typedef std::vector<BoundaryEdge> BoundaryEdgeVector;	
	typedef std::vector<BoundaryLoop> BoundaryLoopVector;

	// for tree search
	struct CADPointParameters{
		double u;
		double v;
		// unsigned int patch_itr;
		Patch::Pointer patch_ptr;
		};

	typedef std::map< NodeType::Pointer, CADPointParameters > PointCloud;
	typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
    typedef std::vector<double>::iterator DistanceIterator;
	typedef Bucket< 3, NodeType, NodeVector, NodeType::Pointer, NodeIterator, DistanceIterator > BucketType;
	typedef Tree< KDTreePartition<BucketType> > tree;
	typedef std::vector<double> DistanceVector;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
	typedef std::list<DataPoint> DataPointsList;

	// For handling of python data
	typedef boost::python::extract<double> extractDouble;
    typedef boost::python::extract<bool> extractBool;
	typedef boost::python::extract<unsigned int> extractUnsignedInt;

	// For convenience
	typedef std::vector <Point <3>> VectorPoint; 

    /// Pointer definition of CADMapper
    KRATOS_CLASS_POINTER_DEFINITION(CADMapper);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CADMapper(ModelPart& fe_model_part, boost::python::dict cad_geometry, boost::python::dict cad_integration_data, CompressedLinearSolverType::Pointer linear_solver)
	: mr_fe_model_part(fe_model_part),
	  mr_cad_geometry(cad_geometry),
	  mr_cad_integration_data(cad_integration_data),
	  m_linear_solver(linear_solver)
    {
        // Set precision for output
        std::cout.precision(12);

		// Initialize reader to read CAD model data from the given python dict (json-format)
		m_cad_reader =  CADModelReader(mr_cad_geometry,mr_cad_integration_data);

		// Read cad geometry data into the loal c++ containers
		m_cad_reader.ReadGeometry(m_patches);

		// Create map to identify position in patch_vector for a given patch_id
		for (unsigned int patch_itr = 0; patch_itr < m_patches.size(); patch_itr++)
			m_patch_position_in_patch_vector[m_patches[patch_itr]->GetId()] = patch_itr;
		
		// Read cad geometry data into the loal c++ containers if available
		if(len(mr_cad_integration_data.keys())>0)
			m_cad_reader.ReadIntegrationData(m_brep_elements);
    }

    /// Destructor.
    virtual ~CADMapper()
    {
    }
	// INTEGRAL APPROACH ########################################################
		// --------------------------------------------------------------------------
		void compute_mapping_matrix(const unsigned int u_resolution, const unsigned int v_resolution)
		{
			std::cout << "\n> Starting computation of mapping matrix..." << std::endl;
			boost::timer function_timer;

			// 1st step: Create for each patch coarse cloud of CAD points in x-y space for neighbor search later 

			// Each correspondingly created point is stored in a list.
			// Further lists in the same order are created to store the respective u & v parameter as well as the patch id
			// As list iterator we use a counter for the number of CAD nodes
			unsigned int cad_node_counter = 0;
			NodeVector list_of_cad_nodes;
			DoubleVector list_of_us_of_cad_nodes;
			DoubleVector list_of_vs_of_cad_nodes;
			IntVector list_of_patch_itrs_of_cad_nodes;

			//Loop over all surface of all patches / faces
			for (unsigned int patch_itr = 0; patch_itr < m_patches.size(); patch_itr++)
			{
				// Get relevant data
				unsigned int patch_id =  m_patches[patch_itr]->GetId();
				DoubleVector& knot_vec_u_i = m_patches[patch_itr]->GetSurface().GetKnotVectorU();
				DoubleVector& knot_vec_v_i = m_patches[patch_itr]->GetSurface().GetKnotVectorV();
				std::cout << "\n> Processing Patch with brep_id " << patch_id << std::endl;

				// Pre-calculations
				unsigned int knot_vector_u_dimension = knot_vec_u_i.size();
				unsigned int knot_vector_v_dimension = knot_vec_v_i.size();

				double u_min = knot_vec_u_i[0];
				double u_max = knot_vec_u_i[knot_vector_u_dimension-1];
				double v_min = knot_vec_v_i[0];
				double v_max = knot_vec_v_i[knot_vector_v_dimension-1];

				double delta_u = (u_max-u_min) / u_resolution;
				double delta_v = (v_max-v_min) / v_resolution;

				// Loop over all u & v according to specified resolution
				for(unsigned int i=1; i<u_resolution; i++)
				{
					// current u-value
					double u_i = u_min + i*delta_u;

					for(unsigned int j=1; j<v_resolution; j++)
					{
						// current v-value
						double v_j = v_min + j*delta_v;

						// Check if u_i and v_j represent a point inside the closed boundary loop
						array_1d<double, 2> point_of_interest;
						point_of_interest[0] = u_i;
						point_of_interest[1] = v_j;
						bool point_is_inside = m_patches[patch_itr]->CheckIfPointIsInside(point_of_interest);

						if(point_is_inside)
						{
							// compute unique point in CAD-model for given u&v
							++cad_node_counter;					
							Point<3> cad_point_coordinates;
							m_patches[patch_itr]->GetSurface().EvaluateSurfacePoint(cad_point_coordinates, u_i, v_j);

							// Add id to point --> node. Add node to list of CAD nodes
							NodeType::Pointer new_cad_node = Node < 3 > ::Pointer(new Node<3>(cad_node_counter, cad_point_coordinates));
							list_of_cad_nodes.push_back(new_cad_node);

							// Store for cad node the corresponding cad information in separate vectors
							list_of_us_of_cad_nodes.push_back(u_i);
							list_of_vs_of_cad_nodes.push_back(v_j);
							list_of_patch_itrs_of_cad_nodes.push_back(patch_itr);
						}
					}
				}
			}

			// 2nd step: Construct KD-Tree with all cad nodes
			std::cout << "\n> Starting construction of search-tree..." << std::endl;
			boost::timer timer;
			typedef Bucket< 3, NodeType, NodeVector, NodeType::Pointer, NodeIterator, DistanceIterator > BucketType;
			typedef Tree< KDTreePartition<BucketType> > tree;
			int bucket_size = 20;
			tree nodes_tree(list_of_cad_nodes.begin(), list_of_cad_nodes.end(), bucket_size);
			std::cout << "> Time needed for constructing search-tree: " << timer.elapsed() << " s" << std::endl;

			// 3rd step: Evaluate nearest CAD nodes for all FEM Gauss points. Flag corresponding control points as relevant for mapping
			NodeVector list_of_nearest_points;
			DoubleVector list_of_u_of_nearest_points;
			DoubleVector list_of_v_of_nearest_points;
			DoubleVector list_of_span_u_of_nearest_points;
			DoubleVector list_of_span_v_of_nearest_points;		
			IntVector list_of_patch_of_nearest_points;

			// Loop over all integration points of fe-model-part and find corresponding closest neighbors of cad-model
			// We assume the surface in the fe-model to be mapped is described by conditions (e.g. ShapeOptimizationConditions)

			std::cout << "\n> Starting to identify neighboring CAD points..." << std::endl;
			boost::timer timer_2;

			for (ModelPart::ConditionsContainerType::iterator cond_i = mr_fe_model_part.ConditionsBegin(); cond_i != mr_fe_model_part.ConditionsEnd(); ++cond_i)
			{
				// Get geometry information of current condition
				Condition::GeometryType& geom_i = cond_i->GetGeometry();

				// Evaluate shape functions of FE model according specified integration methdod
				const Condition::GeometryType::IntegrationPointsArrayType& integration_points = geom_i.IntegrationPoints(m_integration_method);
				const unsigned int number_of_integration_points = integration_points.size();

				for ( unsigned int PointNumber = 0; PointNumber < number_of_integration_points; PointNumber++ )
				{
					// Compute global coordinates of current integration point and get corresponding weight
					NodeType::CoordinatesArrayType ip_coordinates = geom_i.GlobalCoordinates(ip_coordinates, integration_points[PointNumber].Coordinates());
					NodeType::Pointer gauss_point_i = Node < 3 > ::Pointer(new Node<3>(PointNumber, ip_coordinates ));

					// Search nearest cad neighbor of current integration point
					NodeType::Pointer nearest_point = nodes_tree.SearchNearestPoint( *gauss_point_i );

					// Recover CAD information of nearest point
					double u_of_nearest_point = list_of_us_of_cad_nodes[nearest_point->Id()-1];
					double v_of_nearest_point = list_of_vs_of_cad_nodes[nearest_point->Id()-1];
					int patch_itr_of_nearest_point = list_of_patch_itrs_of_cad_nodes[nearest_point->Id()-1];

					// Perform Newton-Raphson for detailed search

					// Initialize P: point on the mesh
					Vector P = ZeroVector(3);
					P(0) = ip_coordinates[0];
					P(1) = ip_coordinates[1];
					P(2) = ip_coordinates[2];
					// Initialize Q_k: point on the CAD surface
					Vector Q_k = ZeroVector(3);
					Q_k(0) = nearest_point->X();
					Q_k(1) = nearest_point->Y();
					Q_k(2) = nearest_point->Z();
					// Initialize what's needed in the Newton-Raphson iteration				
					Vector Q_minus_P = ZeroVector(3); // Distance between current Q_k and P
					Matrix myHessian = ZeroMatrix(2,2);
					Vector myGradient = ZeroVector(2);
					double det_H = 0;
					Matrix InvH = ZeroMatrix(2,2);				
					double u_k = u_of_nearest_point;
					double v_k = v_of_nearest_point;
					Point<3> newtonRaphPoint;

					double norm_deltau = 100000000;
					unsigned int k = 0;
					unsigned int max_itr = 50;
					while (norm_deltau > 1e-5)
					{
						// The distance between Q (on the CAD surface) and P (on the FE-mesh) is evaluated
						Q_minus_P(0) = Q_k(0) - P(0);
						Q_minus_P(1) = Q_k(1) - P(1);
						Q_minus_P(2) = Q_k(2) - P(2);

						// The distance is used to compute Hessian and gradient
						m_patches[patch_itr_of_nearest_point]->GetSurface().EvaluateGradientsForClosestPointSearch(Q_minus_P, myHessian, myGradient , u_k, v_k);

						// u_k and v_k are updated
						MathUtils<double>::InvertMatrix( myHessian, InvH, det_H );
						Vector deltau = prod(InvH,myGradient);
						u_k -= deltau(0);
						v_k -= deltau(1);

						// Q is updated
						m_patches[patch_itr_of_nearest_point]->GetSurface().EvaluateSurfacePoint(newtonRaphPoint, u_k, v_k);
						Q_k(0) = newtonRaphPoint[0];
						Q_k(1) = newtonRaphPoint[1];
						Q_k(2) = newtonRaphPoint[2];
						norm_deltau = norm_2(deltau);

						k++;

						if(k>max_itr)
						{
							std::cout << "WARNING!!! Newton-Raphson to find closest point did not converge in the following number of iterations: " << k-1 << std::endl;
							KRATOS_WATCH(Q_k);
							KRATOS_WATCH(P);
						}
					}

					// Update nearest point
					u_of_nearest_point = u_k;
					v_of_nearest_point = v_k;
					nearest_point->X() = Q_k(0);
					nearest_point->Y() = Q_k(1);
					nearest_point->Z() = Q_k(2);

					// Compute and store span of each parameter to avoid redundant computations later
					IntVector knot_span_nearest_point = m_patches[patch_itr_of_nearest_point]->GetSurface().GetKnotSpan(u_of_nearest_point, v_of_nearest_point);

					// Set flag to mark control point as relevant for mapping
					int span_u_of_np = knot_span_nearest_point[0];
					int span_v_of_np = knot_span_nearest_point[1];
					m_patches[patch_itr_of_nearest_point]->GetSurface().FlagControlPointsForMapping(span_u_of_np, span_v_of_np, u_of_nearest_point, v_of_nearest_point);

					// Store information about nearest point in vector for recovery in the same loop later when the mapping matrix is constructed
					list_of_nearest_points.push_back(nearest_point);
					list_of_u_of_nearest_points.push_back(u_of_nearest_point);
					list_of_v_of_nearest_points.push_back(v_of_nearest_point);
					list_of_span_u_of_nearest_points.push_back(span_u_of_np);
					list_of_span_v_of_nearest_points.push_back(span_v_of_np);				
					list_of_patch_of_nearest_points.push_back(patch_itr_of_nearest_point);
				}
			}
			std::cout << "> Time needed for identify neighboring CAD points: " << timer_2.elapsed() << " s" << std::endl;

			// Then we identify mapping relevant control points required from the specified boundary conditions
			// Accordingly we check all Gauss points of all brep elements for their control points
			for (BREPElementVector::iterator brep_elem_i = m_brep_elements.begin(); brep_elem_i != m_brep_elements.end(); ++brep_elem_i)
			{
				// Get Gauss points of current brep element
				BREPGaussPointVector brep_gps = brep_elem_i->GetGaussPoints();

				// Loop over all Gauss points of current brep element 
				for (BREPGaussPointVector::iterator brep_gp_i = brep_gps.begin(); brep_gp_i != brep_gps.end(); ++brep_gp_i)
				{
					// Flag control points on master patch
					unsigned int master_patch_id = brep_gp_i->GetPatchId();
					Vector location = brep_gp_i->GetLocation();
					m_patches[m_patch_position_in_patch_vector[master_patch_id]]->GetSurface().FlagControlPointsForMapping(-1, -1, location[0], location[1]);

					// Flag control points on slave patch if brep element is a coupling element
					if(brep_elem_i->HasCouplingCondition())
					{
						unsigned int slave_patch_id = brep_gp_i->GetSlavePatchId();
						location = brep_gp_i->GetSlaveLocation();
						m_patches[m_patch_position_in_patch_vector[slave_patch_id]]->GetSurface().FlagControlPointsForMapping(-1, -1, location[0], location[1]);
					}
				}
			}

			// Count relevant control points and assign each a unique mapping matrix Id (iterator over points)

			// First we identify relevant control points affecting the Gauss points on the surface
			m_n_control_points = 0;
			m_n_relevant_control_points = 0;
			unsigned int mapping_matrix_id = 0;
			for (PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
			{
				for (ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
				{
					if(cp_i->IsRelevantForMapping())
					{
						cp_i->SetMappingMatrixId(mapping_matrix_id);
						++m_n_relevant_control_points;
						++mapping_matrix_id;
					}
					++m_n_control_points;
				}
			}
			std::cout << "\n> Number of control points in total = " << m_n_control_points << "." << std::endl;
			std::cout << "\n> Number of control points relevant for mapping = " << m_n_relevant_control_points << ".\n" << std::endl;

			// Count FE nodes and assign each a unique mapping matrix id (iterator over nodes)
			m_n_relevant_fem_points = 0;
			for (ModelPart::NodesContainerType::iterator node_i = mr_fe_model_part.NodesBegin(); node_i != mr_fe_model_part.NodesEnd(); ++node_i)
			{
				node_i->SetValue(MAPPING_MATRIX_ID,m_n_relevant_fem_points);
				m_n_relevant_fem_points++;
			}

			// Initialize mapping matrix and corresponding mapping rhs-vector
			m_mapping_matrix_CAD_CAD.resize(3*m_n_relevant_control_points,3*m_n_relevant_control_points);
			m_mapping_rhs_vector.resize(3*m_n_relevant_control_points);
			m_mapping_matrix_CAD_FEM.resize(3*m_n_relevant_control_points,3*m_n_relevant_fem_points);
			m_mapping_matrix_CAD_CAD.clear();
			m_mapping_matrix_CAD_FEM.clear();
			m_mapping_rhs_vector.clear();
			
			// Compute mapping matrix
			unsigned int fem_gp_itr = 0;
			for (ModelPart::ConditionsContainerType::iterator cond_i = mr_fe_model_part.ConditionsBegin(); cond_i != mr_fe_model_part.ConditionsEnd(); ++cond_i)
			{
				// Get geometry information of current condition
				Condition::GeometryType& geom_i = cond_i->GetGeometry();
				unsigned int n_fem_nodes = geom_i.size();

				// Get and store mapping matrix ids of nodes of current condition
				Vector mapping_matrix_ids_fem = ZeroVector(n_fem_nodes);
				for(unsigned int i=0; i<n_fem_nodes;i++)
					mapping_matrix_ids_fem[i] = geom_i[i].GetValue(MAPPING_MATRIX_ID);

				// Evaluate shape functions of FE model according specified integration methdod
				const Condition::GeometryType::IntegrationPointsArrayType& integration_points = geom_i.IntegrationPoints(m_integration_method);
				const unsigned int number_of_integration_points = integration_points.size();
				const Matrix& N_container = geom_i.ShapeFunctionsValues(m_integration_method);

				for ( unsigned int PointNumber = 0; PointNumber < number_of_integration_points; PointNumber++ )
				{
					// Get weight for integration
					double integration_weight = integration_points[PointNumber].Weight();

					// Get FEM-shape-function-value for current integration point
					Vector N_FEM_GPi = row( N_container, PointNumber);

					// Recover information about nearest CAD point				
					double u_of_nearest_point = list_of_u_of_nearest_points[fem_gp_itr];
					double v_of_nearest_point = list_of_v_of_nearest_points[fem_gp_itr];
					unsigned int span_u_of_nearest_point = list_of_span_u_of_nearest_points[fem_gp_itr];
					unsigned int span_v_of_nearest_point = list_of_span_v_of_nearest_points[fem_gp_itr];
					unsigned int patch_itr_of_nearest_point = list_of_patch_of_nearest_points[fem_gp_itr];
					
					// Get CAD-shape-function-value for all control points affecting the nearest cad point
					matrix<double> R_CAD_Pi;
					m_patches[patch_itr_of_nearest_point]->GetSurface().EvaluateNURBSFunctions( span_u_of_nearest_point,
																							span_v_of_nearest_point,
																							u_of_nearest_point, 
																							v_of_nearest_point,
																							R_CAD_Pi );
					
					// Get the corresponding ids of control points in the mapping matrix
					matrix<unsigned int> mapping_matrix_ids_cad = m_patches[patch_itr_of_nearest_point]->GetSurface().GetMappingMatrixIds( span_u_of_nearest_point, 
																																		span_v_of_nearest_point, 
																																		u_of_nearest_point, 
																																		v_of_nearest_point );

					// Assemble mapping and RHS matrix 
					for(unsigned int i=0; i<mapping_matrix_ids_cad.size2();i++)
					{
						for(unsigned int j=0; j<mapping_matrix_ids_cad.size1();j++)
						{

							unsigned int R_row_id = mapping_matrix_ids_cad(j,i);
							double R_row = R_CAD_Pi(j,i);

							// First we assemble CAD-FEM matrix
							for (unsigned int k=0; k<n_fem_nodes;k++)
							{
								unsigned int N_id = mapping_matrix_ids_fem[k];
								double N = N_FEM_GPi[k];

								m_mapping_matrix_CAD_FEM(3*R_row_id+0,3*N_id+0) += integration_weight * R_row * N;
								m_mapping_matrix_CAD_FEM(3*R_row_id+1,3*N_id+1) += integration_weight * R_row * N;
								m_mapping_matrix_CAD_FEM(3*R_row_id+2,3*N_id+2) += integration_weight * R_row * N;
							}

							// Then we assemble CAD-CAD matrix
							for(unsigned int k=0; k<mapping_matrix_ids_cad.size2();k++)
							{
								for(unsigned int l=0; l<mapping_matrix_ids_cad.size1();l++)
								{
									unsigned int R_coll_id = mapping_matrix_ids_cad(l,k);
									double R_coll = R_CAD_Pi(l,k);

									m_mapping_matrix_CAD_CAD(3*R_row_id+0,3*R_coll_id+0) += integration_weight * R_row * R_coll;
									m_mapping_matrix_CAD_CAD(3*R_row_id+1,3*R_coll_id+1) += integration_weight * R_row * R_coll;
									m_mapping_matrix_CAD_CAD(3*R_row_id+2,3*R_coll_id+2) += integration_weight * R_row * R_coll;
								}
							}
						}
					}
					// Update iterator to count the FEM Gauss points
					fem_gp_itr++;
				}
			}
			std::cout << "\n> Finished computation of mapping matrix in " << function_timer.elapsed() << " s." << std::endl;
		}

		// --------------------------------------------------------------------------
		void apply_boundary_conditions( double penalty_factor_disp, 
										double penalty_factor_rot, 
										double penalty_factor_dirichlet, 
										boost::python::list& edges_with_specific_dirichlet_conditions, 
										boost::python::list& edges_with_enforced_tangent_continuity )
		{
			std::cout << "\n> Starting to apply boundary conditions..." << std::endl;
			boost::timer function_timer;

			// Loop over all brep elements specifying boundary conditions 
			for (BREPElementVector::iterator brep_elem_i = m_brep_elements.begin(); brep_elem_i != m_brep_elements.end(); ++brep_elem_i)
			{
				// Check if brep_elem_i is a element used for coupling or for Dirichlet boundary conditions
				if(brep_elem_i->HasCouplingCondition())
					apply_coupling_condition(brep_elem_i, penalty_factor_disp, penalty_factor_rot, edges_with_enforced_tangent_continuity);
				else if(brep_elem_i->HasDirichletCondition())
					apply_dirichlet_condition(brep_elem_i, penalty_factor_dirichlet, edges_with_specific_dirichlet_conditions);
			}

			std::cout << "\n> Finished applying coupling boundary conditions in " << function_timer.elapsed() << " s." << std::endl;
		}

		// --------------------------------------------------------------------------
		void apply_coupling_condition( BREPElementVector::iterator &brep_elem_i, 
									double penalty_factor_disp, 
									double penalty_factor_rot, 
									boost::python::list& edges_with_enforced_tangent_continuity )
		{
			// Get Gauss points of current brep element
			BREPGaussPointVector brep_gps = brep_elem_i->GetGaussPoints();

			// Check if for current element some continuity is to be enforced
			bool tangent_continuity_to_be_enforced = false;
			double penalty_factor_tangent_continuity = 0.0;
			for (unsigned int i = 0; i < boost::python::len(edges_with_enforced_tangent_continuity); ++i)
			{
				unsigned int listed_edge_id = extractUnsignedInt(edges_with_enforced_tangent_continuity[i][0]);
				if(brep_elem_i->GetEdgeId() == listed_edge_id)
				{
					tangent_continuity_to_be_enforced = true;
					double extracted_factor = extractDouble(edges_with_enforced_tangent_continuity[i][1]);
					penalty_factor_tangent_continuity = extracted_factor;
				}
			}

			// Loop over all Gauss points of current brep element 
			for (BREPGaussPointVector::iterator brep_gp_i = brep_gps.begin(); brep_gp_i != brep_gps.end(); ++brep_gp_i)
			{
				// Read information from Gauss point
				unsigned int master_patch_id = brep_gp_i->GetPatchId();
				unsigned int slave_patch_id = brep_gp_i->GetSlavePatchId();
				Patch::Pointer master_patch = m_patches[m_patch_position_in_patch_vector[master_patch_id]];
				Patch::Pointer slave_patch = m_patches[m_patch_position_in_patch_vector[slave_patch_id]];
				double gp_i_weight = brep_gp_i->GetWeight();
				Vector location_on_master_patch = brep_gp_i->GetLocation();
				Vector location_on_slave_patch = brep_gp_i->GetSlaveLocation();
				Vector tangent_on_master_patch = brep_gp_i->GetTangent();
				Vector tangent_on_slave_patch = brep_gp_i->GetSlaveTangent();

				// Evaluate NURBS basis function for Gauss point on both patches and get the corresponding ids of control points in the mapping matrix
				matrix<double> R_gpi_master;
				double u_m = location_on_master_patch(0);
				double v_m = location_on_master_patch(1);
				master_patch->GetSurface().EvaluateNURBSFunctions(-1,-1,u_m, v_m, R_gpi_master);
				matrix<unsigned int> mapping_matrix_ids_gpi_master = master_patch->GetSurface().GetMappingMatrixIds(-1,-1,u_m, v_m);

				matrix<double> R_gpi_slave;
				double u_s = location_on_slave_patch(0);
				double v_s = location_on_slave_patch(1);
				slave_patch->GetSurface().EvaluateNURBSFunctions(-1,-1,u_s, v_s, R_gpi_slave);	
				matrix<unsigned int> mapping_matrix_ids_gpi_slave = slave_patch->GetSurface().GetMappingMatrixIds(-1,-1,u_s, v_s);							

				// Compute Jacobian J1
				matrix<double> g_master = master_patch->GetSurface().GetBaseVectors(-1,-1,u_m,v_m);
				Vector g1 = ZeroVector(3);
				g1(0) = g_master(0,0);
				g1(1) = g_master(1,0);
				g1(2) = g_master(2,0);
				Vector g2 = ZeroVector(3);
				g2(0) = g_master(0,1);
				g2(1) = g_master(1,1);
				g2(2) = g_master(2,1);
				double J1 = norm_2( g1* tangent_on_master_patch(0) + g2* tangent_on_master_patch(1) );

				// First we introduce coupling of displacements
				apply_displacement_coupling( mapping_matrix_ids_gpi_master, 
											mapping_matrix_ids_gpi_slave, 
											R_gpi_master, 
											R_gpi_slave, 
											J1,
											gp_i_weight,
											penalty_factor_disp );

				// Then check if for current element, tangent continuity is to be enforced. If yes, we enforce the tangent continuity..
				if(tangent_continuity_to_be_enforced)
				{
					enforce_tangent_continuity( master_patch, 
												slave_patch,
												u_m, v_m,
												u_s, v_s,
												tangent_on_master_patch, 
												tangent_on_slave_patch,
												mapping_matrix_ids_gpi_master, 
												mapping_matrix_ids_gpi_slave,
												J1,
												gp_i_weight,
												penalty_factor_tangent_continuity );
				}
				// ...if no, we introduce coupling of rotations
				else
					apply_rotation_coupling( master_patch, 
											slave_patch,
											u_m, v_m,
											u_s, v_s,
											tangent_on_master_patch, 
											tangent_on_slave_patch,
											mapping_matrix_ids_gpi_master, 
											mapping_matrix_ids_gpi_slave,
											J1,
											gp_i_weight,
											penalty_factor_rot );
			}
		}

		// --------------------------------------------------------------------------
		void apply_displacement_coupling( matrix<unsigned int> &mapping_matrix_ids_gpi_master,
										matrix<unsigned int> &mapping_matrix_ids_gpi_slave,
										matrix<double> &R_gpi_master, 
										matrix<double> &R_gpi_slave, 
										double J1,
										double gp_i_weight,
										double penalty_factor_disp )
		{	
			// First we consider the relation Master-Master ( MM )
			for(unsigned int i=0; i<mapping_matrix_ids_gpi_master.size2();i++)
			{
				for(unsigned int j=0; j<mapping_matrix_ids_gpi_master.size1();j++)
				{
					unsigned int R_row_id = mapping_matrix_ids_gpi_master(j,i);
					double R_row = R_gpi_master(j,i);

					for(unsigned int k=0; k<mapping_matrix_ids_gpi_master.size2();k++)
					{
						for(unsigned int l=0; l<mapping_matrix_ids_gpi_master.size1();l++)
						{
							unsigned int R_coll_id = mapping_matrix_ids_gpi_master(l,k);
							double R_coll = R_gpi_master(l,k);

							m_mapping_matrix_CAD_CAD(3*R_row_id+0,3*R_coll_id+0) += penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
							m_mapping_matrix_CAD_CAD(3*R_row_id+1,3*R_coll_id+1) += penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
							m_mapping_matrix_CAD_CAD(3*R_row_id+2,3*R_coll_id+2) += penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
						}
					}
				}
			}

			// Then we consider the relation Slave-Slave ( SS )
			for(unsigned int i=0; i<mapping_matrix_ids_gpi_slave.size2();i++)
			{
				for(unsigned int j=0; j<mapping_matrix_ids_gpi_slave.size1();j++)
				{
					unsigned int R_row_id = mapping_matrix_ids_gpi_slave(j,i);
					double R_row = R_gpi_slave(j,i);

					for(unsigned int k=0; k<mapping_matrix_ids_gpi_slave.size2();k++)
					{
						for(unsigned int l=0; l<mapping_matrix_ids_gpi_slave.size1();l++)
						{
							unsigned int R_coll_id = mapping_matrix_ids_gpi_slave(l,k);
							double R_coll = R_gpi_slave(l,k);

							m_mapping_matrix_CAD_CAD(3*R_row_id+0,3*R_coll_id+0) += penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
							m_mapping_matrix_CAD_CAD(3*R_row_id+1,3*R_coll_id+1) += penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
							m_mapping_matrix_CAD_CAD(3*R_row_id+2,3*R_coll_id+2) += penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
						}
					}
				}
			}			

			// Then we consider the Master-Slave relation ( MS & SM )
			for(unsigned int i=0; i<mapping_matrix_ids_gpi_master.size2();i++)
			{
				for(unsigned int j=0; j<mapping_matrix_ids_gpi_master.size1();j++)
				{
					unsigned int R_row_id = mapping_matrix_ids_gpi_master(j,i);
					double R_row = R_gpi_master(j,i);

					for(unsigned int k=0; k<mapping_matrix_ids_gpi_slave.size2();k++)
					{
						for(unsigned int l=0; l<mapping_matrix_ids_gpi_slave.size1();l++)
						{
							unsigned int R_coll_id = mapping_matrix_ids_gpi_slave(l,k);
							double R_coll = R_gpi_slave(l,k);

							// MS 
							m_mapping_matrix_CAD_CAD(3*R_row_id+0,3*R_coll_id+0) -= penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
							m_mapping_matrix_CAD_CAD(3*R_row_id+1,3*R_coll_id+1) -= penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
							m_mapping_matrix_CAD_CAD(3*R_row_id+2,3*R_coll_id+2) -= penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;		

							// SM
							m_mapping_matrix_CAD_CAD(3*R_coll_id+0,3*R_row_id+0) -= penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
							m_mapping_matrix_CAD_CAD(3*R_coll_id+1,3*R_row_id+1) -= penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
							m_mapping_matrix_CAD_CAD(3*R_coll_id+2,3*R_row_id+2) -= penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;							
						}
					}
				}
			}
		}

		// --------------------------------------------------------------------------
		void apply_rotation_coupling( Patch::Pointer master_patch,
									Patch::Pointer slave_patch,
									double u_m, double v_m,
									double u_s, double v_s,
									Vector &tangent_on_master_patch,
									Vector &tangent_on_slave_patch,
									matrix<unsigned int> &mapping_matrix_ids_gpi_master,
									matrix<unsigned int> &mapping_matrix_ids_gpi_slave,
									double J1,
									double gp_i_weight,
									double penalty_factor_rot )
		{		
			// Variables needed later
			Vector T1_m, T1_s, T2_m, T2_s, T3_m, T3_s;
			std::vector<Vector> t1r_m, t1r_s, t2r_m, t2r_s, t3r_m, t3r_s;				

			// Compute geometric quantities
			master_patch->GetSurface().ComputeVariationOfLocalCSY( u_m, v_m, tangent_on_master_patch, T1_m, T2_m, T3_m, t1r_m, t2r_m, t3r_m );
			slave_patch->GetSurface().ComputeVariationOfLocalCSY( u_s, v_s, tangent_on_slave_patch, T1_s, T2_s, T3_s, t1r_s, t2r_s, t3r_s );

			// Check if master and slave tangent point in same direction. If yes, we have to subtract in the following.
			int sign_factor = 1;
			if( inner_prod(T2_m,T2_s) > 0 )
				sign_factor = -1;

			// Merge boundary conditions into mapping matrix

			// First we consider the relation Master-Master ( MM )
			unsigned int k_coll = 0;
			for(unsigned int i=0; i<mapping_matrix_ids_gpi_master.size2();i++)
			{
				for(unsigned int j=0; j<mapping_matrix_ids_gpi_master.size1();j++)
				{
					unsigned int R_row_id = mapping_matrix_ids_gpi_master(j,i);
					Vector omega_mx_coll = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_coll+0]);
					Vector omega_my_coll = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_coll+1]);
					Vector omega_mz_coll = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_coll+2]);
					double omega_T2_mx_coll = inner_prod(omega_mx_coll,T2_m);
					double omega_T2_my_coll = inner_prod(omega_my_coll,T2_m);
					double omega_T2_mz_coll = inner_prod(omega_mz_coll,T2_m);

					unsigned int k_row = 0;
					for(unsigned int k=0; k<mapping_matrix_ids_gpi_master.size2();k++)
					{
						for(unsigned int l=0; l<mapping_matrix_ids_gpi_master.size1();l++)
						{
							unsigned int R_coll_id = mapping_matrix_ids_gpi_master(l,k);
							Vector omega_mx_row = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_row+0]);
							Vector omega_my_row = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_row+1]);
							Vector omega_mz_row = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_row+2]);
							double omega_T2_mx_row = inner_prod(omega_mx_row,T2_m);
							double omega_T2_my_row = inner_prod(omega_my_row,T2_m);
							double omega_T2_mz_row = inner_prod(omega_mz_row,T2_m);

							m_mapping_matrix_CAD_CAD(3*R_row_id+0,3*R_coll_id+0) += penalty_factor_rot * gp_i_weight * J1 * omega_T2_mx_row * omega_T2_mx_coll;
							m_mapping_matrix_CAD_CAD(3*R_row_id+1,3*R_coll_id+1) += penalty_factor_rot * gp_i_weight * J1 * omega_T2_my_row * omega_T2_my_coll;
							m_mapping_matrix_CAD_CAD(3*R_row_id+2,3*R_coll_id+2) += penalty_factor_rot * gp_i_weight * J1 * omega_T2_mz_row * omega_T2_mz_coll;
							
							k_row++;
						}
					}
					k_coll++;
				}
			}

			// Then we consider the relation Slave-Slave ( SS )
			k_coll = 0;
			for(unsigned int i=0; i<mapping_matrix_ids_gpi_slave.size2();i++)
			{
				for(unsigned int j=0; j<mapping_matrix_ids_gpi_slave.size1();j++)
				{
					unsigned int R_row_id = mapping_matrix_ids_gpi_slave(j,i);
					Vector omega_sx_coll = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_coll+0]);
					Vector omega_sy_coll = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_coll+1]);
					Vector omega_sz_coll = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_coll+2]);
					double omega_T2_sx_coll = inner_prod(omega_sx_coll,T2_s);
					double omega_T2_sy_coll = inner_prod(omega_sy_coll,T2_s);
					double omega_T2_sz_coll = inner_prod(omega_sz_coll,T2_s);

					unsigned int k_row = 0;
					for(unsigned int k=0; k<mapping_matrix_ids_gpi_slave.size2();k++)
					{
						for(unsigned int l=0; l<mapping_matrix_ids_gpi_slave.size1();l++)
						{
							unsigned int R_coll_id = mapping_matrix_ids_gpi_slave(l,k);
							Vector omega_sx_row = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_row+0]);
							Vector omega_sy_row = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_row+1]);
							Vector omega_sz_row = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_row+2]);
							double omega_T2_sx_row = inner_prod(omega_sx_row,T2_s);
							double omega_T2_sy_row = inner_prod(omega_sy_row,T2_s);
							double omega_T2_sz_row = inner_prod(omega_sz_row,T2_s);

							m_mapping_matrix_CAD_CAD(3*R_row_id+0,3*R_coll_id+0) += penalty_factor_rot * gp_i_weight * J1 * omega_T2_sx_row * omega_T2_sx_coll;
							m_mapping_matrix_CAD_CAD(3*R_row_id+1,3*R_coll_id+1) += penalty_factor_rot * gp_i_weight * J1 * omega_T2_sy_row * omega_T2_sy_coll;
							m_mapping_matrix_CAD_CAD(3*R_row_id+2,3*R_coll_id+2) += penalty_factor_rot * gp_i_weight * J1 * omega_T2_sz_row * omega_T2_sz_coll;
							
							k_row++;
						}
					}
					k_coll++;
				}
			}			

			// Then we consider the Master-slave relation ( MS & SM )
			unsigned int k_m = 0;
			for(unsigned int i=0; i<mapping_matrix_ids_gpi_master.size2();i++)
			{
				for(unsigned int j=0; j<mapping_matrix_ids_gpi_master.size1();j++)
				{
					unsigned int R_m_id = mapping_matrix_ids_gpi_master(j,i);
					Vector omega_mx = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_m+0]);
					Vector omega_my = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_m+1]);
					Vector omega_mz = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_m+2]);
					double omega_T2_mx = inner_prod(omega_mx,T2_m);
					double omega_T2_my = inner_prod(omega_my,T2_m);
					double omega_T2_mz = inner_prod(omega_mz,T2_m);

					unsigned int k_s = 0;
					for(unsigned int k=0; k<mapping_matrix_ids_gpi_slave.size2();k++)
					{
						for(unsigned int l=0; l<mapping_matrix_ids_gpi_slave.size1();l++)
						{
							unsigned int R_s_id = mapping_matrix_ids_gpi_slave(l,k);
							Vector omega_sx = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_s+0]);
							Vector omega_sy = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_s+1]);
							Vector omega_sz = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_s+2]);
							double omega_T2_sx = inner_prod(omega_sx,T2_s);
							double omega_T2_sy = inner_prod(omega_sy,T2_s);
							double omega_T2_sz = inner_prod(omega_sz,T2_s);

							// MS
							m_mapping_matrix_CAD_CAD(3*R_m_id+0,3*R_s_id+0) += sign_factor * penalty_factor_rot * gp_i_weight * J1 * omega_T2_mx * omega_T2_sx;
							m_mapping_matrix_CAD_CAD(3*R_m_id+1,3*R_s_id+1) += sign_factor * penalty_factor_rot * gp_i_weight * J1 * omega_T2_my * omega_T2_sy;
							m_mapping_matrix_CAD_CAD(3*R_m_id+2,3*R_s_id+2) += sign_factor * penalty_factor_rot * gp_i_weight * J1 * omega_T2_mz * omega_T2_sz;

							// SM
							m_mapping_matrix_CAD_CAD(3*R_s_id+0,3*R_m_id+0) += sign_factor * penalty_factor_rot * gp_i_weight * J1 * omega_T2_mx * omega_T2_sx;
							m_mapping_matrix_CAD_CAD(3*R_s_id+1,3*R_m_id+1) += sign_factor * penalty_factor_rot * gp_i_weight * J1 * omega_T2_my * omega_T2_sy;
							m_mapping_matrix_CAD_CAD(3*R_s_id+2,3*R_m_id+2) += sign_factor * penalty_factor_rot * gp_i_weight * J1 * omega_T2_mz * omega_T2_sz;						

							k_s++;								
						}
					}
					k_m++;
				}
			}
		}

		// --------------------------------------------------------------------------
		void enforce_tangent_continuity( Patch::Pointer master_patch,
										Patch::Pointer slave_patch,
										double u_m, double v_m,
										double u_s, double v_s,
										Vector& tangent_on_master_patch,
										Vector& tangent_on_slave_patch,
										matrix<unsigned int>& mapping_matrix_ids_gpi_master,
										matrix<unsigned int>& mapping_matrix_ids_gpi_slave,
										double J1,
										double gp_i_weight,
										double penalty_factor_tangent_continuity )
		{
			// Variables needed later
			Vector T1_m, T1_s, T2_m, T2_s, T3_m, T3_s;
			Vector T1_der_m, T1_der_s, T2_der_m, T2_der_s, T3_der_m, T3_der_s;
			std::vector<Vector> t1r_m, t1r_s, t2r_m, t2r_s, t3r_m, t3r_s;
			std::vector<Vector> t1_der_r_m, t1_der_r_s, t2_der_r_m, t2_der_r_s, t3_der_r_m, t3_der_r_s;					
			std::vector<std::vector<Vector>> t1rs_m, t1rs_s, t2rs_m, t2rs_s, t3rs_m, t3rs_s;
			std::vector<std::vector<Vector>> t1_der_rs_m, t1_der_rs_s, t2_der_rs_m, t2_der_rs_s, t3_der_rs_m, t3_der_rs_s;

			std::cout << "Called: cad_mapper::enforce_tangent_continuity()" << std::endl;

			// Compute geometric quantities
			master_patch->GetSurface().ComputeSecondVariationOfLocalCSY( u_m, v_m, 
																		tangent_on_master_patch, 
																		T1_m, T2_m, T3_m, 
																		T1_der_m, T2_der_m, T3_der_m,
																		t1r_m, t2r_m, t3r_m,
																		t1_der_r_m, t2_der_r_m, t3_der_r_m,
																		t1rs_m, t2rs_m, t3rs_m,
																		t1_der_rs_m, t2_der_rs_m, t3_der_rs_m );
			slave_patch->GetSurface().ComputeSecondVariationOfLocalCSY( u_s, v_s, 
																	tangent_on_slave_patch, 
																	T1_s, T2_s, T3_s, 
																	T1_der_s, T2_der_s, T3_der_s,
																	t1r_s, t2r_s, t3r_s,
																	t1_der_r_s, t2_der_r_s, t3_der_r_s,
																	t1rs_s, t2rs_s, t3rs_s,
																	t1_der_rs_s, t2_der_rs_s, t3_der_rs_s );	

			
			double fac = inner_prod(T3_m, T1_s);
			KRATOS_WATCH(fac);

			// First we consider contribution to the m_mapping_rhs_vector

			// Master-Master-relation ( MM )
			unsigned int k_row = 0;
			for(unsigned int i=0; i<mapping_matrix_ids_gpi_master.size2();i++)
			{
				for(unsigned int j=0; j<mapping_matrix_ids_gpi_master.size1();j++)
				{
					unsigned int R_row_id = mapping_matrix_ids_gpi_master(j,i);

					m_mapping_rhs_vector(3*R_row_id+0) += penalty_factor_tangent_continuity * gp_i_weight * J1 * fac * inner_prod(t3r_m[3*k_row+0],T1_s);
					m_mapping_rhs_vector(3*R_row_id+1) += penalty_factor_tangent_continuity * gp_i_weight * J1 * fac * inner_prod(t3r_m[3*k_row+1],T1_s);
					m_mapping_rhs_vector(3*R_row_id+2) += penalty_factor_tangent_continuity * gp_i_weight * J1 * fac * inner_prod(t3r_m[3*k_row+2],T1_s);

					k_row++;
				}
			}

			// Slave-Slave-relation ( SS )
			k_row = 0;
			for(unsigned int i=0; i<mapping_matrix_ids_gpi_slave.size2();i++)
			{
				for(unsigned int j=0; j<mapping_matrix_ids_gpi_slave.size1();j++)
				{
					unsigned int R_row_id = mapping_matrix_ids_gpi_slave(j,i);

					m_mapping_rhs_vector(3*R_row_id+0) += penalty_factor_tangent_continuity * gp_i_weight * J1 * fac * inner_prod(t1r_s[3*k_row+0],T3_m);
					m_mapping_rhs_vector(3*R_row_id+1) += penalty_factor_tangent_continuity * gp_i_weight * J1 * fac * inner_prod(t1r_s[3*k_row+1],T3_m);
					m_mapping_rhs_vector(3*R_row_id+2) += penalty_factor_tangent_continuity * gp_i_weight * J1 * fac * inner_prod(t1r_s[3*k_row+2],T3_m);

					k_row++;
				}
			}	

			// Then we consider the contribution to the m_mapping_matrix_CAD_CAD

			// Master-Master-relation ( MM )
			k_row = 0;
			for(unsigned int i=0; i<mapping_matrix_ids_gpi_master.size2();i++)
			{
				for(unsigned int j=0; j<mapping_matrix_ids_gpi_master.size1();j++)
				{
					unsigned int R_row_id = mapping_matrix_ids_gpi_master(j,i);

					unsigned int k_coll = 0;
					for(unsigned int k=0; k<mapping_matrix_ids_gpi_master.size2();k++)
					{
						for(unsigned int l=0; l<mapping_matrix_ids_gpi_master.size1();l++)
						{
							unsigned int R_coll_id = mapping_matrix_ids_gpi_master(l,k);
							
							double term_1_x = inner_prod(t3r_m[3*k_coll+0],T1_s) * inner_prod(t3r_m[3*k_row+0],T1_s);
							double term_1_y = inner_prod(t3r_m[3*k_coll+1],T1_s) * inner_prod(t3r_m[3*k_row+1],T1_s);
							double term_1_z = inner_prod(t3r_m[3*k_coll+2],T1_s) * inner_prod(t3r_m[3*k_row+2],T1_s);

							double term_2_x = fac * inner_prod(t3rs_m[3*k_row+0][3*k_coll+0],T1_s);
							double term_2_y = fac * inner_prod(t3rs_m[3*k_row+1][3*k_coll+1],T1_s);
							double term_2_z = fac * inner_prod(t3rs_m[3*k_row+2][3*k_coll+2],T1_s);

							m_mapping_matrix_CAD_CAD(3*R_row_id+0,3*R_coll_id+0) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_x + term_2_x );
							m_mapping_matrix_CAD_CAD(3*R_row_id+1,3*R_coll_id+1) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_y + term_2_y );
							m_mapping_matrix_CAD_CAD(3*R_row_id+2,3*R_coll_id+2) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_z + term_2_z );					

							k_coll++;								
						}
					}
					k_row++;
				}
			}

			// Slave-Slave-relation ( SS )
			k_row = 0;
			for(unsigned int i=0; i<mapping_matrix_ids_gpi_slave.size2();i++)
			{
				for(unsigned int j=0; j<mapping_matrix_ids_gpi_slave.size1();j++)
				{
					unsigned int R_row_id = mapping_matrix_ids_gpi_slave(j,i);

					unsigned int k_coll = 0;
					for(unsigned int k=0; k<mapping_matrix_ids_gpi_slave.size2();k++)
					{
						for(unsigned int l=0; l<mapping_matrix_ids_gpi_slave.size1();l++)
						{
							unsigned int R_coll_id = mapping_matrix_ids_gpi_slave(l,k);
							
							double term_1_x = inner_prod(t1r_s[3*k_coll+0],T3_m) * inner_prod(T3_m, t1r_s[3*k_row+0]);
							double term_1_y = inner_prod(t1r_s[3*k_coll+1],T3_m) * inner_prod(T3_m, t1r_s[3*k_row+1]);
							double term_1_z = inner_prod(t1r_s[3*k_coll+2],T3_m) * inner_prod(T3_m, t1r_s[3*k_row+2]);

							double term_2_x = fac * inner_prod(T3_m, t1rs_s[3*k_row+0][3*k_coll+0]);
							double term_2_y = fac * inner_prod(T3_m, t1rs_s[3*k_row+1][3*k_coll+1]);
							double term_2_z = fac * inner_prod(T3_m, t1rs_s[3*k_row+2][3*k_coll+2]);

							m_mapping_matrix_CAD_CAD(3*R_row_id+0,3*R_coll_id+0) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_x + term_2_x );
							m_mapping_matrix_CAD_CAD(3*R_row_id+1,3*R_coll_id+1) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_y + term_2_y );
							m_mapping_matrix_CAD_CAD(3*R_row_id+2,3*R_coll_id+2) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_z + term_2_z );					

							k_coll++;								
						}
					}
					k_row++;
				}
			}	

			// Master-slave-relation ( MS )
			k_row = 0;
			for(unsigned int i=0; i<mapping_matrix_ids_gpi_master.size2();i++)
			{
				for(unsigned int j=0; j<mapping_matrix_ids_gpi_master.size1();j++)
				{
					unsigned int R_row_id = mapping_matrix_ids_gpi_master(j,i);

					unsigned int k_coll = 0;
					for(unsigned int k=0; k<mapping_matrix_ids_gpi_slave.size2();k++)
					{
						for(unsigned int l=0; l<mapping_matrix_ids_gpi_slave.size1();l++)
						{
							unsigned int R_coll_id = mapping_matrix_ids_gpi_slave(l,k);

							double term_1_x = inner_prod(T3_m,t1r_s[3*k_coll+0]) * inner_prod(t3r_m[3*k_row+0],T1_s);
							double term_1_y = inner_prod(T3_m,t1r_s[3*k_coll+1]) * inner_prod(t3r_m[3*k_row+1],T1_s);
							double term_1_z = inner_prod(T3_m,t1r_s[3*k_coll+2]) * inner_prod(t3r_m[3*k_row+2],T1_s);

							double term_2_x = fac * inner_prod(t3r_m[3*k_row+0],t1r_s[3*k_coll+0]);
							double term_2_y = fac * inner_prod(t3r_m[3*k_row+1],t1r_s[3*k_coll+1]);
							double term_2_z = fac * inner_prod(t3r_m[3*k_row+2],t1r_s[3*k_coll+2]);

							m_mapping_matrix_CAD_CAD(3*R_row_id+0,3*R_coll_id+0) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_x + term_2_x );
							m_mapping_matrix_CAD_CAD(3*R_row_id+1,3*R_coll_id+1) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_y + term_2_y );
							m_mapping_matrix_CAD_CAD(3*R_row_id+2,3*R_coll_id+2) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_z + term_2_z );											

							k_coll++;								
						}
					}
					k_row++;
				}
			}

			// Master-slave-relation ( SM )
			k_row = 0;
			for(unsigned int i=0; i<mapping_matrix_ids_gpi_slave.size2();i++)
			{
				for(unsigned int j=0; j<mapping_matrix_ids_gpi_slave.size1();j++)
				{
					unsigned int R_row_id = mapping_matrix_ids_gpi_slave(j,i);

					unsigned int k_coll = 0;
					for(unsigned int k=0; k<mapping_matrix_ids_gpi_master.size2();k++)
					{
						for(unsigned int l=0; l<mapping_matrix_ids_gpi_master.size1();l++)
						{
							unsigned int R_coll_id = mapping_matrix_ids_gpi_master(l,k);

							double term_1_y = inner_prod(t3r_m[3*k_coll+1], T1_s) * inner_prod(T3_m, t1r_s[3*k_row+1]);
							double term_1_z = inner_prod(t3r_m[3*k_coll+2], T1_s) * inner_prod(T3_m, t1r_s[3*k_row+2]);
							double term_1_x = inner_prod(t3r_m[3*k_coll+0], T1_s) * inner_prod(T3_m, t1r_s[3*k_row+0]);

							double term_2_x = fac * inner_prod(t3r_m[3*k_coll+0], t1r_s[3*k_row+0]);
							double term_2_y = fac * inner_prod(t3r_m[3*k_coll+1], t1r_s[3*k_row+1]);
							double term_2_z = fac * inner_prod(t3r_m[3*k_coll+2], t1r_s[3*k_row+2]);

							m_mapping_matrix_CAD_CAD(3*R_row_id+0,3*R_coll_id+0) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_x + term_2_x );
							m_mapping_matrix_CAD_CAD(3*R_row_id+1,3*R_coll_id+1) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_y + term_2_y );
							m_mapping_matrix_CAD_CAD(3*R_row_id+2,3*R_coll_id+2) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_z + term_2_z );											

							k_coll++;								
						}
					}
					k_row++;
				}
			}		
		}

		// --------------------------------------------------------------------------
		void apply_dirichlet_condition( BREPElementVector::iterator brep_elem_i, 
										double penalty_factor_dirichlet,
										boost::python::list& edges_with_specific_dirichlet_conditions )
		{
			// Get Gauss points of current brep element
			BREPGaussPointVector brep_gps = brep_elem_i->GetGaussPoints();

			// Check if for current element some specific dirichlet condition is to be defined
			bool fix_x = true;
			bool fix_y = true;
			bool fix_z = true;
			for (unsigned int i = 0; i < boost::python::len(edges_with_specific_dirichlet_conditions); ++i)
			{
				unsigned int listed_edge_id = extractUnsignedInt(edges_with_specific_dirichlet_conditions[i][0]);
				if(brep_elem_i->GetEdgeId() == listed_edge_id)
				{
					fix_x = extractBool(edges_with_specific_dirichlet_conditions[i][1][0]);
					fix_y = extractBool(edges_with_specific_dirichlet_conditions[i][1][1]);
					fix_z = extractBool(edges_with_specific_dirichlet_conditions[i][1][2]);
				}
			}

			// Loop over all Gauss points of current brep element 
			for (BREPGaussPointVector::iterator brep_gp_i = brep_gps.begin(); brep_gp_i != brep_gps.end(); ++brep_gp_i)
			{
				// Read information from Gauss point
				unsigned int master_patch_id = brep_gp_i->GetPatchId();
				Patch::Pointer master_patch = m_patches[m_patch_position_in_patch_vector[master_patch_id]];
				double gp_i_weight = brep_gp_i->GetWeight();
				Vector location_on_master_patch = brep_gp_i->GetLocation();
				Vector tangent_on_master_patch = brep_gp_i->GetTangent();

				// Evaluate NURBS basis function for Gauss point on both patches and get the corresponding ids of control points in the mapping matrix
				matrix<double> R_gpi_master;
				double u_m = location_on_master_patch(0);
				double v_m = location_on_master_patch(1);
				master_patch->GetSurface().EvaluateNURBSFunctions(-1,-1,u_m, v_m, R_gpi_master);
				matrix<unsigned int> mapping_matrix_ids_gpi_master = master_patch->GetSurface().GetMappingMatrixIds(-1,-1,u_m, v_m);						

				// Compute Jacobian J1
				matrix<double> g_master = master_patch->GetSurface().GetBaseVectors(-1,-1,u_m,v_m);
				Vector g1 = ZeroVector(3);
				g1(0) = g_master(0,0);
				g1(1) = g_master(1,0);
				g1(2) = g_master(2,0);
				Vector g2 = ZeroVector(3);
				g2(0) = g_master(0,1);
				g2(1) = g_master(1,1);
				g2(2) = g_master(2,1);
				double J1 = norm_2( g1* tangent_on_master_patch(0) + g2* tangent_on_master_patch(1) );

				// Merge boundary condition into mapping matrix ( Note we have only a Master-Master relation, so N_m*N_m)
				for(unsigned int i=0; i<mapping_matrix_ids_gpi_master.size1();i++)
				{
					for(unsigned int j=0; j<mapping_matrix_ids_gpi_master.size2();j++)
					{
						unsigned int R_row_id = mapping_matrix_ids_gpi_master(i,j);
						double R_row = R_gpi_master(i,j);

						for(unsigned int k=0; k<mapping_matrix_ids_gpi_master.size1();k++)
						{
							for(unsigned int l=0; l<mapping_matrix_ids_gpi_master.size2();l++)
							{
								unsigned int R_coll_id = mapping_matrix_ids_gpi_master(k,l);
								double R_coll = R_gpi_master(k,l);

								if(fix_x)
									m_mapping_matrix_CAD_CAD(3*R_row_id+0,3*R_coll_id+0) += penalty_factor_dirichlet * gp_i_weight * J1 * R_row * R_coll;
								if(fix_y)								
									m_mapping_matrix_CAD_CAD(3*R_row_id+1,3*R_coll_id+1) += penalty_factor_dirichlet * gp_i_weight * J1 * R_row * R_coll;
								if(fix_z)
									m_mapping_matrix_CAD_CAD(3*R_row_id+2,3*R_coll_id+2) += penalty_factor_dirichlet * gp_i_weight * J1 * R_row * R_coll;
							}
						}
					}
				}
			}		
		}

		// --------------------------------------------------------------------------
		void map_to_cad_space()
		{
			std::cout << "\n> Starting to map to CAD space..." << std::endl;
			boost::timer function_timer;

			// Check for validity of mapping matrix
			for(unsigned int i=0; i<m_mapping_matrix_CAD_CAD.size1();i++)
				if(std::abs(m_mapping_matrix_CAD_CAD(i,i))<1e-10)
				{
					std::cout << "\nWARNING,small value on main diagonal of mapping matrix !!!! " <<std::endl;
					std::cout << "Value = " << m_mapping_matrix_CAD_CAD(i,i) << std::endl;
					std::cout << "Iterator i = " << i << std::endl;
					m_mapping_matrix_CAD_CAD(i,i) = 1e-3;
					// KRATOS_THROW_ERROR(std::runtime_error, "Zero on the main diagonal of the mapping matrix!!!!!!!", m_mapping_matrix_CAD_CAD(i,i));	
				}

			// Initialize vectors needed later
			Vector dx = ZeroVector(3*m_n_relevant_fem_points);
			Vector ds = ZeroVector(3*m_n_relevant_control_points);

			// Prepare RHS vector of mapping system of equation
			for (ModelPart::NodesContainerType::iterator node_i = mr_fe_model_part.NodesBegin(); node_i != mr_fe_model_part.NodesEnd(); ++node_i)
			{
				unsigned int mapping_id = node_i->GetValue(MAPPING_MATRIX_ID);
				dx[3*mapping_id+0] = node_i->GetValue(SHAPE_CHANGE_ABSOLUTE_X);
				dx[3*mapping_id+1] = node_i->GetValue(SHAPE_CHANGE_ABSOLUTE_Y);
				dx[3*mapping_id+2] = node_i->GetValue(SHAPE_CHANGE_ABSOLUTE_Z);
			}
			noalias(m_mapping_rhs_vector) += prod(m_mapping_matrix_CAD_FEM,dx);

			// Assign sparse matrix to compressed matrix required by linear solver
			CompressedMatrixType mapping_matrix_CAD_CAD = m_mapping_matrix_CAD_CAD;

			// Solve linear systems to obtain mapped quantities each in X,Y,Z-direction separately
			// Note that an alternative would be to solve a big block structured matrix at once
			m_linear_solver->Solve(mapping_matrix_CAD_CAD, ds, m_mapping_rhs_vector);
			
			// Update solution (displacement of control points) in cad data set (both in c++ and python)
			m_cad_reader.UpdateControlPoints(m_patches, ds);

			// Test solution
			// KRATOS_WATCH(ds);
			Vector rhs_test = ZeroVector(3*m_n_relevant_control_points);
			noalias(rhs_test) = prod(m_mapping_matrix_CAD_CAD,ds);
			Vector rhs_difference = m_mapping_rhs_vector - rhs_test;
			double normalized_difference_in_rhs = norm_2(rhs_difference);
			std::cout << "\n> Solution of linear system leads to a difference in the RHS of: normalized_difference_in_rhs = " << normalized_difference_in_rhs << std::endl;

			std::cout << "\n> Mapping to CAD space finished in " << function_timer.elapsed() << " s." << std::endl;
		}	

		// --------------------------------------------------------------------------
		void check_c0_continuity( VectorPoint myMaster, VectorPoint mySlave, double &average, double &max )
		{
			if( myMaster.size() != mySlave.size() )
			{
				KRATOS_WATCH(" Size different ");
				average = NAN;
				max = NAN;
			}else
			{
				DoubleVector distance;
				for( size_t i = 0; i < myMaster.size(); i++)
				{
					double X = myMaster[i].X() - mySlave[i].X();
					double Y = myMaster[i].Y() - mySlave[i].Y();
					double Z = myMaster[i].Z() - mySlave[i].Z();
					double eulerian_distance = sqrt( X*X + Y*Y + Z*Z );

					distance.push_back( eulerian_distance );
				}

				average = 0;
				max = distance[0];
				for( size_t i = 0; i < distance.size( ); i++)
				{
					average = average + distance[i];
					if( distance[i] > max ) 
					{
						max = distance[i];
					}
				}

				average = average / distance.size();
			}
		}
	////////////////////////////////////////////////////////

	// Output functions #########################################################
		// --------------------------------------------------------------------------
		void output_gauss_points(std::string output_filename)
		{
			std::cout << "\n> Starting writing gauss points of given FEM mesh..." << std::endl;
		
			unsigned int gauss_point_counter = 0;
			NodeVector list_of_gauss_points;

			// Loop over all integration points of fe-model-part
			for (ModelPart::ConditionsContainerType::iterator cond_i = mr_fe_model_part.ConditionsBegin(); cond_i != mr_fe_model_part.ConditionsEnd(); ++cond_i)
			{
				// Get geometry information of current condition
				Condition::GeometryType& geom_i = cond_i->GetGeometry();

				// Evaluate shape functions of FE model according specified integration methdod
				const Condition::GeometryType::IntegrationPointsArrayType& integration_points = geom_i.IntegrationPoints(m_integration_method);
				const unsigned int number_of_integration_points = integration_points.size();

				for ( unsigned int PointNumber = 0; PointNumber < number_of_integration_points; PointNumber++ )
				{
					// Compute global coordinates of current integration point and get corresponding weight
					NodeType::CoordinatesArrayType ip_coordinates = geom_i.GlobalCoordinates(ip_coordinates, integration_points[PointNumber].Coordinates());
					NodeType::Pointer gauss_point_i = Node < 3 > ::Pointer(new Node<3>(++gauss_point_counter, ip_coordinates ));

					list_of_gauss_points.push_back(gauss_point_i);
				}
			}

			// Write points specified in list above to the given file
			std::ofstream temp_file(output_filename);
			for(NodeVector::iterator it =  list_of_gauss_points.begin(); it!=list_of_gauss_points.end(); it++)
			{
				NodeType::Pointer gp_i = *it;
				temp_file << gp_i->X() << " " << gp_i->Y() << " " << gp_i->Z() << std::endl;
			}
			temp_file.close();

			std::cout << "\n> Finished writing gauss points of given FEM mesh..." << std::endl;
		}		

		// --------------------------------------------------------------------------
		void output_surface_points(std::string output_filename, const unsigned int u_resolution, const unsigned int v_resolution, int specific_patch )
		{
			std::cout << "\n> Writing to file some surface points..." << std::endl;
		
			// Set a max value for the coordinte output to avoid clutter through points that are plotted outside the trimmed surface
			double max_coordinate = 1000;

			// Open file to write all points in
			std::ofstream file_to_write(output_filename);

			// If no specific patch defined, output all
			if(specific_patch<0)
			{
				//Loop over all patches
				for (unsigned int patch_itr = 0; patch_itr < m_patches.size(); patch_itr++)
				{
					// Get relevant data
					DoubleVector& knot_vec_u_i = m_patches[patch_itr]->GetSurface().GetKnotVectorU();
					DoubleVector& knot_vec_v_i = m_patches[patch_itr]->GetSurface().GetKnotVectorV();

					// Pre-calculations
					unsigned int knot_vector_u_dimension = knot_vec_u_i.size();
					unsigned int knot_vector_v_dimension = knot_vec_v_i.size();

					double u_min = knot_vec_u_i[0];
					double u_max = knot_vec_u_i[knot_vector_u_dimension-1];
					double v_min = knot_vec_v_i[0];
					double v_max = knot_vec_v_i[knot_vector_v_dimension-1];

					double delta_u = (u_max-u_min) / u_resolution;
					double delta_v = (v_max-v_min) / v_resolution;

					// Loop over all u & v according to specified resolution
					for(unsigned int i=0; i<=u_resolution; i++)
					{
						// current u-value
						double u_i = u_min + i*delta_u;

						for(unsigned int j=0; j<=v_resolution; j++)
						{
							// current v-value
							double v_j = v_min + j*delta_v;

							// Check if u_i and v_j represent a point inside the closed boundary loop
							array_1d<double, 2> point_of_interest;
							point_of_interest[0] = u_i;
							point_of_interest[1] = v_j;
							bool point_is_inside = m_patches[patch_itr]->CheckIfPointIsInside(point_of_interest);

							if(point_is_inside)
							{
								// compute point in CAD-model for given u&v				
								Point<3> cad_point;
								m_patches[patch_itr]->GetSurface().EvaluateSurfacePoint(cad_point, u_i, v_j);

								// Check and set a max value for irrelevant points outside the trimmed surface
								if(std::abs(cad_point.X())>max_coordinate)
									cad_point.X() = MathUtils<int>::Sign(cad_point.X()) * max_coordinate;
								if(std::abs(cad_point.Y())>max_coordinate)
									cad_point.Y() = MathUtils<int>::Sign(cad_point.Y()) * max_coordinate;
								if(std::abs(cad_point.Z())>max_coordinate)
									cad_point.Z() = MathUtils<int>::Sign(cad_point.Z()) * max_coordinate;														

								// Output point
								file_to_write << cad_point.X() << " " << cad_point.Y() << " " << cad_point.Z() << std::endl;
							}
						}
					}
				}
			}
			else // if specific patch defined
			{
				// Get relevant data
				DoubleVector& knot_vec_u_i = m_patches[specific_patch]->GetSurface().GetKnotVectorU();
				DoubleVector& knot_vec_v_i = m_patches[specific_patch]->GetSurface().GetKnotVectorV();

				// Pre-calculations
				unsigned int knot_vector_u_dimension = knot_vec_u_i.size();
				unsigned int knot_vector_v_dimension = knot_vec_v_i.size();

				double u_min = knot_vec_u_i[0];
				double u_max = knot_vec_u_i[knot_vector_u_dimension-1];
				double v_min = knot_vec_v_i[0];
				double v_max = knot_vec_v_i[knot_vector_v_dimension-1];

				double delta_u = (u_max-u_min) / u_resolution;
				double delta_v = (v_max-v_min) / v_resolution;

				// Loop over all u & v according to specified resolution
				for(unsigned int i=0; i<=u_resolution; i++)
				{
					// current u-value
					double u_i = u_min + i*delta_u;

					for(unsigned int j=0; j<=v_resolution; j++)
					{
						// current v-value
						double v_j = v_min + j*delta_v;

						// Check if u_i and v_j represent a point inside the closed boundary loop
						array_1d<double, 2> point_of_interest;
						point_of_interest[0] = u_i;
						point_of_interest[1] = v_j;
						bool point_is_inside = m_patches[specific_patch]->CheckIfPointIsInside(point_of_interest);

						if(point_is_inside)
						{
							// compute point in CAD-model for given u&v				
							Point<3> cad_point;
							m_patches[specific_patch]->GetSurface().EvaluateSurfacePoint(cad_point, u_i, v_j);

							// Check and set a max value for irrelevant points outside the trimmed surface
							if(std::abs(cad_point.X())>max_coordinate)
								cad_point.X() = MathUtils<int>::Sign(cad_point.X()) * max_coordinate;
							if(std::abs(cad_point.Y())>max_coordinate)
								cad_point.Y() = MathUtils<int>::Sign(cad_point.Y()) * max_coordinate;
							if(std::abs(cad_point.Z())>max_coordinate)
								cad_point.Z() = MathUtils<int>::Sign(cad_point.Z()) * max_coordinate;														

							// Output point
							file_to_write << cad_point.X() << " " << cad_point.Y() << " " << cad_point.Z() << std::endl;
						}
					}
				}
			}

			// Close file
			file_to_write.close();

			std::cout << "\t\t\t DONE" << std::endl;
		}	

		// --------------------------------------------------------------------------
		void output_boundary_loop_points(std::string output_filename, const unsigned int u_resolution, int specific_patch )
		{
			std::cout << "\n> Starting writing points on boundary loop of given CAD geometry to file..." << std::endl;
		
			// Open file to write all points in
			std::ofstream file_to_write(output_filename);

			// If no specific patch defined, output boundary loop of all patches
			if(specific_patch<0)
			{
				//Loop over all patches
				for (unsigned int patch_itr = 0; patch_itr < m_patches.size(); patch_itr++)
				{
					BoundaryLoopVector boundary_loops = m_patches[patch_itr]->GetBoundaryLoops();

					// Loop over all boundary loops of current patch
					for(BoundaryLoopVector::iterator loop_i =  boundary_loops.begin(); loop_i!=boundary_loops.end(); loop_i++)
					{
						BoundaryEdgeVector boundary_edges = loop_i->GetBoundaryEdges();

						// Loop over all edges
						for (BoundaryEdgeVector::iterator edge_i = boundary_edges.begin(); edge_i!=boundary_edges.end(); edge_i++)
						{
							DoubleVector& knot_vec_u_i = edge_i->GetKnotVectorU();
							unsigned int knot_vector_u_dimension = knot_vec_u_i.size();

							double u_min = knot_vec_u_i[0];
							double u_max = knot_vec_u_i[knot_vector_u_dimension-1];

							double delta_u = (u_max-u_min) / u_resolution;

							for(unsigned int i=0; i<=u_resolution; i++)
							{
								double u_i = u_min + i*delta_u;

								// Evaluate point
								Point<3> edge_point;
								edge_i->EvaluateCurvePoint(edge_point, u_i);

								// Output point
								file_to_write << edge_point[0] << " " << edge_point[1] << " " << edge_point[2] << std::endl;
							}
						}
					}
				}
			}
			else // if specific_patch defined
			{
				BoundaryLoopVector boundary_loops = m_patches[specific_patch]->GetBoundaryLoops();

				// Loop over all boundary loops of current patch
				for(BoundaryLoopVector::iterator loop_i =  boundary_loops.begin(); loop_i!=boundary_loops.end(); loop_i++)
				{
					BoundaryEdgeVector boundary_edges = loop_i->GetBoundaryEdges();

					// Loop over all edges
					for (BoundaryEdgeVector::iterator edge_i = boundary_edges.begin(); edge_i!=boundary_edges.end(); edge_i++)
					{
						DoubleVector& knot_vec_u_i = edge_i->GetKnotVectorU();
						unsigned int knot_vector_u_dimension = knot_vec_u_i.size();

						double u_min = knot_vec_u_i[0];
						double u_max = knot_vec_u_i[knot_vector_u_dimension-1];

						double delta_u = (u_max-u_min) / u_resolution;

						for(unsigned int i=0; i<=u_resolution; i++)
						{
							double u_i = u_min + i*delta_u;

							// Evaluate point
							Point<3> edge_point;
							edge_i->EvaluateCurvePoint(edge_point, u_i);

							// Output point
							file_to_write << edge_point[0] << " " << edge_point[1] << " " << edge_point[2] << std::endl;
						}
					}
				}
			}

			// Close file
			file_to_write.close();

			std::cout << "\n> Finished writing points on boundary loop of given CAD geometry to file..." << std::endl;
		}	
		
		// --------------------------------------------------------------------------
		void output_control_point_displacements(std::string output_filename)
		{
			// Outputs the displacements of the control points in a format that may be read by Gid

			std::cout << "\n> Writing control points displacements..." << std::endl;
			std::ofstream output_file(output_filename);

			output_file << "Rhino Post Results File 1.0" << std::endl;
			output_file << "Result \"Displacement\" \"Load Case\" 0 Vector OnNodes" << std::endl;
			output_file << "Values" << std::endl;

			unsigned int cp_itr = 0;
			for (PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
			{
				for (ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
				{
					// It is important to iterate outside to stick to the carat settings
					++cp_itr;
					
					if(cp_i->IsRelevantForMapping())
						output_file << cp_itr << " " << cp_i->getdX() << " " << cp_i->getdY() << " " << cp_i->getdZ() << std::endl;
				}
			}

			output_file << "End Values" << std::endl;
			output_file.close();

			std::cout << "\t\t\t DONE" << std::endl;			
		}

		// --------------------------------------------------------------------------
		void write_updated_georhino_file(std::string new_filename, std::string original_georhino_filename)
		{
			std::cout << "\n> Writing control points displacements..." << std::endl;
			// write control point displacements
				std::string output_filename = new_filename + ".post.res";
				std::ofstream output_file(output_filename);
				output_file << "Rhino Post Results File 1.0" << std::endl;
				output_file << "Result \"Displacement\" \"Load Case\" 0 Vector OnNodes" << std::endl;
				output_file << "Values" << std::endl;
				unsigned int cp_itr = 0;
				for (PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
				{
					for (ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
					{
						// It is important to iterate outside to stick to the carat settings
						++cp_itr;
						
						if(cp_i->IsRelevantForMapping())
							output_file << cp_itr << " " << cp_i->getdX() << " " << cp_i->getdY() << " " << cp_i->getdZ() << std::endl;
					}
				}
				output_file << "End Values" << std::endl;
				output_file.close();
			std::cout << "\t\t\t DONE" << std::endl;			
		
			std::cout << "\n> Writing georhino file..." << std::endl;
			// write georhino file
				std::ifstream input_file(original_georhino_filename);
				output_filename = new_filename + ".georhino.txt";
				output_file.open(output_filename);
				std::string line; // last line to be read
				// copy all lines from the old file before the node section
					while(std::getline(input_file, line))
					{
						// if the line contains the word "NODE"
						if(find_substring(line, "NODE"))
							break;
						// else copy line
							output_file << line << std::endl;
					}
				// read and ignore node section
					while(std::getline(input_file, line))
					{
						// stop as soon as the line doesn't contain the word "NODE" anymore
						if(!find_substring(line, "NODE"))
							break;
					}
				// write node section
					int id = 1;
					for (PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
					{
						for (ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
						{
							// It is important to iterate outside to stick to the carat settings							
							if(cp_i->IsRelevantForMapping())
								// write cp
								// "  NODE  'id'  X  'x_coord'  Y  'y_coord'  Z  'z_coord'"
								output_file << "  NODE  " << id << "  X " << cp_i->getX0() << "  Y "<< cp_i->getY0() << "  Z " << cp_i->getZ0() << std::endl;
							id++;
						}
					}
					output_file << line << std::endl;
				// substitute all the lines with "NODE_ID", ignore those with "GP_POINT_GEO", copy the others
					std::vector<ControlPoint*> control_points_vector;
					for (PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
					{
						for (ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
						{
							control_points_vector.push_back(&(*cp_i));
						}
					}
					id = 1;
					while(std::getline(input_file, line))
					{
						if(find_substring(line, "NODE_ID"))
						{
							ControlPoint cp = *(control_points_vector[id-1]);
							if(cp.IsRelevantForMapping())
							{
								// "  NODE_ID  'id'  W  'weight'"
								output_file << "  NODE_ID  " << id << "  W  " << cp.getWeight() << std::endl;
							}
							else
							{
								// "  NODE_ID  0  X  'x_coord'  Y  'y_coord'  Z  'z_coord'  W  'weight'"
								output_file << "  NODE_ID  0  X  " << cp.getX0() << "  Y  " << cp.getY0() << "  Z  " << cp.getZ0() << "  W  " << cp.getWeight() << std::endl;
							}
							id++;
						}
						else if(find_substring(line, "GP_POINT_GEO"))
						{
							;// ignore line
						}
						else
						{
							// copy
								output_file << line << std::endl;
						}
					}
			std::cout << "\t\t\t DONE" << std::endl;			
		}
		// {
		// 	// write control point displacements
		// 		std::cout << "\n> Writing control points displacements..." << std::endl;

		// 		std::string output_filename = new_filename + ".post.res";
		// 		std::ofstream output_file(output_filename);
		// 		output_file << "Rhino Post Results File 1.0" << std::endl;
		// 		output_file << "Result \"Displacement\" \"Load Case\" 0 Vector OnNodes" << std::endl;
		// 		output_file << "Values" << std::endl;
		// 		unsigned int cp_itr = 0;
		// 		for (PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
		// 		{
		// 			for (ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
		// 			{
		// 				// It is important to iterate outside to stick to the carat settings
		// 				++cp_itr;
		// 				output_file << cp_itr << " " << cp_i->getdX() << " " << cp_i->getdY() << " " << cp_i->getdZ() << std::endl;
		// 			}
		// 		}
		// 		output_file << "End Values" << std::endl;
		// 		output_file.close();

		// 		std::cout << "\t\t\t DONE" << std::endl;			
		
		// 	std::string word; unsigned int node_id;
		// 	// parse original georhino file
		// 		std::list<std::string> lines_to_add_in_the_new_georhino_file;
		// 		int id = 1;
		// 		std::size_t start_pos, end_pos;

		// 		std::ifstream input_file(original_georhino_filename);

		// 		std::string line;
		// 		while(std::getline(input_file, line))
		// 		{
		// 			std::size_t found = line.find("NODE_ID");
		// 			if(found != std::string::npos)
		// 			{
		// 				// line has either the format:
		// 				// "  NODE_ID  0  X  'x_coord'  Y  'y_coord'  Z  'z_coord'  W  'weight'"
		// 				// or the format:
		// 				// "  NODE_ID  'id'  W  'weight'"
		// 				std::istringstream iss(line);
		// 				iss >> word >> node_id;
		// 				if(node_id == 0)
		// 				{
		// 					// create a line with the format:
		// 					// "  NODE  'id'  X  'x_coord'  Y  'y_coord'  Z  'z_coord'"
		// 						start_pos = line.find("X");
		// 						end_pos = line.find("W");
		// 						std::string line_to_save = "  NODE  " + std::to_string(id) + "  ";
		// 						std::string substring = line.substr(start_pos, end_pos-start_pos);
		// 						line_to_save += substring;
		// 					lines_to_add_in_the_new_georhino_file.push_back(line_to_save);
		// 				}
		// 				id++;
		// 			}
		// 		}

		// 	// write a new georhino fil
		// 		output_filename = new_filename + ".georhino.txt";
		// 		output_file.open(output_filename);

		// 		// go to beginning of input file
		// 			// input_file.seekg(0, std::ios::beg);
		// 			input_file.close();
		// 			input_file.open(original_georhino_filename);
		// 		// copy line by line the first section of the file
		// 			bool first_occurrence = false;
		// 			bool last_occurrence = false;
		// 			while(!last_occurrence && std::getline(input_file, line))
		// 			{
		// 				std::size_t found = line.find("NODE");
						
		// 				// before reaching the section
		// 				if(!first_occurrence)
		// 				{
		// 					// copy line by line
		// 					output_file << line << std::endl; // std::endl??

		// 					if(found != std::string::npos)
		// 					{
		// 						first_occurrence = true;
		// 					}
		// 				}
		// 				// while inside the section
		// 				else
		// 				{
		// 					if(found == std::string::npos)
		// 					{
		// 						last_occurrence = true;
		// 					}
		// 					else
		// 					{
		// 						// copy line by line
		// 						output_file << line << std::endl; // std::endl??
		// 					}
		// 				}
		// 			}
		// 		// add saved lines
		// 			while(!lines_to_add_in_the_new_georhino_file.empty())
		// 			{
						
		// 				output_file << lines_to_add_in_the_new_georhino_file.front() << std::endl; // std::endl??
		// 				lines_to_add_in_the_new_georhino_file.pop_front();
		// 			}
		// 			output_file << line << std::endl;
		// 		// copy line by line rest of the file, ignoring Gauss points at the end and
		// 		// changing lines like:
		// 		// "  NODE_ID  0  X  'x_coord'  Y  'y_coord'  Z  'z_coord'  W  'weight'"
		// 		// with lines of type:
		// 		// "  NODE_ID  'id'  W  'weight'"
		// 			id = 1;
		// 			while(std::getline(input_file, line))
		// 			{
		// 				std::size_t found_1 = line.find("NODE_ID");
		// 				std::size_t found_2 = line.find("GP_POINT_GEO");
		// 				if(found_1 != std::string::npos) // the line contains NODE_ID
		// 				{
		// 					std::istringstream iss(line);
		// 					std::string word; unsigned int node_id;
		// 					iss >> word >> node_id;
		// 					if(node_id == 0)
		// 					{
		// 						std::size_t pos = line.find("W");
		// 						output_file << " NODE_ID   " << id << "  " << line.substr(pos) << std::endl;
		// 					}
		// 					else
		// 					{
		// 						output_file << line << std::endl;
		// 					}
		// 					id++;
		// 				}
		// 				else if(found_2 == std::string::npos) // the line DOES NOT contain GP_POINT_GEO
		// 				{
		// 					output_file << line << std::endl;
		// 				}
		// 			}
		// 		input_file.close();
		// 		output_file.close();

		// }		

		// --------------------------------------------------------------------------
		bool find_substring(std::string str, std::string substring)
		{
			std::size_t found = str.find(substring);
			if(found != std::string::npos) return true;
			return false;
		}

		// --------------------------------------------------------------------------
		void output_surface_border_points(std::string output_filename, const unsigned int u_resolution, int specific_patch)
		{
			std::cout << "\n> Starting writing points on surface border of given CAD geometry to file..." << std::endl;
		
			// Set a max value for the coordinte output to avoid clutter through points that are plotted outside the trimmed surface
			double max_coordinate = 1000;

			// Open file to write all points in
			std::ofstream file_to_write(output_filename);

			// If no specific patch defined, output boundary loop of all patches
			if(specific_patch<0)
			{
				//Loop over all patches
				for (unsigned int patch_itr = 0; patch_itr < m_patches.size(); patch_itr++)
				{
					BoundaryLoopVector boundary_loops = m_patches[patch_itr]->GetBoundaryLoops();

					// Loop over all boundary loops of current patch
					for(BoundaryLoopVector::iterator loop_i =  boundary_loops.begin(); loop_i!=boundary_loops.end(); loop_i++)
					{
						BoundaryEdgeVector boundary_edges = loop_i->GetBoundaryEdges();

						// Loop over all edges
						for (BoundaryEdgeVector::iterator edge_i = boundary_edges.begin(); edge_i!=boundary_edges.end(); edge_i++)
						{
							DoubleVector& knot_vec_u_i = edge_i->GetKnotVectorU();
							unsigned int knot_vector_u_dimension = knot_vec_u_i.size();

							double u_min = knot_vec_u_i[0];
							double u_max = knot_vec_u_i[knot_vector_u_dimension-1];

							double delta_u = (u_max-u_min) / u_resolution;

							for(unsigned int i=0; i<=u_resolution; i++)
							{
								double u_i = u_min + i*delta_u;

								// Evaluate point in the parameter space
								Point<3> edge_point;
								edge_i->EvaluateCurvePoint(edge_point, u_i);

								// Check if u_i and v_j represent a point inside the closed boundary loop
								array_1d<double, 2> point_of_interest;
								point_of_interest[0] = edge_point[0];
								point_of_interest[1] = edge_point[1];
								bool point_is_inside = m_patches[patch_itr]->CheckIfPointIsInside(point_of_interest);

								if(point_is_inside)
								{
									// compute point in CAD-model for given u&v				
									Point<3> cad_point;
									m_patches[patch_itr]->GetSurface().EvaluateSurfacePoint(cad_point, edge_point[0], edge_point[1]);

									// Check and set a max value for irrelevant points outside the trimmed surface
									if(std::abs(cad_point.X())>max_coordinate)
										cad_point.X() = MathUtils<int>::Sign(cad_point.X()) * max_coordinate;
									if(std::abs(cad_point.Y())>max_coordinate)
										cad_point.Y() = MathUtils<int>::Sign(cad_point.Y()) * max_coordinate;
									if(std::abs(cad_point.Z())>max_coordinate)
										cad_point.Z() = MathUtils<int>::Sign(cad_point.Z()) * max_coordinate;														

									// Output point
									file_to_write << cad_point.X() << " " << cad_point.Y() << " " << cad_point.Z() << std::endl;
								}

							}
						}
					}
				}
			}
			else // if specific_patch defined
			{}
			// {
			// 	BoundaryLoopVector boundary_loops = m_patches[specific_patch]->GetBoundaryLoops();

			// 	// Loop over all boundary loops of current patch
			// 	for(BoundaryLoopVector::iterator loop_i =  boundary_loops.begin(); loop_i!=boundary_loops.end(); loop_i++)
			// 	{
			// 		BoundaryEdgeVector boundary_edges = loop_i->GetBoundaryEdges();

			// 		// Loop over all edges
			// 		for (BoundaryEdgeVector::iterator edge_i = boundary_edges.begin(); edge_i!=boundary_edges.end(); edge_i++)
			// 		{
			// 			DoubleVector& knot_vec_u_i = edge_i->GetKnotVectorU();
			// 			unsigned int knot_vector_u_dimension = knot_vec_u_i.size();

			// 			double u_min = knot_vec_u_i[0];
			// 			double u_max = knot_vec_u_i[knot_vector_u_dimension-1];

			// 			double delta_u = (u_max-u_min) / u_resolution;

			// 			for(unsigned int i=0; i<=u_resolution; i++)
			// 			{
			// 				double u_i = u_min + i*delta_u;

			// 				// Evaluate point
			// 				Point<3> edge_point;
			// 				edge_i->EvaluateCurvePoint(edge_point, u_i);

			// 				// Output point
			// 				file_to_write << edge_point[0] << " " << edge_point[1] << " " << edge_point[2] << std::endl;
			// 			}
			// 		}
			// 	}
			// }

			// Close file
			file_to_write.close();

			std::cout << "\n> Finished writing points on surface border of given CAD geometry to file..." << std::endl;
		}	

		// --------------------------------------------------------------------------
		void output_surface_border_points_two( std::string output_filename )
		{
			std::ofstream file_to_write(output_filename);
			int patch = 0;
			VectorPoint SlavePointVector;
			VectorPoint MasterPointVector; 
			DoubleVector CosineVector;
		
			for (BREPElementVector::iterator brep_elem_i = m_brep_elements.begin(); brep_elem_i != m_brep_elements.end(); ++brep_elem_i)
			{
				

				if(brep_elem_i->HasCouplingCondition() || brep_elem_i->HasDirichletCondition())
				{
					patch++;
					// Get Gauss points of current brep element
					BREPGaussPointVector brep_gps = brep_elem_i->GetGaussPoints();

					// Loop over all Gauss points of current brep element 
					for (BREPGaussPointVector::iterator brep_gp_i = brep_gps.begin(); brep_gp_i != brep_gps.end(); ++brep_gp_i)
					{
						// Read information from Gauss point
						unsigned int master_patch_id = brep_gp_i->GetPatchId();
						Patch::Pointer master_patch = m_patches[m_patch_position_in_patch_vector[master_patch_id]];
						Vector location_on_master_patch = brep_gp_i->GetLocation();
						matrix<double> R_gpi_master;
						double u_m = location_on_master_patch(0);
						double v_m = location_on_master_patch(1);
						Point<3> cad_point_master;
						master_patch->GetSurface().EvaluateSurfacePoint(cad_point_master, u_m, v_m);
						file_to_write << cad_point_master.X() << " " << cad_point_master.Y() << " " << cad_point_master.Z() << std::endl;

						if(brep_elem_i->HasCouplingCondition())
						{
							unsigned int slave_patch_id = brep_gp_i->GetSlavePatchId();
							Patch::Pointer slave_patch = m_patches[m_patch_position_in_patch_vector[slave_patch_id]];
							Vector location_on_slave_patch = brep_gp_i->GetSlaveLocation();
							double u_s = location_on_slave_patch(0);
							double v_s = location_on_slave_patch(1);
							Point<3> cad_point_slave;
							slave_patch->GetSurface().EvaluateSurfacePoint(cad_point_slave, u_s, v_s);
							file_to_write << cad_point_slave.X() << " " << cad_point_slave.Y() << " " << cad_point_slave.Z() << std::endl;

							// store information needed to evaluate C0-continuity
							MasterPointVector.push_back( cad_point_master );
							SlavePointVector.push_back( cad_point_slave );

							// evaluate C1-continuity
							matrix<double> g_master = master_patch->GetSurface().GetBaseVectors(-1,-1,u_m,v_m);
							Vector g1_m = ZeroVector(3);
							g1_m(0) = g_master(0,0);
							g1_m(1) = g_master(1,0);
							g1_m(2) = g_master(2,0);
							Vector g2_m = ZeroVector(3);
							g2_m(0) = g_master(0,1);
							g2_m(1) = g_master(1,1);
							g2_m(2) = g_master(2,1);

							matrix<double> g_slave = slave_patch->GetSurface().GetBaseVectors(-1,-1,u_s,v_s);
							Vector g1_s = ZeroVector(3);
							g1_s(0) = g_slave(0,0);
							g1_s(1) = g_slave(1,0);
							g1_s(2) = g_slave(2,0);
							Vector g2_s = ZeroVector(3);
							g2_s(0) = g_slave(0,1);
							g2_s(1) = g_slave(1,1);
							g2_s(2) = g_slave(2,1);

							auto normal_m = MathUtils<double>::CrossProduct(g1_m, g2_m);
							auto normal_s = MathUtils<double>::CrossProduct(g1_s, g2_s);

							auto inner_ms = inner_prod( normal_m, normal_s);

							auto cosine_theta = inner_ms/ ( norm_2(normal_m) * norm_2(normal_s) );

							CosineVector.push_back( std::abs(cosine_theta) ); // std::abs!!!!!

						}
						
					}
				}

			}
			double average, max;
			check_c0_continuity( MasterPointVector, SlavePointVector, average, max);

			KRATOS_WATCH("C_Zero Continuity");
			KRATOS_WATCH( max );
			KRATOS_WATCH( average );

			// -----------------------------------------------------------
			average = std::accumulate( CosineVector.begin(), CosineVector.end(), 0.0)/CosineVector.size();
			auto min = CosineVector[0];
			for(size_t i = 0; i < CosineVector.size(); i++)
			{
				if(CosineVector[i] < min)
				{
					min = CosineVector[i];
				}
			}
			KRATOS_WATCH("C_One Continuity");
			std::cout<< "Min: " << min << std::endl;
			KRATOS_WATCH( average );
			// -----------------------------------------------------------

		}
	////////////////////////////////////////////////////////

	// LEAST SQUARES APPROACH ###################################################
		// --------------------------------------------------------------------------
		void compute_nearest_points(const unsigned int u_resolution, const unsigned int v_resolution)
		{
			// Compute CAD points closer to FE-nodes
			// 1: Create for each patch coarse cloud of CAD points in x-y space for neighbor search later 
			// Each correspondingly created point is stored in a list.
			// Further lists in the same order are created to store the respective u & v parameter as well as the patch id
			// As list iterator we use a counter for the number of CAD nodes
			std::cout << "\n> Starting computation of nearest points..." << std::endl;
			unsigned int cad_node_counter = 0;
			NodeVector list_of_cad_nodes;
			DoubleVector list_of_us_of_cad_nodes;
			DoubleVector list_of_vs_of_cad_nodes;
			IntVector list_of_patch_itrs_of_cad_nodes;

			//Loop over all surface of all patches / faces
			for (unsigned int patch_itr = 0; patch_itr < m_patches.size(); patch_itr++)
			{
				// Get relevant data
				unsigned int patch_id =  m_patches[patch_itr]->GetId();
				DoubleVector& knot_vec_u_i = m_patches[patch_itr]->GetSurface().GetKnotVectorU();
				DoubleVector& knot_vec_v_i = m_patches[patch_itr]->GetSurface().GetKnotVectorV();
				std::cout << "\n> Processing Patch with brep_id " << patch_id << std::endl;

				// Pre-calculations
				unsigned int knot_vector_u_dimension = knot_vec_u_i.size();
				unsigned int knot_vector_v_dimension = knot_vec_v_i.size();

				double u_min = knot_vec_u_i[0];
				double u_max = knot_vec_u_i[knot_vector_u_dimension-1];
				double v_min = knot_vec_v_i[0];
				double v_max = knot_vec_v_i[knot_vector_v_dimension-1];

				double delta_u = (u_max-u_min) / u_resolution;
				double delta_v = (v_max-v_min) / v_resolution;

				// Loop over all u & v according to specified resolution
				for(unsigned int i=1; i<u_resolution; i++)
				{
					// current u-value
					double u_i = u_min + i*delta_u;

					for(unsigned int j=1; j<v_resolution; j++)
					{
						// current v-value
						double v_j = v_min + j*delta_v;

						// Check if u_i and v_j represent a point inside the closed boundary loop
						array_1d<double, 2> point_of_interest;
						point_of_interest[0] = u_i;
						point_of_interest[1] = v_j;
						bool point_is_inside = m_patches[patch_itr]->CheckIfPointIsInside(point_of_interest);

						if(point_is_inside)
						{
							// compute unique point in CAD-model for given u&v
							++cad_node_counter;					
							Point<3> cad_point_coordinates;
							m_patches[patch_itr]->GetSurface().EvaluateSurfacePoint(cad_point_coordinates, u_i, v_j);

							// Add id to point --> node. Add node to list of CAD nodes
							NodeType::Pointer new_cad_node = Node < 3 > ::Pointer(new Node<3>(cad_node_counter, cad_point_coordinates));
							list_of_cad_nodes.push_back(new_cad_node);

							// Store for cad node the corresponding cad information in separate vectors
							list_of_us_of_cad_nodes.push_back(u_i);
							list_of_vs_of_cad_nodes.push_back(v_j);
							list_of_patch_itrs_of_cad_nodes.push_back(patch_itr);
						}
					}
				}
			}
			// 2: Construct KD-Tree with all cad nodes
			std::cout << "\n> Starting construction of search-tree..." << std::endl;
			boost::timer timer;
			typedef Bucket< 3, NodeType, NodeVector, NodeType::Pointer, NodeIterator, DistanceIterator > BucketType;
			typedef Tree< KDTreePartition<BucketType> > tree;
			int bucket_size = 20;
			tree nodes_tree(list_of_cad_nodes.begin(), list_of_cad_nodes.end(), bucket_size);
			std::cout << "> Time needed for constructing search-tree: " << timer.elapsed() << " s" << std::endl;
			// 3: Evaluate nearest CAD nodes for all FEM nodes. Flag corresponding control points as relevant for mapping
			m_list_of_nearest_points.clear();
			m_list_of_u_of_nearest_points.clear();
			m_list_of_v_of_nearest_points.clear();
			m_list_of_span_u_of_nearest_points.clear();
			m_list_of_span_v_of_nearest_points.clear();		
			m_list_of_patch_of_nearest_points.clear();

			m_list_of_neighbour_points.clear();
			m_list_of_u_of_neighbour_points.clear();
			m_list_of_v_of_neighbour_points.clear();
			m_list_of_span_u_of_neighbour_points.clear();
			m_list_of_span_v_of_neighbour_points.clear();				
			m_list_of_patch_of_neighbour_points.clear();
			// Loop over all nodes of fe-model-part and find corresponding closest neighbors of cad-model
			std::cout << "\n> Starting to identify neighboring CAD points..." << std::endl;
			boost::timer timer_2;

			for(ModelPart::NodesContainerType::iterator node_itr = mr_fe_model_part.NodesBegin(); node_itr!=mr_fe_model_part.NodesEnd(); node_itr++)
			{	
				// Get node information
				ModelPart::NodeType& node_i = *node_itr;
				array_1d<double,3>  i_coord = node_i.Coordinates();

				// Search nearest cad neighbor of current integration point
				NodeType::Pointer nearest_point = nodes_tree.SearchNearestPoint( node_i ); 

				// Recover CAD information of nearest point
				double u_of_nearest_point = list_of_us_of_cad_nodes[nearest_point->Id()-1];
				double v_of_nearest_point = list_of_vs_of_cad_nodes[nearest_point->Id()-1];
				int patch_itr_of_nearest_point = list_of_patch_itrs_of_cad_nodes[nearest_point->Id()-1];
				IntVector knot_span_nearest_point = m_patches[patch_itr_of_nearest_point]->GetSurface().GetKnotSpan(u_of_nearest_point, v_of_nearest_point);
				int span_u_of_np = knot_span_nearest_point[0];
				int span_v_of_np = knot_span_nearest_point[1];

				// Store CAD information of nearest point
				m_list_of_neighbour_points.push_back(nearest_point);
				m_list_of_u_of_neighbour_points.push_back(u_of_nearest_point);
				m_list_of_v_of_neighbour_points.push_back(v_of_nearest_point);
				m_list_of_span_u_of_neighbour_points.push_back(span_u_of_np);
				m_list_of_span_v_of_neighbour_points.push_back(span_v_of_np);				
				m_list_of_patch_of_neighbour_points.push_back(patch_itr_of_nearest_point);

				// // Perform Newton-Raphson for detailed search
				// 	// Initialize P: point on the mesh
				// 	Vector P = ZeroVector(3);
				// 	P(0) = i_coord[0]; 
				// 	P(1) = i_coord[1];
				// 	P(2) = i_coord[2];
				// 	// Initialize Q_k: point on the CAD surface
				// 	Vector Q_k = ZeroVector(3);
				// 	Q_k(0) = nearest_point->X();
				// 	Q_k(1) = nearest_point->Y();
				// 	Q_k(2) = nearest_point->Z();
				// 	// Initialize what's needed in the Newton-Raphson iteration				
				// 	Vector Q_minus_P = ZeroVector(3); // Distance between current Q_k and P
				// 	Matrix myHessian = ZeroMatrix(2,2);
				// 	Vector myGradient = ZeroVector(2);
				// 	double det_H = 0;
				// 	Matrix InvH = ZeroMatrix(2,2);				
				// 	double u_k = u_of_nearest_point;
				// 	double v_k = v_of_nearest_point;
				// 	Point<3> newtonRaphPoint;

				// 	double norm_deltau = 100000000;
				// 	unsigned int k = 0;
				// 	unsigned int max_itr = 50;
				// 	while (norm_deltau > 1e-5)
				// 	{
				// 		// The distance between Q (on the CAD surface) and P (on the FE-mesh) is evaluated
				// 		Q_minus_P(0) = Q_k(0) - P(0);
				// 		Q_minus_P(1) = Q_k(1) - P(1);
				// 		Q_minus_P(2) = Q_k(2) - P(2);

				// 		// The distance is used to compute Hessian and gradient
				// 		m_patches[patch_itr_of_nearest_point]->GetSurface().EvaluateGradientsForClosestPointSearch(Q_minus_P, myHessian, myGradient , u_k, v_k);

				// 		// u_k and v_k are updated
				// 		MathUtils<double>::InvertMatrix( myHessian, InvH, det_H );
				// 		Vector deltau = prod(InvH,myGradient);
				// 		u_k -= deltau(0);
				// 		v_k -= deltau(1);

				// 		// Q is updated
				// 		m_patches[patch_itr_of_nearest_point]->GetSurface().EvaluateSurfacePoint(newtonRaphPoint, u_k, v_k);
				// 		Q_k(0) = newtonRaphPoint[0];
				// 		Q_k(1) = newtonRaphPoint[1];
				// 		Q_k(2) = newtonRaphPoint[2];
				// 		norm_deltau = norm_2(deltau);

				// 		k++;

				// 		if(k>max_itr)
				// 		{
				// 			std::cout << "WARNING!!! Newton-Raphson to find closest point did not converge in the following number of iterations: " << k-1 << std::endl;
				// 			KRATOS_WATCH(Q_k);
				// 			KRATOS_WATCH(P);
				// 		}
				// 	}

				// 	// Update nearest point
				// 	u_of_nearest_point = u_k;
				// 	v_of_nearest_point = v_k;
				// 	nearest_point->X() = Q_k(0);
				// 	nearest_point->Y() = Q_k(1);
				// 	nearest_point->Z() = Q_k(2);

				// Compute and store span of each parameter to avoid redundant computations later
				knot_span_nearest_point = m_patches[patch_itr_of_nearest_point]->GetSurface().GetKnotSpan(u_of_nearest_point, v_of_nearest_point);

				// Set flag to mark control point as relevant for mapping
				span_u_of_np = knot_span_nearest_point[0];
				span_v_of_np = knot_span_nearest_point[1];
				m_patches[patch_itr_of_nearest_point]->GetSurface().FlagControlPointsForMapping(span_u_of_np, span_v_of_np, u_of_nearest_point, v_of_nearest_point);

				// Store information about nearest point in vector for recovery in the same loop later when the mapping matrix is constructed
				m_list_of_nearest_points.push_back(nearest_point);
				m_list_of_u_of_nearest_points.push_back(u_of_nearest_point);
				m_list_of_v_of_nearest_points.push_back(v_of_nearest_point);
				m_list_of_span_u_of_nearest_points.push_back(span_u_of_np);
				m_list_of_span_v_of_nearest_points.push_back(span_v_of_np);				
				m_list_of_patch_of_nearest_points.push_back(patch_itr_of_nearest_point);
			}
			
			// Then we identify mapping relevant control points required from the specified boundary conditions
			// Accordingly we check all Gauss points of all brep elements for their control points
			for (BREPElementVector::iterator brep_elem_i = m_brep_elements.begin(); brep_elem_i != m_brep_elements.end(); ++brep_elem_i)
			{
				// Get Gauss points of current brep element
				BREPGaussPointVector brep_gps = brep_elem_i->GetGaussPoints();

				// Loop over all Gauss points of current brep element 
				for (BREPGaussPointVector::iterator brep_gp_i = brep_gps.begin(); brep_gp_i != brep_gps.end(); ++brep_gp_i)
				{
					// Flag control points on master patch
					unsigned int master_patch_id = brep_gp_i->GetPatchId();
					Vector location = brep_gp_i->GetLocation();
					m_patches[m_patch_position_in_patch_vector[master_patch_id]]->GetSurface().FlagControlPointsForMapping(-1, -1, location[0], location[1]);

					// Flag control points on slave patch if brep element is a coupling element
					if(brep_elem_i->HasCouplingCondition())
					{
						unsigned int slave_patch_id = brep_gp_i->GetSlavePatchId();
						location = brep_gp_i->GetSlaveLocation();
						m_patches[m_patch_position_in_patch_vector[slave_patch_id]]->GetSurface().FlagControlPointsForMapping(-1, -1, location[0], location[1]);
					}
				}
			}

			std::cout << "> Time needed for identify neighboring CAD points: " << timer_2.elapsed() << " s" << std::endl;
		}

		void print_nearest_points(std::string output_filename)
		{
			std::cout << "\n> Starting writing results of compute_nearest_points() to file..." << std::endl;
			std::cout << "> WARNING: x y z coordinates of nearest neighbour are wrong!!!" << std::endl;

			// Open file to write all points in
			std::ofstream file_to_write(output_filename);
			// write header
			file_to_write	<< "FE-node | nearest neighbour| Newton-Raphson point" << std::endl;
			file_to_write 	<< "x" << " "
							<< "y" << " "
							<< "z" << " "
							<< "|" << " "
							<< "patch-id" << " "
							<< "u" << " "
							<< "v" << " "
							<< "u_min" << " "
							<< "u_max" << " "
							<< "v_min" << " "
							<< "v_max" << " "
							<< "x" << " "
							<< "y" << " "
							<< "z" << " "
							<< "|" << " "
							<< "patch-id" << " "
							<< "u" << " "
							<< "v" << " "
							<< "x" << " "
							<< "y" << " "
							<< "z" << std::endl;
			// write points
			int i = 0, patch_itr;
			double u_i, v_i;
			for(ModelPart::NodesContainerType::iterator node_itr = mr_fe_model_part.NodesBegin(); node_itr!=mr_fe_model_part.NodesEnd(); node_itr++)
			{
				ModelPart::NodeType& node = *node_itr;
				
				// compute updated coordinates of Newton Raphson point
				u_i = m_list_of_u_of_nearest_points[i];
				v_i = m_list_of_v_of_nearest_points[i];
				patch_itr = m_list_of_patch_of_nearest_points[i];
				Point<3> cad_point;
				m_patches[patch_itr]->GetSurface().EvaluateSurfacePoint(cad_point, u_i, v_i);
				

				// write to file
				if(node.Has(SHAPE_CHANGE_ABSOLUTE))
				{
								// print FE-node
					file_to_write 	<< node.Coordinate(1) + node.GetValue(SHAPE_CHANGE_ABSOLUTE_X) << " "
									<< node.Coordinate(2) + node.GetValue(SHAPE_CHANGE_ABSOLUTE_Y) << " "
									<< node.Coordinate(3) + node.GetValue(SHAPE_CHANGE_ABSOLUTE_Z) << " "
								// print neighbour
									<< m_list_of_patch_of_neighbour_points[i] << " "

									<< m_list_of_u_of_neighbour_points[i] << " "
									<< m_list_of_v_of_neighbour_points[i] << " "

									<< m_patches[m_list_of_patch_of_neighbour_points[i]]->GetSurface().GetKnotVectorU().front() << " "
									<< m_patches[m_list_of_patch_of_neighbour_points[i]]->GetSurface().GetKnotVectorU().back() << " "
									<< m_patches[m_list_of_patch_of_neighbour_points[i]]->GetSurface().GetKnotVectorV().front() << " "
									<< m_patches[m_list_of_patch_of_neighbour_points[i]]->GetSurface().GetKnotVectorV().back() << " "

									<< m_list_of_neighbour_points[i]->X() << " "
									<< m_list_of_neighbour_points[i]->Y() << " "
									<< m_list_of_neighbour_points[i]->Z() << " "
								// print nearest point
									<< m_list_of_patch_of_nearest_points[i] << " "

									<< m_list_of_u_of_nearest_points[i] << " "
									<< m_list_of_v_of_nearest_points[i] << " "

									<< cad_point.X() << " "
									<< cad_point.Y() << " "
									<< cad_point.Z() << std::endl;
				}
				else
				{
								// print FE-node
					file_to_write 	<< node.Coordinate(1) << " "
									<< node.Coordinate(2) << " "
									<< node.Coordinate(3) << " "
								// print neighbour
									<< m_list_of_patch_of_neighbour_points[i] << " "

									<< m_list_of_u_of_neighbour_points[i] << " "
									<< m_list_of_v_of_neighbour_points[i] << " "

									<< m_patches[m_list_of_patch_of_neighbour_points[i]]->GetSurface().GetKnotVectorU().front() << " "
									<< m_patches[m_list_of_patch_of_neighbour_points[i]]->GetSurface().GetKnotVectorU().back() << " "
									<< m_patches[m_list_of_patch_of_neighbour_points[i]]->GetSurface().GetKnotVectorV().front() << " "
									<< m_patches[m_list_of_patch_of_neighbour_points[i]]->GetSurface().GetKnotVectorV().back() << " "

									<< m_list_of_neighbour_points[i]->X() << " "
									<< m_list_of_neighbour_points[i]->Y() << " "
									<< m_list_of_neighbour_points[i]->Z() << " "
								// print nearest point
									<< m_list_of_patch_of_nearest_points[i] << " "

									<< m_list_of_u_of_nearest_points[i] << " "
									<< m_list_of_v_of_nearest_points[i] << " "

									<< cad_point.X() << " "
									<< cad_point.Y() << " "
									<< cad_point.Z() << std::endl;
				}
				i++;
			}

			// Close file
			file_to_write.close();
			std::cout << "\n> Finished writing results" << std::endl;
			
		}

		void print_nearest_points_2(std::string output_filename) // modular code version of print_nearest_points()
		{
			std::cout << "\n> Writing data for debugging..." << std::endl;

			// Open file to write all points in
			std::ofstream file_to_write(output_filename);
			// write header
			file_to_write	<< "data point | nearest neighbour| Newton-Raphson point" << std::endl;
			file_to_write 	<< "x" << " "
							<< "y" << " "
							<< "z" << " "
							<< "|" << " "
							<< "patch-id" << " "
							<< "u" << " "
							<< "v" << " "
							<< "u_min" << " "
							<< "u_max" << " "
							<< "v_min" << " "
							<< "v_max" << " "
							<< "x" << " "
							<< "y" << " "
							<< "z" << " "
							<< "|" << " "
							<< "patch-id" << " "
							<< "u" << " "
							<< "v" << " "
							<< "x" << " "
							<< "y" << " "
							<< "z" << std::endl;
			// write points
			int i = 0;
			for(DataPointsList::iterator data_point_i = m_data_points.begin(); data_point_i != m_data_points.end(); data_point_i++)
			{
				// write to file
								// print FE-node
					file_to_write 	<< data_point_i->getX() << " "
									<< data_point_i->getY() << " "
									<< data_point_i->getZ() << " "
								// print neighbour
									<< data_point_i->getPatch()->GetId() - 1<< " "

									<< m_list_of_u_of_neighbour_points[i] << " "
									<< m_list_of_v_of_neighbour_points[i] << " "

									<< data_point_i->getPatch()->GetSurface().GetKnotVectorU().front() << " "
									<< data_point_i->getPatch()->GetSurface().GetKnotVectorU().back() << " "
									<< data_point_i->getPatch()->GetSurface().GetKnotVectorV().front() << " "
									<< data_point_i->getPatch()->GetSurface().GetKnotVectorV().back() << " "

									<< m_list_of_neighbour_points[i]->X() << " "
									<< m_list_of_neighbour_points[i]->Y() << " "
									<< m_list_of_neighbour_points[i]->Z() << " "
								// print nearest point
									<< data_point_i->getPatch()->GetId() - 1<< " "

									<< data_point_i->getU() << " "
									<< data_point_i->getV() << " "

									<< data_point_i->getCADPoint().X() << " "
									<< data_point_i->getCADPoint().Y() << " "
									<< data_point_i->getCADPoint().Z() << std::endl;
				i++;
			}

			// Close file
			file_to_write.close();
			
			std::cout << "\t\t\t DONE" << std::endl;
		}
		// LHS indirect, RHS indirect
			void compute_a_matrix()  // toBeChecked
			{
				std::cout << "\n> Starting computation of a matrix..." << std::endl;
				boost::timer a_matrix_timer;
				// Count relevant control points and assign each a unique mapping matrix Id (iterator over points)		
				// First we identify relevant control points affecting the cad points on the surface
				m_n_control_points = 0; 
				m_n_relevant_control_points = 0;
				unsigned int mapping_matrix_id = 0;
				for (PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
				{
					for (ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
					{
						if(cp_i->IsRelevantForMapping())
						{
							cp_i->SetMappingMatrixId(mapping_matrix_id);
							++m_n_relevant_control_points;
							++mapping_matrix_id;
						}
						++m_n_control_points;
					}
				}
				std::cout << "\n> Number of control points in total = " << m_n_control_points << "." << std::endl;
				std::cout << "\n> Number of control points relevant for mapping = " << m_n_relevant_control_points << ".\n" << std::endl;

				// Count FE nodes and assign each a unique mapping matrix id (iterator over nodes)
				m_n_relevant_fem_points = 0;
				for (ModelPart::NodesContainerType::iterator node_i = mr_fe_model_part.NodesBegin(); node_i != mr_fe_model_part.NodesEnd(); ++node_i)
				{
					node_i->SetValue(MAPPING_MATRIX_ID,m_n_relevant_fem_points); // ??? A_MAPPING_MATRIX
					m_n_relevant_fem_points++;
				}
				
				// Initialize a matrix
				double number_of_FE_dofs = 3 * m_n_relevant_fem_points;
				double number_of_CAD_dofs = 3 * m_n_relevant_control_points;
				m_a_matrix.resize(number_of_FE_dofs, number_of_CAD_dofs);
				m_a_matrix.clear();

				// Compute a matrix
				for(ModelPart::NodesContainerType::iterator node_i = mr_fe_model_part.NodesBegin(); node_i!=mr_fe_model_part.NodesEnd(); node_i++)
				{
						// Get the corresponding id of FE-node in the matrix			
						int row_id = node_i->GetValue(MAPPING_MATRIX_ID); // ??? A_MAPPING_MATRIX
						
						// Recover information about nearest CAD point				
						double u_of_nearest_point = m_list_of_u_of_nearest_points[row_id];
						double v_of_nearest_point = m_list_of_v_of_nearest_points[row_id];
						unsigned int span_u_of_nearest_point = m_list_of_span_u_of_nearest_points[row_id];
						unsigned int span_v_of_nearest_point = m_list_of_span_v_of_nearest_points[row_id];
						unsigned int patch_itr_of_nearest_point = m_list_of_patch_of_nearest_points[row_id];
						
						// Get CAD-shape-function-value for all control points affecting the nearest cad point
						matrix<double> R_CAD_Pi;
						m_patches[patch_itr_of_nearest_point]->GetSurface().EvaluateNURBSFunctions( span_u_of_nearest_point,
																								span_v_of_nearest_point,
																								u_of_nearest_point, 
																								v_of_nearest_point,
																								R_CAD_Pi );
						
						// Get the corresponding id of control points in the matrix
						matrix<unsigned int> mapping_matrix_ids_cad = m_patches[patch_itr_of_nearest_point]->GetSurface().GetMappingMatrixIds( span_u_of_nearest_point, span_v_of_nearest_point, u_of_nearest_point, v_of_nearest_point);

						// Assemble a matrix
						for(unsigned int i=0; i<mapping_matrix_ids_cad.size2();i++)
						{
							for(unsigned int j=0; j<mapping_matrix_ids_cad.size1();j++)
							{
								unsigned int col_id = mapping_matrix_ids_cad(j,i);
								double R_i = R_CAD_Pi(j,i);
								
								for(unsigned int m=0; m<3; m++)
								{
									m_a_matrix(3*row_id + m, 3*col_id + m) = R_i;
								}
							}
						}
				}
				std::cout << "\n> Computation of a matrix finished in " << a_matrix_timer.elapsed() << " s." << std::endl;			
			}
			
			void test_rectangular_matrix()
			{
				// 1: define relevant CPs vector
				Vector x_CP = ZeroVector(3*m_n_relevant_control_points);
				for (PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
				{
					for (ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
					{
						if(cp_i->IsRelevantForMapping())
						{
							int cp_id = cp_i->GetMappingMatrixId(); //??? line 1982:	cp_i->SetMappingMatrixId(mapping_matrix_id);
							x_CP[3*cp_id] = cp_i->getX();
							x_CP[3*cp_id+1] = cp_i->getY();
							x_CP[3*cp_id+2] = cp_i->getZ();
							// std::cout << std::endl; //???
							// KRATOS_WATCH(cp_i->getX()); //???
							// KRATOS_WATCH(cp_i->getY()); //???
							// KRATOS_WATCH(cp_i->getZ()); //???
						}
					}
				}
				// KRATOS_WATCH(x_CP); //???	

				// 2: compute product result
				Vector product_result = ZeroVector(3*m_n_relevant_fem_points);
				noalias(product_result) = prod(m_a_matrix,x_CP);
				
				// 3: define CAD points vector (i.e. expected result)
				Vector x_CAD = ZeroVector(3*m_n_relevant_fem_points);
				for(ModelPart::NodesContainerType::iterator node_i = mr_fe_model_part.NodesBegin(); node_i!=mr_fe_model_part.NodesEnd(); node_i++)
				{
					// std::cout << std::endl; //???
					// ModelPart::NodeType& node = *node_i; //???
					// KRATOS_WATCH(node.Coordinates()); //???
					int i_id = node_i->GetValue(MAPPING_MATRIX_ID); // ??? A_MAPPING_MATRIX	 //??? GetValue() could be avoided as the id are in ascending order		
					// KRATOS_WATCH(m_list_of_nearest_points[i_id]->X()); //???
					// KRATOS_WATCH(m_list_of_nearest_points[i_id]->Y()); //???
					// KRATOS_WATCH(m_list_of_nearest_points[i_id]->Z()); //???
					x_CAD[3*i_id] = m_list_of_nearest_points[i_id]->X(); //??? bla? store m_list_of_nearest_points as attribute or take it as input?
					x_CAD[3*i_id+1] = m_list_of_nearest_points[i_id]->Y();		// bla == i_id, since m_list_of_nearest_points is filled looping over nodes
					x_CAD[3*i_id+2] = m_list_of_nearest_points[i_id]->Z();
				}
				// KRATOS_WATCH(x_CAD); //???
				// KRATOS_WATCH(product_result); //???
				// KRATOS_WATCH(m_a_matrix); //???

				// 4: compare x_CAD with product_result
				Vector difference = x_CAD - product_result;
				double euclidean_norm = norm_2(difference);
				std::cout << "\n> test_rectangular_matrix(): the euclidean norm is " << euclidean_norm << std::endl;
			}

			void map_to_cad_space_2()//std::string output_filename)
			{
				std::cout << "\n> Starting to map to CAD space..." << std::endl;
				boost::timer function_timer;
				// 1: compute transpose
				boost::timer lhs_timer;
				std::cout << "\t> Computing transpose... \n";
				KRATOS_WATCH(m_a_matrix);
				Matrix transpose = trans(m_a_matrix); 
				KRATOS_WATCH(transpose);
				std::cout << "\t\t\t DONE" << std::endl;
				// 2: compute left hand side
				std::cout << "\t> Computing LHS...\n ";		
				Matrix lhs;
				std::cout << "\t\t[" << 3*m_n_relevant_control_points
						<< "x" << 3*m_n_relevant_fem_points
						<< "] [" << 3*m_n_relevant_fem_points
						<< "x" << 3*m_n_relevant_control_points
						<< "]\n";
				lhs.resize(3*m_n_relevant_control_points, 3*m_n_relevant_control_points);
				// noalias(lhs) = prod(transpose, m_a_matrix);
				prod(transpose, m_a_matrix, lhs);
				KRATOS_WATCH(lhs);
				std::cout << "\t\t\t DONE" << std::endl;
				std::cout << "\n> Computing  transpose and LHS finished in " << lhs_timer.elapsed() << " s." << std::endl;

				// 3: compute right hand side
				std::cout << "\t> Computing RHS... \n";		
				Vector x_FEM = ZeroVector(3*m_n_relevant_fem_points);
				for(ModelPart::NodesContainerType::iterator node_i = mr_fe_model_part.NodesBegin(); node_i!=mr_fe_model_part.NodesEnd(); node_i++)
				{
					// std::cout << std::endl; //???
					// ModelPart::NodeType& node = *node_i; //???
					// KRATOS_WATCH(node.Coordinates()); //???
					int i_id = node_i->GetValue(MAPPING_MATRIX_ID); // ??? A_MAPPING_MATRIX	 //??? GetValue() could be avoided as the id are in ascending order		
					// KRATOS_WATCH(m_list_of_nearest_points[i_id]->X()); //???
					// KRATOS_WATCH(m_list_of_nearest_points[i_id]->Y()); //???
					// KRATOS_WATCH(m_list_of_nearest_points[i_id]->Z()); //???
					
					ModelPart::NodeType& node = *node_i;
					array_1d<double,3>  node_coord = node.Coordinates();

					x_FEM[3*i_id] = node_coord[0] + node_i->GetValue(SHAPE_CHANGE_ABSOLUTE_X);
					x_FEM[3*i_id+1] = node_coord[1] + node_i->GetValue(SHAPE_CHANGE_ABSOLUTE_Y);
					x_FEM[3*i_id+2] = node_coord[2] + node_i->GetValue(SHAPE_CHANGE_ABSOLUTE_Z);
				}
				Vector rhs = prod(transpose,x_FEM);		
				KRATOS_WATCH(rhs);		
				std::cout << "\t\t\t DONE" << std::endl;
				// 4: solve for x_CP: updated position of CP
				std::cout << "\t> Solving... \n";	
				Vector x_CP = ZeroVector(3*m_n_relevant_control_points);
				CompressedMatrixType compressed_lhs = lhs;
				m_linear_solver->Solve(compressed_lhs, x_CP, rhs);		
				std::cout << "\t\t\t DONE" << std::endl;

				// 5: update
				std::cout << "\t> Updating control points positions... ";
				Vector x_CP_old = ZeroVector(3*m_n_relevant_control_points);
				for (PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
				{
					for (ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
					{
						if(cp_i->IsRelevantForMapping())
						{
							int cp_id = cp_i->GetMappingMatrixId(); //??? line 1982:	cp_i->SetMappingMatrixId(mapping_matrix_id);
							x_CP_old[3*cp_id] = cp_i->getX();
							x_CP_old[3*cp_id+1] = cp_i->getY();
							x_CP_old[3*cp_id+2] = cp_i->getZ();
							// std::cout << std::endl; //???
							// KRATOS_WATCH(cp_i->getX()); //???
							// KRATOS_WATCH(cp_i->getY()); //???
							// KRATOS_WATCH(cp_i->getZ()); //???
						}
					}
				}
				Vector ds = x_CP - x_CP_old; // CP position update
				m_cad_reader.UpdateControlPoints(m_patches, ds);

				std::cout << "\n> Mapping to CAD space finished in " << function_timer.elapsed() << " s." << std::endl;	

				// 6: testing purpose
				Vector x_CAD = ZeroVector(3*m_n_relevant_fem_points);
				for(ModelPart::NodesContainerType::iterator node_i = mr_fe_model_part.NodesBegin(); node_i!=mr_fe_model_part.NodesEnd(); node_i++)
				{
					// std::cout << std::endl; //???
					// ModelPart::NodeType& node = *node_i; //???
					// KRATOS_WATCH(node.Coordinates()); //???
					int i_id = node_i->GetValue(MAPPING_MATRIX_ID); // ??? A_MAPPING_MATRIX	 //??? GetValue() could be avoided as the id are in ascending order		
					// KRATOS_WATCH(m_list_of_nearest_points[i_id]->X()); //???
					// KRATOS_WATCH(m_list_of_nearest_points[i_id]->Y()); //???
					// KRATOS_WATCH(m_list_of_nearest_points[i_id]->Z()); //???
					x_CAD[3*i_id] = m_list_of_nearest_points[i_id]->X(); //??? bla? store list_of_nearest_points as attribute or take it as input?
					x_CAD[3*i_id+1] = m_list_of_nearest_points[i_id]->Y();		// bla == i_id, since list_of_nearest_points is filled looping over nodes
					x_CAD[3*i_id+2] = m_list_of_nearest_points[i_id]->Z();
				}

				// std::ofstream file_to_write(output_filename);
				// file_to_write << "lhs:" << std::endl;
				// for(unsigned int i = 0; i<3*m_n_relevant_control_points; i++)
				// {
				// 	for(unsigned int j = 0; j<3*m_n_relevant_control_points; j++)
				// 	{
				// 		file_to_write << lhs(i,j) << "\t";
				// 	}
				// 	file_to_write << std::endl;
				// }
				// file_to_write << "rhs:" << std::endl;
				// for(unsigned int i = 0; i<3*m_n_relevant_control_points; i++)
				// {
				// 	file_to_write << rhs[i] << std::endl;
				// }
				// file_to_write.close();
			}
		////////////////////////////////////////////////////////
		
		// LHS direct, RHS indirect
			void compute_lhs_matrix() // toBeChecked
			{
				std::cout << "\n> Starting direct computation of LHS..." << std::endl;
				boost::timer lhs_timer;
				// Count relevant control points and assign each a unique mapping matrix Id (iterator over points)		
				// First we identify relevant control points affecting the cad points on the surface
				m_n_control_points = 0; 
				m_n_relevant_control_points = 0;
				unsigned int mapping_matrix_id = 0;
				for (PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
				{
					for (ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
					{
						if(cp_i->IsRelevantForMapping())
						{
							cp_i->SetMappingMatrixId(mapping_matrix_id);
							++m_n_relevant_control_points;
							++mapping_matrix_id;
						}
						++m_n_control_points;
					}
				}
				std::cout << "\n> Number of control points in total = " << m_n_control_points << "." << std::endl;
				std::cout << "\n> Number of control points relevant for mapping = " << m_n_relevant_control_points << "." << std::endl;

				// Count FE nodes and assign each a unique mapping matrix id (iterator over nodes)
				m_n_relevant_fem_points = 0;
				for (ModelPart::NodesContainerType::iterator node_i = mr_fe_model_part.NodesBegin(); node_i != mr_fe_model_part.NodesEnd(); ++node_i)
				{
					node_i->SetValue(MAPPING_MATRIX_ID,m_n_relevant_fem_points); // ??? A_MAPPING_MATRIX
					m_n_relevant_fem_points++;
				}
				
				// Initialize a matrix
				double number_of_CAD_dofs = 3 * m_n_relevant_control_points;
				m_lhs_matrix.resize(number_of_CAD_dofs, number_of_CAD_dofs);
				m_lhs_matrix.clear();

				// Compute a matrix
				for(ModelPart::NodesContainerType::iterator node_i = mr_fe_model_part.NodesBegin(); node_i!=mr_fe_model_part.NodesEnd(); node_i++)
				{
						// Get the corresponding id of FE-node in the matrix			
						int row_id = node_i->GetValue(MAPPING_MATRIX_ID); // ??? A_MAPPING_MATRIX
						
						// Recover information about nearest CAD point				
						double u_of_nearest_point = m_list_of_u_of_nearest_points[row_id];
						double v_of_nearest_point = m_list_of_v_of_nearest_points[row_id];
						unsigned int span_u_of_nearest_point = m_list_of_span_u_of_nearest_points[row_id];
						unsigned int span_v_of_nearest_point = m_list_of_span_v_of_nearest_points[row_id];
						unsigned int patch_itr_of_nearest_point = m_list_of_patch_of_nearest_points[row_id];
						
						// Get CAD-shape-function-value for all control points affecting the nearest cad point
						matrix<double> R_CAD_Pi;
						m_patches[patch_itr_of_nearest_point]->GetSurface().EvaluateNURBSFunctions( span_u_of_nearest_point,
																								span_v_of_nearest_point,
																								u_of_nearest_point, 
																								v_of_nearest_point,
																								R_CAD_Pi );
						
						// Get the corresponding id of control points in the matrix
						matrix<unsigned int> mapping_matrix_ids_cad = m_patches[patch_itr_of_nearest_point]->GetSurface().GetMappingMatrixIds( span_u_of_nearest_point, span_v_of_nearest_point, u_of_nearest_point, v_of_nearest_point);

						// Assemble a matrix
						for(unsigned int i=0; i<mapping_matrix_ids_cad.size1(); i++)
						{
							for(unsigned int j=0; j<mapping_matrix_ids_cad.size2(); j++)
							{
								for(unsigned int k=0; k<mapping_matrix_ids_cad.size1(); k++)
								{
									for(unsigned int l=0; l<mapping_matrix_ids_cad.size2(); l++)
									{
										unsigned int row_id = mapping_matrix_ids_cad(i,j);
										unsigned int col_id = mapping_matrix_ids_cad(k,l);
										double R_ij = R_CAD_Pi(i,j);
										double R_kl = R_CAD_Pi(k,l);
										for(unsigned int m=0; m<3; m++)
										{
											// add contribution of 
											m_lhs_matrix(3*row_id+m, 3*col_id+m) += R_ij*R_kl;
										}
									}
								}
							}
						}
				}
				std::cout << "\n> Direct computation of LHS finished in " << lhs_timer.elapsed() << " s." << std::endl;
			}
			
			void compute_rhs_vector(Vector& rhs) // toBeChecked
			{
				// 3: compute right hand side
				std::cout << "\t> Computing RHS... \n";		
				Vector x_FEM = ZeroVector(3*m_n_relevant_fem_points);
				for(ModelPart::NodesContainerType::iterator node_i = mr_fe_model_part.NodesBegin(); node_i!=mr_fe_model_part.NodesEnd(); node_i++)
				{
					int i_id = node_i->GetValue(MAPPING_MATRIX_ID); // ??? A_MAPPING_MATRIX	 //??? GetValue() could be avoided as the id are in ascending order		
					
					ModelPart::NodeType& node = *node_i;
					array_1d<double,3>  node_coord = node.Coordinates();

					x_FEM[3*i_id] = node_coord[0] + node_i->GetValue(SHAPE_CHANGE_ABSOLUTE_X);
					x_FEM[3*i_id+1] = node_coord[1] + node_i->GetValue(SHAPE_CHANGE_ABSOLUTE_Y);
					x_FEM[3*i_id+2] = node_coord[2] + node_i->GetValue(SHAPE_CHANGE_ABSOLUTE_Z);
				}
				Matrix transpose = trans(m_a_matrix);
				rhs.clear(); 
				rhs = prod(transpose,x_FEM);		
				std::cout << "\t\t\t DONE" << std::endl;
			}

			void map_to_cad_space_3()//) // toBeChecked
			{
				std::cout << "\n> Starting to map to CAD space..." << std::endl;
				boost::timer function_timer;

				// 1: compute lhs and rhs
				// m_lhs_matrix = compute_lhs_matrix();
				KRATOS_WATCH(m_a_matrix);
				compute_lhs_matrix();
				KRATOS_WATCH(m_lhs_matrix);
				Vector rhs;
				compute_rhs_vector(rhs);
				KRATOS_WATCH(rhs);
			
				// 2: solve for x_CP: updated position of CP
				std::cout << "\t> Solving... \n";	
				Vector x_CP = ZeroVector(3*m_n_relevant_control_points);
				CompressedMatrixType compressed_lhs = m_lhs_matrix; //???
				m_linear_solver->Solve(compressed_lhs, x_CP, rhs);		
				std::cout << "\t\t\t DONE" << std::endl;

				// 3: update
				std::cout << "\t> Updating control points positions... ";
				Vector x_CP_old = ZeroVector(3*m_n_relevant_control_points);
				for (PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
				{
					for (ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
					{
						if(cp_i->IsRelevantForMapping())
						{
							int cp_id = cp_i->GetMappingMatrixId(); //??? line 1982:	cp_i->SetMappingMatrixId(mapping_matrix_id);
							x_CP_old[3*cp_id] = cp_i->getX();
							x_CP_old[3*cp_id+1] = cp_i->getY();
							x_CP_old[3*cp_id+2] = cp_i->getZ();
						}
					}
				}
				Vector ds = x_CP - x_CP_old; // CP position update
				m_cad_reader.UpdateControlPoints(m_patches, ds);

				std::cout << "\n> Mapping to CAD space finished in " << function_timer.elapsed() << " s." << std::endl;	

				// std::ofstream file_to_write(output_filename);
				// file_to_write << "lhs:" << std::endl;
				// for(unsigned int i = 0; i<3*m_n_relevant_control_points; i++)
				// {
				// 	for(unsigned int j = 0; j<3*m_n_relevant_control_points; j++)
				// 	{
				// 		file_to_write << m_lhs_matrix(i,j) << "\t";
				// 	}
				// 	file_to_write << std::endl;
				// }
				// file_to_write << "rhs:" << std::endl;
				// for(unsigned int i = 0; i<3*m_n_relevant_control_points; i++)
				// {
				// 	file_to_write << rhs[i] << std::endl;
				// }
				// file_to_write.close();

				// 4: Test solution
				Vector rhs_test = ZeroVector(3*m_n_relevant_control_points);
				noalias(rhs_test) = prod(compressed_lhs,x_CP);
				Vector rhs_difference = rhs - rhs_test;
				double normalized_difference_in_rhs = norm_2(rhs_difference);
				std::cout << "\n> Solution of linear system leads to a difference in the RHS of: normalized_difference_in_rhs = " << normalized_difference_in_rhs << std::endl;

			}

		////////////////////////////////////////////////////////

		/////////////// LHS direct, RHS direct + ///////////////
		/////////////// SOLVING ONE COORDINATE AT A TIME: m_lhs_x, m_lhs_y, m_lhs_z, m_rhs_x, m_rhs_y, m_rhs_z ///////////////
			void compute_lhs_small()
				{
					std::cout << "\n> Starting direct computation of LHS..." << std::endl;
					boost::timer lhs_timer;
					// Count relevant control points and assign each a unique mapping matrix Id (iterator over points)		
						// First we identify relevant control points affecting the cad points on the surface
						m_n_control_points = 0; 
						m_n_relevant_control_points = 0;
						unsigned int mapping_matrix_id = 0;
						for (PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
						{
							for (ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
							{
								if(cp_i->IsRelevantForMapping())
								{
									cp_i->SetMappingMatrixId(mapping_matrix_id);
									++m_n_relevant_control_points;
									++mapping_matrix_id;
								}
								++m_n_control_points;
							}
						}
						std::cout << "\n> Number of control points in total = " << m_n_control_points << "." << std::endl;
						std::cout << "\n> Number of control points relevant for mapping = " << m_n_relevant_control_points << "." << std::endl;

						// Count FE nodes and assign each a unique mapping matrix id (iterator over nodes)
						m_n_relevant_fem_points = 0;
						for (ModelPart::NodesContainerType::iterator node_i = mr_fe_model_part.NodesBegin(); node_i != mr_fe_model_part.NodesEnd(); ++node_i)
						{
							node_i->SetValue(MAPPING_MATRIX_ID,m_n_relevant_fem_points); // ??? A_MAPPING_MATRIX
							m_n_relevant_fem_points++;
						}
						
					// Initialize a matrix
						m_lhs_x.resize(m_n_relevant_control_points, m_n_relevant_control_points);
						m_lhs_x.clear();
						m_lhs_y.resize(m_n_relevant_control_points, m_n_relevant_control_points);
						m_lhs_y.clear();
						m_lhs_z.resize(m_n_relevant_control_points, m_n_relevant_control_points);
						m_lhs_z.clear();

					// Compute a matrix
						for(ModelPart::NodesContainerType::iterator node_i = mr_fe_model_part.NodesBegin(); node_i!=mr_fe_model_part.NodesEnd(); node_i++)
						{
							// Get the corresponding id of FE-node in the matrix			
							int row_id = node_i->GetValue(MAPPING_MATRIX_ID); // ??? A_MAPPING_MATRIX
							
							// Recover information about nearest CAD point				
							double u_of_nearest_point = m_list_of_u_of_nearest_points[row_id];
							double v_of_nearest_point = m_list_of_v_of_nearest_points[row_id];
							unsigned int span_u_of_nearest_point = m_list_of_span_u_of_nearest_points[row_id];
							unsigned int span_v_of_nearest_point = m_list_of_span_v_of_nearest_points[row_id];
							unsigned int patch_itr_of_nearest_point = m_list_of_patch_of_nearest_points[row_id];
							
							// Get CAD-shape-function-value for all control points affecting the nearest cad point
							matrix<double> R_CAD_Pi;
							m_patches[patch_itr_of_nearest_point]->GetSurface().EvaluateNURBSFunctions( span_u_of_nearest_point,
																									span_v_of_nearest_point,
																									u_of_nearest_point, 
																									v_of_nearest_point,
																									R_CAD_Pi );
							
							// Get the corresponding id of control points in the matrix
							matrix<unsigned int> mapping_matrix_ids_cad = m_patches[patch_itr_of_nearest_point]->GetSurface().GetMappingMatrixIds( span_u_of_nearest_point, span_v_of_nearest_point, u_of_nearest_point, v_of_nearest_point);

							// Assemble a matrix
							for(unsigned int i=0; i<mapping_matrix_ids_cad.size1(); i++)
							{
								for(unsigned int j=0; j<mapping_matrix_ids_cad.size2(); j++)
								{
									for(unsigned int k=0; k<mapping_matrix_ids_cad.size1(); k++)
									{
										for(unsigned int l=0; l<mapping_matrix_ids_cad.size2(); l++)
										{
											unsigned int row_id = mapping_matrix_ids_cad(i,j);
											unsigned int col_id = mapping_matrix_ids_cad(k,l);
											double R_ij = R_CAD_Pi(i,j);
											double R_kl = R_CAD_Pi(k,l);
											// add contribution of 
											m_lhs_x(row_id, col_id) += R_ij*R_kl;
											m_lhs_y(row_id, col_id) += R_ij*R_kl;
											m_lhs_z(row_id, col_id) += R_ij*R_kl;
										}
									}
								}
							}
						}
					std::cout << "\n> Direct computation of LHS finished in " << lhs_timer.elapsed() << " s." << std::endl;

				}

			void compute_rhs_small()
				{
					std::cout << "\t> Computing RHS... \n";
					// Initialize rhs vectors
						m_rhs_x.resize(m_n_relevant_control_points);	
						m_rhs_y.resize(m_n_relevant_control_points);	
						m_rhs_z.resize(m_n_relevant_control_points);	
						m_rhs_x.clear();
						m_rhs_y.clear();
						m_rhs_z.clear();

					// Compute rhs vectors
						for(ModelPart::NodesContainerType::iterator node_i = mr_fe_model_part.NodesBegin(); node_i!=mr_fe_model_part.NodesEnd(); node_i++)
						{
							// Get the corresponding id of FE-node in the matrix			
							int row_id = node_i->GetValue(MAPPING_MATRIX_ID); // ??? A_MAPPING_MATRIX

							ModelPart::NodeType& node = *node_i;
							array_1d<double,3>  node_coord = node.Coordinates();
							
							// Recover information about nearest CAD point				
							double u_of_nearest_point = m_list_of_u_of_nearest_points[row_id];
							double v_of_nearest_point = m_list_of_v_of_nearest_points[row_id];
							unsigned int span_u_of_nearest_point = m_list_of_span_u_of_nearest_points[row_id];
							unsigned int span_v_of_nearest_point = m_list_of_span_v_of_nearest_points[row_id];
							unsigned int patch_itr_of_nearest_point = m_list_of_patch_of_nearest_points[row_id];
							
							// Get CAD-shape-function-value for all control points affecting the nearest cad point
							matrix<double> R_CAD_Pi;
							m_patches[patch_itr_of_nearest_point]->GetSurface().EvaluateNURBSFunctions( span_u_of_nearest_point,
																									span_v_of_nearest_point,
																									u_of_nearest_point, 
																									v_of_nearest_point,
																									R_CAD_Pi );
							
							// Get the corresponding id of control points in the vector
							matrix<unsigned int> mapping_matrix_ids_cad = m_patches[patch_itr_of_nearest_point]->GetSurface().GetMappingMatrixIds( span_u_of_nearest_point, span_v_of_nearest_point, u_of_nearest_point, v_of_nearest_point);

							// Assemble a matrix
							for(unsigned int i=0; i<mapping_matrix_ids_cad.size1(); i++)
							{
								for(unsigned int j=0; j<mapping_matrix_ids_cad.size2(); j++)
								{
									unsigned int row_id = mapping_matrix_ids_cad(i,j);								
									double R_ij = R_CAD_Pi(i,j);
									m_rhs_x(row_id) += R_ij * (node_coord[0] + node_i->GetValue(SHAPE_CHANGE_ABSOLUTE_X));
									m_rhs_y(row_id) += R_ij * (node_coord[1] + node_i->GetValue(SHAPE_CHANGE_ABSOLUTE_Y));
									m_rhs_z(row_id) += R_ij * (node_coord[2] + node_i->GetValue(SHAPE_CHANGE_ABSOLUTE_Z));
								}
							}					
					}	
					std::cout << "\t\t\t DONE" << std::endl;
				}

			void map_to_cad_space_4(double alpha, double beta, double delta,
									double penalty_factor_disp,
									double penalty_factor_rot,
									boost::python::list& edges_with_enforced_tangent_continuity)
				{
					std::cout << "\n> Starting to map to CAD space..." << std::endl;
					boost::timer function_timer;

					// 1: compute lhs and rhs
						compute_lhs_small(); // computes m_lhs_x, m_lhs_y, m_lhs_z
						compute_rhs_small(); // computes m_rhs_x, m_rhs_y, m_rhs_z
						print_mapping_ids(true); ///////////////////////////////////////////////////////////////////////////////////////////
					// 2: apply regularization
						if(alpha > 0) alpha_criterion(alpha);
						if(beta  > 0)  beta_criterion(beta );
						if(delta > 0) delta_criterion(delta);
					// 2bis: apply B.C.s
						apply_boundary_conditions_small(penalty_factor_disp, penalty_factor_rot, edges_with_enforced_tangent_continuity);
					// 3: solve for x_CP: updated position of CP
						CompressedMatrixType compressed_lhs_x = m_lhs_x; //???
						CompressedMatrixType compressed_lhs_y = m_lhs_y; //???
						CompressedMatrixType compressed_lhs_z = m_lhs_z; //???

						std::cout << "\t> Solving x... \n";	
						Vector x_CP = ZeroVector(m_n_relevant_control_points);
						m_linear_solver->Solve(compressed_lhs_x, x_CP, m_rhs_x);		
						std::cout << "\t\t\t DONE" << std::endl;

						std::cout << "\t> Solving y... \n";	
						Vector y_CP = ZeroVector(m_n_relevant_control_points);
						m_linear_solver->Solve(compressed_lhs_y, y_CP, m_rhs_y);		
						std::cout << "\t\t\t DONE" << std::endl;

						std::cout << "\t> Solving z... \n";	
						Vector z_CP = ZeroVector(m_n_relevant_control_points);
						m_linear_solver->Solve(compressed_lhs_z, z_CP, m_rhs_z);		
						std::cout << "\t\t\t DONE" << std::endl;

					// 4: update
						std::cout << "\t> Updating control points positions... ";
						Vector ds = ZeroVector(3*m_n_relevant_control_points);
						for (PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
						{
							for (ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
							{
								if(cp_i->IsRelevantForMapping())
								{
									int cp_id = cp_i->GetMappingMatrixId(); //??? line 1982:	cp_i->SetMappingMatrixId(mapping_matrix_id);
									ds[3*cp_id]   = x_CP[cp_id] - cp_i->getX();
									ds[3*cp_id+1] = y_CP[cp_id] - cp_i->getY();
									ds[3*cp_id+2] = z_CP[cp_id] - cp_i->getZ();
								}
							}
						}
						m_cad_reader.UpdateControlPoints(m_patches, ds);

					std::cout << "\n> Mapping to CAD space finished in " << function_timer.elapsed() << " s." << std::endl;	
					// 5: Test solution
						Vector rhs_test = ZeroVector(m_n_relevant_control_points);
						noalias(rhs_test) = prod(compressed_lhs_x,x_CP);
						Vector rhs_difference = m_rhs_x - rhs_test;
						double normalized_difference_in_rhs = norm_2(rhs_difference);
						std::cout << "\n> Solution of linear system leads to a difference in the RHS_x of: normalized_difference_in_rhs = " << normalized_difference_in_rhs << std::endl;
			
						rhs_test = ZeroVector(m_n_relevant_control_points);
						noalias(rhs_test) = prod(compressed_lhs_y,y_CP);
						rhs_difference = m_rhs_y - rhs_test;
						normalized_difference_in_rhs = norm_2(rhs_difference);
						std::cout << "\n> Solution of linear system leads to a difference in the RHS_y of: normalized_difference_in_rhs = " << normalized_difference_in_rhs << std::endl;
			
						rhs_test = ZeroVector(m_n_relevant_control_points);
						noalias(rhs_test) = prod(compressed_lhs_z,z_CP);
						rhs_difference = m_rhs_z - rhs_test;
						normalized_difference_in_rhs = norm_2(rhs_difference);
						std::cout << "\n> Solution of linear system leads to a difference in the RHS_z of: normalized_difference_in_rhs = " << normalized_difference_in_rhs << std::endl;	
				}

			// regularization techniques
				void alpha_criterion(double alpha)
						{

					unsigned int u_resolution = 200;
					unsigned int v_resolution = 200;


					// Each correspondingly created point is stored in a list.
					// Further lists in the same order are created to store the respective u & v parameter as well as the patch id
					// As list iterator we use a counter for the number of CAD nodes
					unsigned int cad_node_counter = 0;
					NodeVector list_of_cad_nodes;
					DoubleVector list_of_us_of_cad_nodes;
					DoubleVector list_of_vs_of_cad_nodes;
					IntVector list_of_patch_itrs_of_cad_nodes;

					//Loop over all surface of all patches / faces
					for (unsigned int patch_itr = 0; patch_itr < m_patches.size(); patch_itr++)
					{
						// Get relevant data
						unsigned int patch_id =  m_patches[patch_itr]->GetId();
						DoubleVector& knot_vec_u_i = m_patches[patch_itr]->GetSurface().GetKnotVectorU();
						DoubleVector& knot_vec_v_i = m_patches[patch_itr]->GetSurface().GetKnotVectorV();
						std::cout << "\n> Processing Patch with brep_id " << patch_id << std::endl;

						// Pre-calculations
						unsigned int knot_vector_u_dimension = knot_vec_u_i.size();
						unsigned int knot_vector_v_dimension = knot_vec_v_i.size();

						double u_min = knot_vec_u_i[0];
						double u_max = knot_vec_u_i[knot_vector_u_dimension-1];
						double v_min = knot_vec_v_i[0];
						double v_max = knot_vec_v_i[knot_vector_v_dimension-1];

						double delta_u = (u_max-u_min) / u_resolution;
						double delta_v = (v_max-v_min) / v_resolution;

						// Loop over all u & v according to specified resolution
						for(unsigned int i=1; i<u_resolution; i++)
						{
							// current u-value
							double u_i = u_min + i*delta_u;

							for(unsigned int j=1; j<v_resolution; j++)
							{
								// current v-value
								double v_j = v_min + j*delta_v;

								// Check if u_i and v_j represent a point inside the closed boundary loop
								array_1d<double, 2> point_of_interest;
								point_of_interest[0] = u_i;
								point_of_interest[1] = v_j;
								bool point_is_inside = m_patches[patch_itr]->CheckIfPointIsInside(point_of_interest);

								if(point_is_inside)
								{
									// compute unique point in CAD-model for given u&v
									++cad_node_counter;					
									Point<3> cad_point_coordinates;
									m_patches[patch_itr]->GetSurface().EvaluateSurfacePoint(cad_point_coordinates, u_i, v_j);

									// Add id to point --> node. Add node to list of CAD nodes
									NodeType::Pointer new_cad_node = Node < 3 > ::Pointer(new Node<3>(cad_node_counter, cad_point_coordinates));
									list_of_cad_nodes.push_back(new_cad_node);

									// Store for cad node the corresponding cad information in separate vectors
									list_of_us_of_cad_nodes.push_back(u_i);
									list_of_vs_of_cad_nodes.push_back(v_j);
								}
							}
						}
					}

					// 2nd step: Construct KD-Tree with all cad nodes
					std::cout << "\n> Starting construction of search-tree..." << std::endl;
					boost::timer timer;
					typedef Bucket< 3, NodeType, NodeVector, NodeType::Pointer, NodeIterator, DistanceIterator > BucketType;
					typedef Tree< KDTreePartition<BucketType> > tree;
					int bucket_size = 20;
					tree nodes_tree(list_of_cad_nodes.begin(), list_of_cad_nodes.end(), bucket_size);
					std::cout << "> Time needed for constructing search-tree: " << timer.elapsed() << " s" << std::endl;





























						// loop over patches
						for (unsigned int patch_cp = 0; patch_cp < m_patches.size(); patch_cp++)
						{
							NURBSSurface surface = m_patches[patch_cp]->GetSurface();
							ControlPointVector CP = surface.GetControlPoints();

							// loop over relevant control points
							for(unsigned int cp_index=0; cp_index<CP.size(); cp_index++)
							{
								ControlPoint cp = CP[cp_index];
								if(cp.IsRelevantForMapping())
								{
									unsigned int x = cp.GetMappingMatrixId();
									// // compute coordinates of the point on surface associated to cp
									// surface.SetGrevilleAbscissae(cp_index, u_x, v_x);


									Node<3> closest_point(0,cp.getX(),cp.getY(),cp.getZ());
									// NodeType::CoordinatesArrayType closest_point_coords = ;
									// NodeType::Pointer &closest_point;// = Node < 3 > ::Pointer(new Node<3>(1, ((closest_point_coords)) ));


									// Search nearest cad neighbor of current integration point
									NodeType::Pointer nearest_point = nodes_tree.SearchNearestPoint( closest_point );

									// Recover CAD information of nearest point
									double u_x = list_of_us_of_cad_nodes[nearest_point->Id()-1];
									double v_x = list_of_vs_of_cad_nodes[nearest_point->Id()-1];










									matrix<double> R; surface.EvaluateNURBSFunctions( -1, -1, u_x, v_x, R);						
									matrix<int> control_points_ids = surface.GetRelevantControlPointsIndexes( -1, -1, u_x, v_x);

									// LHS contribution 1
									m_lhs_x(x, x) += alpha;
									m_lhs_y(x, x) += alpha;
									m_lhs_z(x, x) += alpha;

									for(unsigned int k=0; k<control_points_ids.size1(); k++)
									{
										for(unsigned int l=0; l<control_points_ids.size2(); l++)
										{
											ControlPoint cp_j = CP[control_points_ids(k,l)];
											double R_j = R(k,l);

											if(cp_j.IsRelevantForMapping())
											{
												unsigned int j = cp_j.GetMappingMatrixId();
												// LHS contribution 2
												m_lhs_x(x, j) -= alpha * R_j;
												m_lhs_y(x, j) -= alpha * R_j;
												m_lhs_z(x, j) -= alpha * R_j;
												// LHS contribution 3
												m_lhs_x(j, x) -= alpha * R_j;
												m_lhs_y(j, x) -= alpha * R_j;
												m_lhs_z(j, x) -= alpha * R_j;

												for(unsigned int m=0; m<control_points_ids.size1(); m++)
												{
													for(unsigned int n=0; n<control_points_ids.size2(); n++)
													{
														ControlPoint cp_i = CP[control_points_ids(m,n)];
														double R_i = R(m,n);
														if(cp_i.IsRelevantForMapping())
														{
															unsigned int i = cp_i.GetMappingMatrixId();
															// LHS contribution 4
															m_lhs_x(j, i) += alpha * R_j*R_i;
															m_lhs_y(j, i) += alpha * R_j*R_i;
															m_lhs_z(j, i) += alpha * R_j*R_i;
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				// previous version
					// {
					// 	double u_x, v_x;
					// 	// loop over patches
					// 	for (unsigned int patch_cp = 0; patch_cp < m_patches.size(); patch_cp++)
					// 	{
					// 		NURBSSurface surface = m_patches[patch_cp]->GetSurface();
					// 		ControlPointVector CP = surface.GetControlPoints();

					// 		// loop over relevant control points
					// 		for(unsigned int cp_index=0; cp_index<CP.size(); cp_index++)
					// 		{
					// 			ControlPoint cp = CP[cp_index];
					// 			if(cp.IsRelevantForMapping())
					// 			{
					// 				unsigned int x = cp.GetMappingMatrixId();
					// 				// compute coordinates of the point on surface associated to cp
					// 				surface.SetGrevilleAbscissae(cp_index, u_x, v_x);

					// 				matrix<double> R; surface.EvaluateNURBSFunctions( -1, -1, u_x, v_x, R);						
					// 				matrix<int> control_points_ids = surface.GetRelevantControlPointsIndexes( -1, -1, u_x, v_x);

					// 				// LHS contribution 1
					// 				m_lhs(x, x) += alpha;

					// 				for(unsigned int k=0; k<control_points_ids.size1(); k++)
					// 				{
					// 					for(unsigned int l=0; l<control_points_ids.size2(); l++)
					// 					{
					// 						ControlPoint cp_j = CP[control_points_ids(k,l)];
					// 						double R_j = R(k,l);

					// 						if(cp_j.IsRelevantForMapping())
					// 						{
					// 							unsigned int j = cp_j.GetMappingMatrixId();
					// 							// LHS contribution 2
					// 							m_lhs(x, j) -= alpha * R_j;
					// 							// LHS contribution 3
					// 							m_lhs(j, x) -= alpha * R_j;

					// 							for(unsigned int m=0; m<control_points_ids.size1(); m++)
					// 							{
					// 								for(unsigned int n=0; n<control_points_ids.size2(); n++)
					// 								{
					// 									ControlPoint cp_i = CP[control_points_ids(m,n)];
					// 									double R_i = R(m,n);
					// 									if(cp_i.IsRelevantForMapping())
					// 									{
					// 										unsigned int i = cp_i.GetMappingMatrixId();
					// 										// LHS contribution 4
					// 										m_lhs(j, i) += alpha * R_j*R_i;
					// 									}
					// 									// else
					// 									// {
					// 									// 	// RHS contribution 2
					// 									// 	m_rhs_x(j) += alpha * R_j*R_i*cp_i.getX();
					// 									// 	m_rhs_y(j) += alpha * R_j*R_i*cp_i.getY();
					// 									// 	m_rhs_z(j) += alpha * R_j*R_i*cp_i.getZ();
					// 									// }
					// 								}
					// 							}
					// 						}

					// 						// else
					// 						// {
					// 						// 	// RHS contribution 1									
					// 						// 	m_rhs_x(x) += alpha * R_j*cp_j.getX();
					// 						// 	m_rhs_y(x) += alpha * R_j*cp_j.getY();
					// 						// 	m_rhs_z(x) += alpha * R_j*cp_j.getZ();
					// 						// }
					// 					}
					// 				}
					// 			}
					// 		}
					// 	}
					// }			

				void beta_criterion(double beta)
				{
					// add contribution to lhs
						for(unsigned int i=0; i<m_n_relevant_control_points; i++)
						{
							m_lhs_x(i,i) += beta;
							m_lhs_y(i,i) += beta;
							m_lhs_z(i,i) += beta;
						}
					// add contribution to rhs
						int cp_id;
						for (PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
						{
							for (ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
							{
								if(cp_i->IsRelevantForMapping())
								{
									cp_id = cp_i->GetMappingMatrixId();
									m_rhs_x[cp_id] += beta * cp_i->getX();
									m_rhs_y[cp_id] += beta * cp_i->getY();
									m_rhs_z[cp_id] += beta * cp_i->getZ();
								}
							}
						}
				}

				void delta_criterion(double delta)
				{
					// add contribution to lhs
						for(unsigned int i=0; i<m_n_relevant_control_points; i++)
						{
							m_lhs_x(i,i) += delta;
							m_lhs_y(i,i) += delta;
							m_lhs_z(i,i) += delta;
						}
				}
			
			// boundary conditions
				void apply_boundary_conditions_small( double penalty_factor_disp, 
													double penalty_factor_rot, 
													//double penalty_factor_dirichlet, 
													//boost::python::list& edges_with_specific_dirichlet_conditions, 
													boost::python::list& edges_with_enforced_tangent_continuity
													)
					{
						std::cout << "\n> Starting to apply boundary conditions..." << std::endl;
						boost::timer function_timer;
						
						// Loop over all brep elements specifying boundary conditions 
						for (BREPElementVector::iterator brep_elem_i = m_brep_elements.begin(); brep_elem_i != m_brep_elements.end(); ++brep_elem_i)
						{
							// Check if brep_elem_i is a element used for coupling or for Dirichlet boundary conditions
							if(brep_elem_i->HasCouplingCondition())
								apply_coupling_condition_small(brep_elem_i, penalty_factor_disp, penalty_factor_rot, edges_with_enforced_tangent_continuity);
							// else if(brep_elem_i->HasDirichletCondition())
								// apply_dirichlet_condition(brep_elem_i, penalty_factor_dirichlet, edges_with_specific_dirichlet_conditions);
						}
					
						std::cout << "\n> Finished applying coupling boundary conditions in " << function_timer.elapsed() << " s." << std::endl;
					}

				void apply_coupling_condition_small( BREPElementVector::iterator &brep_elem_i, 
													double penalty_factor_disp, 
													double penalty_factor_rot, 
													boost::python::list& edges_with_enforced_tangent_continuity )
					{
						// Get Gauss points of current brep element
							BREPGaussPointVector brep_gps = brep_elem_i->GetGaussPoints();

						// Check if for current element some continuity is to be enforced
							bool tangent_continuity_to_be_enforced = false;
							double penalty_factor_tangent_continuity = 0.0;
							for (unsigned int i = 0; i < boost::python::len(edges_with_enforced_tangent_continuity); ++i)
							{
								unsigned int listed_edge_id = extractUnsignedInt(edges_with_enforced_tangent_continuity[i][0]);
								if(brep_elem_i->GetEdgeId() == listed_edge_id)
								{
									tangent_continuity_to_be_enforced = true;
									double extracted_factor = extractDouble(edges_with_enforced_tangent_continuity[i][1]);
									penalty_factor_tangent_continuity = extracted_factor;
									std::cout << "tangent continuity will been enforced on " << brep_elem_i->GetEdgeId() << std::endl;
									
								}
							}

						// Loop over all Gauss points of current brep element 
						for (BREPGaussPointVector::iterator brep_gp_i = brep_gps.begin(); brep_gp_i != brep_gps.end(); ++brep_gp_i)
						{
							// Read information from Gauss point
							unsigned int master_patch_id = brep_gp_i->GetPatchId();
							unsigned int slave_patch_id = brep_gp_i->GetSlavePatchId();
							Patch::Pointer master_patch = m_patches[m_patch_position_in_patch_vector[master_patch_id]];
							Patch::Pointer slave_patch = m_patches[m_patch_position_in_patch_vector[slave_patch_id]];
							double gp_i_weight = brep_gp_i->GetWeight();
							Vector location_on_master_patch = brep_gp_i->GetLocation();
							Vector location_on_slave_patch = brep_gp_i->GetSlaveLocation();
							Vector tangent_on_master_patch = brep_gp_i->GetTangent();
							Vector tangent_on_slave_patch = brep_gp_i->GetSlaveTangent();

							// Evaluate NURBS basis function for Gauss point on both patches and get the corresponding ids of control points in the mapping matrix
							matrix<double> R_gpi_master;
							double u_m = location_on_master_patch(0);
							double v_m = location_on_master_patch(1);
							master_patch->GetSurface().EvaluateNURBSFunctions(-1,-1,u_m, v_m, R_gpi_master);
							matrix<unsigned int> mapping_matrix_ids_gpi_master = master_patch->GetSurface().GetMappingMatrixIds(-1,-1,u_m, v_m);

							matrix<double> R_gpi_slave;
							double u_s = location_on_slave_patch(0);
							double v_s = location_on_slave_patch(1);
							slave_patch->GetSurface().EvaluateNURBSFunctions(-1,-1,u_s, v_s, R_gpi_slave);	
							matrix<unsigned int> mapping_matrix_ids_gpi_slave = slave_patch->GetSurface().GetMappingMatrixIds(-1,-1,u_s, v_s);							

							// Compute Jacobian J1
							matrix<double> g_master = master_patch->GetSurface().GetBaseVectors(-1,-1,u_m,v_m);
							Vector g1 = ZeroVector(3);
							g1(0) = g_master(0,0);
							g1(1) = g_master(1,0);
							g1(2) = g_master(2,0);
							Vector g2 = ZeroVector(3);
							g2(0) = g_master(0,1);
							g2(1) = g_master(1,1);
							g2(2) = g_master(2,1);
							double J1 = norm_2( g1* tangent_on_master_patch(0) + g2* tangent_on_master_patch(1) );

							// First we introduce coupling of displacements
							apply_displacement_coupling_small( mapping_matrix_ids_gpi_master, 
															mapping_matrix_ids_gpi_slave, 
															R_gpi_master, 
															R_gpi_slave, 
															J1,
															gp_i_weight,
															penalty_factor_disp );

							// Then check if for current element, tangent continuity is to be enforced. If yes, we enforce the tangent continuity..
							if(tangent_continuity_to_be_enforced)
							{
								enforce_tangent_continuity_small( master_patch, 
																slave_patch,
																u_m, v_m,
																u_s, v_s,
																tangent_on_master_patch, 
																tangent_on_slave_patch,
																mapping_matrix_ids_gpi_master, 
																mapping_matrix_ids_gpi_slave,
																J1,
																gp_i_weight,
																penalty_factor_tangent_continuity );
								std::cout << "tangent continuity has been enforced on " << brep_elem_i->GetEdgeId() << std::endl;
							}
							// ...if no, we introduce coupling of rotations
							else
								apply_rotation_coupling_small( master_patch, 
															slave_patch,
															u_m, v_m,
															u_s, v_s,
															tangent_on_master_patch, 
															tangent_on_slave_patch,
															mapping_matrix_ids_gpi_master, 
															mapping_matrix_ids_gpi_slave,
															J1,
															gp_i_weight,
															penalty_factor_rot );
						}
					}
					
				void apply_displacement_coupling_small( matrix<unsigned int> &mapping_matrix_ids_gpi_master,
														matrix<unsigned int> &mapping_matrix_ids_gpi_slave,
														matrix<double> &R_gpi_master, 
														matrix<double> &R_gpi_slave, 
														double J1,
														double gp_i_weight,
														double penalty_factor_disp )
					{	
						// First we consider the relation Master-Master ( MM )
						for(unsigned int i=0; i<mapping_matrix_ids_gpi_master.size2();i++)
						{
							for(unsigned int j=0; j<mapping_matrix_ids_gpi_master.size1();j++)
							{
								unsigned int R_row_id = mapping_matrix_ids_gpi_master(j,i);
								double R_row = R_gpi_master(j,i);

								for(unsigned int k=0; k<mapping_matrix_ids_gpi_master.size2();k++)
								{
									for(unsigned int l=0; l<mapping_matrix_ids_gpi_master.size1();l++)
									{
										unsigned int R_coll_id = mapping_matrix_ids_gpi_master(l,k);
										double R_coll = R_gpi_master(l,k);

										m_lhs_x(R_row_id, R_coll_id) += penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
										m_lhs_y(R_row_id, R_coll_id) += penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
										m_lhs_z(R_row_id, R_coll_id) += penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
									}
								}
							}
						}

						// Then we consider the relation Slave-Slave ( SS )
						for(unsigned int i=0; i<mapping_matrix_ids_gpi_slave.size2();i++)
						{
							for(unsigned int j=0; j<mapping_matrix_ids_gpi_slave.size1();j++)
							{
								unsigned int R_row_id = mapping_matrix_ids_gpi_slave(j,i);
								double R_row = R_gpi_slave(j,i);

								for(unsigned int k=0; k<mapping_matrix_ids_gpi_slave.size2();k++)
								{
									for(unsigned int l=0; l<mapping_matrix_ids_gpi_slave.size1();l++)
									{
										unsigned int R_coll_id = mapping_matrix_ids_gpi_slave(l,k);
										double R_coll = R_gpi_slave(l,k);

										m_lhs_x(R_row_id, R_coll_id) += penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
										m_lhs_y(R_row_id, R_coll_id) += penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
										m_lhs_z(R_row_id, R_coll_id) += penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
									}
								}
							}
						}			

						// Then we consider the Master-Slave relation ( MS & SM )
						for(unsigned int i=0; i<mapping_matrix_ids_gpi_master.size2();i++)
						{
							for(unsigned int j=0; j<mapping_matrix_ids_gpi_master.size1();j++)
							{
								unsigned int R_row_id = mapping_matrix_ids_gpi_master(j,i);
								double R_row = R_gpi_master(j,i);

								for(unsigned int k=0; k<mapping_matrix_ids_gpi_slave.size2();k++)
								{
									for(unsigned int l=0; l<mapping_matrix_ids_gpi_slave.size1();l++)
									{
										unsigned int R_coll_id = mapping_matrix_ids_gpi_slave(l,k);
										double R_coll = R_gpi_slave(l,k);

										// MS 
										m_lhs_x(R_row_id, R_coll_id) -= penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
										m_lhs_y(R_row_id, R_coll_id) -= penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
										m_lhs_z(R_row_id, R_coll_id) -= penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;

										// SM
										m_lhs_x(R_coll_id, R_row_id) -= penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
										m_lhs_y(R_coll_id, R_row_id) -= penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
										m_lhs_z(R_coll_id, R_row_id) -= penalty_factor_disp * gp_i_weight * J1 * R_row * R_coll;
									}
								}
							}
						}
					}
				void apply_rotation_coupling_small( Patch::Pointer master_patch,
													Patch::Pointer slave_patch,
													double u_m, double v_m,
													double u_s, double v_s,
													Vector &tangent_on_master_patch,
													Vector &tangent_on_slave_patch,
													matrix<unsigned int> &mapping_matrix_ids_gpi_master,
													matrix<unsigned int> &mapping_matrix_ids_gpi_slave,
													double J1,
													double gp_i_weight,
													double penalty_factor_rot )
					{		
						// Variables needed later
							Vector T1_m, T1_s, T2_m, T2_s, T3_m, T3_s;
							std::vector<Vector> t1r_m, t1r_s, t2r_m, t2r_s, t3r_m, t3r_s;				

						// Compute geometric quantities
							master_patch->GetSurface().ComputeVariationOfLocalCSY( u_m, v_m, tangent_on_master_patch, T1_m, T2_m, T3_m, t1r_m, t2r_m, t3r_m );
							slave_patch->GetSurface().ComputeVariationOfLocalCSY( u_s, v_s, tangent_on_slave_patch, T1_s, T2_s, T3_s, t1r_s, t2r_s, t3r_s );

						// Check if master and slave tangent point in same direction. If yes, we have to subtract in the following.
							int sign_factor = 1;
							if( inner_prod(T2_m,T2_s) > 0 )
								sign_factor = -1;

						// Merge boundary conditions into mapping matrix

						// First we consider the relation Master-Master ( MM )
							unsigned int k_coll = 0;
							for(unsigned int i=0; i<mapping_matrix_ids_gpi_master.size2();i++)
							{
								for(unsigned int j=0; j<mapping_matrix_ids_gpi_master.size1();j++)
								{
									unsigned int R_row_id = mapping_matrix_ids_gpi_master(j,i);
									Vector omega_mx_coll = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_coll+0]);
									Vector omega_my_coll = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_coll+1]);
									Vector omega_mz_coll = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_coll+2]);
									double omega_T2_mx_coll = inner_prod(omega_mx_coll,T2_m);
									double omega_T2_my_coll = inner_prod(omega_my_coll,T2_m);
									double omega_T2_mz_coll = inner_prod(omega_mz_coll,T2_m);

									unsigned int k_row = 0;
									for(unsigned int k=0; k<mapping_matrix_ids_gpi_master.size2();k++)
									{
										for(unsigned int l=0; l<mapping_matrix_ids_gpi_master.size1();l++)
										{
											unsigned int R_coll_id = mapping_matrix_ids_gpi_master(l,k);
											Vector omega_mx_row = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_row+0]);
											Vector omega_my_row = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_row+1]);
											Vector omega_mz_row = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_row+2]);
											double omega_T2_mx_row = inner_prod(omega_mx_row,T2_m);
											double omega_T2_my_row = inner_prod(omega_my_row,T2_m);
											double omega_T2_mz_row = inner_prod(omega_mz_row,T2_m);

											m_lhs_x(R_row_id, R_coll_id) += penalty_factor_rot * gp_i_weight * J1 * omega_T2_mx_row * omega_T2_mx_coll;
											m_lhs_y(R_row_id, R_coll_id) += penalty_factor_rot * gp_i_weight * J1 * omega_T2_my_row * omega_T2_my_coll;
											m_lhs_z(R_row_id, R_coll_id) += penalty_factor_rot * gp_i_weight * J1 * omega_T2_mz_row * omega_T2_mz_coll;
											
											k_row++;
										}
									}
									k_coll++;
								}
							}

						// Then we consider the relation Slave-Slave ( SS )
							k_coll = 0;
							for(unsigned int i=0; i<mapping_matrix_ids_gpi_slave.size2();i++)
							{
								for(unsigned int j=0; j<mapping_matrix_ids_gpi_slave.size1();j++)
								{
									unsigned int R_row_id = mapping_matrix_ids_gpi_slave(j,i);
									Vector omega_sx_coll = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_coll+0]);
									Vector omega_sy_coll = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_coll+1]);
									Vector omega_sz_coll = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_coll+2]);
									double omega_T2_sx_coll = inner_prod(omega_sx_coll,T2_s);
									double omega_T2_sy_coll = inner_prod(omega_sy_coll,T2_s);
									double omega_T2_sz_coll = inner_prod(omega_sz_coll,T2_s);

									unsigned int k_row = 0;
									for(unsigned int k=0; k<mapping_matrix_ids_gpi_slave.size2();k++)
									{
										for(unsigned int l=0; l<mapping_matrix_ids_gpi_slave.size1();l++)
										{
											unsigned int R_coll_id = mapping_matrix_ids_gpi_slave(l,k);
											Vector omega_sx_row = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_row+0]);
											Vector omega_sy_row = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_row+1]);
											Vector omega_sz_row = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_row+2]);
											double omega_T2_sx_row = inner_prod(omega_sx_row,T2_s);
											double omega_T2_sy_row = inner_prod(omega_sy_row,T2_s);
											double omega_T2_sz_row = inner_prod(omega_sz_row,T2_s);

											m_lhs_x(R_row_id, R_coll_id) += penalty_factor_rot * gp_i_weight * J1 * omega_T2_sx_row * omega_T2_sx_coll;
											m_lhs_y(R_row_id, R_coll_id) += penalty_factor_rot * gp_i_weight * J1 * omega_T2_sy_row * omega_T2_sy_coll;
											m_lhs_z(R_row_id, R_coll_id) += penalty_factor_rot * gp_i_weight * J1 * omega_T2_sz_row * omega_T2_sz_coll;
											
											k_row++;
										}
									}
									k_coll++;
								}
							}			

						// Then we consider the Master-slave relation ( MS & SM )
							unsigned int k_m = 0;
							for(unsigned int i=0; i<mapping_matrix_ids_gpi_master.size2();i++)
							{
								for(unsigned int j=0; j<mapping_matrix_ids_gpi_master.size1();j++)
								{
									unsigned int R_m_id = mapping_matrix_ids_gpi_master(j,i);
									Vector omega_mx = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_m+0]);
									Vector omega_my = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_m+1]);
									Vector omega_mz = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_m+2]);
									double omega_T2_mx = inner_prod(omega_mx,T2_m);
									double omega_T2_my = inner_prod(omega_my,T2_m);
									double omega_T2_mz = inner_prod(omega_mz,T2_m);

									unsigned int k_s = 0;
									for(unsigned int k=0; k<mapping_matrix_ids_gpi_slave.size2();k++)
									{
										for(unsigned int l=0; l<mapping_matrix_ids_gpi_slave.size1();l++)
										{
											unsigned int R_s_id = mapping_matrix_ids_gpi_slave(l,k);
											Vector omega_sx = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_s+0]);
											Vector omega_sy = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_s+1]);
											Vector omega_sz = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_s+2]);
											double omega_T2_sx = inner_prod(omega_sx,T2_s);
											double omega_T2_sy = inner_prod(omega_sy,T2_s);
											double omega_T2_sz = inner_prod(omega_sz,T2_s);

											// MS
											m_lhs_x(R_m_id, R_s_id) += sign_factor * penalty_factor_rot * gp_i_weight * J1 * omega_T2_mx * omega_T2_sx;
											m_lhs_y(R_m_id, R_s_id) += sign_factor * penalty_factor_rot * gp_i_weight * J1 * omega_T2_my * omega_T2_sy;
											m_lhs_z(R_m_id, R_s_id) += sign_factor * penalty_factor_rot * gp_i_weight * J1 * omega_T2_mz * omega_T2_sz;

											// SM
											m_lhs_x(R_s_id, R_m_id) += sign_factor * penalty_factor_rot * gp_i_weight * J1 * omega_T2_mx * omega_T2_sx;
											m_lhs_y(R_s_id, R_m_id) += sign_factor * penalty_factor_rot * gp_i_weight * J1 * omega_T2_my * omega_T2_sy;
											m_lhs_z(R_s_id, R_m_id) += sign_factor * penalty_factor_rot * gp_i_weight * J1 * omega_T2_mz * omega_T2_sz;						

											k_s++;								
										}
									}
									k_m++;
								}
							}
					}
				void enforce_tangent_continuity_small( Patch::Pointer master_patch,
												Patch::Pointer slave_patch,
												double u_m, double v_m,
												double u_s, double v_s,
												Vector& tangent_on_master_patch,
												Vector& tangent_on_slave_patch,
												matrix<unsigned int>& mapping_matrix_ids_gpi_master,
												matrix<unsigned int>& mapping_matrix_ids_gpi_slave,
												double J1,
												double gp_i_weight,
												double penalty_factor_tangent_continuity )
					{
						// Variables needed later
							Vector T1_m, T1_s, T2_m, T2_s, T3_m, T3_s;
							Vector T1_der_m, T1_der_s, T2_der_m, T2_der_s, T3_der_m, T3_der_s;
							std::vector<Vector> t1r_m, t1r_s, t2r_m, t2r_s, t3r_m, t3r_s;
							std::vector<Vector> t1_der_r_m, t1_der_r_s, t2_der_r_m, t2_der_r_s, t3_der_r_m, t3_der_r_s;					
							std::vector<std::vector<Vector>> t1rs_m, t1rs_s, t2rs_m, t2rs_s, t3rs_m, t3rs_s;
							std::vector<std::vector<Vector>> t1_der_rs_m, t1_der_rs_s, t2_der_rs_m, t2_der_rs_s, t3_der_rs_m, t3_der_rs_s;

						std::cout << "Called: cad_mapper::enforce_tangent_continuity()" << std::endl;

						// Compute geometric quantities
							master_patch->GetSurface().ComputeSecondVariationOfLocalCSY( u_m, v_m, 
																						tangent_on_master_patch, 
																						T1_m, T2_m, T3_m, 
																						T1_der_m, T2_der_m, T3_der_m,
																						t1r_m, t2r_m, t3r_m,
																						t1_der_r_m, t2_der_r_m, t3_der_r_m,
																						t1rs_m, t2rs_m, t3rs_m,
																						t1_der_rs_m, t2_der_rs_m, t3_der_rs_m );
							slave_patch->GetSurface().ComputeSecondVariationOfLocalCSY( u_s, v_s, 
																					tangent_on_slave_patch, 
																					T1_s, T2_s, T3_s, 
																					T1_der_s, T2_der_s, T3_der_s,
																					t1r_s, t2r_s, t3r_s,
																					t1_der_r_s, t2_der_r_s, t3_der_r_s,
																					t1rs_s, t2rs_s, t3rs_s,
																					t1_der_rs_s, t2_der_rs_s, t3_der_rs_s );	

							
							double fac = - inner_prod(T3_m, T1_s);
							KRATOS_WATCH(fac);

						// First we consider contribution to the m_mapping_rhs_vector

						// Master-Master-relation ( MM )
							unsigned int k_row = 0;
							for(unsigned int i=0; i<mapping_matrix_ids_gpi_master.size2();i++)
							{
								for(unsigned int j=0; j<mapping_matrix_ids_gpi_master.size1();j++)
								{
									unsigned int R_row_id = mapping_matrix_ids_gpi_master(j,i);

									m_rhs_x(R_row_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * fac * inner_prod(t3r_m[3*k_row+0],T1_s);
									m_rhs_y(R_row_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * fac * inner_prod(t3r_m[3*k_row+1],T1_s);
									m_rhs_z(R_row_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * fac * inner_prod(t3r_m[3*k_row+2],T1_s);

									k_row++;
								}
							}

						// Slave-Slave-relation ( SS )
							k_row = 0;
							for(unsigned int i=0; i<mapping_matrix_ids_gpi_slave.size2();i++)
							{
								for(unsigned int j=0; j<mapping_matrix_ids_gpi_slave.size1();j++)
								{
									unsigned int R_row_id = mapping_matrix_ids_gpi_slave(j,i);

									m_rhs_x(R_row_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * fac * inner_prod(t1r_s[3*k_row+0],T3_m);
									m_rhs_y(R_row_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * fac * inner_prod(t1r_s[3*k_row+1],T3_m);
									m_rhs_z(R_row_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * fac * inner_prod(t1r_s[3*k_row+2],T3_m);

									k_row++;
								}
							}	

						// Then we consider the contribution to the m_mapping_matrix_CAD_CAD

						// Master-Master-relation ( MM )
							k_row = 0;
							for(unsigned int i=0; i<mapping_matrix_ids_gpi_master.size2();i++)
							{
								for(unsigned int j=0; j<mapping_matrix_ids_gpi_master.size1();j++)
								{
									unsigned int R_row_id = mapping_matrix_ids_gpi_master(j,i);

									unsigned int k_coll = 0;
									for(unsigned int k=0; k<mapping_matrix_ids_gpi_master.size2();k++)
									{
										for(unsigned int l=0; l<mapping_matrix_ids_gpi_master.size1();l++)
										{
											unsigned int R_coll_id = mapping_matrix_ids_gpi_master(l,k);
											
											double term_1_x = inner_prod(t3r_m[3*k_coll+0],T1_s) * inner_prod(t3r_m[3*k_row+0],T1_s);
											double term_1_y = inner_prod(t3r_m[3*k_coll+1],T1_s) * inner_prod(t3r_m[3*k_row+1],T1_s);
											double term_1_z = inner_prod(t3r_m[3*k_coll+2],T1_s) * inner_prod(t3r_m[3*k_row+2],T1_s);

											double term_2_x = fac * inner_prod(t3rs_m[3*k_row+0][3*k_coll+0],T1_s);
											double term_2_y = fac * inner_prod(t3rs_m[3*k_row+1][3*k_coll+1],T1_s);
											double term_2_z = fac * inner_prod(t3rs_m[3*k_row+2][3*k_coll+2],T1_s);

											m_lhs_x(R_row_id, R_coll_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_x + term_2_x );
											m_lhs_y(R_row_id, R_coll_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_y + term_2_y );
											m_lhs_z(R_row_id, R_coll_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_z + term_2_z );					

											k_coll++;								
										}
									}
									k_row++;
								}
							}

						// Slave-Slave-relation ( SS )
							k_row = 0;
							for(unsigned int i=0; i<mapping_matrix_ids_gpi_slave.size2();i++)
							{
								for(unsigned int j=0; j<mapping_matrix_ids_gpi_slave.size1();j++)
								{
									unsigned int R_row_id = mapping_matrix_ids_gpi_slave(j,i);

									unsigned int k_coll = 0;
									for(unsigned int k=0; k<mapping_matrix_ids_gpi_slave.size2();k++)
									{
										for(unsigned int l=0; l<mapping_matrix_ids_gpi_slave.size1();l++)
										{
											unsigned int R_coll_id = mapping_matrix_ids_gpi_slave(l,k);
											
											double term_1_x = inner_prod(t1r_s[3*k_coll+0],T3_m) * inner_prod(T3_m, t1r_s[3*k_row+0]);
											double term_1_y = inner_prod(t1r_s[3*k_coll+1],T3_m) * inner_prod(T3_m, t1r_s[3*k_row+1]);
											double term_1_z = inner_prod(t1r_s[3*k_coll+2],T3_m) * inner_prod(T3_m, t1r_s[3*k_row+2]);

											double term_2_x = fac * inner_prod(T3_m, t1rs_s[3*k_row+0][3*k_coll+0]);
											double term_2_y = fac * inner_prod(T3_m, t1rs_s[3*k_row+1][3*k_coll+1]);
											double term_2_z = fac * inner_prod(T3_m, t1rs_s[3*k_row+2][3*k_coll+2]);

											m_lhs_x(R_row_id, R_coll_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_x + term_2_x );
											m_lhs_y(R_row_id, R_coll_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_y + term_2_y );
											m_lhs_z(R_row_id, R_coll_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_z + term_2_z );					

											k_coll++;								
										}
									}
									k_row++;
								}
							}	

						// Master-slave-relation ( MS )
							k_row = 0;
							for(unsigned int i=0; i<mapping_matrix_ids_gpi_master.size2();i++)
							{
								for(unsigned int j=0; j<mapping_matrix_ids_gpi_master.size1();j++)
								{
									unsigned int R_row_id = mapping_matrix_ids_gpi_master(j,i);

									unsigned int k_coll = 0;
									for(unsigned int k=0; k<mapping_matrix_ids_gpi_slave.size2();k++)
									{
										for(unsigned int l=0; l<mapping_matrix_ids_gpi_slave.size1();l++)
										{
											unsigned int R_coll_id = mapping_matrix_ids_gpi_slave(l,k);

											double term_1_x = inner_prod(T3_m,t1r_s[3*k_coll+0]) * inner_prod(t3r_m[3*k_row+0],T1_s);
											double term_1_y = inner_prod(T3_m,t1r_s[3*k_coll+1]) * inner_prod(t3r_m[3*k_row+1],T1_s);
											double term_1_z = inner_prod(T3_m,t1r_s[3*k_coll+2]) * inner_prod(t3r_m[3*k_row+2],T1_s);

											double term_2_x = fac * inner_prod(t3r_m[3*k_row+0],t1r_s[3*k_coll+0]);
											double term_2_y = fac * inner_prod(t3r_m[3*k_row+1],t1r_s[3*k_coll+1]);
											double term_2_z = fac * inner_prod(t3r_m[3*k_row+2],t1r_s[3*k_coll+2]);

											m_lhs_x(R_row_id, R_coll_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_x + term_2_x );
											m_lhs_y(R_row_id, R_coll_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_y + term_2_y );
											m_lhs_z(R_row_id, R_coll_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_z + term_2_z );											

											k_coll++;								
										}
									}
									k_row++;
								}
							}

						// Master-slave-relation ( SM )
							k_row = 0;
							for(unsigned int i=0; i<mapping_matrix_ids_gpi_slave.size2();i++)
							{
								for(unsigned int j=0; j<mapping_matrix_ids_gpi_slave.size1();j++)
								{
									unsigned int R_row_id = mapping_matrix_ids_gpi_slave(j,i);

									unsigned int k_coll = 0;
									for(unsigned int k=0; k<mapping_matrix_ids_gpi_master.size2();k++)
									{
										for(unsigned int l=0; l<mapping_matrix_ids_gpi_master.size1();l++)
										{
											unsigned int R_coll_id = mapping_matrix_ids_gpi_master(l,k);

											double term_1_y = inner_prod(t3r_m[3*k_coll+1], T1_s) * inner_prod(T3_m, t1r_s[3*k_row+1]);
											double term_1_z = inner_prod(t3r_m[3*k_coll+2], T1_s) * inner_prod(T3_m, t1r_s[3*k_row+2]);
											double term_1_x = inner_prod(t3r_m[3*k_coll+0], T1_s) * inner_prod(T3_m, t1r_s[3*k_row+0]);

											double term_2_x = fac * inner_prod(t3r_m[3*k_coll+0], t1r_s[3*k_row+0]);
											double term_2_y = fac * inner_prod(t3r_m[3*k_coll+1], t1r_s[3*k_row+1]);
											double term_2_z = fac * inner_prod(t3r_m[3*k_coll+2], t1r_s[3*k_row+2]);

											m_lhs_x(R_row_id, R_coll_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_x + term_2_x );
											m_lhs_y(R_row_id, R_coll_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_y + term_2_y );
											m_lhs_z(R_row_id, R_coll_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * ( term_1_z + term_2_z );											

											k_coll++;								
										}
									}
									k_row++;
								}
							}		
					}

				void enforce_tangent_continuity_small_2( Patch::Pointer master_patch,
												Patch::Pointer slave_patch,
												double u_m, double v_m,
												double u_s, double v_s,
												Vector& tangent_on_master_patch,
												Vector& tangent_on_slave_patch,
												matrix<unsigned int>& mapping_matrix_ids_gpi_master,
												matrix<unsigned int>& mapping_matrix_ids_gpi_slave,
												double J1,
												double gp_i_weight,
												double penalty_factor_tangent_continuity )
					{
						// Variables needed later
							Vector T1_m, T1_s, T2_m, T2_s, T3_m, T3_s;
							Vector T1_der_m, T1_der_s, T2_der_m, T2_der_s, T3_der_m, T3_der_s;
							std::vector<Vector> t1r_m, t1r_s, t2r_m, t2r_s, t3r_m, t3r_s;
							std::vector<Vector> t1_der_r_m, t1_der_r_s, t2_der_r_m, t2_der_r_s, t3_der_r_m, t3_der_r_s;					
							std::vector<std::vector<Vector>> t1rs_m, t1rs_s, t2rs_m, t2rs_s, t3rs_m, t3rs_s;
							std::vector<std::vector<Vector>> t1_der_rs_m, t1_der_rs_s, t2_der_rs_m, t2_der_rs_s, t3_der_rs_m, t3_der_rs_s;

						std::cout << "Called: cad_mapper::enforce_tangent_continuity()" << std::endl;

						// Compute geometric quantities
							master_patch->GetSurface().ComputeSecondVariationOfLocalCSY( u_m, v_m, 
																						tangent_on_master_patch, 
																						T1_m, T2_m, T3_m, 
																						T1_der_m, T2_der_m, T3_der_m,
																						t1r_m, t2r_m, t3r_m,
																						t1_der_r_m, t2_der_r_m, t3_der_r_m,
																						t1rs_m, t2rs_m, t3rs_m,
																						t1_der_rs_m, t2_der_rs_m, t3_der_rs_m );
							slave_patch->GetSurface().ComputeSecondVariationOfLocalCSY( u_s, v_s, 
																					tangent_on_slave_patch, 
																					T1_s, T2_s, T3_s, 
																					T1_der_s, T2_der_s, T3_der_s,
																					t1r_s, t2r_s, t3r_s,
																					t1_der_r_s, t2_der_r_s, t3_der_r_s,
																					t1rs_s, t2rs_s, t3rs_s,
																					t1_der_rs_s, t2_der_rs_s, t3_der_rs_s );	

							
							double cosine = inner_prod(T3_m, T3_s);
							double delta_omega = - acos(cosine);
							KRATOS_WATCH(delta_omega);
							// double delta_omega = 0.2;
							int sign_factor = 1;
							if( inner_prod(T2_m,T2_s) > 0 )
								sign_factor = -1;

						// First we consider contribution to the m_mapping_rhs_vector

						// Master-Master-relation ( MM )
							unsigned int k_row = 0;
							for(unsigned int k=0; k<mapping_matrix_ids_gpi_master.size2();k++)
							{
								for(unsigned int l=0; l<mapping_matrix_ids_gpi_master.size1();l++)
								{
									unsigned int R_row_id = mapping_matrix_ids_gpi_master(l,k);
									Vector omega_mx_row = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_row+0]);
									Vector omega_my_row = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_row+1]);
									Vector omega_mz_row = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_row+2]);
									double omega_T2_mx_row = inner_prod(omega_mx_row,T2_m);
									double omega_T2_my_row = inner_prod(omega_my_row,T2_m);
									double omega_T2_mz_row = inner_prod(omega_mz_row,T2_m);

									m_rhs_x(R_row_id) -= penalty_factor_tangent_continuity * gp_i_weight * J1 * delta_omega * omega_T2_mx_row;
									m_rhs_y(R_row_id) -= penalty_factor_tangent_continuity * gp_i_weight * J1 * delta_omega * omega_T2_my_row;
									m_rhs_z(R_row_id) -= penalty_factor_tangent_continuity * gp_i_weight * J1 * delta_omega * omega_T2_mz_row;
									
									k_row++;
								}
							}

						// Slave-Slave-relation ( SS )
							k_row = 0;
							for(unsigned int k=0; k<mapping_matrix_ids_gpi_slave.size2();k++)
							{
								for(unsigned int l=0; l<mapping_matrix_ids_gpi_slave.size1();l++)
								{
									unsigned int R_row_id = mapping_matrix_ids_gpi_slave(l,k);
									Vector omega_sx_row = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_row+0]);
									Vector omega_sy_row = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_row+1]);
									Vector omega_sz_row = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_row+2]);
									double omega_T2_sx_row = inner_prod(omega_sx_row,T2_s);
									double omega_T2_sy_row = inner_prod(omega_sy_row,T2_s);
									double omega_T2_sz_row = inner_prod(omega_sz_row,T2_s);

									m_rhs_x(R_row_id) -= sign_factor * penalty_factor_tangent_continuity * gp_i_weight * J1 * delta_omega * omega_T2_sx_row;
									m_rhs_y(R_row_id) -= sign_factor * penalty_factor_tangent_continuity * gp_i_weight * J1 * delta_omega * omega_T2_sy_row;
									m_rhs_z(R_row_id) -= sign_factor * penalty_factor_tangent_continuity * gp_i_weight * J1 * delta_omega * omega_T2_sz_row;
									
									k_row++;
								}
							}

						// Then we consider the contribution to the m_mapping_matrix_CAD_CAD
						// First we consider the relation Master-Master ( MM )
							unsigned int k_coll = 0;
							for(unsigned int i=0; i<mapping_matrix_ids_gpi_master.size2();i++)
							{
								for(unsigned int j=0; j<mapping_matrix_ids_gpi_master.size1();j++)
								{
									unsigned int R_row_id = mapping_matrix_ids_gpi_master(j,i);
									Vector omega_mx_coll = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_coll+0]);
									Vector omega_my_coll = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_coll+1]);
									Vector omega_mz_coll = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_coll+2]);
									double omega_T2_mx_coll = inner_prod(omega_mx_coll,T2_m);
									double omega_T2_my_coll = inner_prod(omega_my_coll,T2_m);
									double omega_T2_mz_coll = inner_prod(omega_mz_coll,T2_m);

									unsigned int k_row = 0;
									for(unsigned int k=0; k<mapping_matrix_ids_gpi_master.size2();k++)
									{
										for(unsigned int l=0; l<mapping_matrix_ids_gpi_master.size1();l++)
										{
											unsigned int R_coll_id = mapping_matrix_ids_gpi_master(l,k);
											Vector omega_mx_row = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_row+0]);
											Vector omega_my_row = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_row+1]);
											Vector omega_mz_row = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_row+2]);
											double omega_T2_mx_row = inner_prod(omega_mx_row,T2_m);
											double omega_T2_my_row = inner_prod(omega_my_row,T2_m);
											double omega_T2_mz_row = inner_prod(omega_mz_row,T2_m);

											Vector omega_mx_row_col = MathUtils<double>::CrossProduct(T3_m,t3rs_m[3*k_row+0][3*k_coll+0]);
											Vector omega_my_row_col = MathUtils<double>::CrossProduct(T3_m,t3rs_m[3*k_row+1][3*k_coll+1]);
											Vector omega_mz_row_col = MathUtils<double>::CrossProduct(T3_m,t3rs_m[3*k_row+2][3*k_coll+2]);
											double omega_T2_mx_row_col = inner_prod(omega_mx_row_col,T2_m);
											double omega_T2_my_row_col = inner_prod(omega_my_row_col,T2_m);
											double omega_T2_mz_row_col = inner_prod(omega_mz_row_col,T2_m);

											m_lhs_x(R_row_id, R_coll_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * (omega_T2_mx_row * omega_T2_mx_coll - delta_omega * omega_T2_mx_row_col);
											m_lhs_y(R_row_id, R_coll_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * (omega_T2_my_row * omega_T2_my_coll - delta_omega * omega_T2_my_row_col);
											m_lhs_z(R_row_id, R_coll_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * (omega_T2_mz_row * omega_T2_mz_coll - delta_omega * omega_T2_mz_row_col);
											
											k_row++;
										}
									}
									k_coll++;
								}
							}

						// Then we consider the relation Slave-Slave ( SS )
							k_coll = 0;
							for(unsigned int i=0; i<mapping_matrix_ids_gpi_slave.size2();i++)
							{
								for(unsigned int j=0; j<mapping_matrix_ids_gpi_slave.size1();j++)
								{
									unsigned int R_row_id = mapping_matrix_ids_gpi_slave(j,i);
									Vector omega_sx_coll = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_coll+0]);
									Vector omega_sy_coll = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_coll+1]);
									Vector omega_sz_coll = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_coll+2]);
									double omega_T2_sx_coll = inner_prod(omega_sx_coll,T2_s);
									double omega_T2_sy_coll = inner_prod(omega_sy_coll,T2_s);
									double omega_T2_sz_coll = inner_prod(omega_sz_coll,T2_s);

									unsigned int k_row = 0;
									for(unsigned int k=0; k<mapping_matrix_ids_gpi_slave.size2();k++)
									{
										for(unsigned int l=0; l<mapping_matrix_ids_gpi_slave.size1();l++)
										{
											unsigned int R_coll_id = mapping_matrix_ids_gpi_slave(l,k);
											Vector omega_sx_row = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_row+0]);
											Vector omega_sy_row = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_row+1]);
											Vector omega_sz_row = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_row+2]);
											double omega_T2_sx_row = inner_prod(omega_sx_row,T2_s);
											double omega_T2_sy_row = inner_prod(omega_sy_row,T2_s);
											double omega_T2_sz_row = inner_prod(omega_sz_row,T2_s);

											Vector omega_sx_row_col = MathUtils<double>::CrossProduct(T3_s,t3rs_s[3*k_row+0][3*k_coll+0]);
											Vector omega_sy_row_col = MathUtils<double>::CrossProduct(T3_s,t3rs_s[3*k_row+1][3*k_coll+1]);
											Vector omega_sz_row_col = MathUtils<double>::CrossProduct(T3_s,t3rs_s[3*k_row+2][3*k_coll+2]);
											double omega_T2_sx_row_col = inner_prod(omega_sx_row_col,T2_s);
											double omega_T2_sy_row_col = inner_prod(omega_sy_row_col,T2_s);
											double omega_T2_sz_row_col = inner_prod(omega_sz_row_col,T2_s);

											m_lhs_x(R_row_id, R_coll_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * (omega_T2_sx_row * omega_T2_sx_coll - sign_factor * delta_omega * omega_T2_sx_row_col);
											m_lhs_y(R_row_id, R_coll_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * (omega_T2_sy_row * omega_T2_sy_coll - sign_factor * delta_omega * omega_T2_sy_row_col);
											m_lhs_z(R_row_id, R_coll_id) += penalty_factor_tangent_continuity * gp_i_weight * J1 * (omega_T2_sz_row * omega_T2_sz_coll - sign_factor * delta_omega * omega_T2_sz_row_col);
											
											k_row++;
										}
									}
									k_coll++;
								}
							}			

						// Then we consider the Master-slave relation ( MS & SM )
							unsigned int k_m = 0;
							for(unsigned int i=0; i<mapping_matrix_ids_gpi_master.size2();i++)
							{
								for(unsigned int j=0; j<mapping_matrix_ids_gpi_master.size1();j++)
								{
									unsigned int R_m_id = mapping_matrix_ids_gpi_master(j,i);
									Vector omega_mx = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_m+0]);
									Vector omega_my = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_m+1]);
									Vector omega_mz = MathUtils<double>::CrossProduct(T3_m,t3r_m[3*k_m+2]);
									double omega_T2_mx = inner_prod(omega_mx,T2_m);
									double omega_T2_my = inner_prod(omega_my,T2_m);
									double omega_T2_mz = inner_prod(omega_mz,T2_m);

									unsigned int k_s = 0;
									for(unsigned int k=0; k<mapping_matrix_ids_gpi_slave.size2();k++)
									{
										for(unsigned int l=0; l<mapping_matrix_ids_gpi_slave.size1();l++)
										{
											unsigned int R_s_id = mapping_matrix_ids_gpi_slave(l,k);
											Vector omega_sx = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_s+0]);
											Vector omega_sy = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_s+1]);
											Vector omega_sz = MathUtils<double>::CrossProduct(T3_s,t3r_s[3*k_s+2]);
											double omega_T2_sx = inner_prod(omega_sx,T2_s);
											double omega_T2_sy = inner_prod(omega_sy,T2_s);
											double omega_T2_sz = inner_prod(omega_sz,T2_s);

											// MS
											m_lhs_x(R_m_id, R_s_id) += sign_factor * penalty_factor_tangent_continuity * gp_i_weight * J1 * omega_T2_mx * omega_T2_sx;
											m_lhs_y(R_m_id, R_s_id) += sign_factor * penalty_factor_tangent_continuity * gp_i_weight * J1 * omega_T2_my * omega_T2_sy;
											m_lhs_z(R_m_id, R_s_id) += sign_factor * penalty_factor_tangent_continuity * gp_i_weight * J1 * omega_T2_mz * omega_T2_sz;

											// SM
											m_lhs_x(R_s_id, R_m_id) += sign_factor * penalty_factor_tangent_continuity * gp_i_weight * J1 * omega_T2_mx * omega_T2_sx;
											m_lhs_y(R_s_id, R_m_id) += sign_factor * penalty_factor_tangent_continuity * gp_i_weight * J1 * omega_T2_my * omega_T2_sy;
											m_lhs_z(R_s_id, R_m_id) += sign_factor * penalty_factor_tangent_continuity * gp_i_weight * J1 * omega_T2_mz * omega_T2_sz;						

											k_s++;								
										}
									}
									k_m++;
								}
							}
						}					
			// measuring quality of the results
				double measure_g0_continuity(boost::python::list& edges_with_enforced_tangent_continuity)
					{
						double integral_d = 0;
						double length = 0;
						for (BREPElementVector::iterator brep_elem_i = m_brep_elements.begin(); brep_elem_i != m_brep_elements.end(); ++brep_elem_i)
						{
							// Check if brep_elem_i is a element used for coupling or for Dirichlet boundary conditions
							if(brep_elem_i->HasCouplingCondition())
							{
								// Get Gauss points of current brep element
									BREPGaussPointVector brep_gps = brep_elem_i->GetGaussPoints();
								
								// Loop over all Gauss points of current brep element 
								for (BREPGaussPointVector::iterator brep_gp_i = brep_gps.begin(); brep_gp_i != brep_gps.end(); ++brep_gp_i)
								{
									// Read information from Gauss point
										unsigned int master_patch_id = brep_gp_i->GetPatchId();
										unsigned int slave_patch_id = brep_gp_i->GetSlavePatchId();
										Patch::Pointer master_patch = m_patches[m_patch_position_in_patch_vector[master_patch_id]];
										Patch::Pointer slave_patch = m_patches[m_patch_position_in_patch_vector[slave_patch_id]];
										double gp_i_weight = brep_gp_i->GetWeight();
										Vector location_on_master_patch = brep_gp_i->GetLocation();
										Vector location_on_slave_patch = brep_gp_i->GetSlaveLocation();
										Vector tangent_on_master_patch = brep_gp_i->GetTangent();
										Vector tangent_on_slave_patch = brep_gp_i->GetSlaveTangent();

									// Evaluate NURBS basis function for Gauss point on both patches and get the corresponding ids of control points in the mapping matrix
										matrix<double> R_gpi_master;
										double u_m = location_on_master_patch(0);
										double v_m = location_on_master_patch(1);
										master_patch->GetSurface().EvaluateNURBSFunctions(-1,-1,u_m, v_m, R_gpi_master);
										ControlPointVector cps_master = master_patch->GetSurface().GetControlPoints();
										matrix<unsigned int> cp_ids_gpi_master = master_patch->GetSurface().
																				GetRelevantControlPointsIndexes(-1,-1,u_m, v_m);

										matrix<double> R_gpi_slave;
										double u_s = location_on_slave_patch(0);
										double v_s = location_on_slave_patch(1);
										slave_patch->GetSurface().EvaluateNURBSFunctions(-1,-1,u_s, v_s, R_gpi_slave);	
										ControlPointVector cps_slave = slave_patch->GetSurface().GetControlPoints();
										matrix<unsigned int> cp_ids_gpi_slave = slave_patch->GetSurface().
																				GetRelevantControlPointsIndexes(-1,-1,u_s, v_s);							

									// Compute Jacobian J1
										matrix<double> g_master = master_patch->GetSurface().GetBaseVectors(-1,-1,u_m,v_m);
										Vector g1 = ZeroVector(3);
										g1(0) = g_master(0,0);
										g1(1) = g_master(1,0);
										g1(2) = g_master(2,0);
										Vector g2 = ZeroVector(3);
										g2(0) = g_master(0,1);
										g2(1) = g_master(1,1);
										g2(2) = g_master(2,1);
										double J1 = norm_2( g1* tangent_on_master_patch(0) + g2* tangent_on_master_patch(1) );

									// Compute distance between master and slave Gauss point
										Vector master_point = ZeroVector(3);
										Vector slave_point = ZeroVector(3);
										for(unsigned int i=0; i<cp_ids_gpi_master.size1(); i++)
										{
											for(unsigned int j=0; j<cp_ids_gpi_master.size2(); j++)
											{
												master_point[0] += R_gpi_master(i,j) * cps_master[cp_ids_gpi_master(i,j)].getX();
												master_point[1] += R_gpi_master(i,j) * cps_master[cp_ids_gpi_master(i,j)].getY();
												master_point[2] += R_gpi_master(i,j) * cps_master[cp_ids_gpi_master(i,j)].getZ();
											}
										}
										for(unsigned int i=0; i<cp_ids_gpi_slave.size1(); i++)
										{
											for(unsigned int j=0; j<cp_ids_gpi_slave.size2(); j++)
											{
												slave_point[0] += R_gpi_slave(i,j) * cps_slave[cp_ids_gpi_slave(i,j)].getX();
												slave_point[1] += R_gpi_slave(i,j) * cps_slave[cp_ids_gpi_slave(i,j)].getY();
												slave_point[2] += R_gpi_slave(i,j) * cps_slave[cp_ids_gpi_slave(i,j)].getZ();
											}
										}
										double distance = norm_2(master_point - slave_point);
									// Sum contribution
										integral_d += distance * J1 * gp_i_weight;
										length += J1 * gp_i_weight;
								}
							}
						}
						return integral_d/length;
					}
				double measure_g1_continuity(boost::python::list& edges_with_enforced_tangent_continuity)
					{
						double integral_angle = 0;
						double length = 0;
						for (BREPElementVector::iterator brep_elem_i = m_brep_elements.begin(); brep_elem_i != m_brep_elements.end(); ++brep_elem_i)
						{
							// Check if brep_elem_i is a element used for coupling or for Dirichlet boundary conditions
							if(brep_elem_i->HasCouplingCondition())
							{
								// Get Gauss points of current brep element
									BREPGaussPointVector brep_gps = brep_elem_i->GetGaussPoints();

								// Check if for current element some continuity is to be enforced
									bool tangent_continuity_to_be_enforced = false;
									for (unsigned int i = 0; i < boost::python::len(edges_with_enforced_tangent_continuity); ++i)
									{
										unsigned int listed_edge_id = extractUnsignedInt(edges_with_enforced_tangent_continuity[i][0]);
										if(brep_elem_i->GetEdgeId() == listed_edge_id)
										{
											tangent_continuity_to_be_enforced = true;
										}
									}
								if(!tangent_continuity_to_be_enforced)
								{
									// Loop over all Gauss points of current brep element 
									for (BREPGaussPointVector::iterator brep_gp_i = brep_gps.begin(); brep_gp_i != brep_gps.end(); ++brep_gp_i)
									{
										// Read information from Gauss point
											// master u, v
											unsigned int master_patch_id = brep_gp_i->GetPatchId();
											Patch::Pointer master_patch = m_patches[m_patch_position_in_patch_vector[master_patch_id]];
											Vector location_on_master_patch = brep_gp_i->GetLocation();
											double u_m = location_on_master_patch(0);
											double v_m = location_on_master_patch(1);
											// slave u, v
											unsigned int slave_patch_id = brep_gp_i->GetSlavePatchId();
											Patch::Pointer slave_patch = m_patches[m_patch_position_in_patch_vector[slave_patch_id]];
											Vector location_on_slave_patch = brep_gp_i->GetSlaveLocation();
											double u_s = location_on_slave_patch(0);
											double v_s = location_on_slave_patch(1);
											// tangent on master
											Vector tangent_on_master_patch = brep_gp_i->GetTangent();
											// gp weight
											double gp_i_weight = brep_gp_i->GetWeight();

										// Compute base vectors
											matrix<double> g_master = master_patch->GetSurface().GetBaseVectors(-1,-1,u_m,v_m);
											Vector g1 = ZeroVector(3);
											g1(0) = g_master(0,0);
											g1(1) = g_master(1,0);
											g1(2) = g_master(2,0);
											Vector g2 = ZeroVector(3);
											g2(0) = g_master(0,1);
											g2(1) = g_master(1,1);
											g2(2) = g_master(2,1);
											Vector g3 = ZeroVector(3);
											g3(0) = g_master(0,2);
											g3(1) = g_master(1,2);
											g3(2) = g_master(2,2);
											matrix<double> g_slave = slave_patch->GetSurface().GetBaseVectors(-1,-1,u_s,v_s);
											Vector g3_s = ZeroVector(3);
											g3_s(0) = g_slave(0,2);
											g3_s(1) = g_slave(1,2);
											g3_s(2) = g_slave(2,2);
											
										// Compute Jacobian J1
											double J1 = norm_2( g1* tangent_on_master_patch(0) + g2* tangent_on_master_patch(1) );

										// Compute angle between normals on master and slave Gauss point
											double cosine = inner_prod(g3, g3_s);
											double angle_rad = acos(std::abs(cosine)); // abs() => ( 0 < angle_rad < PI/2 )
											double angle_deg = angle_rad * 180 / 3.1415926535;

										// Sum contribution
											integral_angle += angle_deg * J1 * gp_i_weight;
											length += J1 * gp_i_weight;
									}
								}
							}
						}
						return integral_angle/length;
					}

		////////////////////////////////////////////////////////


		///////////////////////// MODULAR CODE ///////////////////////////////
		/*
			parametrisation + reconstruction:
			parametrisation
				input	- FE-nodes: (X0, Y0, Z0, X, Y, Z)
						- user-defined nodes: (X0, Y0, Z0, X, Y, Z, (patch, u_approx, v_approx))
						- data_points: (patch, u, v, X, Y, Z)
				output	- data_points: (patch, u, v, X, Y, Z)

				only FE-nodes need to be parametrised
				user-defined nodes -> create data point -> optimise_parametrisation
			
			reconstruction
				input	- data_points

				setup	- activate correct CPs
				setup	- without BCs - BCs penalty - BCs penalty 2-step - BCs augmented Lagrange - BCs augmented Lagrange 2-step

				output	- updated CPs position
		*/
		void print_mapping_ids(bool old_version=false)
		{
			if(old_version)
			{
				std::cout << "write cps of version 4" << std::endl;
				std::ofstream file_to_write("cps_4.txt");
				for(PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
				{
					file_to_write << "patch: " << (*patch_i)->GetId() << std::endl;
					for(ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
					{
						if(cp_i->IsRelevantForMapping())
						{
							file_to_write << "\tcp " << cp_i->getGlobalId() << ": " << cp_i->GetMappingMatrixId() << std::endl;
						}
						else
						{
							file_to_write << "\tcp " << cp_i->getGlobalId() << ":" << std::endl;
						}
					}
				}				
				file_to_write.close();
			}
			else
			{
				std::cout << "write cps of version 5" << std::endl;				
				std::ofstream file_to_write("cps_5.txt");
				for(PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
				{
					file_to_write << "patch: " << (*patch_i)->GetId() << std::endl;
					for(ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
					{
						if(cp_i->IsActive())
						{
							file_to_write << "\tcp " << cp_i->getGlobalId() << ": " << cp_i->GetMappingMatrixId() << std::endl;
						}
						else
						{
							file_to_write << "\tcp " << cp_i->getGlobalId() << ":" << std::endl;
						}
					}
				}
				file_to_write.close();
			}
		}
		void print_lhs(bool old_version=false)
		{
			if(old_version)
			{
				{				std::cout << "write lhs of version 4" << std::endl;
								std::ofstream file_to_write("lhs_4x.txt");
								for(unsigned int i=0; i<m_lhs_x.size1(); i++)
								{
									for(unsigned int j=0; j<m_lhs_x.size2(); j++)
									{
										file_to_write << m_lhs_x(i,j) << " ";
									}
									file_to_write << std::endl;
								}		
								file_to_write.close();
				}
				{				std::cout << "write lhs of version 4" << std::endl;
								std::ofstream file_to_write("lhs_4y.txt");
								for(unsigned int i=0; i<m_lhs_x.size1(); i++)
								{
									for(unsigned int j=0; j<m_lhs_x.size2(); j++)
									{
										file_to_write << m_lhs_y(i,j) << " ";
									}
									file_to_write << std::endl;
								}		
								file_to_write.close();
				}
				{				std::cout << "write lhs of version 4" << std::endl;
								std::ofstream file_to_write("lhs_4z.txt");
								for(unsigned int i=0; i<m_lhs_x.size1(); i++)
								{
									for(unsigned int j=0; j<m_lhs_x.size2(); j++)
									{
										file_to_write << m_lhs_z(i,j) << " ";
									}
									file_to_write << std::endl;
								}		
								file_to_write.close();
				}
				{				std::cout << "write rhs of version 4" << std::endl;
								std::ofstream file_to_write("rhs_4x.txt");
								for(unsigned int i=0; i<m_rhs_x.size(); i++)
								{
									file_to_write << m_rhs_x(i) << " ";
								}		
								file_to_write.close();
				}				
				{				std::cout << "write rhs of version 4" << std::endl;
								std::ofstream file_to_write("rhs_4y.txt");
								for(unsigned int i=0; i<m_rhs_y.size(); i++)
								{
									file_to_write << m_rhs_y(i) << " ";
								}		
								file_to_write.close();
				}				
				{				std::cout << "write rhs of version 4" << std::endl;
								std::ofstream file_to_write("rhs_4z.txt");
								for(unsigned int i=0; i<m_rhs_z.size(); i++)
								{
									file_to_write << m_rhs_z(i) << " ";
								}		
								file_to_write.close();
				}				
			}
			else
			{
				{				std::cout << "write lhs of version 5" << std::endl;				
								std::ofstream file_to_write("lhs_5x.txt");
								for(unsigned int i=0; i<m_lhs_x.size1(); i++)
								{
									for(unsigned int j=0; j<m_lhs_x.size2(); j++)
									{
										file_to_write << m_lhs_x(i,j) << " ";
									}
									file_to_write << std::endl;
								}	
								file_to_write.close();
				}
				{				std::cout << "write lhs of version 5" << std::endl;				
								std::ofstream file_to_write("lhs_5y.txt");
								for(unsigned int i=0; i<m_lhs_x.size1(); i++)
								{
									for(unsigned int j=0; j<m_lhs_x.size2(); j++)
									{
										file_to_write << m_lhs_y(i,j) << " ";
									}
									file_to_write << std::endl;
								}	
								file_to_write.close();
				}
				{				std::cout << "write lhs of version 5" << std::endl;				
								std::ofstream file_to_write("lhs_5z.txt");
								for(unsigned int i=0; i<m_lhs_x.size1(); i++)
								{
									for(unsigned int j=0; j<m_lhs_x.size2(); j++)
									{
										file_to_write << m_lhs_z(i,j) << " ";
									}
									file_to_write << std::endl;
								}	
								file_to_write.close();
				}
				{				std::cout << "write rhs of version 5" << std::endl;
								std::ofstream file_to_write("rhs_5x.txt");
								for(unsigned int i=0; i<m_rhs_x.size(); i++)
								{
									file_to_write << m_rhs_x(i) << " ";
								}		
								file_to_write.close();
				}				
				{				std::cout << "write rhs of version 5" << std::endl;
								std::ofstream file_to_write("rhs_5y.txt");
								for(unsigned int i=0; i<m_rhs_y.size(); i++)
								{
									file_to_write << m_rhs_y(i) << " ";
								}		
								file_to_write.close();
				}				
				{				std::cout << "write rhs of version 5" << std::endl;
								std::ofstream file_to_write("rhs_5z.txt");
								for(unsigned int i=0; i<m_rhs_z.size(); i++)
								{
									file_to_write << m_rhs_z(i) << " ";
								}		
								file_to_write.close();
				}
			}
		}
			// define data points
				void use_all_FE_nodes_as_data_points()
				{
					for(ModelPart::NodesContainerType::iterator node_i = mr_fe_model_part.NodesBegin(); node_i!=mr_fe_model_part.NodesEnd(); node_i++)
					{
						
						NodeType::Pointer node_ptr = mr_fe_model_part.pGetNode( node_i->GetId() );
						DataPoint data_point(node_ptr);
						m_data_points.push_back(data_point);
					}

				}
			// parametrisation
				void parametrisation(const unsigned int u_resolution, const unsigned int v_resolution)
				{
					std::cout << "\n\n\n ****************** PARAMETRISATION ******************" << std::endl;
					boost::timer overall;
		
					NodeVector list_of_cad_nodes;
					create_point_cloud(u_resolution, v_resolution, list_of_cad_nodes);

					// create KD-Tree
					std::cout << "\n\t> Constructing search-tree..." << std::endl;
					boost::timer timer;
					int bucket_size = 20;
					tree nodes_tree(list_of_cad_nodes.begin(), list_of_cad_nodes.end(), bucket_size);
					std::cout << "\t\t\t DONE (" << timer.elapsed() << " s)" << std::endl;

					compute_nearest_neighbours(nodes_tree);
					compute_nearest_points_Newton_Raphson();

					std::cout << "\n ****************** (" << overall.elapsed() << " s) ********************" << std::endl;
				}

				// auxiliary functions
					void create_point_cloud(const unsigned int u_resolution, const unsigned int v_resolution,
											NodeVector& list_of_cad_nodes)
					{
						//
							// input:
							//		u_resolution, v_resolution
							// output:
							//		list_of_cad_nodes, list_of_us_of_cad_nodes, list_of_vs_of_cad_nodes, list_of_patch_itrs_of_cad_nodes
							
							// for each patch
							// 1:	Create a grid of points, equidistant in the parameter space.
							//		All the points belong to the inner of the patch (i.e. no points lying on the border)
							// for each point
							// 	if the point is inside
							// 2:	its coordinates in the geometrical space are computed
							// 3:	all the information is stored in the output lists:	- u,
							//															- v,
							//															- patch_id,
							//															- CAD node (containing X, Y, Z coordinates)
							
							// N.B. CAD nodes ids start with 1

						std::cout << "\n\t> Creating point cloud..." << std::endl;
						boost::timer timer;		

						unsigned int cad_node_counter = 0;

						//Loop over all surface of all patches
						for (unsigned int patch_itr = 0; patch_itr < m_patches.size(); patch_itr++)
						{
							Patch::Pointer patch_ptr = m_patches[patch_itr];
							NURBSSurface& surface = patch_ptr->GetSurface();

							// Get relevant data
								double u_min = surface.GetKnotVectorU().front();
								double u_max = surface.GetKnotVectorU().back();
								double v_min = surface.GetKnotVectorV().front();
								double v_max = surface.GetKnotVectorV().back();
								double delta_u = (u_max-u_min) / u_resolution;
								double delta_v = (v_max-v_min) / v_resolution;

							// Loop over all u & v according to specified resolution
							for(unsigned int i=1; i<u_resolution; i++)
							{
								double u_i = u_min + i*delta_u;

								for(unsigned int j=1; j<v_resolution; j++)
								{
									double v_j = v_min + j*delta_v;

									// Check if u_i and v_j represent a point inside the closed boundary loop
									if(patch_ptr->CheckIfPointIsInside(u_i, v_j))
									{
										++cad_node_counter;					
										// compute unique point in CAD-model for given u&v
										Point<3> cad_point;
										surface.EvaluateSurfacePoint(cad_point, u_i, v_j);

										// Add id to point --> node. Add node to list of CAD nodes
										NodeType::Pointer cad_node_ptr = Node<3>::Pointer(new Node<3>(cad_node_counter, cad_point));

										// // Store for cad node the corresponding cad information in separate vectors
										list_of_cad_nodes.push_back(cad_node_ptr);
										m_point_cloud[cad_node_ptr] = {u_i, v_j, patch_ptr};
									}
								}
							}
						}						
			
						std::cout << "\t\t\t DONE (" << timer.elapsed() << " s)" << std::endl;
					}

					void compute_nearest_neighbours(tree& nodes_tree)
					{

						std::cout << "\n\t> Computing nearest neighbours..." << std::endl;
						boost::timer timer;

						// loop over all data points (knowing XYZ in undeformed configuration, find an approximation of patch, u, v of corresponding CAD point)
						for(DataPointsList::iterator data_point_i = m_data_points.begin(); data_point_i != m_data_points.end(); data_point_i++)
						{
							NodeType::Pointer data_point_node_ptr = data_point_i->getNodePtr();
							// Search nearest cad neighbor of current data point
							NodeType::Pointer neighbour = nodes_tree.SearchNearestPoint( *data_point_node_ptr );

							// Store CAD information of neighbour
							data_point_i->setPatch(m_point_cloud[neighbour].patch_ptr);
							// data_point_i->setU(m_point_cloud[neighbour].u);
							// data_point_i->setV(m_point_cloud[neighbour].v);
							data_point_i->updateUAndV(m_point_cloud[neighbour].u, m_point_cloud[neighbour].v);

							// only for debugging
								m_list_of_neighbour_points.push_back(neighbour);
								m_list_of_u_of_neighbour_points.push_back(m_point_cloud[neighbour].u);
								m_list_of_v_of_neighbour_points.push_back(m_point_cloud[neighbour].v);
								// m_list_of_span_u_of_neighbour_points.push_back(span_u_of_np);
								// m_list_of_span_v_of_neighbour_points.push_back(span_v_of_np);				
								// m_list_of_patch_of_neighbour_points.push_back(patch_itr_of_nearest_point);
						}
						
						std::cout << "\t\t\t DONE (" << timer.elapsed() << " s)" << std::endl;
					}

					void compute_nearest_points_Newton_Raphson()
					{
						std::cout << "\n\t> Optimising coordinates with Newton-Raphson..." << std::endl;
						boost::timer timer;

						// Specify a tolerance and a maximum number of iteration for the Newton-Raphson optimizer
						double tol = 1e-5;
						unsigned int max_itr = 50;
						// Loop over data points and find for each the closest cad point
						for(DataPointsList::iterator data_point_i = m_data_points.begin(); data_point_i != m_data_points.end(); data_point_i++)
						{
							data_point_i->optimize_parametrisation(tol, max_itr);
						}

						std::cout << "\t\t\t DONE (" << timer.elapsed() << " s)" << std::endl;
					}
				//
			// reconstruction
				// flag control points
					void activate_control_point() // one cp is active, cp's neighbours are relevant
						{
							// cp active
							// cp relevant
							// neighbourhood relevant
						}

					void activate_patch(Patch::Pointer patch) // all relevant cps are active
						{
							double u, v;
							// loop over data points
							for(DataPointsList::iterator data_point_i = m_data_points.begin(); data_point_i != m_data_points.end(); data_point_i++)
							{
								u = data_point_i->getU();
								v = data_point_i->getV();
								if( data_point_i->getPatch()->GetId() == patch->GetId() &&
									patch->CheckIfPointIsInside(u, v)
								  )
								{
									// flag cps as relevant and active
									data_point_i->flagControlPointsAsRelevantAndActive();
								}
							}	
						}

					void activate_knotspan() // some cps are relevant, inner cps are active (if any)
						{}



					void deactivate_all()
						{
							for(PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
							{
								for(ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
								{
									cp_i->Reset();
								}	
							}
							m_n_active_control_points = 0;
							m_n_relevant_control_points = 0;
						}
					
					void activate_borders()
						{
							std::cout << "WARNING: in cad_mapper.h, void activate_borders() to be cleaned." << std::endl;
							// Then we identify mapping relevant control points required from the specified boundary conditions
							// Accordingly we check all Gauss points of all brep elements for their control points
							for (BREPElementVector::iterator brep_elem_i = m_brep_elements.begin(); brep_elem_i != m_brep_elements.end(); ++brep_elem_i)
							{
								// Get Gauss points of current brep element
								BREPGaussPointVector brep_gps = brep_elem_i->GetGaussPoints();

								// Loop over all Gauss points of current brep element 
								for (BREPGaussPointVector::iterator brep_gp_i = brep_gps.begin(); brep_gp_i != brep_gps.end(); ++brep_gp_i)
								{
									// Flag control points on master patch
									unsigned int master_patch_id = brep_gp_i->GetPatchId();
									Vector location = brep_gp_i->GetLocation();
									m_patches[m_patch_position_in_patch_vector[master_patch_id]]->GetSurface().FlagControlPointsAsRelevantAndActive(-1, -1, location[0], location[1]);

									// Flag control points on slave patch if brep element is a coupling element
									if(brep_elem_i->HasCouplingCondition())
									{
										unsigned int slave_patch_id = brep_gp_i->GetSlavePatchId();
										location = brep_gp_i->GetSlaveLocation();
										m_patches[m_patch_position_in_patch_vector[slave_patch_id]]->GetSurface().FlagControlPointsAsRelevantAndActive(-1, -1, location[0], location[1]);
									}
								}
							}
						}
					
					void activate_all()
						{


							for(PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
							{
								Patch::Pointer patch = (*patch_i);
								activate_patch(patch);
							}
						}

				// assign ids
					void assign_mapping_matrix_ids()
						{
							m_n_active_control_points = 0;
							m_n_relevant_control_points = 0;

							for(PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
							{
								for(ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
								{
									if(cp_i->IsActive())
									{
										cp_i->SetMappingMatrixId(m_n_active_control_points);
										m_n_active_control_points++;
									}
									if(cp_i->IsRelevantForMapping())
									{
										m_n_relevant_control_points++;
									}
								}	
							}
							
						}
				// setup system
					void initialize_lse() // parameters: type of norm, if one cp at a time
						{
							// Initialize lhs matrices
								m_lhs_x.resize(m_n_active_control_points, m_n_active_control_points);
								m_lhs_y.resize(m_n_active_control_points, m_n_active_control_points);
								m_lhs_z.resize(m_n_active_control_points, m_n_active_control_points);
								m_lhs_x.clear();
								m_lhs_y.clear();
								m_lhs_z.clear();

							// Initialize rhs vectors
								m_rhs_x.resize(m_n_active_control_points);	
								m_rhs_y.resize(m_n_active_control_points);	
								m_rhs_z.resize(m_n_active_control_points);	
								m_rhs_x.clear();
								m_rhs_y.clear();
								m_rhs_z.clear();
						}
				// 
					void least_square_minimization_of_data_points()
						{
							// Loop over data points
							for(DataPointsList::iterator data_point_i = m_data_points.begin(); data_point_i != m_data_points.end(); data_point_i++)
							{
								// HOW TO DIRECTLY LOOP OVER DATA POINTS WHICH HAVE AT LEAST AN ACTIVE CP???
								//
								// get relevant control points and correspoinding NURBS function values
								matrix<double> R_data_point = data_point_i->EvaluateNURBSFunctions();
								matrix<unsigned int> cp_ids = data_point_i->getRelavantControlPointsIds();
								ControlPointVector CP = data_point_i->getControlPoints();

								for(unsigned int i=0; i<cp_ids.size1(); i++)
								{
									for(unsigned int j=0; j<cp_ids.size2(); j++)
									{
										ControlPoint cp_ij = CP[cp_ids(i,j)];

										// modify the LSE iff at least one CP is active  
										if(cp_ij.IsActive())
										{
											double R_ij = R_data_point(i,j);
											unsigned int mapping_matrix_id_ij = cp_ij.GetMappingMatrixId();

											for(unsigned int k=0; k<cp_ids.size1(); k++)
											{
												for(unsigned int l=0; l<cp_ids.size2(); l++)
												{
													ControlPoint cp_kl = CP[cp_ids(k,l)];
													double R_kl = R_data_point(k,l);
													if(cp_kl.IsActive())
													{
														unsigned int mapping_matrix_id_kl = cp_kl.GetMappingMatrixId();
														// contribution to LHS
														m_lhs_x(mapping_matrix_id_ij, mapping_matrix_id_kl) += R_ij*R_kl;
														m_lhs_y(mapping_matrix_id_ij, mapping_matrix_id_kl) += R_ij*R_kl;
														m_lhs_z(mapping_matrix_id_ij, mapping_matrix_id_kl) += R_ij*R_kl;
													}
													else
													{
														// dirichlet BC contribution
														m_rhs_x(mapping_matrix_id_ij) -= R_ij*R_kl*cp_kl.getX();
														m_rhs_y(mapping_matrix_id_ij) -= R_ij*R_kl*cp_kl.getY();
														m_rhs_z(mapping_matrix_id_ij) -= R_ij*R_kl*cp_kl.getZ();
													}
												}
											}

											// contribution to RHS
											m_rhs_x(mapping_matrix_id_ij) += R_ij * data_point_i->getX();
											m_rhs_y(mapping_matrix_id_ij) += R_ij * data_point_i->getY();
											m_rhs_z(mapping_matrix_id_ij) += R_ij * data_point_i->getZ();
										}
									}
								}
							}
						}
				// solve
					void solve_update_and_test()
						{
							// solve for updated position of CP: x_CP, y_CP, z_CP
								CompressedMatrixType compressed_lhs_x = m_lhs_x; //???
								CompressedMatrixType compressed_lhs_y = m_lhs_y; //???
								CompressedMatrixType compressed_lhs_z = m_lhs_z; //???

								std::cout << "\n\t> Solving x... \n";	
								Vector x_CP = ZeroVector(m_n_active_control_points);
								m_linear_solver->Solve(compressed_lhs_x, x_CP, m_rhs_x);		
								std::cout << "\t\t\t DONE" << std::endl;

								std::cout << "\t> Solving y... \n";	
								Vector y_CP = ZeroVector(m_n_active_control_points);
								m_linear_solver->Solve(compressed_lhs_y, y_CP, m_rhs_y);		
								std::cout << "\t\t\t DONE" << std::endl;

								std::cout << "\t> Solving z... \n";	
								Vector z_CP = ZeroVector(m_n_active_control_points);
								m_linear_solver->Solve(compressed_lhs_z, z_CP, m_rhs_z);		
								std::cout << "\t\t\t DONE" << std::endl;

							// update
								std::cout << "\n\t> Updating control points positions... " << std::endl;
								m_cad_reader.UpdateControlPointsPositions(m_patches, x_CP, y_CP, z_CP);
								std::cout << "\t\t\t DONE" << std::endl;

							// test computing residuals
								Vector rhs_test = ZeroVector(m_n_active_control_points);
								noalias(rhs_test) = prod(compressed_lhs_x,x_CP);
								Vector rhs_difference = m_rhs_x - rhs_test;
								double normalized_difference_in_rhs = norm_2(rhs_difference);
								std::cout << "\n\t> Residual when solving for x (Euclidean norm): " << normalized_difference_in_rhs << std::endl;
					
								rhs_test = ZeroVector(m_n_active_control_points);
								noalias(rhs_test) = prod(compressed_lhs_y,y_CP);
								rhs_difference = m_rhs_y - rhs_test;
								normalized_difference_in_rhs = norm_2(rhs_difference);
								std::cout << "\t> Residual when solving for y (Euclidean norm): " << normalized_difference_in_rhs << std::endl;
					
								rhs_test = ZeroVector(m_n_active_control_points);
								noalias(rhs_test) = prod(compressed_lhs_z,z_CP);
								rhs_difference = m_rhs_z - rhs_test;
								normalized_difference_in_rhs = norm_2(rhs_difference);
								std::cout << "\t> Residual when solving for z (Euclidean norm): " << normalized_difference_in_rhs << std::endl;	
							
						}
				// regularization techniques
					void alpha_criterion()
						{
							double u_x, v_x;
							// loop over patches
							for (unsigned int patch_cp = 0; patch_cp < m_patches.size(); patch_cp++)
							{
								NURBSSurface surface = m_patches[patch_cp]->GetSurface();
								ControlPointVector CP = surface.GetControlPoints();

								// loop over relevant control points
								for(unsigned int cp_index=0; cp_index<CP.size(); cp_index++)
								{
									ControlPoint cp = CP[cp_index];
									if(cp.IsActive())
									{
										unsigned int x = cp.GetMappingMatrixId();
										// compute coordinates of the point on surface associated to cp
										surface.SetGrevilleAbscissae(cp_index, u_x, v_x);

										matrix<double> R; surface.EvaluateNURBSFunctions( -1, -1, u_x, v_x, R);						
										matrix<int> control_points_ids = surface.GetRelevantControlPointsIndexes( -1, -1, u_x, v_x);
										// LHS contribution 1
										m_lhs_x(x, x) += m_alpha;
										m_lhs_y(x, x) += m_alpha;
										m_lhs_z(x, x) += m_alpha;

										for(unsigned int k=0; k<control_points_ids.size1(); k++)
										{
											for(unsigned int l=0; l<control_points_ids.size2(); l++)
											{
												ControlPoint cp_j = CP[control_points_ids(k,l)];
												double R_j = R(k,l);

												if(cp_j.IsActive())
												{
													unsigned int j = cp_j.GetMappingMatrixId();
													// LHS contribution 2
													m_lhs_x(x, j) -= m_alpha * R_j;
													m_lhs_y(x, j) -= m_alpha * R_j;
													m_lhs_z(x, j) -= m_alpha * R_j;
													// LHS contribution 3
													m_lhs_x(j, x) -= m_alpha * R_j;
													m_lhs_y(j, x) -= m_alpha * R_j;
													m_lhs_z(j, x) -= m_alpha * R_j;

													for(unsigned int m=0; m<control_points_ids.size1(); m++)
													{
														for(unsigned int n=0; n<control_points_ids.size2(); n++)
														{
															ControlPoint cp_i = CP[control_points_ids(m,n)];
															double R_i = R(m,n);
															if(cp_i.IsActive())
															{
																unsigned int i = cp_i.GetMappingMatrixId();
																// LHS contribution 4
																m_lhs_x(j, i) += m_alpha * R_j*R_i;
																m_lhs_y(j, i) += m_alpha * R_j*R_i;
																m_lhs_z(j, i) += m_alpha * R_j*R_i;
															}
															else
															{
																// RHS contribution 2
																m_rhs_x(j) -= m_alpha * R_j*R_i*cp_i.getX();
																m_rhs_y(j) -= m_alpha * R_j*R_i*cp_i.getY();
																m_rhs_z(j) -= m_alpha * R_j*R_i*cp_i.getZ();
															}
														}
													}
												}

												else
												{
													// RHS contribution 1									
													m_rhs_x(x) += m_alpha * R_j*cp_j.getX();
													m_rhs_y(x) += m_alpha * R_j*cp_j.getY();
													m_rhs_z(x) += m_alpha * R_j*cp_j.getZ();
												}
											}
										}
									}
								}
							}
						}
							
					void alpha_criterion_NR()
						{
							double u_x, v_x;
							// loop over patches
							for (unsigned int patch_cp = 0; patch_cp < m_patches.size(); patch_cp++)
							{
								NURBSSurface surface = m_patches[patch_cp]->GetSurface();
								ControlPointVector CP = surface.GetControlPoints();

								// loop over relevant control points
								for(unsigned int cp_index=0; cp_index<CP.size(); cp_index++)
								{
									ControlPoint cp = CP[cp_index];
									if(cp.IsActive())
									{
										unsigned int x = cp.GetMappingMatrixId();
										// compute coordinates of the point on surface associated to cp
										surface.SetGrevilleAbscissae(cp_index, u_x, v_x);
										std::cout << "\n\n\nold coords:" << std::endl;
										KRATOS_WATCH(u_x); KRATOS_WATCH(v_x);
										DataPoint data_point(cp.getX(), cp.getY(), cp.getZ());
										data_point.setPatch(m_patches[patch_cp]);
										data_point.updateUAndV(u_x, v_x);
										data_point.optimize_parametrisation(1e-5, 50);
										u_x = data_point.getU();
										v_x = data_point.getV();
										std::cout << "new coords:" << std::endl;
										KRATOS_WATCH(u_x); KRATOS_WATCH(v_x);

										matrix<double> R; surface.EvaluateNURBSFunctions( -1, -1, u_x, v_x, R);						
										matrix<int> control_points_ids = surface.GetRelevantControlPointsIndexes( -1, -1, u_x, v_x);
										// LHS contribution 1
										m_lhs_x(x, x) += m_alpha;
										m_lhs_y(x, x) += m_alpha;
										m_lhs_z(x, x) += m_alpha;

										for(unsigned int k=0; k<control_points_ids.size1(); k++)
										{
											for(unsigned int l=0; l<control_points_ids.size2(); l++)
											{
												ControlPoint cp_j = CP[control_points_ids(k,l)];
												double R_j = R(k,l);

												if(cp_j.IsActive())
												{
													unsigned int j = cp_j.GetMappingMatrixId();
													// LHS contribution 2
													m_lhs_x(x, j) -= m_alpha * R_j;
													m_lhs_y(x, j) -= m_alpha * R_j;
													m_lhs_z(x, j) -= m_alpha * R_j;
													// LHS contribution 3
													m_lhs_x(j, x) -= m_alpha * R_j;
													m_lhs_y(j, x) -= m_alpha * R_j;
													m_lhs_z(j, x) -= m_alpha * R_j;

													for(unsigned int m=0; m<control_points_ids.size1(); m++)
													{
														for(unsigned int n=0; n<control_points_ids.size2(); n++)
														{
															ControlPoint cp_i = CP[control_points_ids(m,n)];
															double R_i = R(m,n);
															if(cp_i.IsActive())
															{
																unsigned int i = cp_i.GetMappingMatrixId();
																// LHS contribution 4
																m_lhs_x(j, i) += m_alpha * R_j*R_i;
																m_lhs_y(j, i) += m_alpha * R_j*R_i;
																m_lhs_z(j, i) += m_alpha * R_j*R_i;
															}
															else
															{
																// RHS contribution 2
																m_rhs_x(j) -= m_alpha * R_j*R_i*cp_i.getX();
																m_rhs_y(j) -= m_alpha * R_j*R_i*cp_i.getY();
																m_rhs_z(j) -= m_alpha * R_j*R_i*cp_i.getZ();
															}
														}
													}
												}

												else
												{
													// RHS contribution 1									
													m_rhs_x(x) += m_alpha * R_j*cp_j.getX();
													m_rhs_y(x) += m_alpha * R_j*cp_j.getY();
													m_rhs_z(x) += m_alpha * R_j*cp_j.getZ();
												}
											}
										}
									}
								}
							}
						}

					void beta_criterion()
						{
							// add contribution to lhs
								for(unsigned int i=0; i<m_n_active_control_points; i++)
								{
									m_lhs_x(i,i) += m_beta;
									m_lhs_y(i,i) += m_beta;
									m_lhs_z(i,i) += m_beta;
								}
							// add contribution to rhs
								int cp_id;
								for (PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
								{
									for (ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
									{
										if(cp_i->IsActive())
										{
											cp_id = cp_i->GetMappingMatrixId();
											m_rhs_x[cp_id] += m_beta * cp_i->getX();
											m_rhs_y[cp_id] += m_beta * cp_i->getY();
											m_rhs_z[cp_id] += m_beta * cp_i->getZ();
										}
									}
								}
						}
				// boundary conditions
					void apply_boundary_conditions_5(double penalty_factor_disp, 
													double penalty_factor_rot
													//double penalty_factor_dirichlet, 
													//boost::python::list& edges_with_specific_dirichlet_conditions, 
													// boost::python::list& edges_with_enforced_tangent_continuity
													)
						{
							// Loop over all brep elements specifying boundary conditions 
							for (BREPElementVector::iterator brep_elem_i = m_brep_elements.begin(); brep_elem_i != m_brep_elements.end(); ++brep_elem_i)
							{
								// Check if brep_elem_i is a element used for coupling or for Dirichlet boundary conditions
								if(brep_elem_i->HasCouplingCondition())
									apply_coupling_condition_small(brep_elem_i, penalty_factor_disp, penalty_factor_rot, m_edges_with_enforced_tangent_continuity);
								// else if(brep_elem_i->HasDirichletCondition())
									// apply_dirichlet_condition(brep_elem_i, penalty_factor_dirichlet, edges_with_specific_dirichlet_conditions);
							}
												}
				// functions exposed to user (different strategies are implemented here)
					void user_exposed_function_general_structure()
						{
							// activate cps
							// assign ids
							// initialize lse
							// add contributions
								// 1.
								// 2.
								// ...
							// solve and update
							// deactivate all
						}
					void map_all_patches()
						{
							std::cout << "\n\n\n ****************** RECONSTRUCTION *******************" << std::endl;
							boost::timer overall;
							// activate cps
														std::cout << "\n\t> Activating control points..." << std::endl;
														boost::timer timer1;
							activate_all();
							if(m_penalty_factor_disp > 0 || m_penalty_factor_rot > 0) activate_borders();

														std::cout << "\t\t\t DONE (" << timer1.elapsed() << " s)" << std::endl;

							// assign ids
														std::cout << "\n\t> Assigning mapping matrix ids..." << std::endl;
														boost::timer timer2;
							assign_mapping_matrix_ids();
														std::cout << "\t\t\t DONE (" << timer2.elapsed() << " s)" << std::endl;

							// initialize lse
														std::cout << "\n\t> Setting up the linear system of equations..." << std::endl;
														boost::timer timer3;
							initialize_lse();
							// add contributions
								least_square_minimization_of_data_points();
								if(m_alpha > 0) alpha_criterion_NR();
								if(m_beta > 0)  beta_criterion();
								// if boundaries
								if(m_penalty_factor_disp > 0 || m_penalty_factor_rot > 0) apply_boundary_conditions_5(m_penalty_factor_disp, m_penalty_factor_rot);

														std::cout << "\t\t\t DONE (" << timer3.elapsed() << " s)" << std::endl;
							// solve and update
							solve_update_and_test();
							// deactivate all
							// deactivate_all();
							std::cout << "\n ****************** (" << overall.elapsed() << " s) ********************" << std::endl;
						}
					void map_patch_by_patch()
						{
							std::cout << "\n\n\n ****************** RECONSTRUCTION *******************" << std::endl;
							boost::timer overall;
							for(PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
							{
								std::cout << "\n ****************** PATCH " << (*patch_i)->GetId() << " **************************" << std::endl;
								// activate cps
														std::cout << "\n\t> Activating control points..." << std::endl;
														boost::timer timer1;
								activate_patch(*patch_i);
														std::cout << "\t\t\t DONE (" << timer1.elapsed() << " s)" << std::endl;
								// assign ids
														std::cout << "\n\t> Assigning mapping matrix ids..." << std::endl;
														boost::timer timer2;
								assign_mapping_matrix_ids();
														std::cout << "\t\t\t DONE (" << timer2.elapsed() << " s)" << std::endl;
								// initialize lse
														std::cout << "\n\t> Setting up the linear system of equations..." << std::endl;
														boost::timer timer3;
								initialize_lse();
								// add contributions
									least_square_minimization_of_data_points();
									if(m_beta > 0)  beta_criterion();
									// if boundaries
														std::cout << "\t\t\t DONE (" << timer3.elapsed() << " s)" << std::endl;
								// solve and update
								solve_update_and_test();
								// deactivate all
								deactivate_all();
								
								output_surface_points(std::to_string((*patch_i)->GetId()) + ".txt", 150, 150, -1); ////////////////////////////////////////////////////
							}
							std::cout << "\n ****************** (" << overall.elapsed() << " s) ********************" << std::endl;
						}
					void map_CP_by_CP()
						{

						}

					void map_boundary_conditions()
						{
							// activate cps
							activate_borders();
							// assign ids
							assign_mapping_matrix_ids();
							// initialize lse
							initialize_lse();
							// add contributions
								least_square_minimization_of_data_points();
								if(m_beta > 0)  beta_criterion();
								if(m_penalty_factor_disp > 0 || m_penalty_factor_rot > 0) apply_boundary_conditions_5(m_penalty_factor_disp, m_penalty_factor_rot);
							// solve and update
							solve_update_and_test();
							// deactivate all
							deactivate_all();
							activate_all();
							activate_borders();
						}
					void crazy_step_back()
						{
							// activate cps
							activate_borders();
							// assign ids
							assign_mapping_matrix_ids();
							// initialize lse

							// add contributions

							// solve and update
								// set CP back to original position
									Vector x_CP = ZeroVector(m_n_active_control_points);
									Vector y_CP = ZeroVector(m_n_active_control_points);
									Vector z_CP = ZeroVector(m_n_active_control_points);
									for (PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
									{
										for (ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
										{
											if(cp_i->IsActive())
											{
												// Updating c++ data base
												unsigned int cp_mapping_matrix_id = cp_i->GetMappingMatrixId();
												x_CP[cp_mapping_matrix_id] = cp_i->getX0();
												y_CP[cp_mapping_matrix_id] = cp_i->getY0();
												z_CP[cp_mapping_matrix_id] = cp_i->getZ0();
											}
										}
									}
									m_cad_reader.UpdateControlPointsPositions(m_patches, x_CP, y_CP, z_CP);
							// deactivate all
							deactivate_all();
						}
					void map_boundary_conditions_augmented_Lagrange(boost::python::list& edges_with_enforced_tangent_continuity,
																	double r1,
																	double r2)
						{
							// activate cps
							activate_borders();
							// assign ids
							assign_mapping_matrix_ids();
							double lambda_1 = 0;
							double lambda_2 = 0;
							double up_1 = 1e5;
							double up_2 = 1e5;
							while(up_1 > 1e-2 || up_2 > 1e-2)
							{
								undeform();
								//
									// initialize lse
									initialize_lse();
									// add contributions
										least_square_minimization_of_data_points();
										if(m_beta > 0)  beta_criterion(); 					// depends on active CP position !!!
										apply_boundary_conditions_5(lambda_1, lambda_2);	// depends on active CP position !!!
									// solve and update
									solve_update_and_test();
								//
								up_1 = measure_g0_continuity(edges_with_enforced_tangent_continuity);
								up_2 = measure_g1_continuity(edges_with_enforced_tangent_continuity);
								lambda_1 = lambda_1 + r1 * up_1;
								lambda_2 = lambda_2 + r2 * up_2;
								KRATOS_WATCH(up_1);
								KRATOS_WATCH(up_2);
								KRATOS_WATCH(lambda_1);
								KRATOS_WATCH(lambda_2);
								double a;
								std::cin >> a;
								if(a==2) break;
							}
							// deactivate all							
							deactivate_all();
						}
					void map_all_patches_augmented_Lagrange(boost::python::list& edges_with_enforced_tangent_continuity,
																	double r1,
																	double r2)
						{
							std::cout << "\n\n\n ****************** RECONSTRUCTION *******************" << std::endl;
							boost::timer overall;
							// activate cps
														std::cout << "\n\t> Activating control points..." << std::endl;
														boost::timer timer1;
							activate_all();
							activate_borders();

														std::cout << "\t\t\t DONE (" << timer1.elapsed() << " s)" << std::endl;

							// assign ids
														std::cout << "\n\t> Assigning mapping matrix ids..." << std::endl;
														boost::timer timer2;
							assign_mapping_matrix_ids();
														std::cout << "\t\t\t DONE (" << timer2.elapsed() << " s)" << std::endl;

							double lambda_1 = 0;
							double lambda_2 = 0;
							double up_1 = 1e5;
							double up_2 = 1e5;
							while(up_1 > 1e-2 || up_2 > 1e-2)
							{
								undeform();
									// initialize lse
																std::cout << "\n\t> Setting up the linear system of equations..." << std::endl;
																boost::timer timer3;
									initialize_lse();
									// add contributions
										least_square_minimization_of_data_points();
										if(m_beta > 0)  beta_criterion();
										// if boundaries
										apply_boundary_conditions_5(lambda_1, lambda_2);

																std::cout << "\t\t\t DONE (" << timer3.elapsed() << " s)" << std::endl;
									// solve and update
									solve_update_and_test();
								up_1 = measure_g0_continuity(edges_with_enforced_tangent_continuity);
								up_2 = measure_g1_continuity(edges_with_enforced_tangent_continuity);
								lambda_1 = lambda_1 + r1 * up_1;
								lambda_2 = lambda_2 + r2 * up_2;
								KRATOS_WATCH(up_1);
								KRATOS_WATCH(up_2);
								KRATOS_WATCH(lambda_1);
								KRATOS_WATCH(lambda_2);
								double a;
								std::cin >> a;
								if(a==2) break;
							}
							// deactivate all							
							deactivate_all();
							std::cout << "\n ****************** (" << overall.elapsed() << " s) ********************" << std::endl;
						}
												
				// setup
					void apply_regularization_schemes(double alpha, double beta, double delta)
						{
							m_alpha = alpha;
							m_beta = beta;
							m_delta = delta;
						}
					void apply_penalty_factors(
											//    double penalty_factor_dirichlet,
											   double penalty_factor_disp,
											   double penalty_factor_rot,
											//    double penalty_factor_tangent_continuity,
											   boost::python::list edges_with_enforced_tangent_continuity
											   )
						{
							m_penalty_factor_disp = penalty_factor_disp;
							m_penalty_factor_rot = penalty_factor_rot;
							m_edges_with_enforced_tangent_continuity = edges_with_enforced_tangent_continuity;
							KRATOS_WATCH("apply enforcement");
						}
				// to clean later
					void undeform()
						{
							for(PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
							{
								for(ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
								{
									cp_i->setdX(0);
									cp_i->setdY(0);
									cp_i->setdZ(0);
								}	
							}
						}
		////////////////////////////////////////////////////////
		
		/////////////// ONE CONTROL POINT AT A TIME: exploitation of NURBS locality ///////////////
			// cp: control point to move
			// Q: data points affected
			// P: control points affecting Q

			// // get_surrounding_cps()
			// given cp, compute P
			// cp -> i_u_cp, j_v_cp
			// cp_(i_u_cp - p)(j_v_cp - q), ..., cp(i_u_cp + p)(j_v_cp + q) belong to p

			// given cp, compute Q
			// u_min = u_i of i_u_cp; u_max = u_i of i_u_cp + p + 1 (or u_i+1 of i_u_cp + p);
			// v_min = v_i of i_v_cp; v_max = v_i of i_v_cp + q + 1 (or v_i+1 of i_v_cp + q);
			// loop over all data points 
			// 	if u_min < u_i < u_max and v_min < v_i < v_max
			// 		point belongs to Q

			// given cp, compute R
			// given cp, compute R*
			// given cp, compute R
			// given cp, compute R
		////////////////////////////////////////////////////////
		
		// EXTERNAL: separating FE-mesh data from computation //
			void external_compute_lhs()
			{
				std::cout << "\n\t\t> Starting direct computation of LHS..." << std::endl;
				boost::timer lhs_timer;
				// 1: Identify relevant control points affecting the cad points on the surface and assign them an id
					// Flag relevant control points		
					for(int i = 0; i < m_number_of_points_external; i++)
					{
						double u = m_list_of_u_external[i];
						double v = m_list_of_v_external[i];
						unsigned int patch_id = m_list_of_patch_ids_external[i];
						m_patches[patch_id]->GetSurface().FlagControlPointsForMapping(-1, 
																					-1,
																					u,
																					v);
					}
					std::cout << "\n> Control points flagged for mapping." << std::endl;
					// Count them and assign each a unique mapping matrix Id	
					m_n_control_points = 0; 
					m_n_relevant_control_points = 0;
					unsigned int mapping_matrix_id = 0;
					for (PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
					{
						for (ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
						{
							if(cp_i->IsRelevantForMapping())
							{
								cp_i->SetMappingMatrixId(mapping_matrix_id);
								++m_n_relevant_control_points;
								++mapping_matrix_id;
							}
							++m_n_control_points;
						}
					}
					std::cout << "\n> Number of control points in total = " << m_n_control_points << "." << std::endl;
					std::cout << "\n> Number of control points relevant for mapping = " << m_n_relevant_control_points << ".\n" << std::endl;

				// 2: Compute lhs matrix
					// Initialize lhs matrix
					double number_of_CAD_dofs = 3 * m_n_relevant_control_points;
					m_lhs_external.resize(number_of_CAD_dofs, number_of_CAD_dofs);
					m_lhs_external.clear();
					// loop over points to map and add their contribution
					for(int i = 0; i < m_number_of_points_external; i++)
					{
						double u = m_list_of_u_external[i];
						double v = m_list_of_v_external[i];
						unsigned int patch_id = m_list_of_patch_ids_external[i];

						// Get CAD-shape-function-value for all control points affecting the nearest cad point
						matrix<double> R_CAD_Pi;
						m_patches[patch_id]->GetSurface().EvaluateNURBSFunctions( -1,
																				-1,
																				u, 
																				v,
																				R_CAD_Pi );
																							
						// Get the corresponding id of control points in the matrix
						matrix<unsigned int> mapping_matrix_ids_cad = m_patches[patch_id]->GetSurface().GetMappingMatrixIds( -1,
																															-1,
																															u,
																															v);

						// Assemble a matrix
						for(unsigned int i=0; i<mapping_matrix_ids_cad.size1(); i++)
						{
							for(unsigned int j=0; j<mapping_matrix_ids_cad.size2(); j++)
							{
								for(unsigned int k=0; k<mapping_matrix_ids_cad.size1(); k++)
								{
									for(unsigned int l=0; l<mapping_matrix_ids_cad.size2(); l++)
									{
										unsigned int row_id = mapping_matrix_ids_cad(i,j);
										unsigned int col_id = mapping_matrix_ids_cad(k,l);
										double R_ij = R_CAD_Pi(i,j);
										double R_kl = R_CAD_Pi(k,l);
										for(unsigned int m=0; m<3; m++)
										{
											// add contribution of 
											m_lhs_external(3*row_id+m, 3*col_id+m) += R_ij*R_kl;
										}
									}
								}
							}
						}
					}
				std::cout << "\t\t> Direct computation of LHS finished in " << lhs_timer.elapsed() << " s." << std::endl;
				KRATOS_WATCH(m_lhs_external);			
			}

			void external_compute_rhs()
			{
				std::cout << "\n\t\t> Starting direct computation of RHS..." << std::endl;
				boost::timer rhs_timer;
				// 1: Compute transposed rectangular matrix
					// Initialize a matrix
					double number_of_points_dofs = 3 * m_number_of_points_external;
					double number_of_CAD_dofs = 3 * m_n_relevant_control_points;
					m_a_matrix.resize(number_of_CAD_dofs, number_of_points_dofs);
					m_a_matrix.clear();

					// Compute a matrix
					for(int i = 0; i < m_number_of_points_external; i++)
					{
						// Get the corresponding id of FE-node in the matrix			
						int row_id = i;
						double u = m_list_of_u_external[i];
						double v = m_list_of_v_external[i];
						unsigned int patch_id = m_list_of_patch_ids_external[i];

						// Get CAD-shape-function-value for all control points affecting the nearest cad point
						matrix<double> R_CAD_Pi;
						m_patches[patch_id]->GetSurface().EvaluateNURBSFunctions( -1,
																				-1,
																				u, 
																				v,
																				R_CAD_Pi );
					
						// Get the corresponding id of control points in the matrix
						matrix<unsigned int> mapping_matrix_ids_cad = m_patches[patch_id]->GetSurface().GetMappingMatrixIds( -1,
																															-1,
																															u,
																															v);
						// Assemble a matrix
						for(unsigned int i=0; i<mapping_matrix_ids_cad.size1();i++)
						{
							for(unsigned int j=0; j<mapping_matrix_ids_cad.size2();j++)
							{
								unsigned int col_id = mapping_matrix_ids_cad(i,j);
								double R_i = R_CAD_Pi(i,j);
							
								for(unsigned int m=0; m<3; m++)
								{
									m_a_matrix(3*col_id + m, 3*row_id + m) = R_i;
								}
							}
						}
					}
				// 2: Compute points coordinates vector
					Vector x_points = ZeroVector(3*m_number_of_points_external);
					for(int i = 0; i < m_number_of_points_external; i++)
					{
						x_points[3*i] = m_list_of_updated_x_external[i];
						x_points[3*i+1] = m_list_of_updated_y_external[i];
						x_points[3*i+2] = m_list_of_updated_z_external[i];
					}
				// 3: Compute matrix vector multiplication
					m_rhs_external.clear(); 
					m_rhs_external = prod(m_a_matrix,x_points);
				
				std::cout << "\t\t> Direct computation of RHS finished in " << rhs_timer.elapsed() << " s." << std::endl;
				KRATOS_WATCH(m_rhs_external);					
			}
			
			void external_map_to_cad_space()
			{
				std::cout << "\n> Mapping to CAD space..." << std::endl;			
				boost::timer map_timer;
				// 1: Solve system
					external_compute_lhs();
					external_compute_rhs();
					// CompressedMatrixType compressed_lhs = m_lhs_matrix; //???
					Vector x_CP = ZeroVector(3*m_n_relevant_control_points);
					m_linear_solver->Solve(m_lhs_external, x_CP, m_rhs_external);		
				// 2: Update control points position //??? CHANGE FUNCTION UpdateControlPoints
					Vector x_CP_old = ZeroVector(3*m_n_relevant_control_points);
					for (PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
					{
						for (ControlPointVector::iterator cp_i = (*patch_i)->GetSurface().GetControlPoints().begin(); cp_i != (*patch_i)->GetSurface().GetControlPoints().end(); ++cp_i)
						{
							if(cp_i->IsRelevantForMapping())
							{
								int cp_id = cp_i->GetMappingMatrixId(); //??? line 1982:	cp_i->SetMappingMatrixId(mapping_matrix_id);
								x_CP_old[3*cp_id] = cp_i->getX();
								x_CP_old[3*cp_id+1] = cp_i->getY();
								x_CP_old[3*cp_id+2] = cp_i->getZ();
							}
						}
					}
					Vector ds = x_CP - x_CP_old; // CP position update
					m_cad_reader.UpdateControlPoints(m_patches, ds);
				std::cout << "\n> Mapping finished in " << map_timer.elapsed() << " s." << std::endl;			
			}

			void set_point(unsigned int patch_id, const double u, const double v,
						const double updated_x, const double updated_y, const double updated_z)
			{
					m_number_of_points_external++;
					m_list_of_patch_ids_external.push_back(patch_id);
					m_list_of_u_external.push_back(u);
					m_list_of_v_external.push_back(v);
					m_list_of_updated_x_external.push_back(updated_x);
					m_list_of_updated_y_external.push_back(updated_y);
					m_list_of_updated_z_external.push_back(updated_z);
			}
			double compute_objective()
			{
				unsigned int patch_id;
				double u, v, x, y, z, dx, dy, dz;
				double obj = 0;
				for(int i = 0; i < m_number_of_points_external; i++)
				{
					patch_id = m_list_of_patch_ids_external[i];
					u = m_list_of_u_external[i];
					v = m_list_of_v_external[i];
					Point<3> cad_point_coordinates;
					m_patches[patch_id]->GetSurface().EvaluateSurfacePoint(cad_point_coordinates, u, v);

					x = m_list_of_updated_x_external[i];
					y = m_list_of_updated_y_external[i];
					z = m_list_of_updated_z_external[i];
					dx = x - cad_point_coordinates[0];
					dy = y - cad_point_coordinates[1];
					dz = z - cad_point_coordinates[2];
					obj += pow(dx,2) + pow(dy,2) + pow(dz,2);
				}
				return obj;
			}

		////////////////////////////////////////////////////////

		// old functions for testing
			void compare_lhs() // compare lhs obtained by map_to_cad_space_3 and map_to_cad_space_4
			{
				std::cout << "\nComparing lhs: ##############################################################" << std::endl;
				// compute lhs with both methods
					compute_lhs_matrix();
					compute_lhs_small();

				// compute difference in lhs
					Matrix lhs_difference;
					lhs_difference.resize(m_n_relevant_control_points, m_n_relevant_control_points);
					lhs_difference.clear();

					double max=0;
					double frob=0;
					for(unsigned int i = 0; i < m_lhs.size1(); i++)
					{
						for(unsigned int j=0; j<m_lhs.size2(); j++)
						{
							lhs_difference(i,j) = m_lhs(i,j) - m_lhs_matrix(3*i, 3*j);
							if( std::abs(lhs_difference(i,j)) > max) max = std::abs(lhs_difference(i,j)); // std::abs!!!!!
							frob += lhs_difference(i,j)*lhs_difference(i,j);

						}
					}
					frob = sqrt(frob);

				// quantify difference in lhs: norm, ...
					std::cout << "\t\tinfinity norm is " << max << std::endl;
					std::cout << "\t\tFrobenius norm is " << frob << std::endl;
			}

			void compare_rhs()
			{
				std::cout << "\nComparing rhs: ##############################################################" << std::endl;
				// compute rhs with both methods
					compute_a_matrix();
					Vector rhs;
					compute_rhs_vector(rhs);
					compute_rhs_small();

				// compute difference in rhs
					Vector rhs_difference;
					rhs_difference.resize(3*m_n_relevant_control_points);
					rhs_difference.clear();
					for(unsigned int i = 0; i<m_n_relevant_control_points; i++)
					{
						rhs_difference[3*i] = rhs[3*i] - m_rhs_x[i];
						rhs_difference[3*i+1] = rhs[3*i+1] - m_rhs_y[i];
						rhs_difference[3*i+2] = rhs[3*i+2] - m_rhs_z[i];
					}
					
				// quantify difference in rhs
					std::cout << "\t\tinfinity norm is " << norm_inf(rhs_difference) << std::endl;
					std::cout << "\t\teuclidean norm is " << norm_2(rhs_difference) << std::endl;
			}
		////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////

	double compute_real_length(unsigned int patch_itr, const double u_a, const double v_a, const double u_b, const double v_b)
	{
		/*
		Given are two points in the parameter space: A = (u_a, v_a); B = (u_b, v_b).
		The points on the line segment from A to B define a 3D-curve, when mapped in the real space.
		This function returns the length of the mapped curve, assuming an Euclidean norm.
		
		To ease the calculation the curve is defined as:
			gamma: [0, 1] -> R^3,
						s |-> (x,y,z), 
		through the linear mapping functions:
			u: [0, 1] -> KnotVectorU,
					s |-> u = u_a + S * (u_b-u_a)
			v: [0, 1] -> KnotVectorV,
					s |-> v = v_a + S * (v_b-v_a).
		Thus the length is:
						   /
			L_2(gamma) =  | 	|| d(gamma)/ds ||_2 ds  ("integral over [0,1] of the Euclidean norm of the first derivative of gamma")
						 / [0,1]

		N.B.: the first derivative of gamma is the tangent vector in the real space.
		*/
		std::cout << "\n> Starting computation of real length..." << std::endl;
		
		double u_s, v_s, s, R_m_du, R_m_dv, R_m_ds,weight_s; 
		double delta_u = u_b - u_a;
		double delta_v = v_b - v_a;
		matrix<double> R_s;
		std::vector<Matrix> dR_s;
		matrix<int> CP_s;
		double length = 0;
		Vector gamma_prime_s;
		double norm_2_gamma_prime_s = 0;

		// compute gauss points
		NURBSSurface surface = m_patches[patch_itr]->GetSurface();
		int p = surface.GetDegP();
		int q = surface.GetDegQ();
		int deg = p*q-1; //2->3 //3->8
		int number_of_Gauss_points = (deg+1) / 2; //2->2 //3->5
		KRATOS_WATCH(number_of_Gauss_points);
		std::cout << "> Warning: Gauss quadrature coordinates and weights currently hardcoded!!!" << std::endl;

		///////HARD CODING GAUSS POINTS COORDINATES AND WEIGHTS////////////////////////////////////////////////////
			DoubleVector x, w;
			/* n = 2 */
					if(number_of_Gauss_points == 2)
			{
			x.push_back(-0.5773502691896257645091488);
			x.push_back(0.5773502691896257645091488);
			w.push_back(1.0000000000000000000000000);
			w.push_back(1.0000000000000000000000000);
			}
			// /* n = 3 */
			// 		if(number_of_Gauss_points == 3)
			// {
			// double x[3] = {-0.7745966692414833770358531,0.0000000000000000000000000,0.7745966692414833770358531};
			// double w[3] = {0.5555555555555555555555556,0.8888888888888888888888889,0.5555555555555555555555556};
			// }
			// /* n = 4 */
			// 		if(number_of_Gauss_points == 4)
			// {
			// double x[4] = {-0.8611363115940525752239465,-0.3399810435848562648026658,0.3399810435848562648026658,0.8611363115940525752239465};
			// double w[4] = {0.3478548451374538573730639,0.6521451548625461426269361,0.6521451548625461426269361,0.3478548451374538573730639};
			// }

			/* n = 5 */
					if(number_of_Gauss_points == 5)
			{
			x.push_back(-0.9061798459386639927976269);
			x.push_back(-0.5384693101056830910363144);
			x.push_back(0.0000000000000000000000000);
			x.push_back(0.5384693101056830910363144);
			x.push_back(0.9061798459386639927976269);
			w.push_back(0.2369268850561890875142640);
			w.push_back(0.4786286704993664680412915);
			w.push_back(0.5688888888888888888888889);
			w.push_back(0.4786286704993664680412915);
			w.push_back(0.2369268850561890875142640);
			}
			// /* n = 6 */
			// 		if(number_of_Gauss_points == 6)
			// {
			// double x[6] = {-0.9324695142031520278123016,-0.6612093864662645136613996,-0.2386191860831969086305017,0.2386191860831969086305017,0.6612093864662645136613996,0.9324695142031520278123016};
			// double w[6] = {0.1713244923791703450402961,0.3607615730481386075698335,0.4679139345726910473898703,0.4679139345726910473898703,0.3607615730481386075698335,0.1713244923791703450402961};
			// }
			// /* n = 7 */
			// 		if(number_of_Gauss_points == 7)
			// {
			// double x[7] = {-0.9491079123427585245261897,-0.7415311855993944398638648,-0.4058451513773971669066064,0.0000000000000000000000000,0.4058451513773971669066064,0.7415311855993944398638648,0.9491079123427585245261897};
			// double w[7] = {0.1294849661688696932706114,0.2797053914892766679014678,0.3818300505051189449503698,0.4179591836734693877551020,0.3818300505051189449503698,0.2797053914892766679014678,0.1294849661688696932706114};
			// }
			// /* n = 8 */
			// 		if(number_of_Gauss_points == 8)
			// {
			// double x[8] = {-0.9602898564975362316835609,-0.7966664774136267395915539,-0.5255324099163289858177390,-0.1834346424956498049394761,0.1834346424956498049394761,0.5255324099163289858177390,0.7966664774136267395915539,0.9602898564975362316835609};
			// double w[8] = {0.1012285362903762591525314,0.2223810344533744705443560,0.3137066458778872873379622,0.3626837833783619829651504,0.3626837833783619829651504,0.3137066458778872873379622,0.2223810344533744705443560,0.1012285362903762591525314	};
			// }
			// /* n = 9 */
			// 		if(number_of_Gauss_points == 9)
			// {
			// double x[9] = {-0.9681602395076260898355762,-0.8360311073266357942994298,-0.6133714327005903973087020,-0.3242534234038089290385380,0.0000000000000000000000000,0.3242534234038089290385380,0.6133714327005903973087020,0.8360311073266357942994298,0.9681602395076260898355762};
			// double w[9] = {0.0812743883615744119718922,0.1806481606948574040584720,0.2606106964029354623187429,0.3123470770400028400686304,0.3302393550012597631645251,0.3123470770400028400686304,0.2606106964029354623187429,0.1806481606948574040584720,0.0812743883615744119718922};
			// }
			// /* n = 10 */
			// 		if(number_of_Gauss_points == 10)
			// {
			// double x[10] = {-0.9739065285171717200779640,-0.8650633666889845107320967,-0.6794095682990244062343274,-0.4333953941292471907992659,-0.1488743389816312108848260,0.1488743389816312108848260,0.4333953941292471907992659,0.6794095682990244062343274,0.8650633666889845107320967,0.9739065285171717200779640};
			// double w[10] = {0.0666713443086881375935688,0.1494513491505805931457763,0.2190863625159820439955349,0.2692667193099963550912269,0.2955242247147528701738930,0.2955242247147528701738930,0.2692667193099963550912269,0.2190863625159820439955349,0.1494513491505805931457763,0.0666713443086881375935688};
			// }
			// /* n = 11 */
			// 		if(number_of_Gauss_points == 11)
			// {
			// double x[11] = {-0.9782286581460569928039380,-0.8870625997680952990751578,-0.7301520055740493240934163,-0.5190961292068118159257257,-0.2695431559523449723315320,0.0000000000000000000000000,0.2695431559523449723315320,0.5190961292068118159257257,0.7301520055740493240934163,0.8870625997680952990751578,0.9782286581460569928039380};
			// double w[11] = {0.0556685671161736664827537,0.1255803694649046246346943,0.1862902109277342514260976,0.2331937645919904799185237,0.2628045445102466621806889,0.2729250867779006307144835,0.2628045445102466621806889,0.2331937645919904799185237,0.1862902109277342514260976,0.1255803694649046246346943,0.0556685671161736664827537};
			// }
			// /* n = 12 */
			// 		if(number_of_Gauss_points == 12)
			// {
			// double x[12] = {-0.9815606342467192506905491,-0.9041172563704748566784659,-0.7699026741943046870368938,-0.5873179542866174472967024,-0.3678314989981801937526915,-0.1252334085114689154724414,0.1252334085114689154724414,0.3678314989981801937526915,0.5873179542866174472967024,0.7699026741943046870368938,0.9041172563704748566784659,0.9815606342467192506905491};
			// double w[12] = {0.0471753363865118271946160,0.1069393259953184309602547,0.1600783285433462263346525,0.2031674267230659217490645,0.2334925365383548087608499,0.2491470458134027850005624,0.2491470458134027850005624,0.2334925365383548087608499,0.2031674267230659217490645,0.1600783285433462263346525,0.1069393259953184309602547,0.0471753363865118271946160};
			// }
		///////////////////////////////////////////////////////////////////////////////////////////////////////////

		// for s in GP: //???
		if(number_of_Gauss_points !=2 and number_of_Gauss_points != 5) return -1;
		for(int s_it = 0; s_it < number_of_Gauss_points; s_it++)
		{
			s = (x[s_it]+1)/2; // mapping: x goes belongs to [-1, 1], s belongs to [0, 1]
			weight_s = w[s_it];

			// 1: compute derivative of gamma with respect to s at s (tangent vector)
			gamma_prime_s = ZeroVector(3);
			u_s = u_a + s * delta_u;
			v_s = v_a + s * delta_v;
			surface.EvaluateNURBSFunctionsAndDerivative( -1,
													     -1,
													     u_s,
													     v_s,
													     R_s,
													     dR_s);
			CP_s = surface.GetRelevantControlPointsIndexes( -1,
													 		-1,
													 		u_s,
													 		v_s);
			// loop over relevant control points: CP_m
			for(int i=0; i < q+1; i++)
			{
				for(int j=0; j < p+1; j++)
				{
					ControlPointVector CP  = surface.GetControlPoints();

					// compute derivative of NURBS functions with respect to s
					R_m_du = dR_s[0](i,j);
					R_m_dv = dR_s[1](i,j);
					R_m_ds = R_m_du * delta_u + R_m_dv * delta_v; // d(Nu*Nv)/ds = d(Nu*Nv)/du * du/ds + d(Nu*Nv)/dv * dv/ds

					// add contribution to derivative of gamma with respect to s
					gamma_prime_s[0] += R_m_ds * CP[CP_s(i,j)].getX(); // 3d vector 
					gamma_prime_s[1] += R_m_ds * CP[CP_s(i,j)].getY(); // 3d vector 
					gamma_prime_s[2] += R_m_ds * CP[CP_s(i,j)].getZ(); // 3d vectoZ 

				}
			}

			// 2: compute euclidean norm of tangent vector
			norm_2_gamma_prime_s = norm_2(gamma_prime_s);
			
			// 3: add contribution to real length of the curve gamma
			length += norm_2_gamma_prime_s * weight_s;
		}

		std::cout << "> Finished computing real length." << std::endl;
		return length/2; // <-- Jacobian is 1/2
	}
    // ==============================================================================

    /// Turn back information as a string.
    virtual std::string Info() const
    {
		return "CADMapper";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
		rOStream << "CADMapper";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }

  private:
    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    ModelPart &mr_fe_model_part;
	CADModelReader m_cad_reader;
    boost::python::dict mr_cad_geometry;
	boost::python::dict mr_cad_integration_data;
    PatchVector m_patches;
	BREPElementVector m_brep_elements;
	unsigned int m_n_control_points;
	unsigned int m_n_relevant_control_points;
	std::map<unsigned int, unsigned int> m_patch_position_in_patch_vector;

    // ==============================================================================
    // General working arrays
    // ==============================================================================
	unsigned int m_n_relevant_fem_points;
    SparseMatrixType m_mapping_matrix_CAD_CAD;
	SparseMatrixType m_mapping_matrix_CAD_FEM;
	Vector m_mapping_rhs_vector;
	Matrix m_a_matrix; // rectangular NURBS matrix
	Matrix m_lhs_matrix;
	// arrays for debugging 1 //???
	NodeVector m_list_of_nearest_points;	
	DoubleVector m_list_of_u_of_nearest_points;
	DoubleVector m_list_of_v_of_nearest_points;
	DoubleVector m_list_of_span_u_of_nearest_points;
	DoubleVector m_list_of_span_v_of_nearest_points;		
	IntVector m_list_of_patch_of_nearest_points;
	// arrays for debugging 2 //???
	NodeVector m_list_of_neighbour_points;	
	DoubleVector m_list_of_u_of_neighbour_points;
	DoubleVector m_list_of_v_of_neighbour_points;
	DoubleVector m_list_of_span_u_of_neighbour_points;
	DoubleVector m_list_of_span_v_of_neighbour_points;		
	IntVector m_list_of_patch_of_neighbour_points;
	// ONE COORDINATE AT A TIME //
	CompressedMatrixType m_lhs;	
	CompressedMatrixType m_lhs_x;
	CompressedMatrixType m_lhs_y;	
	CompressedMatrixType m_lhs_z;	
	Vector m_rhs_x;
	Vector m_rhs_y;
	Vector m_rhs_z;
	// MODULAR //
	unsigned int m_n_active_control_points;
	double m_alpha = 0;
	double m_beta = 0;
	double m_delta = 0;
	double m_penalty_factor_disp;
	double m_penalty_factor_rot;
	boost::python::list m_edges_with_enforced_tangent_continuity;
	// EXTERNAL: separating FE-mesh data from computation //
	double m_number_of_points_external = 0;
	IntVector m_list_of_patch_ids_external;
	DoubleVector m_list_of_u_external;
	DoubleVector m_list_of_v_external;
	DoubleVector m_list_of_updated_x_external;
	DoubleVector m_list_of_updated_y_external;
	DoubleVector m_list_of_updated_z_external;
	CompressedMatrixType m_lhs_external;
	Vector m_rhs_external;
	// clean:
	PointCloud m_point_cloud;
	DataPointsList m_data_points;
	
	const Condition::GeometryType::IntegrationMethod m_integration_method = GeometryData::GI_GAUSS_5;

	// ==============================================================================
    // Solver and strategies
    // ==============================================================================
	CompressedLinearSolverType::Pointer m_linear_solver;

    /// Assignment operator.
    //      CADMapper& operator=(CADMapper const& rOther);

    /// Copy constructor.
    //      CADMapper(CADMapper const& rOther);

}; // Class CADMapper
} // namespace Kratos.

#endif // CAD_MAPPER_H
