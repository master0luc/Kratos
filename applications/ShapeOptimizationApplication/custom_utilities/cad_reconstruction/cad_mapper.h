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
    typedef CompressedSpaceType::MatrixType CompressedMatrixType;
	typedef SparseSpaceType::MatrixType SparseMatrixType;
    typedef std::vector<double> DoubleVector;
    typedef std::vector<int> IntVector;
    typedef std::vector<ControlPoint> ControlPointVector;
	typedef std::vector<Patch> PatchVector;
	typedef std::vector<BREPElement> BREPElementVector;
	typedef std::vector<BREPGaussPoint> BREPGaussPointVector;
	typedef Node<3> NodeType;
    typedef std::vector<NodeType::Pointer> NodeVector;
    typedef LinearSolver<CompressedSpaceType, DenseSpaceType > CompressedLinearSolverType;
	typedef std::vector<BoundaryEdge> BoundaryEdgeVector;	
	typedef std::vector<BoundaryLoop> BoundaryLoopVector;

	// for tree search
	typedef std::vector<double> DistanceVector;
    typedef std::vector<double>::iterator DistanceIterator;
	typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

	// For handling of python data
	typedef boost::python::extract<double> extractDouble;
    typedef boost::python::extract<bool> extractBool;
	typedef boost::python::extract<unsigned int> extractUnsignedInt;

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
			m_patch_position_in_patch_vector[m_patches[patch_itr].GetId()] = patch_itr;
		
		// Read cad geometry data into the loal c++ containers if available
		if(len(mr_cad_integration_data.keys())>0)
			m_cad_reader.ReadIntegrationData(m_brep_elements);
    }

    /// Destructor.
    virtual ~CADMapper()
    {
    }

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
			unsigned int patch_id =  m_patches[patch_itr].GetId();
			DoubleVector& knot_vec_u_i = m_patches[patch_itr].GetSurface().GetKnotVectorU();
			DoubleVector& knot_vec_v_i = m_patches[patch_itr].GetSurface().GetKnotVectorV();
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
					bool point_is_inside = m_patches[patch_itr].CheckIfPointIsInside(point_of_interest);

					if(point_is_inside)
					{
						// compute unique point in CAD-model for given u&v
						++cad_node_counter;					
						Point<3> cad_point_coordinates;
						m_patches[patch_itr].GetSurface().EvaluateSurfacePoint(cad_point_coordinates, u_i, v_j);

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
        		NodeType resulting_nearest_point;
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
					m_patches[patch_itr_of_nearest_point].GetSurface().EvaluateGradientsForClosestPointSearch(Q_minus_P, myHessian, myGradient , u_k, v_k);

					// u_k and v_k are updated
					MathUtils<double>::InvertMatrix( myHessian, InvH, det_H );
					Vector deltau = prod(InvH,myGradient);
					u_k -= deltau(0);
					v_k -= deltau(1);

					// Q is updated
					m_patches[patch_itr_of_nearest_point].GetSurface().EvaluateSurfacePoint(newtonRaphPoint, u_k, v_k);
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
				IntVector knot_span_nearest_point = m_patches[patch_itr_of_nearest_point].GetSurface().GetKnotSpan(u_of_nearest_point, v_of_nearest_point);

				// Set flag to mark control point as relevant for mapping
				int span_u_of_np = knot_span_nearest_point[0];
				int span_v_of_np = knot_span_nearest_point[1];
				m_patches[patch_itr_of_nearest_point].GetSurface().FlagControlPointsForMapping(span_u_of_np, span_v_of_np, u_of_nearest_point, v_of_nearest_point);

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
				m_patches[m_patch_position_in_patch_vector[master_patch_id]].GetSurface().FlagControlPointsForMapping(-1, -1, location[0], location[1]);

				// Flag control points on slave patch if brep element is a coupling element
				if(brep_elem_i->HasCouplingCondition())
				{
					unsigned int slave_patch_id = brep_gp_i->GetSlavePatchId();
					location = brep_gp_i->GetSlaveLocation();
					m_patches[m_patch_position_in_patch_vector[slave_patch_id]].GetSurface().FlagControlPointsForMapping(-1, -1, location[0], location[1]);
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
			for (ControlPointVector::iterator cp_i = patch_i->GetSurface().GetControlPoints().begin(); cp_i != patch_i->GetSurface().GetControlPoints().end(); ++cp_i)
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
				m_patches[patch_itr_of_nearest_point].GetSurface().EvaluateNURBSFunctions( span_u_of_nearest_point,
																						   span_v_of_nearest_point,
																						   u_of_nearest_point, 
																						   v_of_nearest_point,
																						   R_CAD_Pi );
				
				// Get the corresponding ids of control points in the mapping matrix
				matrix<unsigned int> mapping_matrix_ids_cad = m_patches[patch_itr_of_nearest_point].GetSurface().GetMappingMatrixIds( span_u_of_nearest_point, 
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
			Patch& master_patch = m_patches[m_patch_position_in_patch_vector[master_patch_id]];
			Patch& slave_patch = m_patches[m_patch_position_in_patch_vector[slave_patch_id]];
			double gp_i_weight = brep_gp_i->GetWeight();
			Vector location_on_master_patch = brep_gp_i->GetLocation();
			Vector location_on_slave_patch = brep_gp_i->GetSlaveLocation();
			Vector tangent_on_master_patch = brep_gp_i->GetTangent();
			Vector tangent_on_slave_patch = brep_gp_i->GetSlaveTangent();

			// Evaluate NURBS basis function for Gauss point on both patches and get the corresponding ids of control points in the mapping matrix
			matrix<double> R_gpi_master;
			double u_m = location_on_master_patch(0);
			double v_m = location_on_master_patch(1);
			master_patch.GetSurface().EvaluateNURBSFunctions(-1,-1,u_m, v_m, R_gpi_master);
			matrix<unsigned int> mapping_matrix_ids_gpi_master = master_patch.GetSurface().GetMappingMatrixIds(-1,-1,u_m, v_m);

			matrix<double> R_gpi_slave;
			double u_s = location_on_slave_patch(0);
			double v_s = location_on_slave_patch(1);
			slave_patch.GetSurface().EvaluateNURBSFunctions(-1,-1,u_s, v_s, R_gpi_slave);	
			matrix<unsigned int> mapping_matrix_ids_gpi_slave = slave_patch.GetSurface().GetMappingMatrixIds(-1,-1,u_s, v_s);							

			// Compute Jacobian J1
			matrix<double> g_master = master_patch.GetSurface().GetBaseVectors(-1,-1,u_m,v_m);
			Vector g1 = ZeroVector(3);
			g1(0) = g_master(0,0);
			g1(1) = g_master(1,0);
			g1(2) = g_master(2,0);
			Vector g2 = ZeroVector(3);
			g2(0) = g_master(0,1);
			g2(1) = g_master(1,1);
			g2(2) = g_master(2,1);
			double J1 = norm_2( g1* tangent_on_master_patch(0) + g2* tangent_on_master_patch(1) );

			// if(brep_elem_i->GetEdgeId()==1005)
			// 	m_length_circle += gp_i_weight * J1;

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
	void apply_rotation_coupling( Patch &master_patch,
							      Patch &slave_patch,
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
		master_patch.GetSurface().ComputeVariationOfLocalCSY( u_m, v_m, tangent_on_master_patch, T1_m, T2_m, T3_m, t1r_m, t2r_m, t3r_m );
		slave_patch.GetSurface().ComputeVariationOfLocalCSY( u_s, v_s, tangent_on_slave_patch, T1_s, T2_s, T3_s, t1r_s, t2r_s, t3r_s );

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
	void enforce_tangent_continuity( Patch& master_patch,
								     Patch& slave_patch,
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
		master_patch.GetSurface().ComputeSecondVariationOfLocalCSY( u_m, v_m, 
																	tangent_on_master_patch, 
																	T1_m, T2_m, T3_m, 
																	T1_der_m, T2_der_m, T3_der_m,
																	t1r_m, t2r_m, t3r_m,
																	t1_der_r_m, t2_der_r_m, t3_der_r_m,
																	t1rs_m, t2rs_m, t3rs_m,
																	t1_der_rs_m, t2_der_rs_m, t3_der_rs_m );
		slave_patch.GetSurface().ComputeSecondVariationOfLocalCSY( u_s, v_s, 
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
			Patch& master_patch = m_patches[m_patch_position_in_patch_vector[master_patch_id]];
			double gp_i_weight = brep_gp_i->GetWeight();
			Vector location_on_master_patch = brep_gp_i->GetLocation();
			Vector tangent_on_master_patch = brep_gp_i->GetTangent();

			// Evaluate NURBS basis function for Gauss point on both patches and get the corresponding ids of control points in the mapping matrix
			matrix<double> R_gpi_master;
			double u_m = location_on_master_patch(0);
			double v_m = location_on_master_patch(1);
			master_patch.GetSurface().EvaluateNURBSFunctions(-1,-1,u_m, v_m, R_gpi_master);
			matrix<unsigned int> mapping_matrix_ids_gpi_master = master_patch.GetSurface().GetMappingMatrixIds(-1,-1,u_m, v_m);						

			// Compute Jacobian J1
			matrix<double> g_master = master_patch.GetSurface().GetBaseVectors(-1,-1,u_m,v_m);
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
		std::cout << "\n> Starting writing surface points of given CAD geometry to file..." << std::endl;
	
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
				DoubleVector& knot_vec_u_i = m_patches[patch_itr].GetSurface().GetKnotVectorU();
				DoubleVector& knot_vec_v_i = m_patches[patch_itr].GetSurface().GetKnotVectorV();

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
						bool point_is_inside = m_patches[patch_itr].CheckIfPointIsInside(point_of_interest);

						if(point_is_inside)
						{
							// compute point in CAD-model for given u&v				
							Point<3> cad_point;
							m_patches[patch_itr].GetSurface().EvaluateSurfacePoint(cad_point, u_i, v_j);

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
			DoubleVector& knot_vec_u_i = m_patches[specific_patch].GetSurface().GetKnotVectorU();
			DoubleVector& knot_vec_v_i = m_patches[specific_patch].GetSurface().GetKnotVectorV();

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
					bool point_is_inside = m_patches[specific_patch].CheckIfPointIsInside(point_of_interest);

					if(point_is_inside)
					{
						// compute point in CAD-model for given u&v				
						Point<3> cad_point;
						m_patches[specific_patch].GetSurface().EvaluateSurfacePoint(cad_point, u_i, v_j);

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

		std::cout << "\n> Finished writing surface points of given CAD geometry to file..." << std::endl;
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
				BoundaryLoopVector boundary_loops = m_patches[patch_itr].GetBoundaryLoops();

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
			BoundaryLoopVector boundary_loops = m_patches[specific_patch].GetBoundaryLoops();

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
    void output_control_point_displacements()
    {
		// Outputs the displacements of the control points in a format that may be read by Gid

		std::cout << "\n> Starting to write displacement of control points..." << std::endl;
		std::ofstream output_file("control_point_displacements.post.res");

		output_file << "Rhino Post Results File 1.0" << std::endl;
		output_file << "Result \"Displacement\" \"Load Case\" 0 Vector OnNodes" << std::endl;
		output_file << "Values" << std::endl;

		unsigned int cp_itr = 0;
		for (PatchVector::iterator patch_i = m_patches.begin(); patch_i != m_patches.end(); ++patch_i)
		{
			for (ControlPointVector::iterator cp_i = patch_i->GetSurface().GetControlPoints().begin(); cp_i != patch_i->GetSurface().GetControlPoints().end(); ++cp_i)
			{
				// It is important to iterate outside to stick to the carat settings
				++cp_itr;
				
				if(cp_i->IsRelevantForMapping())
					output_file << cp_itr << " " << cp_i->getdX() << " " << cp_i->getdY() << " " << cp_i->getdZ() << std::endl;
			}
		}

		output_file << "End Values" << std::endl;

		output_file.close();

		std::cout << "\n> Fished writing displacements of control points..." << std::endl;
	}

    // --------------------------------------------------------------------------
    void evaluate_geometric_deviations(const unsigned int u_resolution, const unsigned int v_resolution, boost::python::dict deviations)
    {
		std::cout << "\n> Starting to evaluate L2 norm of deviation vector..." << std::endl;
		boost::timer function_timer;

		// 1st step: Create for each patch coarse cloud of CAD points in x-y space for neighbor search later 
		unsigned int cad_node_counter = 0;
		NodeVector list_of_cad_nodes;
		DoubleVector list_of_us_of_cad_nodes;
		DoubleVector list_of_vs_of_cad_nodes;
		IntVector list_of_patch_itrs_of_cad_nodes;

		//Loop over all surface of all patches / faces
		for (unsigned int patch_itr = 0; patch_itr < m_patches.size(); patch_itr++)
		{
			// Get relevant data
			unsigned int patch_id =  m_patches[patch_itr].GetId();
			DoubleVector& knot_vec_u_i = m_patches[patch_itr].GetSurface().GetKnotVectorU();
			DoubleVector& knot_vec_v_i = m_patches[patch_itr].GetSurface().GetKnotVectorV();
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
					bool point_is_inside = m_patches[patch_itr].CheckIfPointIsInside(point_of_interest);

					if(point_is_inside)
					{
						// compute unique point in CAD-model for given u&v
						++cad_node_counter;					
						Point<3> cad_point_coordinates;
						m_patches[patch_itr].GetSurface().EvaluateSurfacePoint(cad_point_coordinates, u_i, v_j);

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

		// 3rd step: Evaluate nearest CAD nodes for all FEM nodes and add mutual distance to L2 norm
		double l2_norm_of_geometric_deviation = 0.0;
		double max_nodal_deviation = -100000000;
		double average_nodal_deviation = 0.0;
		for(auto& node_i : mr_fe_model_part.Nodes())
		{
			// Search nearest cad neighbor of current integration point
			NodeType resulting_nearest_point;
			NodeType::Pointer nearest_point = nodes_tree.SearchNearestPoint( node_i );

			// Recover CAD information of nearest point
			double u_of_nearest_point = list_of_us_of_cad_nodes[nearest_point->Id()-1];
			double v_of_nearest_point = list_of_vs_of_cad_nodes[nearest_point->Id()-1];
			int patch_itr_of_nearest_point = list_of_patch_itrs_of_cad_nodes[nearest_point->Id()-1];

			// Perform Newton-Raphson for detailed search

			// Initialize P: point on the mesh
			Vector P = ZeroVector(3);
			P(0) = node_i.X();
			P(1) = node_i.Y();
			P(2) = node_i.Z();
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
				m_patches[patch_itr_of_nearest_point].GetSurface().EvaluateGradientsForClosestPointSearch(Q_minus_P, myHessian, myGradient , u_k, v_k);

				// u_k and v_k are updated
				MathUtils<double>::InvertMatrix( myHessian, InvH, det_H );
				Vector deltau = prod(InvH,myGradient);
				u_k -= deltau(0);
				v_k -= deltau(1);

				// Q is updated
				m_patches[patch_itr_of_nearest_point].GetSurface().EvaluateSurfacePoint(newtonRaphPoint, u_k, v_k);
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
			nearest_point->X() = Q_k(0);
			nearest_point->Y() = Q_k(1);
			nearest_point->Z() = Q_k(2);

			// Distance between FEM node and nearest CAD node
			double dist_x = node_i.X() - nearest_point->X();
			double dist_y = node_i.Y() - nearest_point->Y();
			double dist_z = node_i.Z() - nearest_point->Z();

			double squared_nodal_deviation = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
			double nodal_deviation = sqrt(squared_nodal_deviation);

			if(nodal_deviation>max_nodal_deviation)
				max_nodal_deviation = nodal_deviation;
			l2_norm_of_geometric_deviation += squared_nodal_deviation;
			average_nodal_deviation += nodal_deviation;
		}
		std::cout << "\n> Finished evaluating L2 norm of deviation vector in: " << function_timer.elapsed() << " s." << std::endl;

		// Compute geometric quality values and add to output argument
		l2_norm_of_geometric_deviation = sqrt(l2_norm_of_geometric_deviation);
		deviations["l2_norm"] = l2_norm_of_geometric_deviation;
		deviations["nodal_max"] = max_nodal_deviation;
		deviations["nodal_average"] = average_nodal_deviation / mr_fe_model_part.Nodes().size();
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
