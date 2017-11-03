// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "spaces/ublas_space.h"
#include "includes/properties.h"
#include "includes/model_part.h"

/* GAUSS-LEGENDRE */

/* Utilities */
#include "utilities/mortar_utilities.h"
#include "utilities/exact_mortar_segmentation_utility.h"
#include "custom_utilities/derivatives_utilities.h"

namespace Kratos 
{
    namespace Testing 
    {

        typedef Point                                                    PointType;
        typedef Node<3>                                                   NodeType;
        typedef Geometry<NodeType>                                GeometryNodeType;
        typedef Geometry<PointType>                              GeometryPointType;
        typedef Vector                                                  VectorType;
        typedef Matrix                                                  MatrixType;
        
        ///Type definition for integration methods
        typedef GeometryData::IntegrationMethod                   IntegrationMethod;
        typedef IntegrationPoint<2>                            IntegrationPointType;
        typedef GeometryNodeType::IntegrationPointsArrayType integration_pointsType;

        /** 
         * Checks if the derivatives of the shape functions work as expected
         * Case of the Triangle3D3
         */
    
        KRATOS_TEST_CASE_IN_SUITE(TestDerivativesShapeFunctionTriangle, ContactStructuralApplicationFastSuite)
        {
            // Type definitions
            typedef PointBelong<3> PointBelongType;
            typedef array_1d<PointBelongType, 3> ConditionArrayType;
            typedef typename std::vector<ConditionArrayType> ConditionArrayListType;
            typedef Triangle3D3<PointType> TriangleType;
            typedef DerivativesUtilities<3, 3, false> DerivativesUtilitiesType;
    
            ModelPart model_part("Main");
            model_part.SetBufferSize(2);
            
            Properties::Pointer p_cond_prop = model_part.pGetProperties(0);
            
            // Variables addition
            model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            
            PointType aux_point;
            aux_point.Coordinates() = ZeroVector(3);
            
            // First we create the nodes 
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1,-0.2,0.1,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0,0.1,0.0);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 0.2,1.2,0.0);
            
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, 0.6,0.8,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(5, 1.0,0.1,1.0e-3);
            NodeType::Pointer p_node_5 = model_part.CreateNewNode(6,-0.2,0.1,1.0e-3);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            
            Triangle3D3 <Node<3>> triangle_0( condition_nodes_0 );
            const array_1d<double, 3>& normal_0 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 1, triangle_0, p_cond_prop);
            p_node_1->SetValue(NORMAL, normal_0);
            p_node_2->SetValue(NORMAL, normal_0);
            p_node_3->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            
            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            
            Triangle3D3 <Node<3>> triangle_1( condition_nodes_1 );
            const array_1d<double, 3>& normal_1 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 2, triangle_1, p_cond_prop);
            p_node_4->SetValue(NORMAL, normal_1);
            p_node_5->SetValue(NORMAL, normal_1);
            p_node_6->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            
            // Create and initialize condition variables
            MortarKinematicVariablesWithDerivatives<3, 3> rVariables;
    
            // Create the current contact data
            DerivativeData<3, 3> rDerivativeData;
            rDerivativeData.Initialize(triangle_0, model_part.GetProcessInfo());
            
            // We call the exact integration utility
            ExactMortarIntegrationUtility<3, 3, true>  integration_utility = ExactMortarIntegrationUtility<3, 3, true> (2);
            
            unsigned int number_of_iterations = 3;
            Vector error_vector(number_of_iterations, 0.0);
            for (unsigned int iter = 0; iter < number_of_iterations; ++iter)
            {
                // We add displacement to the node 4
                array_1d<double, 3> aux_delta_disp = ZeroVector(3);
                aux_delta_disp[1] = iter * 0.1;
                array_1d<double, 3>& current_disp = p_node_4->FastGetSolutionStepValue(DISPLACEMENT);
                array_1d<double, 3>& previous_disp = p_node_4->FastGetSolutionStepValue(DISPLACEMENT, 1);
                previous_disp = current_disp;
                current_disp = aux_delta_disp;
                // Finally we move the mesh
                noalias(p_node_4->Coordinates())  = p_node_4->GetInitialPosition().Coordinates();
                noalias(p_node_4->Coordinates()) += p_node_4->FastGetSolutionStepValue(DISPLACEMENT);
                
                // Reading integration points
                ConditionArrayListType conditions_points_slave;
                const bool is_inside = integration_utility.GetExactIntegration(triangle_0, normal_0, triangle_1, normal_1, conditions_points_slave);
                
                IntegrationMethod this_integration_method = GeometryData::GI_GAUSS_2;
                
                if (is_inside == true)
                {
                    // Initialize general variables for the current master element
                    rVariables.Initialize();
                    
                    // Update slave element info
                    rDerivativeData.UpdateMasterPair(p_cond_1);
                    
                    for (unsigned int i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom)
                    {
                        std::vector<PointType::Pointer> points_array (3); // The points are stored as local coordinates, we calculate the global coordinates of this points
                        array_1d<PointBelongsTriangle3D3N, 3> belong_array;
                        for (unsigned int i_node = 0; i_node < 3; ++i_node)
                        {
                            PointType global_point;
                            triangle_0.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                            points_array[i_node] = boost::make_shared<PointType>(global_point);
                            belong_array[i_node] = conditions_points_slave[i_geom][i_node].GetBelong();
                        }
                        
                        TriangleType decomp_geom( points_array );
                        
                        const bool bad_shape = MortarUtilities::HeronCheck(decomp_geom);
                        
                        if (bad_shape == false)
                        {
                            const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );
                            
                            // Integrating the mortar operators
                            for ( unsigned int point_number = 0; point_number < integration_points_slave.size(); ++point_number )
                            {
                                const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                                PointType local_point_parent;
                                PointType gp_global;
                                decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                                triangle_0.PointLocalCoordinates(local_point_parent, gp_global);
                                
                                // SLAVE KINEMATIC COMPUTATIONS
                                triangle_0.ShapeFunctionsValues( rVariables.NSlave, local_point_parent.Coordinates() );
                                triangle_0.ShapeFunctionsLocalGradients( rVariables.DNDeSlave, local_point_parent );
                                
                                rVariables.jSlave = decomp_geom.Jacobian( rVariables.jSlave, local_point_decomp.Coordinates());
                                rVariables.DetjSlave = decomp_geom.DeterminantOfJacobian( local_point_decomp );
        
                                // MASTER KINEMATIC COMPUTATIONS
                                PointType projected_gp_global;
                                const array_1d<double,3> gp_normal = MortarUtilities::GaussPointUnitNormal(rVariables.NSlave, triangle_0);
                                
                                GeometryType::CoordinatesArrayType slave_gp_global;
                                triangle_0.GlobalCoordinates( slave_gp_global, local_point_parent );
                                MortarUtilities::FastProjectDirection( triangle_1, slave_gp_global, projected_gp_global, normal_1, -gp_normal ); // The opposite direction
                                
                                GeometryType::CoordinatesArrayType projected_gp_local;
                                
                                triangle_1.PointLocalCoordinates( projected_gp_local, projected_gp_global.Coordinates( ) ) ;

                                triangle_1.ShapeFunctionsValues( rVariables.NMaster,    projected_gp_local );         
                                triangle_1.ShapeFunctionsLocalGradients( rVariables.DNDeMaster, projected_gp_local );
                                rVariables.jMaster = triangle_1.Jacobian( rVariables.jMaster, projected_gp_local);
                        
                                // Now we compute the derivatives
                                const bool consider_normal_variation = false;
                                DerivativesUtilitiesType::CalculateDeltaCellVertex(rVariables, rDerivativeData, belong_array, consider_normal_variation, triangle_0, triangle_1, normal_0);
    
                                // Update the derivatives of the shape functions and the gap
                                DerivativesUtilitiesType::CalculateDeltaN(rVariables, rDerivativeData, triangle_0, triangle_1, normal_0, normal_1, decomp_geom, local_point_decomp, local_point_parent);
                            }
                        }
                    }
                }
                else
                {
                    KRATOS_ERROR << "WRONG, YOU ARE SUPPOSED TO HAVE AN INTERSECTION" << std::endl;
                }
            }
//             const double tolerance = 1.0e-6;
//             for (unsigned int inode = 0; inode < 3; inode++)
//             {
//                 for (unsigned int jnode = 0; jnode < 3; jnode++)
//                 {
//                     KRATOS_CHECK_NEAR(mass_matrix_0(inode,jnode), mass_matrix_1(inode,jnode), tolerance);
//                 }
//             }
        }
        
    } // namespace Testing
}  // namespace Kratos.
