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
            NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0,0.0,0.0);
            NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0,0.0,0.0);
            NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 0.0,1.0,0.0);
            
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, 0.0,1.0,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(5, 1.0,0.0,1.0e-3);
            NodeType::Pointer p_node_5 = model_part.CreateNewNode(6, 0.0,0.0,1.0e-3);
//             NodeType::Pointer p_node_1 = model_part.CreateNewNode(1,-0.2,0.1,0.0);
//             NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0,0.1,0.0);
//             NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 0.2,1.2,0.0);
//             
//             NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, 0.6,0.8,1.0e-3);
//             NodeType::Pointer p_node_6 = model_part.CreateNewNode(5, 1.0,0.1,1.0e-3);
//             NodeType::Pointer p_node_5 = model_part.CreateNewNode(6,-0.2,0.1,1.0e-3);
            
            NodeType::Pointer p_node0_1 = model_part.CreateNewNode(7, p_node_1->X(), p_node_1->Y(), p_node_1->Z());
            NodeType::Pointer p_node0_2 = model_part.CreateNewNode(8, p_node_2->X(), p_node_2->Y(), p_node_2->Z());
            NodeType::Pointer p_node0_3 = model_part.CreateNewNode(9, p_node_3->X(), p_node_3->Y(), p_node_3->Z());
            
            NodeType::Pointer p_node0_4 = model_part.CreateNewNode(10, p_node_4->X(), p_node_4->Y(), p_node_4->Z());
            NodeType::Pointer p_node0_5 = model_part.CreateNewNode(11, p_node_5->X(), p_node_5->Y(), p_node_5->Z());
            NodeType::Pointer p_node0_6 = model_part.CreateNewNode(12, p_node_6->X(), p_node_6->Y(), p_node_6->Z());
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            Triangle3D3 <Node<3>> triangle_0( condition_nodes_0 );
            
            std::vector<NodeType::Pointer> condition_nodes0_0 (3);
            condition_nodes0_0[0] = p_node0_1;
            condition_nodes0_0[1] = p_node0_2;
            condition_nodes0_0[2] = p_node0_3;
            Triangle3D3 <Node<3>> triangle0_0( condition_nodes0_0 );
            
            const array_1d<double, 3>& normal_0 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 1, triangle_0, p_cond_prop);
            Condition::Pointer p_cond0_0 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 3, triangle0_0, p_cond_prop);
            p_node_1->SetValue(NORMAL, normal_0);
            p_node_2->SetValue(NORMAL, normal_0);
            p_node_3->SetValue(NORMAL, normal_0);
            p_cond_0->SetValue(NORMAL, normal_0);
            
            p_node0_1->SetValue(NORMAL, normal_0);
            p_node0_2->SetValue(NORMAL, normal_0);
            p_node0_3->SetValue(NORMAL, normal_0);
            p_cond0_0->SetValue(NORMAL, normal_0);
            
            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            condition_nodes_1[0] = p_node_4;
            condition_nodes_1[1] = p_node_5;
            condition_nodes_1[2] = p_node_6;
            Triangle3D3 <Node<3>> triangle_1( condition_nodes_1 );
            
            std::vector<NodeType::Pointer> condition_nodes0_1 (3);
            condition_nodes0_1[0] = p_node0_4;
            condition_nodes0_1[1] = p_node0_5;
            condition_nodes0_1[2] = p_node0_6;
            Triangle3D3 <Node<3>> triangle0_1( condition_nodes0_1 );
            
            const array_1d<double, 3>& normal_1 = triangle_0.UnitNormal(aux_point);
            Condition::Pointer p_cond_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 2, triangle_1, p_cond_prop);
            Condition::Pointer p_cond0_1 = model_part.CreateNewCondition("ALMFrictionlessMortarContactCondition3D3N", 4, triangle0_1, p_cond_prop);
            p_node_4->SetValue(NORMAL, normal_1);
            p_node_5->SetValue(NORMAL, normal_1);
            p_node_6->SetValue(NORMAL, normal_1);
            p_cond_1->SetValue(NORMAL, normal_1);
            
            p_node0_4->SetValue(NORMAL, normal_1);
            p_node0_5->SetValue(NORMAL, normal_1);
            p_node0_6->SetValue(NORMAL, normal_1);
            p_cond0_1->SetValue(NORMAL, normal_1);
            
            // Create and initialize condition variables
            MortarKinematicVariablesWithDerivatives<3, 3> rVariables0; // These are the kinematic variables for the initial configuration
            MortarKinematicVariablesWithDerivatives<3, 3> rVariables; // These are the kinematic variables for the current configuration
            
            // Create the current contact data
            DerivativeData<3, 3> rDerivativeData;
            rDerivativeData.Initialize(triangle_0, model_part.GetProcessInfo());
            
            // We call the exact integration utility
            ExactMortarIntegrationUtility<3, 3, true>  integration_utility = ExactMortarIntegrationUtility<3, 3, true> (2);
            
            const unsigned int number_of_iterations = 6;
            Vector error_vector_slave(number_of_iterations, 0.0);
            Vector error_vector_master(number_of_iterations, 0.0);
            for (unsigned int iter = 0; iter < number_of_iterations; ++iter)
            {                
                // We add displacement to the node 4
                array_1d<double, 3> aux_delta_disp = ZeroVector(3);
                aux_delta_disp[1] = - static_cast<double>(iter + 1) * 5.0e-2;
//                 aux_delta_disp[1] = static_cast<double>(iter + 1) * 5.0e-4;
                Node<3>::Pointer node_to_move = p_node_4;
                array_1d<double, 3>& current_disp = node_to_move->FastGetSolutionStepValue(DISPLACEMENT);
                array_1d<double, 3>& previous_disp = node_to_move->FastGetSolutionStepValue(DISPLACEMENT, 1);
                previous_disp = current_disp;
                current_disp = aux_delta_disp;
                // Finally we move the mesh
                noalias(node_to_move->Coordinates())  = node_to_move->GetInitialPosition().Coordinates() + node_to_move->FastGetSolutionStepValue(DISPLACEMENT);
                
                // Reading integration points
                ConditionArrayListType conditions_points_slave0, conditions_points_slave;
                const bool is_inside0 = integration_utility.GetExactIntegration(triangle0_0, normal_0, triangle0_1, normal_1, conditions_points_slave0);
                const bool is_inside = integration_utility.GetExactIntegration(triangle_0, normal_0, triangle_1, normal_1, conditions_points_slave);
                
                IntegrationMethod this_integration_method = GeometryData::GI_GAUSS_2;
                
                if (is_inside && is_inside0)
                {
                    // Initialize general variables for the current master element
                    rVariables0.Initialize();
                    rVariables.Initialize();
                    
                    // Update slave element info
                    rDerivativeData.UpdateMasterPair(p_cond_1);
                    
                    if (conditions_points_slave.size() == conditions_points_slave0.size()) // Just in case we have the "same configuration"
                    {
//                         // Mathematica debug
//                         std::cout << "\nGraphics3D[{EdgeForm[{Thick,Dashed,Red}],FaceForm[],Polygon[{{";
//                 
//                         for (unsigned int i = 0; i < 3; ++i)
//                         {
//                             std::cout << triangle_0[i].X() << "," << triangle_0[i].Y() << "," << triangle_0[i].Z();
//                             
//                             if (i < 3 - 1) std::cout << "},{";
//                         }
//                         std::cout << "}}],Text[Style["<< p_cond_0->Id() <<", Tiny],{"<< triangle_0.Center().X() << "," << triangle_0.Center().Y() << ","<< triangle_0.Center().Z() << "}]}],";// << std::endl;
//                         
//                         std::cout << "\nGraphics3D[{EdgeForm[{Thick,Dashed,Blue}],FaceForm[],Polygon[{{";
//                 
//                         for (unsigned int i = 0; i < 3; ++i)
//                         {
//                             std::cout << triangle_1[i].X() << "," << triangle_1[i].Y() << "," << triangle_1[i].Z();
//                             
//                             if (i < 3 - 1) std::cout << "},{";
//                         }
//                         
//                         std::cout << "}}],Text[Style["<< p_cond_1->Id() <<", Tiny],{"<< triangle_1.Center().X() << "," << triangle_1.Center().Y() << ","<< triangle_1.Center().Z() << "}]}],";// << std::endl;
//                         
//                         for (unsigned int i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom)
//                         {
//                             std::vector<PointType::Pointer> points_array (3); // The points are stored as local coordinates, we calculate the global coordinates of this points
//                             for (unsigned int i_node = 0; i_node < 3; ++i_node)
//                             {
//                                 PointType global_point;
//                                 triangle_0.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
//                                 points_array[i_node] = boost::make_shared<PointType>(global_point);
//                             }
//                             
//                             TriangleType decomp_geom( points_array );
//                             
//                             std::cout << "\nGraphics3D[{Opacity[.3],Triangle[{{"; 
//                             for (unsigned int i = 0; i < 3; ++i)
//                             {
//                                 std::cout << std::setprecision(16) << decomp_geom[i].X() << "," << decomp_geom[i].Y() << "," << decomp_geom[i].Z();
//                                 
//                                 if (i < 2) std::cout << "},{";
//                             }
//                             std::cout << "}}]}],";// << std::endl;
//                         }
//                         
//                         std::cout << std::endl;
                        
                        for (unsigned int i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom)
                        {
                            std::vector<PointType::Pointer> points_array (3); // The points are stored as local coordinates, we calculate the global coordinates of this points
                            std::vector<PointType::Pointer> points_array0 (3);
                            array_1d<PointBelongsTriangle3D3N, 3> belong_array, belong_array0;
                            for (unsigned int i_node = 0; i_node < 3; ++i_node)
                            {
                                PointType global_point;
                                triangle_0.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                                points_array[i_node] = boost::make_shared<PointType>(global_point);
                                belong_array[i_node] = conditions_points_slave[i_geom][i_node].GetBelong();
                                triangle0_0.GlobalCoordinates(global_point, conditions_points_slave0[i_geom][i_node]);
                                points_array0[i_node] = boost::make_shared<PointType>(global_point);
                                belong_array0[i_node] = conditions_points_slave0[i_geom][i_node].GetBelong();
                            }
                            
                            TriangleType decomp_geom( points_array );
                            TriangleType decomp_geom0( points_array0 );
                            
                            if ((MortarUtilities::HeronCheck(decomp_geom) == false) && (MortarUtilities::HeronCheck(decomp_geom0) == false))
                            {
                                const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );
//                                 const GeometryType::IntegrationPointsArrayType& integration_points_slave0 = decomp_geom0.IntegrationPoints( this_integration_method );
                                
                                // Integrating the mortar operators
                                for ( unsigned int point_number = 0; point_number < integration_points_slave.size(); ++point_number )
                                {
                                    // We reset the derivatives
                                    rDerivativeData.ResetDerivatives();
                                    
                                    // We compute the local coordinates 
                                    const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                                    PointType local_point_parent;
                                    PointType gp_global;
                                    
                                    // We compute the initial configuration
                                    decomp_geom0.GlobalCoordinates(gp_global, local_point_decomp);
                                    triangle0_0.PointLocalCoordinates(local_point_parent, gp_global);
                                    
                                    // SLAVE KINEMATIC COMPUTATIONS
                                    triangle0_0.ShapeFunctionsValues( rVariables0.NSlave, local_point_parent.Coordinates() );
                                    triangle0_0.ShapeFunctionsLocalGradients( rVariables0.DNDeSlave, local_point_parent );
                                    
                                    rVariables0.jSlave = decomp_geom.Jacobian( rVariables0.jSlave, local_point_decomp.Coordinates());
                                    rVariables0.DetjSlave = decomp_geom.DeterminantOfJacobian( local_point_decomp );
            
                                    // MASTER KINEMATIC COMPUTATIONS
                                    PointType projected_gp_global;
                                    array_1d<double,3> gp_normal = MortarUtilities::GaussPointUnitNormal(rVariables0.NSlave, triangle0_0);
                                    
                                    GeometryType::CoordinatesArrayType slave_gp_global;
                                    triangle0_0.GlobalCoordinates( slave_gp_global, local_point_parent );
                                    MortarUtilities::FastProjectDirection( triangle0_1, slave_gp_global, projected_gp_global, normal_1, -gp_normal ); // The opposite direction
                                    
                                    GeometryType::CoordinatesArrayType projected_gp_local;
                                    
                                    triangle0_1.PointLocalCoordinates( projected_gp_local, projected_gp_global.Coordinates( ) ) ;

                                    triangle0_1.ShapeFunctionsValues( rVariables0.NMaster,    projected_gp_local );         
                                    triangle0_1.ShapeFunctionsLocalGradients( rVariables0.DNDeMaster, projected_gp_local );
                                    rVariables0.jMaster = triangle0_1.Jacobian( rVariables0.jMaster, projected_gp_local);

                                    // We compute the current configuration
                                    decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                                    triangle_0.PointLocalCoordinates(local_point_parent, gp_global);
                                    
                                    // SLAVE KINEMATIC COMPUTATIONS
                                    triangle_0.ShapeFunctionsValues( rVariables.NSlave, local_point_parent.Coordinates() );
                                    triangle_0.ShapeFunctionsLocalGradients( rVariables.DNDeSlave, local_point_parent );
                                    
                                    rVariables.jSlave = decomp_geom.Jacobian( rVariables.jSlave, local_point_decomp.Coordinates());
                                    rVariables.DetjSlave = decomp_geom.DeterminantOfJacobian( local_point_decomp );
            
                                    // MASTER KINEMATIC COMPUTATIONS
                                    gp_normal = MortarUtilities::GaussPointUnitNormal(rVariables.NSlave, triangle_0);
                                    
                                    triangle_0.GlobalCoordinates( slave_gp_global, local_point_parent );
                                    MortarUtilities::FastProjectDirection( triangle_1, slave_gp_global, projected_gp_global, normal_1, -gp_normal ); // The opposite direction
                                    
                                    triangle_1.PointLocalCoordinates( projected_gp_local, projected_gp_global.Coordinates( ) ) ;

                                    triangle_1.ShapeFunctionsValues( rVariables.NMaster,    projected_gp_local );         
                                    triangle_1.ShapeFunctionsLocalGradients( rVariables.DNDeMaster, projected_gp_local );
                                    rVariables.jMaster = triangle_1.Jacobian( rVariables.jMaster, projected_gp_local);
                                    
                                    // Now we compute the derivatives
                                    const bool consider_normal_variation = false;
                                    DerivativesUtilitiesType::CalculateDeltaCellVertex(rVariables, rDerivativeData, belong_array, consider_normal_variation, triangle_0, triangle_1, normal_0);
        
                                    // Update the derivatives of the shape functions and the gap
                                    DerivativesUtilitiesType::CalculateDeltaN(rVariables, rDerivativeData, triangle_0, triangle_1, normal_0, normal_1, decomp_geom, local_point_decomp, local_point_parent, consider_normal_variation);
                                                                    
                                    // Now we compute the error of the delta N
                                    Vector aux_N_dx_slave  = rVariables0.NSlave;
                                    Vector aux_N_dx_master = rVariables0.NMaster;
                                    for (unsigned int i_node = 0; i_node < 2 * 3; ++i_node)
                                    {
                                        array_1d<double, 3> delta_disp;
                                        if (i_node < 3) delta_disp = triangle_0[i_node].FastGetSolutionStepValue(DISPLACEMENT);
                                        else delta_disp = triangle_1[i_node - 3].FastGetSolutionStepValue(DISPLACEMENT);
                                        for (unsigned int i_dof = 0; i_dof < 3; ++i_dof)
                                        {
                                            const auto& delta_n1 = rDerivativeData.DeltaN1[i_node * 3 + i_dof];
                                            const auto& delta_n2 = rDerivativeData.DeltaN2[i_node * 3 + i_dof];
                                            aux_N_dx_slave  += delta_n1 * delta_disp[i_dof];
                                            aux_N_dx_master += delta_n2 * delta_disp[i_dof];
                                        }
                                    }
                                    
                                    error_vector_slave[iter]  += norm_2(rVariables.NSlave  - aux_N_dx_slave); 
                                    error_vector_master[iter] += norm_2(rVariables.NMaster - aux_N_dx_master); 
                                }
                            }
                        }
                    }
                    else
                    {
                        KRATOS_ERROR << "YOUR INITIAL SPLITTING DOES NOT COINCIDE WITH THE CURRENT ONE" << std::endl;
                    }
                }
                else
                {
                    KRATOS_ERROR << "WRONG, YOU ARE SUPPOSED TO HAVE AN INTERSECTION" << std::endl;
                }
            }
            
            const double tolerance = 1.0e-6;
            for (unsigned int iter = 0; iter < number_of_iterations; ++iter)
            {
                KRATOS_CHECK_LESS_EQUAL(error_vector_slave[iter], tolerance);
                KRATOS_CHECK_LESS_EQUAL(error_vector_master[iter], tolerance);
            }
        }
        
    } // namespace Testing
}  // namespace Kratos.
