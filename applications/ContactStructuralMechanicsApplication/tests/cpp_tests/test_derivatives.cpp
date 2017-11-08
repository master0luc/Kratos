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

        template<unsigned int TDim, unsigned int TNumNodes>
        void TestDerivativesShapeFunction(
            ModelPart& ThisModelPart,
            Condition::Pointer SlaveCondition0,
            Condition::Pointer MasterCondition0,
            Condition::Pointer SlaveCondition1,
            Condition::Pointer MasterCondition1,
            const unsigned int NumberIterations,
            const bool Check = true
            )
        {
            // Type definitions
            typedef PointBelong<TNumNodes> PointBelongType;
            typedef array_1d<PointBelongType, TDim> ConditionArrayType;
            typedef typename std::vector<ConditionArrayType> ConditionArrayListType;
//             typedef Line2D2<PointType> LineType;
            typedef Triangle3D3<PointType> TriangleType;
//             typedef typename std::conditional<TDim == 2, LineType, TriangleType >::type DecompositionType;
            typedef DerivativesUtilities<TDim, TNumNodes, false> DerivativesUtilitiesType;
            
            GeometryType& slave_geometry_0 = SlaveCondition0->GetGeometry();
            GeometryType& master_geometry_0 = MasterCondition0->GetGeometry();
            GeometryType& slave_geometry_1 = SlaveCondition1->GetGeometry();
            GeometryType& master_geometry_1 = MasterCondition1->GetGeometry();
            
            const array_1d<double, 3> normal_0 = SlaveCondition1->GetValue(NORMAL);
            const array_1d<double, 3> normal_1 = MasterCondition1->GetValue(NORMAL);
            
            // Create and initialize condition variables
            MortarKinematicVariablesWithDerivatives<TDim, TNumNodes> rVariables0; // These are the kinematic variables for the initial configuration
            MortarKinematicVariablesWithDerivatives<TDim, TNumNodes> rVariables; // These are the kinematic variables for the current configuration
            
            // Create the current contact data
            DerivativeData<3, 3> rDerivativeData;
            rDerivativeData.Initialize(slave_geometry_1, ThisModelPart.GetProcessInfo());
            
            // We call the exact integration utility
            ExactMortarIntegrationUtility<TDim, TNumNodes, true>  integration_utility = ExactMortarIntegrationUtility<TDim, TNumNodes, true> (2);
            
            Vector error_vector_slave(NumberIterations, 0.0);
            Vector error_vector_master(NumberIterations, 0.0);
            for (unsigned int iter = 0; iter < NumberIterations; ++iter)
            {                
                // We add displacement to the node 4
                array_1d<double, 3> aux_delta_disp = ZeroVector(3);
                aux_delta_disp[1] = - static_cast<double>(iter + 1) * 5.0e-2;
//                 aux_delta_disp[1] = static_cast<double>(iter + 1) * 5.0e-4;
                auto& node_to_move = master_geometry_1[0];
                array_1d<double, 3>& current_disp = node_to_move.FastGetSolutionStepValue(DISPLACEMENT);
                array_1d<double, 3>& previous_disp = node_to_move.FastGetSolutionStepValue(DISPLACEMENT, 1);
                previous_disp = current_disp;
                current_disp = aux_delta_disp;
                // Finally we move the mesh
                noalias(node_to_move.Coordinates())  = node_to_move.GetInitialPosition().Coordinates() + node_to_move.FastGetSolutionStepValue(DISPLACEMENT);
                
                // Reading integration points
                ConditionArrayListType conditions_points_slave0, conditions_points_slave;
                const bool is_inside0 = integration_utility.GetExactIntegration(slave_geometry_0, normal_0, master_geometry_0, normal_1, conditions_points_slave0);
                const bool is_inside = integration_utility.GetExactIntegration(slave_geometry_1, normal_0, master_geometry_1, normal_1, conditions_points_slave);
                
                IntegrationMethod this_integration_method = GeometryData::GI_GAUSS_2;
                
                if (is_inside && is_inside0)
                {
                    // Initialize general variables for the current master element
                    rVariables0.Initialize();
                    rVariables.Initialize();
                    
                    // Update slave element info
                    rDerivativeData.UpdateMasterPair(MasterCondition1);
                    
                    if (conditions_points_slave.size() == conditions_points_slave0.size()) // Just in case we have the "same configuration"
                    {
//                         // Mathematica debug
//                         std::cout << "\nGraphics3D[{EdgeForm[{Thick,Dashed,Red}],FaceForm[],Polygon[{{";
//                 
//                         for (unsigned int i = 0; i < TNumNodes; ++i)
//                         {
//                             std::cout << slave_geometry_1[i].X() << "," << slave_geometry_1[i].Y() << "," << slave_geometry_1[i].Z();
//                             
//                             if (i < TNumNodes - 1) std::cout << "},{";
//                         }
//                         std::cout << "}}],Text[Style["<< SlaveCondition1->Id() <<", Tiny],{"<< slave_geometry_1.Center().X() << "," << slave_geometry_1.Center().Y() << ","<< slave_geometry_1.Center().Z() << "}]}],";// << std::endl;
//                         
//                         std::cout << "\nGraphics3D[{EdgeForm[{Thick,Dashed,Blue}],FaceForm[],Polygon[{{";
//                 
//                         for (unsigned int i = 0; i < TNumNodes; ++i)
//                         {
//                             std::cout << master_geometry_1[i].X() << "," << master_geometry_1[i].Y() << "," << master_geometry_1[i].Z();
//                             
//                             if (i < TNumNodes - 1) std::cout << "},{";
//                         }
//                         
//                         std::cout << "}}],Text[Style["<< MasterCondition1->Id() <<", Tiny],{"<< master_geometry_1.Center().X() << "," << master_geometry_1.Center().Y() << ","<< master_geometry_1.Center().Z() << "}]}],";// << std::endl;
//                         
//                         for (unsigned int i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom)
//                         {
//                             std::vector<PointType::Pointer> points_array (3); // The points are stored as local coordinates, we calculate the global coordinates of this points
//                             for (unsigned int i_node = 0; i_node < 3; ++i_node)
//                             {
//                                 PointType global_point;
//                                 slave_geometry_1.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
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
                            for (unsigned int i_node = 0; i_node < TDim; ++i_node)
                            {
                                PointType global_point;
                                slave_geometry_1.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                                points_array[i_node] = boost::make_shared<PointType>(global_point);
                                belong_array[i_node] = conditions_points_slave[i_geom][i_node].GetBelong();
                                slave_geometry_0.GlobalCoordinates(global_point, conditions_points_slave0[i_geom][i_node]);
                                points_array0[i_node] = boost::make_shared<PointType>(global_point);
                                belong_array0[i_node] = conditions_points_slave0[i_geom][i_node].GetBelong();
                            }
                            
                            if (Check == false) KRATOS_WATCH(belong_array);
                            
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
                                    slave_geometry_0.PointLocalCoordinates(local_point_parent, gp_global);
                                    
                                    // SLAVE KINEMATIC COMPUTATIONS
                                    slave_geometry_0.ShapeFunctionsValues( rVariables0.NSlave, local_point_parent.Coordinates() );
                                    slave_geometry_0.ShapeFunctionsLocalGradients( rVariables0.DNDeSlave, local_point_parent );
                                    
                                    rVariables0.jSlave = decomp_geom.Jacobian( rVariables0.jSlave, local_point_decomp.Coordinates());
                                    rVariables0.DetjSlave = decomp_geom.DeterminantOfJacobian( local_point_decomp );
            
                                    // MASTER KINEMATIC COMPUTATIONS
                                    PointType projected_gp_global;
                                    array_1d<double,3> gp_normal = MortarUtilities::GaussPointUnitNormal(rVariables0.NSlave, slave_geometry_0);
                                    
                                    GeometryType::CoordinatesArrayType slave_gp_global;
                                    slave_geometry_0.GlobalCoordinates( slave_gp_global, local_point_parent );
                                    MortarUtilities::FastProjectDirection( master_geometry_0, slave_gp_global, projected_gp_global, normal_1, -gp_normal ); // The opposite direction
                                    
                                    GeometryType::CoordinatesArrayType projected_gp_local;
                                    
                                    master_geometry_0.PointLocalCoordinates( projected_gp_local, projected_gp_global.Coordinates( ) ) ;

                                    master_geometry_0.ShapeFunctionsValues( rVariables0.NMaster,    projected_gp_local );         
                                    master_geometry_0.ShapeFunctionsLocalGradients( rVariables0.DNDeMaster, projected_gp_local );
                                    rVariables0.jMaster = master_geometry_0.Jacobian( rVariables0.jMaster, projected_gp_local);

                                    // We compute the current configuration
                                    decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                                    slave_geometry_1.PointLocalCoordinates(local_point_parent, gp_global);
                                    
                                    // SLAVE KINEMATIC COMPUTATIONS
                                    slave_geometry_1.ShapeFunctionsValues( rVariables.NSlave, local_point_parent.Coordinates() );
                                    slave_geometry_1.ShapeFunctionsLocalGradients( rVariables.DNDeSlave, local_point_parent );
                                    
                                    rVariables.jSlave = decomp_geom.Jacobian( rVariables.jSlave, local_point_decomp.Coordinates());
                                    rVariables.DetjSlave = decomp_geom.DeterminantOfJacobian( local_point_decomp );
            
                                    // MASTER KINEMATIC COMPUTATIONS
                                    gp_normal = MortarUtilities::GaussPointUnitNormal(rVariables.NSlave, slave_geometry_1);
                                    
                                    slave_geometry_1.GlobalCoordinates( slave_gp_global, local_point_parent );
                                    MortarUtilities::FastProjectDirection( master_geometry_1, slave_gp_global, projected_gp_global, normal_1, -gp_normal ); // The opposite direction
                                    
                                    master_geometry_1.PointLocalCoordinates( projected_gp_local, projected_gp_global.Coordinates( ) ) ;

                                    master_geometry_1.ShapeFunctionsValues( rVariables.NMaster,    projected_gp_local );         
                                    master_geometry_1.ShapeFunctionsLocalGradients( rVariables.DNDeMaster, projected_gp_local );
                                    rVariables.jMaster = master_geometry_1.Jacobian( rVariables.jMaster, projected_gp_local);
                                    
                                    // Now we compute the derivatives
                                    const bool consider_normal_variation = false;
                                    DerivativesUtilitiesType::CalculateDeltaCellVertex(rVariables, rDerivativeData, belong_array, consider_normal_variation, slave_geometry_1, master_geometry_1, normal_0);
        
                                    // Update the derivatives of the shape functions and the gap
                                    DerivativesUtilitiesType::CalculateDeltaN(rVariables, rDerivativeData, slave_geometry_1, master_geometry_1, normal_0, normal_1, decomp_geom, local_point_decomp, local_point_parent, consider_normal_variation);
                                                                    
                                    // Now we compute the error of the delta N
                                    Vector aux_N_dx_slave  = rVariables0.NSlave;
                                    Vector aux_N_dx_master = rVariables0.NMaster;
                                    for (unsigned int i_node = 0; i_node < 2 * TNumNodes; ++i_node)
                                    {
                                        array_1d<double, 3> delta_disp;
                                        if (i_node < 3) delta_disp = slave_geometry_1[i_node].FastGetSolutionStepValue(DISPLACEMENT);
                                        else delta_disp = master_geometry_1[i_node - 3].FastGetSolutionStepValue(DISPLACEMENT);
                                        for (unsigned int i_dof = 0; i_dof < TDim; ++i_dof)
                                        {
                                            const auto& delta_n1 = rDerivativeData.DeltaN1[i_node * TNumNodes + i_dof];
                                            const auto& delta_n2 = rDerivativeData.DeltaN2[i_node * TNumNodes + i_dof];
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
            
            if (Check == true)
            {
                const double tolerance = 1.0e-6;
                for (unsigned int iter = 0; iter < NumberIterations; ++iter)
                {
                    KRATOS_CHECK_LESS_EQUAL(error_vector_slave[iter], tolerance);
                    KRATOS_CHECK_LESS_EQUAL(error_vector_master[iter], tolerance);
                }
            }
            else // DEBUG
            {
                KRATOS_WATCH(error_vector_slave)
                KRATOS_WATCH(error_vector_master)
            }
        }
        
        /** 
         * Checks if the derivatives of the shape functions work as expected
         * Case 1 of the Triangle3D3
         */
    
        KRATOS_TEST_CASE_IN_SUITE(TestDerivativesShapeFunctionTriangle1, ContactStructuralApplicationFastSuite)
        {
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
            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5, 0.0,0.0,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.0,0.0,1.0e-3);
            
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
            
            TestDerivativesShapeFunction<3,3>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, 6, true);
        }
        
        /** 
         * Checks if the derivatives of the shape functions work as expected
         * Case 2 of the Triangle3D3
         */
    
        KRATOS_TEST_CASE_IN_SUITE(TestDerivativesShapeFunctionTriangle2, ContactStructuralApplicationFastSuite)
        {
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
            
            NodeType::Pointer p_node_4 = model_part.CreateNewNode(4,-0.1,1.0,1.0e-3);
            NodeType::Pointer p_node_5 = model_part.CreateNewNode(5, 0.0,0.0,1.0e-3);
            NodeType::Pointer p_node_6 = model_part.CreateNewNode(6, 1.0,0.0,1.0e-3);
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
            
            TestDerivativesShapeFunction<3,3>( model_part, p_cond0_0, p_cond0_1, p_cond_0, p_cond_1, 1, false);
        }
        
    } // namespace Testing
}  // namespace Kratos.
