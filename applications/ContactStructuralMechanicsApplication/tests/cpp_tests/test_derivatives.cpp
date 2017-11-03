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
            ModelPart ModelPart("Main");
            
            // First we create the nodes 
            NodeType::Pointer p_node_1 = ModelPart.CreateNewNode(0,-0.2,0.1,0.0);
            NodeType::Pointer p_node_2 = ModelPart.CreateNewNode(1,1.0,0.1,0.0);
            NodeType::Pointer p_node_3 = ModelPart.CreateNewNode(2,0.2,1.2,0.0);
            NodeType::Pointer p_node_4 = ModelPart.CreateNewNode(3,0.6,0.4,0.0);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> condition_nodes_0 (3);
            
            condition_nodes_0[0] = p_node_1;
            condition_nodes_0[1] = p_node_2;
            condition_nodes_0[2] = p_node_3;
            
            Triangle3D3 <Node<3>> triangle0( condition_nodes_0 );
            
            std::vector<NodeType::Pointer> condition_nodes_1 (3);
            
            condition_nodes_1[0] = p_node_1;
            condition_nodes_1[1] = p_node_2;
            condition_nodes_1[2] = p_node_4;
            
            Triangle3D3 <Node<3>> triangle_1( condition_nodes_1 );
            
            std::vector<NodeType::Pointer> condition_nodes_2 (3);
            
            condition_nodes_2[0] = p_node_2;
            condition_nodes_2[1] = p_node_3;
            condition_nodes_2[2] = p_node_4;
            
            Triangle3D3 <Node<3>> triangle_2( condition_nodes_2 );
            
            std::vector<NodeType::Pointer> condition_nodes_3 (3);
            
            condition_nodes_3[0] = p_node_3;
            condition_nodes_3[1] = p_node_1;
            condition_nodes_3[2] = p_node_4;
            
            Triangle3D3 <Node<3>> triangle_3( condition_nodes_3 );
            
//             // We calculate the integral of the mass matrix (assuming constant density)
//             GeometryNodeType::IntegrationPointsArrayType integration_points = Quadrature<TriangleGaussLegendreIntegrationPoints2, 2, IntegrationPoint<3> >::GenerateIntegrationPoints();
//             
//             bounded_matrix<double, 3, 3> mass_matrix_0 = ZeroMatrix(3, 3);
//             
//             for (unsigned int point_number = 0; point_number < integration_points.size(); point_number++)
//             {
//                 Vector N;
//                 const PointType& local_point = integration_points[point_number].Coordinates();
//                 triangle0.ShapeFunctionsValues( N, local_point );
//                 const double det_j = triangle0.DeterminantOfJacobian( local_point );
//                 const double weight = integration_points[point_number].Weight();
//                 
//                 for (unsigned int inode = 0; inode < 3; inode++)
//                 {
//                     for (unsigned int jnode = 0; jnode < 3; jnode++)
//                     {
//                         mass_matrix_0(inode, jnode) += det_j * weight * N(inode) * N(jnode);
//                     }
//                 }
//             }
//      
//             bounded_matrix<double, 3, 3> mass_matrix_1 = ZeroMatrix(3, 3);
//             
//             for (unsigned int point_number = 0; point_number < integration_points.size(); point_number++)
//             {
//                 Vector N1;
//                 Vector N2;
//                 Vector N3;
//                 
//                 const PointType& local_point_0 = integration_points[point_number].Coordinates();
//                 
//                 PointType gp_global;
//                 
//                 triangle_1.GlobalCoordinates(gp_global, local_point_0);
//                 PointType local_point_1;
//                 triangle0.PointLocalCoordinates(local_point_1, gp_global);
//                 
//                 triangle_2.GlobalCoordinates(gp_global, local_point_0);
//                 PointType local_point_2; 
//                 triangle0.PointLocalCoordinates(local_point_2, gp_global);
//                 
//                 triangle_3.GlobalCoordinates(gp_global, local_point_0);
//                 PointType local_point_3; 
//                 triangle0.PointLocalCoordinates(local_point_3, gp_global);
//                 
//                 triangle0.ShapeFunctionsValues( N1, local_point_1 );
//                 triangle0.ShapeFunctionsValues( N2, local_point_2 );
//                 triangle0.ShapeFunctionsValues( N3, local_point_3 );
//                 
//                 const double det_j_1 = triangle_1.DeterminantOfJacobian( local_point_0 );
//                 const double det_j_2 = triangle_2.DeterminantOfJacobian( local_point_0 );
//                 const double det_j_3 = triangle_3.DeterminantOfJacobian( local_point_0 );
//                 
//                 const double weight = integration_points[point_number].Weight();
//                 
//                 for (unsigned int inode = 0; inode < 3; inode++)
//                 {
//                     for (unsigned int jnode = 0; jnode < 3; jnode++)
//                     {
//                         mass_matrix_1(inode, jnode) += det_j_1 * weight * N1[inode] * N1[jnode] \
//                                                    + det_j_2 * weight * N2[inode] * N2[jnode] \
//                                                    + det_j_3 * weight * N3[inode] * N3[jnode];
//                     }
//                 }
//             }
//             
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
