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
#include "custom_conditions/ALM_frictionless_mortar_contact_condition.h"

namespace Kratos 
{
    namespace Testing 
    {

    typedef ModelPart::IndexType           IndexType;
    typedef ModelPart::NodeIterator NodeIteratorType;

    /** Checks the EmbeddedNavierStokes2D3N element.
    * Checks the LHS and RHS computation.
    */
    KRATOS_TEST_CASE_IN_SUITE(TestALMFrictionlessMortarContactCondition3D3N, ContactStructuralApplicationFastSuite)
    {
//         ModelPart modelPart("Main");
//         modelPart.SetBufferSize(3);
// 
//         // Variables addition
//         modelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
// 
//         // Process info creation
//         double delta_time = 0.1;
//         modelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
// 
// 
//         // Set the element properties
//         Properties::Pointer pElemProp = modelPart.pGetProperties(0);
//         pElemProp->SetValue(DENSITY, 1000.0);
// 
//         // Geometry creation
//         modelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
//         modelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
//         modelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
//         std::vector<ModelPart::IndexType> elemNodes {1, 2, 3};
//         modelPart.CreateNewElement("ALMFrictionlessMortarContactCondition3D3N", 1, elemNodes, pElemProp);
// 
//         Element::Pointer pElement = modelPart.pGetElement(1);
// 
//         // Compute RHS and LHS
//         Vector RHS = ZeroVector(9);
//         Matrix LHS = ZeroMatrix(9,9);
// 
//         pElement->Initialize(); // Initialize the element to initialize the constitutive law
//         pElement->CalculateLocalSystem(LHS, RHS, modelPart.GetProcessInfo());
// 
//         // Check the RHS values (the RHS is computed as the LHS x previous_solution, 
//         // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
// //             KRATOS_CHECK_NEAR(RHS(0), 0.0475309, 1e-7);
// //             KRATOS_CHECK_NEAR(RHS(1), 0.0975309, 1e-7);
// //             KRATOS_CHECK_NEAR(RHS(2), -0.0545696, 1e-7);
// //             KRATOS_CHECK_NEAR(RHS(3), 0.0469136, 1e-7);
// //             KRATOS_CHECK_NEAR(RHS(4), 0.0969136, 1e-7);
// //             KRATOS_CHECK_NEAR(RHS(5), 0.0176796, 1e-7);
// //             KRATOS_CHECK_NEAR(RHS(6), 45172.3, 1e-1);
// //             KRATOS_CHECK_NEAR(RHS(7), 92991.6, 1e-1);
// //             KRATOS_CHECK_NEAR(RHS(8), 0.0202233, 1e-7);
}
        
    } // namespace Testing
}  // namespace Kratos.
