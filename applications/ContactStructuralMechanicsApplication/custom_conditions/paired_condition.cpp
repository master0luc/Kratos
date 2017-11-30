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
#include "custom_conditions/paired_condition.h"

namespace Kratos 
{
/************************************* OPERATIONS **********************************/
/***********************************************************************************/

Condition::Pointer PairedCondition::Create( 
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new PairedCondition( NewId, this->GetGeometry().Create( rThisNodes ), pProperties ));
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer PairedCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer( new PairedCondition( NewId, pGeom, pProperties ));
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer PairedCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pPairedGeom) const
{
    return Condition::Pointer( new PairedCondition( NewId, pGeom, pProperties, pPairedGeom));
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

PairedCondition::~PairedCondition( )
= default;
} // Namespace Kratos
