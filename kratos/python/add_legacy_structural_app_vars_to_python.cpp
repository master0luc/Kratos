//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/legacy_structural_app_vars.h"



namespace Kratos
{

namespace Python
{
using namespace boost::python;

void  AddLegacyStructuralAppVarsToPython()
{
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAMBDAS_T )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( DELTA_LAMBDAS_T )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONTACT_LINK_M )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( AUXILIARY_MATRIX_1 )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( ELASTIC_LEFT_CAUCHY_GREEN_OLD )

        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_OLD )
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_DT )
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_NULL )
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_NULL_DT )
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ACCELERATION_NULL )
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_EINS )
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_EINS_DT )
        KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( ACCELERATION_EINS )

        KRATOS_REGISTER_IN_PYTHON_VARIABLE( PENALTY_T )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( LAMBDAS )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( GAPS )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( DELTA_LAMBDAS )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( STICK )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONTACT_STICK )

        KRATOS_REGISTER_IN_PYTHON_VARIABLE( USE_DISTRIBUTED_PROPERTIES )

        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONTACT_RAMP )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONTACT_SLAVE_INTEGRATION_POINT_INDEX )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CONTACT_DOUBLE_CHECK )
//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( FIRST_TIME_STEP )
//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( QUASI_STATIC_ANALYSIS )


//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( OSS_SWITCH )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( WRINKLING_APPROACH )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( CALCULATE_INSITU_STRESS )
//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( PERIODIC_PAIR_INDEX )
//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( STATIONARY )

//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( PARTITION_MASK )

//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( FACE_HEAT_FLUX )

//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( NODAL_VOLUME )

        KRATOS_REGISTER_IN_PYTHON_VARIABLE( PRESSURE_DT )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( WATER_PRESSURE_DT )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( WATER_PRESSURE_NULL )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( WATER_PRESSURE_NULL_DT )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( WATER_PRESSURE_NULL_ACCELERATION )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( WATER_PRESSURE_EINS )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( WATER_PRESSURE_EINS_DT )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( WATER_PRESSURE_EINS_ACCELERATION )

        KRATOS_REGISTER_IN_PYTHON_VARIABLE( AIR_PRESSURE_DT )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( AIR_PRESSURE_ACCELERATION )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( AIR_PRESSURE_NULL )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( AIR_PRESSURE_NULL_DT )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( AIR_PRESSURE_NULL_ACCELERATION )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( AIR_PRESSURE_EINS )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( AIR_PRESSURE_EINS_DT )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( AIR_PRESSURE_EINS_ACCELERATION )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( SUCTION )
//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( FLAG_VARIABLE )

//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( DP_EPSILON )
//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( DP_ALPHA1 )
//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( DP_K )

//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( EQ_STRAIN_RATE )
//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( RHS_WATER )
//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( RHS_AIR )

        KRATOS_REGISTER_IN_PYTHON_VARIABLE( PERMEABILITY_28_DAYS )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( PERMEABILITY_1_DAY )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( PERMEABILITY_TRANSITION )

//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( TEMPERATURE_OLD_IT )
//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( EFFECTIVE_VISCOSITY )
//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( KINEMATIC_VISCOSITY)
//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( DYNAMIC_VISCOSITY)
//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( WEIGHT_FATHER_NODES )

        KRATOS_REGISTER_IN_PYTHON_VARIABLE( PARENT_ELEMENT_ID )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( INTEGRATION_POINT_INDEX )

        KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_CONTACT_MASTER )
        KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_CONTACT_SLAVE )
//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_BOUNDARY )
//        KRATOS_REGISTER_IN_PYTHON_VARIABLE( IS_VISITED )
}
}  // namespace Python.
} // Namespace Kratos

#undef KRATOS_REGISTER_IN_PYTHON_IN_PYTHON_FLAG
#undef KRATOS_REGISTER_IN_PYTHON_IN_PYTHON_VARIABLE
#undef KRATOS_REGISTER_IN_PYTHON_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS
