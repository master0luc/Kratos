//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

//Utilities
#ifdef INCLUDE_FEAST
#include "custom_utilities/feast_condition_number_utility.h"
#endif

namespace Kratos
{
namespace Python
{
void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;
  
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    
#ifdef INCLUDE_FEAST
    // Feast condition number utility
    class_<FEASTConditionNumberUtility<SparseSpaceType, LocalSpaceType>>("FEASTConditionNumberUtility", init<>())
    .def("GetConditionNumber",&FEASTConditionNumberUtility<SparseSpaceType, LocalSpaceType>::GetConditionNumber)
    ;
#endif
}

}  // namespace Python.

} // Namespace Kratos

