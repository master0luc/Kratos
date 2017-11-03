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
//



// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "python/matrix_python_interface.h"


namespace Kratos
{
namespace Python
{

using namespace boost::python;


void  AddCompressedMatrixToPython()
{
    MatrixPythonInterface<compressed_matrix<double> >::CreateInterface("CompressedMatrix")
    .def(init<compressed_matrix<double>::size_type, compressed_matrix<double>::size_type>())
    ;
}

}  // namespace Python.

} // Namespace Kratos
