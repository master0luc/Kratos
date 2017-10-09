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
//   Date:                $Date:                   December 2016 $
//   Revision:            $Revision:                         0.0 $
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/python.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/optimization_utilities.h"
#include "custom_utilities/geometry_utilities.h"
#include "custom_utilities/vertex_morphing_mapper.h"
#include "custom_utilities/response_functions/strain_energy_response_function.h"
#include "custom_utilities/response_functions/mass_response_function.h"
#include "custom_utilities/cad_reconstruction/cad_mapper.h"
#include "linear_solvers/linear_solver.h"

// ==============================================================================

namespace Kratos
{

namespace Python
{


void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;

    typedef UblasSpace<double, CompressedMatrix, Vector> CompressedSpaceType;
    typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;
    typedef LinearSolver<CompressedSpaceType, DenseSpaceType > SparseLinearSolverType;

    // ================================================================
    // For perfoming the mapping according to Vertex Morphing
    // ================================================================
    class_<VertexMorphingMapper, bases<Process> >("VertexMorphingMapper", init<ModelPart&, std::string, bool, double, bool, boost::python::list>())
            .def("compute_mapping_matrix", &VertexMorphingMapper::compute_mapping_matrix)
            .def("map_sensitivities_to_design_space", &VertexMorphingMapper::map_sensitivities_to_design_space)
            .def("map_design_update_to_geometry_space", &VertexMorphingMapper::map_design_update_to_geometry_space)
            ;

    // ========================================================================
    // For performing individual steps of an optimization algorithm
    // ========================================================================
    class_<OptimizationUtilities, bases<Process> >("OptimizationUtilities", init<ModelPart&, boost::python::dict, boost::python::dict, double, bool>())

            // ----------------------------------------------------------------
            // For running unconstrained descent methods
            // ----------------------------------------------------------------
            .def("compute_search_direction_steepest_descent", &OptimizationUtilities::compute_search_direction_steepest_descent)

            // ----------------------------------------------------------------
            // For running augmented Lagrange method
            // ----------------------------------------------------------------
            .def("initialize_augmented_lagrange", &OptimizationUtilities::initialize_augmented_lagrange)
            .def("compute_search_direction_augmented_lagrange", &OptimizationUtilities::compute_search_direction_augmented_lagrange)
            .def("udpate_augmented_lagrange_parameters", &OptimizationUtilities::udpate_augmented_lagrange_parameters)
            .def("get_penalty_fac", &OptimizationUtilities::get_penalty_fac)
            .def("get_lambda", &OptimizationUtilities::get_lambda)
            .def("get_value_of_augmented_lagrangian", &OptimizationUtilities::get_value_of_augmented_lagrangian)

            // ----------------------------------------------------------------
            // For running penalized projection method
            // ----------------------------------------------------------------
            .def("compute_projected_search_direction", &OptimizationUtilities::compute_projected_search_direction)
        	.def("correct_projected_search_direction", &OptimizationUtilities::correct_projected_search_direction)

            // ----------------------------------------------------------------
            // General optimization operations
            // ----------------------------------------------------------------
            .def("compute_design_update", &OptimizationUtilities::compute_design_update)
            ;

    // ========================================================================
    // For pre- and post-processing of geometry data
    // ========================================================================
    class_<GeometryUtilities, bases<Process> >("GeometryUtilities", init<ModelPart&>())
            .def("compute_unit_surface_normals", &GeometryUtilities::compute_unit_surface_normals)
            .def("project_grad_on_unit_surface_normal", &GeometryUtilities::project_grad_on_unit_surface_normal)
            .def("extract_surface_nodes", &GeometryUtilities::extract_surface_nodes)
            ;

    // ========================================================================
    // For calculations related to response functions
    // ========================================================================
    class_<StrainEnergyResponseFunction, bases<Process> >("StrainEnergyResponseFunction", init<ModelPart&, boost::python::dict>())
            .def("initialize", &StrainEnergyResponseFunction::initialize)
            .def("calculate_value", &StrainEnergyResponseFunction::calculate_value)
            .def("calculate_gradient", &StrainEnergyResponseFunction::calculate_gradient) 
            .def("get_value", &StrainEnergyResponseFunction::get_value)
            .def("get_initial_value", &StrainEnergyResponseFunction::get_initial_value)  
            .def("get_gradient", &StrainEnergyResponseFunction::get_gradient)                              
            ; 
    class_<MassResponseFunction, bases<Process> >("MassResponseFunction", init<ModelPart&, boost::python::dict>())
            .def("initialize", &MassResponseFunction::initialize)
            .def("calculate_value", &MassResponseFunction::calculate_value)
            .def("calculate_gradient", &MassResponseFunction::calculate_gradient)  
            .def("get_value", &MassResponseFunction::get_value)
            .def("get_initial_value", &MassResponseFunction::get_initial_value) 
            .def("get_gradient", &MassResponseFunction::get_gradient)                              
            ;                 

    // ========================================================================
    // For CAD reconstruction
    // ========================================================================
    class_<CADMapper, bases<Process> >("CADMapper", init<ModelPart&, boost::python::dict, boost::python::dict, SparseLinearSolverType::Pointer>())
            .def("compute_mapping_matrix", &CADMapper::compute_mapping_matrix)
            .def("apply_boundary_conditions", &CADMapper::apply_boundary_conditions)
            .def("map_to_cad_space", &CADMapper::map_to_cad_space)
            .def("output_gauss_points", &CADMapper::output_gauss_points)
            .def("output_surface_points", &CADMapper::output_surface_points)
            .def("output_boundary_loop_points", &CADMapper::output_boundary_loop_points)
            .def("output_control_point_displacements", &CADMapper::output_control_point_displacements)
            .def("output_surface_border_points", &CADMapper::output_surface_border_points)
            .def("output_surface_border_points_two", &CADMapper::output_surface_border_points_two)
            .def("compute_nearest_points", &CADMapper::compute_nearest_points)            
            .def("compute_a_matrix", &CADMapper::compute_a_matrix)
            .def("map_to_cad_space_2", &CADMapper::map_to_cad_space_2)            
            .def("print_nearest_points", &CADMapper::print_nearest_points)            
            .def("compute_real_length", &CADMapper::compute_real_length)            
            .def("compute_lhs_matrix", &CADMapper::compute_lhs_matrix)            
            .def("compute_rhs_vector", &CADMapper::compute_rhs_vector)            
            .def("map_to_cad_space_3", &CADMapper::map_to_cad_space_3)            
            // ONE COORDINATE AT A TIME //
            .def("map_to_cad_space_4", &CADMapper::map_to_cad_space_4)  
            .def("apply_boundary_conditions_small", &CADMapper::apply_boundary_conditions_small)  
            .def("measure_g0_continuity", &CADMapper::measure_g0_continuity)  
            .def("measure_g1_continuity", &CADMapper::measure_g1_continuity)
            // MODULAR //
            .def("use_all_FE_nodes_as_data_points", &CADMapper::use_all_FE_nodes_as_data_points)            
            .def("parametrisation", &CADMapper::parametrisation)            
            .def("print_nearest_points_2", &CADMapper::print_nearest_points_2)            
            .def("apply_regularization_schemes", &CADMapper::apply_regularization_schemes)  
            .def("apply_penalty_factors", &CADMapper::apply_penalty_factors)  
            .def("map_all_patches", &CADMapper::map_all_patches)  
            .def("map_patch_by_patch", &CADMapper::map_patch_by_patch)  
            .def("map_boundary_conditions", &CADMapper::map_boundary_conditions)  
            .def("crazy_step_back", &CADMapper::crazy_step_back)  
            .def("map_boundary_conditions_augmented_Lagrange", &CADMapper::map_boundary_conditions_augmented_Lagrange)  
            .def("map_all_patches_augmented_Lagrange", &CADMapper::map_all_patches_augmented_Lagrange)  
          
            .def("write_updated_georhino_file", &CADMapper::write_updated_georhino_file)  
          
            // EXTERNAL: separating FE-mesh data from computation //
            .def("external_map_to_cad_space", &CADMapper::external_map_to_cad_space)            
            .def("set_point", &CADMapper::set_point)            
            .def("compute_objective", &CADMapper::compute_objective)            

            .def("compare_lhs", &CADMapper::compare_lhs)            
            .def("compare_rhs", &CADMapper::compare_rhs)            

            ;                      
}


}  // namespace Python.

} // Namespace Kratos

