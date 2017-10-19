// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
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
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/optimization_utilities.h"
#include "custom_utilities/geometry_utilities.h"
#include "custom_utilities/mapping/mapper_vertex_morphing.h"
#include "custom_utilities/mapping/mapper_vertex_morphing_matrix_free.h"
#include "custom_utilities/mapping/mapper_vertex_morphing_improved_integration.h"
#include "custom_utilities/damping/damping_utilities.h"
#include "custom_utilities/response_functions/strain_energy_response_function.h"
#include "custom_utilities/response_functions/mass_response_function.h"
#include "custom_utilities/input_output/universal_file_io.h"
#include "custom_utilities/input_output/vtk_file_io.h"

#include "linear_solvers/linear_solver.h"
#include "custom_utilities/cad_reconstruction/reconstruction_conditions/reconstruction_condition_container.h"
#include "custom_utilities/cad_reconstruction/cad_reconstruction_solver.h"
#include "custom_utilities/cad_reconstruction/data_management/reconstruction_output_writer.h"
#include "custom_utilities/cad_reconstruction/data_management/quality_evaluation_utility.h"

// ==============================================================================

namespace Kratos
{

namespace Python
{


void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;

    typedef UblasSpace<double, SparseMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, CompressedMatrix, Vector> CompressedSpaceType;
    typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;
    typedef LinearSolver<CompressedSpaceType, DenseSpaceType > CompressedLinearSolverType;

    // ================================================================
    // For perfoming the mapping according to Vertex Morphing
    // ================================================================
    class_<MapperVertexMorphing, bases<Process> >("MapperVertexMorphing", init<ModelPart&, Parameters&>())
        .def("MapToDesignSpace", &MapperVertexMorphing::MapToDesignSpace)
        .def("MapToGeometrySpace", &MapperVertexMorphing::MapToGeometrySpace)
        ;

    class_<MapperVertexMorphingMatrixFree, bases<Process> >("MapperVertexMorphingMatrixFree", init<ModelPart&, Parameters&>())
        .def("MapToDesignSpace", &MapperVertexMorphingMatrixFree::MapToDesignSpace)
        .def("MapToGeometrySpace", &MapperVertexMorphingMatrixFree::MapToGeometrySpace)
        ;

    class_<MapperVertexMorphingImprovedIntegration, bases<Process> >("MapperVertexMorphingImprovedIntegration", init<ModelPart&, Parameters&>())
        .def("MapToDesignSpace", &MapperVertexMorphingImprovedIntegration::MapToDesignSpace)
        .def("MapToGeometrySpace", &MapperVertexMorphingImprovedIntegration::MapToGeometrySpace)
        ;

    // ================================================================
    // For a possible damping of nodal variables
    // ================================================================
    class_<DampingUtilities, bases<Process> >("DampingUtilities", init<ModelPart&, boost::python::dict, Parameters&>())
        .def("DampNodalVariable", &DampingUtilities::DampNodalVariable)
        ;

    // ========================================================================
    // For performing individual steps of an optimization algorithm
    // ========================================================================
    class_<OptimizationUtilities, bases<Process> >("OptimizationUtilities", init<ModelPart&, Parameters::Pointer>())
        // ----------------------------------------------------------------
        // For running unconstrained descent methods
        // ----------------------------------------------------------------
        .def("compute_search_direction_steepest_descent", &OptimizationUtilities::compute_search_direction_steepest_descent)
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
        .def("project_nodal_variable_on_unit_surface_normals", &GeometryUtilities::project_nodal_variable_on_unit_surface_normals)
        .def("update_coordinates_according_to_input_variable", &GeometryUtilities::update_coordinates_according_to_input_variable)
        .def("extract_surface_nodes", &GeometryUtilities::extract_surface_nodes)
        ;

    // ========================================================================
    // For calculations related to response functions
    // ========================================================================
    class_<StrainEnergyResponseFunction, bases<Process> >("StrainEnergyResponseFunction", init<ModelPart&, Parameters&>())
        .def("initialize", &StrainEnergyResponseFunction::initialize)
        .def("calculate_value", &StrainEnergyResponseFunction::calculate_value)
        .def("calculate_gradient", &StrainEnergyResponseFunction::calculate_gradient)
        .def("get_value", &StrainEnergyResponseFunction::get_value)
        .def("get_initial_value", &StrainEnergyResponseFunction::get_initial_value)
        .def("get_gradient", &StrainEnergyResponseFunction::get_gradient)
        ;
    class_<MassResponseFunction, bases<Process> >("MassResponseFunction", init<ModelPart&, Parameters&>())
        .def("initialize", &MassResponseFunction::initialize)
        .def("calculate_value", &MassResponseFunction::calculate_value)
        .def("calculate_gradient", &MassResponseFunction::calculate_gradient)
        .def("get_value", &MassResponseFunction::get_value)
        .def("get_initial_value", &MassResponseFunction::get_initial_value)
        .def("get_gradient", &MassResponseFunction::get_gradient)
        ;

    // ========================================================================
    // For input / output
    // ========================================================================
    class_<UniversalFileIO, bases<Process> >("UniversalFileIO", init<ModelPart&, Parameters&>())
        .def("initializeLogging", &UniversalFileIO::initializeLogging)
        .def("logNodalResults", &UniversalFileIO::logNodalResults)
        ;

    class_<VTKFileIO, bases<Process> >("VTKFileIO", init<ModelPart&, Parameters&>())
        .def("initializeLogging", &VTKFileIO::initializeLogging)
        .def("logNodalResults", &VTKFileIO::logNodalResults)
        ;

    // ========================================================================
    // For CAD reconstruction
    // ========================================================================
    class_<ReconstructionDataBase, bases<Process> >("ReconstructionDataBase", init<ModelPart&, boost::python::dict, boost::python::dict>())
        .def("Create", &ReconstructionDataBase::Create)
        ;
    class_<ReconstructionConditionContainer, bases<Process> >("ReconstructionConditionContainer", init<ReconstructionDataBase&>())
        .def("CreateDisplacementMappingConditions", &ReconstructionConditionContainer::CreateDisplacementMappingConditions)
        .def("CreateDisplacementCouplingConstraintsOnAllCouplingPoints", &ReconstructionConditionContainer::CreateDisplacementCouplingConstraintsOnAllCouplingPoints)    
        .def("CreateRotationCouplingConstraintsOnAllCouplingPoints", &ReconstructionConditionContainer::CreateRotationCouplingConstraintsOnAllCouplingPoints)    
        .def("CreateDirichletConditions", &ReconstructionConditionContainer::CreateDirichletConditions)
        .def("CreateMinimalControlPointDisplacementCondition", &ReconstructionConditionContainer::CreateMinimalControlPointDisplacementCondition)       
        ;           
    class_<CADReconstructionSolver, bases<Process> >("CADReconstructionSolver", init<ReconstructionDataBase&, ReconstructionConditionContainer&, CompressedLinearSolverType::Pointer>())
        .def("InitializeEquationSystem", &CADReconstructionSolver::InitializeEquationSystem) 
        .def("ComputeLHS", &CADReconstructionSolver::ComputeLHS) 
        .def("ComputeRHS", &CADReconstructionSolver::ComputeRHS) 
        .def("SolveEquationSystem", &CADReconstructionSolver::SolveEquationSystem) 
        .def("MultiplyAllPenaltyFactorsByInputFactor", &CADReconstructionSolver::MultiplyAllPenaltyFactorsByInputFactor)         
        .def("UpdateControlPointsAccordingReconstructionStrategy", &CADReconstructionSolver::UpdateControlPointsAccordingReconstructionStrategy)         
        ;   
    class_<ReconstructionOutputWriter, bases<Process> >("ReconstructionOutputWriter", init<ReconstructionDataBase&, std::string, std::string, std::string>())
        .def("OutputCADSurfacePoints", &ReconstructionOutputWriter::OutputCADSurfacePoints)
        .def("OutputGaussPointsOfFEMesh", &ReconstructionOutputWriter::OutputGaussPointsOfFEMesh)
        .def("OutputResultsInRhinoFormat", &ReconstructionOutputWriter::OutputResultsInRhinoFormat)        
        ;                            
    
    class_<QualityEvaluationUtility, bases<Process> >("QualityEvaluationUtility", init<ReconstructionDataBase&, ReconstructionConditionContainer&>())
        .def("EvaluateQualityOfProjection", &QualityEvaluationUtility::EvaluateQualityOfProjection)
        .def("EvaluateGlobalQuality", &QualityEvaluationUtility::EvaluateGlobalQuality)
        .def("EvaluateDisplacementCoupling", &QualityEvaluationUtility::EvaluateDisplacementCoupling)
        .def("EvaluateRotationCoupling", &QualityEvaluationUtility::EvaluateRotationCoupling)        
        ;                            

}

}  // namespace Python.

} // Namespace Kratos

