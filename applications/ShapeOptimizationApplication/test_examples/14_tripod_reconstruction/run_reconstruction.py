# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division 

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Additional imports
import json as json

# ======================================================================================================================================
# Parameters
# ======================================================================================================================================

# Input parameters
fem_input_filename = "tripod"
cad_geometry_input_filename = "tripod_geometry.json" 
cad_integration_input_filename = "tripod_integration_data.json" 

# Output parameters
cad_geometry_output_filename = cad_geometry_input_filename.replace(".json","_updated.json")
rhino_result_file = "tripod.post.res"

# ======================================================================================================================================
# Reading
# ======================================================================================================================================

# Read the FE model
fe_model_part = ModelPart("name_of_empty_mdpa")
fe_model_part.AddNodalSolutionStepVariable(SHAPE_CHANGE_ABSOLUTE)
model_part_io = ModelPartIO(fem_input_filename)
model_part_io.ReadModelPart(fe_model_part)

# Read CAD data
cad_geometry = {}
with open(cad_geometry_input_filename) as cad_data1:
    cad_geometry = json.load(cad_data1)

cad_integration_data = {}
with open(cad_integration_input_filename) as cad_data2:
    cad_integration_data = json.load(cad_data2)    

# Output read fe model part
from gid_output import GiDOutput
fem_output_filename = fem_input_filename+"_as_seen_by_fe_solver"
nodal_results=["SHAPE_CHANGE_ABSOLUTE"]
gauss_points_results=[]
VolumeOutput = True
GiDPostMode = "Binary"
GiDWriteMeshFlag = False
GiDWriteConditionsFlag = True
GiDWriteParticlesFlag = False
GiDMultiFileFlag = "Single"
gig_io = GiDOutput(fem_output_filename, VolumeOutput, GiDPostMode, GiDMultiFileFlag, GiDWriteMeshFlag, GiDWriteConditionsFlag)
gig_io.initialize_results(fe_model_part)
gig_io.write_results(1, fe_model_part, nodal_results, gauss_points_results)
gig_io.finalize_results()

# ======================================================================================================================================
# Mapping
# ======================================================================================================================================    

# Create CAD-mapper
linear_solver = SuperLUSolver()
# DiagPrecond = DiagonalPreconditioner()
# linear_solver =  BICGSTABSolver(1e-9, 5000, DiagPrecond)
# linear_solver = AMGCLSolver(AMGCLSmoother.GAUSS_SEIDEL, AMGCLIterativeSolverType.BICGSTAB, 1e-9, 300, 2, 10)
mapper = CADMapper(fe_model_part,cad_geometry,cad_integration_data,linear_solver)

# Compute mapping matrix
u_resolution = 300
v_resolution = 300
mapper.compute_mapping_matrix(u_resolution,v_resolution)

# Apply boundary conditions
penalty_factor_displacement_coupling = 1e3
penalty_factor_rotation_coupling = 1e3
penalty_factor_dirichlet_condition = 1e3
edges_with_specific_dirichlet_conditions = [ ]
edges_with_enforced_tangent_continuity = [ ]
mapper.apply_boundary_conditions( penalty_factor_displacement_coupling, 
                                  penalty_factor_rotation_coupling, 
                                  penalty_factor_dirichlet_condition,
                                  edges_with_specific_dirichlet_conditions,
                                  edges_with_enforced_tangent_continuity )

# Perform mapping
mapper.map_to_cad_space()

# Output control point update as result file to be used in Rhino
mapper.output_control_point_displacements(rhino_result_file)

# Output json file with updated geometry
with open(cad_geometry_output_filename, 'w') as fp:
    json.dump(cad_geometry, fp)