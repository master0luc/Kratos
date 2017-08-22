# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division 
from math import sin, cos, pi, fabs, sqrt
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
fem_input_filename = "test_1"
# cad_geometry_input_filename = "test_1_geometry.json" 
cad_geometry_input_filename = "test_1_geometry.json" 
cad_integration_input_filename = "Benchmark_halbkreis_32x16_multipatch_integration_data.json" 

# Output parameters
cad_geometry_output_filename = cad_geometry_input_filename
cad_geometry_output_filename = cad_geometry_output_filename.replace(".json","_updated.json")

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

# ======================================================================================================================================
# surface definition
# ======================================================================================================================================    
def surf_1(u,v):
    x_1 = sin(pi/4) *  u - cos(pi/4) * v
    return [u,v,0.5*x_1**2]    
def surf_2(u,v):
    x_1 = sin(pi/4) *  u - cos(pi/4) * v
    if x_1>0:
        return [u,v, -0.5*x_1**2]
    else:
        return [u,v, +0.5*x_1**2]
def surf_3(u,v):
    x_1 = sin(pi/4) *  u - cos(pi/4) * v
    return [u,v,-0.5*x_1**3]
def surf_4(u,v):
    x_1 = sin(pi/4) *  u - cos(pi/4) * v
    return [u,v,fabs(x_1)]
def surf_5(u,v):
    return [u,v,u**2]
def surf_6(u,v):
    return [u,v,0]
def surf_7(u,v):
    if u < 2.5:
        return [u,v,0.5*(u-2.5)**2]
    else:
        return [u,v,-0.5*(u-2.5)**2]
        
# ======================================================================================================================================
# mesh definition
# ======================================================================================================================================    
def mesh_1(surf): # "mesh of 4 nodes"           CORNERS           
    for u, v in [[0,0], [0,5], [5,0], [5,5]]:
        [x,y,z] = surf(u,v)
        mapper.set_point(0, u, v, x, y, z)
def mesh_2(surf, n): # "mesh of n nodes"        MAIN DIAGONAL + CORNERS
    mesh_1(surf) # 4 nodes
    n = n-3
    for i in range(1,n): # n-1 more nodes
        u = i*5/n
        v = u
        [x,y,z] = surf(u,v)
        mapper.set_point(0, u, v, x, y, z)
def mesh_3(surf, n): # "mesh of n nodes"        SECOND DIAGONAL + CORNERS
    mesh_1(surf)
    n = n-3
    for i in range(1,n):
        u = i*5/n
        v = 5 - u
        [x,y,z] = surf(u,v)
        mapper.set_point(0, u, v, x, y, z)
def mesh_4(surf,n): # "mesh of n nodes"         REGULAR GRID n x n
    n = int(sqrt(n))
    for i in range(n):
        for j in range(n):
            u = i * 5/(n-1)
            v = j * 5/(n-1)
            [x,y,z] = surf(u,v)
            mapper.set_point(0, u, v, x, y, z)
def mesh_5(surf): # "mesh of 8 nodes"           CORNERS + 4 INNER CORNERS
    mesh_1(surf)
    for u, v in [[1.5,1.5], [1.5,3.5], [3.5,1.5], [3.5,3.5]]:
        [x,y,z] = surf(u,v)
        mapper.set_point(0, u, v, x, y, z)
def mesh_6(surf,n,m): # "mesh of n x m nodes"   REGULAR GRID n x m
    for i in range(n):
        for j in range(m):
            u = i * 5/(n-1)
            v = j * 5/(m-1)
            [x,y,z] = surf(u,v)
            mapper.set_point(0, u, v, x, y, z)
def mesh_7(surf,n): #                           EDGES' MIDDLE POINTS + SECOND DIAGONAL + CORNERS
    mesh_3(surf, n-4)
    for u, v in [[2.5,0], [0,2.5], [2.5,5], [5,2.5]]:
        [x,y,z] = surf(u,v)
        mapper.set_point(0, u, v, x, y, z)
def mesh_8(surf,n): #                           REGULAR MESH ON EDGES
    n = int(n/4)
    for i in range(n):
        u = i * 5/n
        v = 0
        [x,y,z] = surf(u,v)
        mapper.set_point(0, u, v, x, y, z)
    for i in range(n):
        u = 5
        v = i * 5/n
        [x,y,z] = surf(u,v)
        mapper.set_point(0, u, v, x, y, z)
    for i in range(n):
        u = 5 - i * 5/n
        v = 5
        [x,y,z] = surf(u,v)
        mapper.set_point(0, u, v, x, y, z)
    for i in range(n):
        u = 0
        v = 5 - i * 5/n
        [x,y,z] = surf(u,v)
        mapper.set_point(0, u, v, x, y, z)
def mesh_9(surf,n): # "mesh of 4 x r x r"       4 INNER REGULAR MESHES
    r = int(sqrt(n/4))
    for i in range(r):
        for j in range(r):
            u = (i+1) * 5/(2*(r+1))
            v = (j+1) * 5/(2*(r+1))
            [x,y,z] = surf(u,v)
            mapper.set_point(0, u, v, x, y, z)
            u += 2.5
            [x,y,z] = surf(u,v)
            mapper.set_point(0, u, v, x, y, z)
            v += 2.5
            [x,y,z] = surf(u,v)
            mapper.set_point(0, u, v, x, y, z)
            u -= 2.5
            [x,y,z] = surf(u,v)
            mapper.set_point(0, u, v, x, y, z)
def mesh_10(surf,n):
    n = n/4
    n = int(sqrt(n))
    for k in range(4):
        for i in range(n):
            for j in range(n):
                u = i * 2.5/(n-1)
                v = j * 2.5/(n-1)

                if k == 1 or k == 3:
                    u_abs = u + 2.5
                else:
                    u_abs = u
                if k == 0 or k == 1:
                    v_abs = v + 2.5
                else:
                    v_abs = v
                [x,y,z] = surf(u_abs,v_abs)
                mapper.set_point(k, u, v, x, y, z)
# ======================================================================================================================================
# Mapping
# ======================================================================================================================================    

# Create CAD-mapper
linear_solver = SuperLUSolver()
mapper = CADMapper(fe_model_part,cad_geometry,cad_integration_data,linear_solver)

# Set nearest point + shape update to map
surf = surf_7
mesh = mesh_4
n = 16
m = 3

try:
    mesh(surf)
except TypeError:
    try:
        mesh(surf,n)
    except TypeError:
        mesh(surf,n,m)

# # Perform mapping
mapper.external_map_to_cad_space()

print("Objective =", mapper.compute_objective())
# ======================================================================================================================================
# Writing results
# ======================================================================================================================================

# Output some surface nodes of updated cad geometry
file_to_write = "surface_nodes_of_updated_cad_geometry.txt"
u_resolution = 30
v_resolution = 30
mapper.output_surface_points(file_to_write, u_resolution, v_resolution, -1)

# Output control point update in gid-format 
mapper.output_control_point_displacements()

# Output json file with updated geometry
with open(cad_geometry_output_filename, 'w') as fp:
    json.dump(cad_geometry, fp)
