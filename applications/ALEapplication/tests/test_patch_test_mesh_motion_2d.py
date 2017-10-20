from __future__ import print_function, absolute_import, division
import KratosMultiphysics 

import KratosMultiphysics.ALEApplication as ALEApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

import math

class TestPatchTestALEMeshMotion(KratosUnittest.TestCase):
    def setUp(self):
        pass
    
    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY) 
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_REACTION) 

        #return mp


    def _add_dofs(self,mp):
        for node in mp.Nodes:
            node.AddDof(KratosMultiphysics.MESH_DISPLACEMENT_X, KratosMultiphysics.MESH_REACTION_X)
            node.AddDof(KratosMultiphysics.MESH_DISPLACEMENT_Y, KratosMultiphysics.MESH_REACTION_Y)
            node.AddDof(KratosMultiphysics.MESH_DISPLACEMENT_Z, KratosMultiphysics.MESH_REACTION_Z)

        #return mp
        
    
    def _apply_BCs(self,mp):
        bcs = mp.GetSubModelPart("BoundaryCondtions")
        for node in bcs.Nodes:
            node.Fix(KratosMultiphysics.MESH_DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.MESH_DISPLACEMENT_Y)
            node.Fix(KratosMultiphysics.MESH_DISPLACEMENT_Z)

        #return mp        
            
    def _define_movement(self,dim,mp,time):
        bcmn = mp.GetSubModelPart("PrescribedMotion")
        center = KratosMultiphysics.Vector(2)
        center[0] = 0.75
        center[1] = 0.75

        disp = KratosMultiphysics.Vector(2)

        omega = math.pi * time 

        for node in bcmn.Nodes:
            relative_coordinates = KratosMultiphysics.Vector(2)
            relative_coordinates[0] = node.X0 - center[0]
            relative_coordinates[1] = node.Y0 - center[1]

            disp[0] = math.cos(omega) * relative_coordinates[0] - math.sin(omega) * relative_coordinates[1]
            disp[1] = math.sin(omega) * relative_coordinates[0] + math.cos(omega) * relative_coordinates[1]    

            node.SetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_X,0, disp[0])
            node.SetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_Y,0, disp[1])

            #node.SetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_Y,0, time)
            node.Fix(KratosMultiphysics.MESH_DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.MESH_DISPLACEMENT_Y)


        # #define the applied motion - the idea is that the displacement is defined as u = A*xnode + b
        # #so that the displcement is linear and the exact F = I + A
        # A = KratosMultiphysics.Matrix(3,3)
        # A[0,0] = time*10e6;  A[0,1] = time*10e6; A[0,2] = 0.0
        # A[1,0] = time*10e6;  A[1,1] = time*10e6; A[1,2] = 0.0
        # A[2,1] = 0.0;  A[2,1] = 0.0; A[2,2] = 0.0
                    
        # b = KratosMultiphysics.Vector(3)
        # b[0] = 5000
        # b[1] = 5000  
        # b[2] = 0.0

        # for node in bcmn.Nodes:
        #     xvec = KratosMultiphysics.Vector(3)
        #     xvec[0] = node.X0
        #     xvec[1] = node.Y0
        #     xvec[2] = node.Z0
            
        #     u = KratosMultiphysics.Vector()
        #     u = A*xvec
        #     u += b
            
        #     node.SetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT,0,u)
        
    def _setup_solver(self,mp):
        
        #define a minimal newton raphson solver
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        time_order = 2
        reform_dofs_each_step = False
        compute_reactions = False

        self.strategy = ALEApplication.StructuralMeshMovingStrategy(mp,
                                                             linear_solver,
                                                             time_order,
                                                             reform_dofs_each_step,
                                                             compute_reactions)
        self.strategy.SetEchoLevel(0)
        
        self.strategy.Check()

        return mp
        

    def _solve(self):
         self.strategy.Solve()
        

    def _set_and_fill_buffer(self,mp,buffer_size,delta_time):
        # Set buffer size
        mp.SetBufferSize(buffer_size)
        # Fill buffer
        time = mp.ProcessInfo[KratosMultiphysics.TIME]
        time = time - delta_time * (buffer_size)
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        for size in range(0, buffer_size):
            step = size - (buffer_size -1)
            mp.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            time = time + delta_time
            #delta_time is computed from previous time in process_info
            mp.CloneTimeStep(time)

        mp.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False     


    def _setup_model_part(self,mp):
            
        #create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,0.5,0.0,0.0)
        mp.CreateNewNode(3,1,0.0,0.0)
        mp.CreateNewNode(4,1.5,0.0,0.0)

        mp.CreateNewNode(5,0.0,0.5,0.0)
        mp.CreateNewNode(6,0.5,0.5,0.0)
        mp.CreateNewNode(7,0.75,0.5,0.0)
        mp.CreateNewNode(8,1,0.5,0.0)
        mp.CreateNewNode(9,1.5,0.5,0.0)

        mp.CreateNewNode(10,0.625,0.625,0.0)

        mp.CreateNewNode(11,0.5,0.75,0.0)
        mp.CreateNewNode(12,0.75,0.75,0.0)

        mp.CreateNewNode(13,0.625,0.875,0.0)

        mp.CreateNewNode(14,0.0,1,0.0)
        mp.CreateNewNode(15,0.5,1,0.0)
        mp.CreateNewNode(16,0.75,1.0,0.0)
        mp.CreateNewNode(17,1,1,0.0)
        mp.CreateNewNode(18,1.5,1,0.0)

        mp.CreateNewNode(19,0.0,1.5,0.0)
        mp.CreateNewNode(20,0.5,1.5,0.0)
        mp.CreateNewNode(21,1,1.5,0.0)
        mp.CreateNewNode(22,1.5,1.5,0.0)
        


        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([1,2,3,4,5,14,19,20,21,22,18,9])
        move_mp = mp.CreateSubModelPart("PrescribedMotion")
        move_mp.AddNodes([15,13,11])
            
                
        #create Elements
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 1, [1,2,6], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 2, [2,3,7], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 3, [3,4,8], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 4, [1,6,5], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 5, [2,7,6], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 6, [7,3,8], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 7, [4,9,8], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 8, [5,6,11], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 9, [6,10,11], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 10, [6,7,10], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 11, [7,12,10], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 12, [7,8,12], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 13, [8,9,18], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 14, [5,11,14], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 15, [11,10,12], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 16, [12,8,17], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 17, [8,18,17], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 18, [14,11,15], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 19, [11,13,15], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 20, [11,12,13], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 21, [13,12,16], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 22, [12,17,16], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 23, [15,13,16], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 24, [14,15,20], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 25, [15,16,20], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 26, [16,17,21], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 27, [17,18,21], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 28, [14,20,19], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 29, [20,16,21], mp.GetProperties()[0])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 30, [21,18,22], mp.GetProperties()[0])

        return mp


    def _check_results(self,mp):
            node_6
            disp2 = node.GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_Y,0)





    def test_StructuralMeshMovingElement_2D_triangle(self):
        
        mp = KratosMultiphysics.ModelPart("ale_part")
        self._add_variables(mp)  
        mp = self._setup_model_part(mp)
        self._add_dofs(mp)
        mp = self._setup_solver(mp)
        
        #A,b = self._define_movement(dim)
        
        #time integration parameters
        dt = 10
        time = 0.0
        end_time = 50
        step = 0
        
        self._set_and_fill_buffer(mp,2,dt)

        # Applying boundary conditions
        self._apply_BCs(mp)


        # Solving the system of equations       
        while(time <= end_time):
            time = time + dt
            step = step + 1
            print('############ Time :: ', time, ' ### step ', step)
            mp.CloneTimeStep(time)   
                
            self._define_movement(2,mp,time)

            self._solve()

        # Checking the results
        #self._check_results(mp)
        #self._reset()
        #self._check_results(mp,A,b)
        #self._check_outputs(mp,A,dim)
        
                    
            self.__post_process(mp)
                    
   
   
    def __post_process(self, main_model_part):
        from gid_output_process import GiDOutputProcess
        self.gid_output = GiDOutputProcess(main_model_part,
                                    "gid_output",
                                    KratosMultiphysics.Parameters("""
                                        {
                                            "result_file_configuration" : {
                                                "gidpost_flags": {
                                                    "GiDPostMode": "GiD_PostBinary",
                                                    "WriteDeformedMeshFlag": "WriteUndeformed",
                                                    "WriteConditionsFlag": "WriteConditions",
                                                    "MultiFileFlag": "SingleFile"
                                                },        
                                                "file_label"          : "step",
                                                "output_control_type" : "step",
                                                "output_frequency"    : 1,
                                                "nodal_results"       : ["MESH_DISPLACEMENT", "MESH_VELOCITY"],
                                                "gauss_point_results" : []
                                            }
                                        }
                                        """)
                                    )

        self.gid_output.ExecuteInitialize()
        self.gid_output.ExecuteBeforeSolutionLoop()
        self.gid_output.ExecuteInitializeSolutionStep()
        self.gid_output.PrintOutput()
        self.gid_output.ExecuteFinalizeSolutionStep()
        self.gid_output.ExecuteFinalize()
        
if __name__ == '__main__':
    KratosUnittest.main()
