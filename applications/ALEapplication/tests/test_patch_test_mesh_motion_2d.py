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
        
    
    def _apply_BCs(self,mp):
        for node in mp.Nodes:
            node.Fix(KratosMultiphysics.MESH_DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.MESH_DISPLACEMENT_Y)
            node.Fix(KratosMultiphysics.MESH_DISPLACEMENT_Z)
        
        # for node in mp.Nodes:
        #     xvec = KratosMultiphysics.Vector(3)
        #     xvec[0] = node.X0
        #     xvec[1] = node.Y0
        #     xvec[2] = node.Z0
            
        #     u = KratosMultiphysics.Vector()
        #     u = A*xvec
        #     u += b
            
        #     node.SetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT,0,u)


            
    def _define_movement(self,dim,mp,time):

        center = KratosMultiphysics.Vector(2)
        center[0] = 0.75
        center[1] = 0.75

        disp = KratosMultiphysics.Vector(2)

        omega = time

        print (omega)

        for node in mp.Nodes:

            relative_coordinates = KratosMultiphysics.Vector(2)
            relative_coordinates[0] = node.X0 - center[0]
            relative_coordinates[1] = node.Y0 + center[1]

            disp[0] = math.cos(omega) * relative_coordinates[0] - math.sin(omega) * relative_coordinates[1]
            disp[1] = math.sin(omega) * relative_coordinates[0] - math.cos(omega) * relative_coordinates[1]    

            node.SetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_X,0, disp[0])
            node.SetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT_Y,0, disp[1])

        
    def _solve(self,mp):
        
        #define a minimal newton raphson solver
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        time_order = 2
        reform_dofs_each_step = False
        compute_reactions = False

        strategy = ALEApplication.StructuralMeshMovingStrategy(mp,
                                                             linear_solver,
                                                             time_order,
                                                             reform_dofs_each_step,
                                                             compute_reactions)
        strategy.SetEchoLevel(0)
        
        #strategy.Check()
        strategy.Solve()
        
    
    def _check_results(self,mp,A,b):
        
        ##check that the results are exact on the nodes
        for node in mp.Nodes:
            xvec = KratosMultiphysics.Vector(len(b))
            xvec[0] = node.X0
            xvec[1] = node.Y0
            xvec[2] = node.Z0
            
            u = KratosMultiphysics.Vector(2)
            u = A*xvec
            u += b            
            
            d = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
            self.assertAlmostEqual(d[0], u[0])
            self.assertAlmostEqual(d[1], u[1])
            self.assertAlmostEqual(d[2], u[2])

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

        return mp
            
    #def _check_outputs(self,mp,A,dim):

    def test_StructuralMeshMovingElement_2D_triangle(self):
        dim = 2
        mp = KratosMultiphysics.ModelPart("ale_part")
        self._add_variables(mp)
        
        #create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,0.5,0.0,0.0)
        mp.CreateNewNode(3,1.0,0.0,0.0)
        mp.CreateNewNode(4,1.5,0.0,0.0)

        mp.CreateNewNode(5,0.0,0.5,0.0)
        mp.CreateNewNode(6,0.5,0.5,0.0)
        mp.CreateNewNode(7,1.0,0.5,0.0)
        mp.CreateNewNode(8,1.5,0.5,0.0)

        mp.CreateNewNode(9,0.0,1.0,0.0)
        mp.CreateNewNode(10,0.5,1.0,0.0)
        mp.CreateNewNode(11,1.0,1.0,0.0)
        mp.CreateNewNode(12,1.5,1.0,0.0)

        mp.CreateNewNode(13,0.0,1.5,0.0)
        mp.CreateNewNode(14,0.5,1.5,0.0)
        mp.CreateNewNode(15,1.0,1.5,0.0)
        mp.CreateNewNode(16,1.5,1.5,0.0)
        


        
        for node in mp.Nodes:
            node.AddDof(KratosMultiphysics.MESH_DISPLACEMENT_X, KratosMultiphysics.MESH_REACTION_X)
            node.AddDof(KratosMultiphysics.MESH_DISPLACEMENT_Y, KratosMultiphysics.MESH_REACTION_Y)
            node.AddDof(KratosMultiphysics.MESH_DISPLACEMENT_Z, KratosMultiphysics.MESH_REACTION_Z)
            
        #create a submodelpart for boundary conditions
        bcs = mp.CreateSubModelPart("BoundaryCondtions")
        bcs.AddNodes([1,2,3,4,5,9,13,14,15,16,12,8])
        move_mp = mp.CreateSubModelPart("PrescribedMotion")
        #move_mp.AddNodes([6,7,10,11])
        move_mp.AddNodes([7,11])
                
        #create Element
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 1, [1,2,6], mp.GetProperties()[1])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 2, [2,3,7], mp.GetProperties()[1])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 3, [3,4,8], mp.GetProperties()[1])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 4, [1,5,6], mp.GetProperties()[1])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 5, [2,6,7], mp.GetProperties()[1])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 6, [3,7,8], mp.GetProperties()[1])

        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 7, [5,6,10], mp.GetProperties()[1])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 8, [6,7,11], mp.GetProperties()[1])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 9, [7,8,12], mp.GetProperties()[1])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 10, [5,9,10], mp.GetProperties()[1])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 11, [6,10,11], mp.GetProperties()[1])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 12, [7,11,12], mp.GetProperties()[1])

        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 13, [9,10,14], mp.GetProperties()[1])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 14, [10,11,15], mp.GetProperties()[1])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 15, [11,12,16], mp.GetProperties()[1])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 16, [9,13,14], mp.GetProperties()[1])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 17, [10,14,15], mp.GetProperties()[1])
        mp.CreateNewElement("StructuralMeshMovingElement2D3N", 18, [11,15,16], mp.GetProperties()[1])

        
        #A,b = self._define_movement(dim)
        
        #time integration parameters
        dt = 0.5
        time = 0.0
        end_time = 10
        step = 0
        
        mp = self._set_and_fill_buffer(mp,2,dt)
        # Applying boundary conditions
        self._apply_BCs(bcs)
        # Applying constraints
        #cm = KratosMultiphysics.StructuralMechanicsApplication.ApplyMultipointConstraintsProcess(mp)
        #mp, cm = self._apply_mpc_constraints(mp,cm)
        # Solving the system of equations        
        #self._setup_solver(mp)

        while(time <= end_time):
            time = time + dt
            step = step + 1
            print('############ Time :: ', time, ' ### step ', step)
            mp.CloneTimeStep(time)
            self._define_movement(2,move_mp,time)

            self._solve(mp)
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
