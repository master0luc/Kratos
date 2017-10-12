# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division 

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Additional imports
import json as json
import time
import os

# ======================================================================================================================================
# Class definition
# ======================================================================================================================================

class CADReconstrutionUtilities():

    # --------------------------------------------------------------------------
    def __init__( self, fem_filename, cad_geometry_filename, cad_integration_data_filename ):
        self.FEMInputFilename = fem_filename
        self.CADGeometryFilename = cad_geometry_filename
        self.CADIntegrationDataFilename = cad_integration_data_filename

        # Internal parameters to specify reconstruction method
        
        # Gernal strategy parameters
        self.ReconstructionStrategy = "mapping"    # mapping / least_squares
        self.FEMGaussIntegrationDegree = 5

        # Solution parameters
        self.SolutionIterations = 1
        self.PenaltyFactorForDisplacementCoupling = 1e3
        self.PenaltyFactorForRotationCoupling = 1e3
        self.PenaltyFactorForDirichletConstraints = 1e3
        self.PenaltyMultiplier = 1.0
       
        # Parameters to edit input data       
        self.FERefinementLevel = 1

        # Projection settings
        self.ParameterResolutionForProjection = [ 100, 100 ]
        self.MaxProjectionIterations = 20
        self.ProjectionTolerance = 1e-5

        # Regularization settings
        self.MinimizeControlPointDisplacement = True
        self.Beta = 0.002

        # Linear solver
        self.LinearSolver = SuperLUSolver()
        # DiagPrecond = DiagonalPreconditioner()
        # self.LinearSolver =  BICGSTABSolver(1e-9, 5000, DiagPrecond)
        # self.LinearSolver = AMGCLSolver(AMGCLSmoother.GAUSS_SEIDEL, AMGCLIterativeSolverType.BICGSTAB, 1e-9, 300, 2, 10)        

        # Parameters for output
        self.OutputFolder = "01_Results"
        self.ParameterResolutionForOutputOfSurfacePoints = [ 50, 50 ]
        self.GeorhinoFilename = "tripod.georhino.txt"
        self.RhinoResultsFilename = "tripod.post.res"

    # --------------------------------------------------------------------------
    def Initialize( self ):
        self.__ReadFEData()
        self.__RefineFEModel()
        self.__ReadCADData()
        self.__InitializeConstraints()
        self.__CreateReconstructionDataBase()
        self.__CreateReconstructionOutputWriter()
    # --------------------------------------------------------------------------
    def SetDisplacementCouplingOnAllCouplingPoints( self ):
        self.IsDisplacementCouplingSpecifiedForAllCouplingPoints = True

    # --------------------------------------------------------------------------
    def SetRotationCouplingOnAllCouplingPoints( self ):
        self.IsRotationCouplingSpecifiedForAllCouplingPoints = True
        
    # --------------------------------------------------------------------------
    def SetDirichletConstraints( self, list_of_condition_settings ):
        self.AreDirichletConstraintsSpecified = True
        self.DirichletConditions = list_of_condition_settings

    # --------------------------------------------------------------------------
    def PerformReconstruction( self ):
        self.__CreateReconstructionConditions()
        self.__CreateSolverForReconstruction()
        self.__RunSolutionAlorithm()
    
    # --------------------------------------------------------------------------
    def EvaluateReconstructionQuality( self ):
        quality_evaluator = ReconstructionQualityEvaluationUtility( self.DataBase )
        # evaluate global fitting
        # measure distances vector: very fine point cloud (discard wrong patch)
        #   - norm
        #   - max
        #   - mean
        quality_evaluator.EvaluateGlobalQuality(self.ReconstructionStrategy,
                                                self.ParameterResolutionForProjection,
                                                self.FEMGaussIntegrationDegree,
                                                self.MaxProjectionIterations,
                                                self.ProjectionTolerance)
        
        # evaluate G0
        #   - integral mean over length
        #   - max pointwise
        #   - mean pointwise
        quality_evaluator.EvaluateDisplacementCoupling()
        
        # evaluate G1
        #   - integral mean over length
        #   - max pointwise
        #   - mean pointwise
        quality_evaluator.EvaluateRotationCoupling()

    # --------------------------------------------------------------------------
    def OutputFEData( self ):
        from gid_output import GiDOutput
        fem_output_filename = self.OutputFolder + "/" + self.FEMInputFilename + "_as_used_for_reconstruction"
        nodal_results=["SHAPE_CHANGE_ABSOLUTE"]
        gauss_points_results=[]
        VolumeOutput = True
        GiDPostMode = "Binary"
        GiDWriteMeshFlag = False
        GiDWriteConditionsFlag = True
        GiDWriteParticlesFlag = False
        GiDMultiFileFlag = "Single"
        gig_io = GiDOutput(fem_output_filename, VolumeOutput, GiDPostMode, GiDMultiFileFlag, GiDWriteMeshFlag, GiDWriteConditionsFlag)
        gig_io.initialize_results(self.FEModelPart)
        gig_io.write_results(1, self.FEModelPart, nodal_results, gauss_points_results)
        gig_io.finalize_results()

    # --------------------------------------------------------------------------
    def OutputCADSurfacePoints( self, file_to_write ):
        self.OutputWriter.OutputCADSurfacePoints( file_to_write, self.ParameterResolutionForOutputOfSurfacePoints )

    # --------------------------------------------------------------------------
    def OutputGaussPointsOfFEMesh( self, file_to_write ):
        self.OutputWriter.OutputGaussPointsOfFEMesh( file_to_write, self.FEMGaussIntegrationDegree )

    # --------------------------------------------------------------------------
    def __ReadFEData( self ):
        print("\n> Start importing FE data")
        self.FEModelPart = ModelPart("name_of_empty_mdpa")
        self.FEModelPart.AddNodalSolutionStepVariable(SHAPE_CHANGE_ABSOLUTE)
        model_part_io = ModelPartIO(self.FEMInputFilename)
        model_part_io.ReadModelPart(self.FEModelPart)
        print("> Importing FE data finished.")        

    # --------------------------------------------------------------------------
    def __RefineFEModel( self ):       

        # Assign pseudo material to elements (required by refinement)
        prop_id = 1
        prop = self.FEModelPart.Properties[prop_id]
        mat = LinearElasticPlaneStress2DLaw()
        prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())

        for refinement_level in range(0,self.FERefinementLevel):

            number_of_avg_elems = 10
            number_of_avg_nodes = 10
            nodal_neighbour_search = FindNodalNeighboursProcess(self.FEModelPart, number_of_avg_elems, number_of_avg_nodes)
            neighbour_calculator = FindElementalNeighboursProcess(self.FEModelPart,2,10)
            nodal_neighbour_search.Execute()
            neighbour_calculator.Execute()

            for elem in self.FEModelPart.Elements:
                elem.SetValue(SPLIT_ELEMENT,True)

            refine_on_reference = False
            interpolate_internal_variables = True
            Refine = LocalRefineTriangleMesh(self.FEModelPart)
            Refine.LocalRefineMesh(refine_on_reference, interpolate_internal_variables)

    # --------------------------------------------------------------------------
    def __ReadCADData( self ):
        print("\n> Start importing CAD data")
        self.CADGeometry = {}
        with open(self.CADGeometryFilename) as cad_data1:
            self.CADGeometry = json.load(cad_data1)
        self.CADIntegrationData = {}
        with open(self.CADIntegrationDataFilename) as cad_data2:
            self.CADIntegrationData = json.load(cad_data2)
        print("> Importing CAD data finished.")

    # --------------------------------------------------------------------------
    def __InitializeConstraints( self ):
        self.IsDisplacementCouplingSpecifiedForAllCouplingPoints = False
        self.IsRotationCouplingSpecifiedForAllCouplingPoints = False
        self.AreDirichletConstraintsSpecified = False
    
    # --------------------------------------------------------------------------
    def __CreateReconstructionDataBase( self ):
        self.DataBase = ReconstructionDataBase(self.FEModelPart, self.CADGeometry, self.CADIntegrationData)
        self.DataBase.Create()

    # --------------------------------------------------------------------------
    def __CreateReconstructionOutputWriter( self ):
        if not os.path.exists( self.OutputFolder ):
            os.makedirs( self.OutputFolder )    
        self.OutputWriter = ReconstructionOutputWriter( self.DataBase, self.OutputFolder, self.RhinoResultsFilename, self.GeorhinoFilename )    

    # --------------------------------------------------------------------------
    def __CreateReconstructionConditions( self ):

        # Container to store all conditions (including constraints and reguarlization )
        self.ConditionsContainer = ReconstructionConditionContainer( self.DataBase )

        # Basic reconstruction condition
        if self.ReconstructionStrategy == "mapping":
            self.ConditionsContainer.CreateDisplacementMappingConditions( self.ParameterResolutionForProjection, 
                                                                          self.FEMGaussIntegrationDegree,
                                                                          self.MaxProjectionIterations,
                                                                          self.ProjectionTolerance )
        else:
            raise ValueError( "The following reconstruction strategy does not exist: ", self.ReconstructionStrategy )

        # Reconstruction constraints
        if self.IsDisplacementCouplingSpecifiedForAllCouplingPoints: 
            self.ConditionsContainer.CreateDisplacementCouplingConstraintsOnAllCouplingPoints( self.PenaltyFactorForDisplacementCoupling )
        if self.IsRotationCouplingSpecifiedForAllCouplingPoints: 
            self.ConditionsContainer.CreateRotationCouplingConstraintsOnAllCouplingPoints( self.PenaltyFactorForRotationCoupling )            
        if self.AreDirichletConstraintsSpecified:
            self.ConditionsContainer.CreateDirichletConstraints( self.DirichletConstraints, self.PenaltyFactorForDirichletConstraints )

        # Regularization
        if self.MinimizeControlPointDisplacement: 
            self.ConditionsContainer.CreateMinimalControlPointDisplacementCondition( self.DataBase, self.Beta, self.ReconstructionStrategy )                   

    # --------------------------------------------------------------------------
    def __CreateSolverForReconstruction( self ):
        self.ReconstructionSolver = CADReconstructionSolver( self.DataBase, self.ConditionsContainer, self.LinearSolver )

    # --------------------------------------------------------------------------
    def __RunSolutionAlorithm( self ): 
        
        self.ReconstructionSolver.InitializeEquationSystem()

        for iteration in range(1,self.SolutionIterations+1): 
            print("\n===========================================")
            print("Starting reconstruction iteration ", iteration,"...")
            print("===========================================")            

            if iteration > 1:  
                self.ReconstructionSolver.MultiplyAllPenaltyFactorsByInputFactor( self.PenaltyMultiplier )               

            self.ReconstructionSolver.ComputeLHS()
            self.ReconstructionSolver.ComputeRHS()
            self.ReconstructionSolver.SolveEquationSystem()

            # self.ReconstructionSolver.UpdateControlPointsAccordingReconstructionStrategy( self.ReconstructionStrategy )
            
            self.OutputWriter.OutputResultsInRhinoFormat( iteration ) 

        print("\n===========================================")
        print("Finished reconstruction loop.")
        print("===========================================")                    

# ======================================================================================================================================
# Reconstruction
# ======================================================================================================================================    

print("\n\n========================================================================================================")
print("> Start reconstruction...")
print("========================================================================================================")

# Measure time
start_time = time.time()

# Input parameters
fem_filename = "tripod"
cad_geometry_filename = "tripod_geometry.json" 
cad_integration_data_filename = "tripod_integration_data.json"

# Initialize Reconstruction
CADReconstructionUtility = CADReconstrutionUtilities( fem_filename, cad_geometry_filename, cad_integration_data_filename )
CADReconstructionUtility.Initialize()

# Set Boundary Conditions
CADReconstructionUtility.SetDisplacementCouplingOnAllCouplingPoints()
CADReconstructionUtility.SetRotationCouplingOnAllCouplingPoints()

# Some output before reconstruction
CADReconstructionUtility.OutputCADSurfacePoints( "surface_points_of_cad_geometry.txt" )

# Perform reconstruction
CADReconstructionUtility.PerformReconstruction()

# Measure quality of reconstruction
CADReconstructionUtility.EvaluateReconstructionQuality()

# Some output
CADReconstructionUtility.OutputFEData()
CADReconstructionUtility.OutputCADSurfacePoints( "surface_points_of_updated_cad_geometry.txt" )
CADReconstructionUtility.OutputGaussPointsOfFEMesh( "gauss_points_of_fe_mesh.txt" )

print("\n========================================================================================================")
print("> Finished reconstruction in " ,round( time.time()-start_time, 3 ), " s.")
print("========================================================================================================")