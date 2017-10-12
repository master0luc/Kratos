// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef RECONSTRUCTION_QUALITY_EVALUATION_UTILITY_H
#define RECONSTRUCTION_QUALITY_EVALUATION_UTILITY_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <string>
#include <cmath>
#include <iostream>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "../basic_nurbs_brep_handling/brep_gauss_point.h"
#include "../data_management/cad_projection_utility.h"

#include "reconstruction_data_base.h"

// #include "processes/process.h"
// #include "includes/kratos_flags.h"
// #include "../reconstruction_conditions/reconstruction_condition_displacement_mapping.h"
// #include "../reconstruction_conditions/reconstruction_constraint_displacement_coupling.h"
// #include "../reconstruction_conditions/reconstruction_constraint_rotation_coupling.h"
// #include "../reconstruction_conditions/regularization_condition_min_control_point_displacement.h"

// ==============================================================================

namespace Kratos
{
class ReconstructionQualityEvaluationUtility
{
public:
    ///@name Type Definitions
    ///@{
    typedef Node<3> NodeType;
    typedef std::vector<NodeType::Pointer> NodeVector;   
    typedef std::vector<Patch> PatchVector;
    typedef std::vector<BREPElement> BREPElementVector;
    typedef std::vector<BREPGaussPoint> BREPGaussPointVector;    
    typedef Element::GeometryType::IntegrationMethod IntegrationMethodType;
    
    /// Pointer definition of ReconstructionQualityEvaluationUtility
    KRATOS_CLASS_POINTER_DEFINITION(ReconstructionQualityEvaluationUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ReconstructionQualityEvaluationUtility( ReconstructionDataBase& reconstruction_data_base):
     mrReconstructionDataBase( reconstruction_data_base )
    {      
    }

    /// Destructor.
    virtual ~ReconstructionQualityEvaluationUtility()
    {
    }

    // --------------------------------------------------------------------------
    void EvaluateGlobalQuality(std::string reconstruction_strategy,
                               boost::python::list rParameterResolution,
                               int integration_degree,
                               int max_iterations,
                               double projection_tolerance)
    {
      std::cout << "EvaluateGlobalQuality CALLED" << std::endl;

      // compute distances from surface
      if(reconstruction_strategy.compare("mapping") == 0)
        ComputeDistancesFromGaussPoints(rParameterResolution, integration_degree, max_iterations, projection_tolerance );
      else
        KRATOS_THROW_ERROR(std::invalid_argument, "Reconstruction strategy specified to evaluate global quality is not recognized!","");

      // analyze the global quality of reconstruction
    }

    void ComputeDistancesFromGaussPoints( boost::python::list rParameterResolution, int integration_degree, int max_iterations, double projection_tolerance )
    {
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        PatchVector& patch_vector = mrReconstructionDataBase.GetPatchVector();
        ModelPart& fe_model_part = mrReconstructionDataBase.GetFEModelPart();                       // INITIALIZATION
        DoubleVector distancesVector;
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        IntegrationMethodType fem_integration_method;
        switch(integration_degree)
        {
            case 1 : fem_integration_method = GeometryData::GI_GAUSS_1; break;
            case 2 : fem_integration_method = GeometryData::GI_GAUSS_2; break;
            case 3 : fem_integration_method = GeometryData::GI_GAUSS_3; break;
            case 4 : fem_integration_method = GeometryData::GI_GAUSS_4; break;                       // INITIALIZE PROJECTOR 
            case 5 : fem_integration_method = GeometryData::GI_GAUSS_5; break;
        }
        
        CADProjectionUtility FE2CADProjector( patch_vector, max_iterations, projection_tolerance );
        FE2CADProjector.Initialize( rParameterResolution );
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        int points_counter=0;
        int outside_points_counter=0;

        for (auto & elem_i : fe_model_part.Elements())                                                // COMPUTE POINTS, EVALUATE DISTANCES
        {
            Element::GeometryType& geom_i = elem_i.GetGeometry();
            const Element::GeometryType::IntegrationPointsArrayType& integration_points = geom_i.IntegrationPoints(fem_integration_method);

            for (auto & integration_point_i : integration_points)
            {
                // int integration_point_number = &integration_point_i - &integration_points[0];
                NodeType::CoordinatesArrayType ip_coordinates = geom_i.GlobalCoordinates(ip_coordinates, integration_point_i.Coordinates());
                NodeType::Pointer node_of_interest = Node <3>::Pointer(new Node<3>(1, ip_coordinates));

                array_1d<double,2> parameter_values_of_nearest_point;
                int patch_index_of_nearest_point = -1;

                FE2CADProjector.DetermineNearestCADPoint( node_of_interest, 
                                                          parameter_values_of_nearest_point, 
                                                          patch_index_of_nearest_point );

                Patch patch_of_nearest_point = patch_vector[patch_index_of_nearest_point];
                bool is_inside = patch_of_nearest_point.IsPointInside(parameter_values_of_nearest_point);
                if(is_inside)
                {
                  Point<3> nearest_point_coordinates;
                  patch_of_nearest_point.EvaluateSurfacePoint(parameter_values_of_nearest_point, nearest_point_coordinates);
                  NodeType::Pointer nearest_cad_node = Node <3>::Pointer(new Node<3>(1, nearest_point_coordinates));

                  double distance;
                  EvaluateDistanceBetweenNodes(node_of_interest, nearest_cad_node, distance);
                  distancesVector.push_back(distance);
                }
                else
                  outside_points_counter++;
                points_counter++;

            }
        }
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // norm, max, average, (% of beyond-trimming-boundaries points)                       // STATISTICS
        
        
    }
    // --------------------------------------------------------------------------
    void EvaluateDisplacementCoupling()
    {
      std::cout << "EvaluateDisplacementCoupling CALLED" << std::endl;
    }

    // --------------------------------------------------------------------------
    void EvaluateRotationCoupling()
    {
      std::cout << "EvaluateRotationCoupling CALLED" << std::endl;
    }

    // --------------------------------------------------------------------------
    void EvaluateDistanceBetweenNodes(NodeType::Pointer first, NodeType::Pointer second, double& distance)
    {
      double dx = first->X() - second->X();
      double dy = first->Y() - second->Y();
      double dz = first->Z() - second->Z();
      distance = std::sqrt( std::pow(dx,2) + std::pow(dy,2) + std::pow(dz,2) );
    }
    
    // ==============================================================================

    /// Turn back information as a string.
    virtual std::string Info() const
    {
		return "ReconstructionQualityEvaluationUtility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
		rOStream << "ReconstructionQualityEvaluationUtility";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }

private:

    ReconstructionDataBase& mrReconstructionDataBase;

    /// Assignment operator.
    //      ReconstructionQualityEvaluationUtility& operator=(ReconstructionQualityEvaluationUtility const& rOther);

    /// Copy constructor.
    //      ReconstructionQualityEvaluationUtility(ReconstructionQualityEvaluationUtility const& rOther);

}; // Class ReconstructionQualityEvaluationUtility
} // namespace Kratos.

#endif // RECONSTRUCTION_QUALITY_EVALUATION_UTILITY_H
