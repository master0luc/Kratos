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
      DoubleVector distances;
      if(reconstruction_strategy.compare("mapping") == 0)
        ComputeDistancesFromGaussPoints(rParameterResolution, integration_degree, max_iterations, projection_tolerance, distances );
      else
        KRATOS_THROW_ERROR(std::invalid_argument, "Reconstruction strategy specified to evaluate global quality is not recognized!","");

      // analyze the global quality of reconstruction
        Vector distancesVector = ZeroVector(distances.size());
        for (unsigned int i=0; i<distances.size(); ++i)
        {
          distancesVector(i) = distances[i];
        }
        std::cout << "> Max value = " << norm_inf(distancesVector) << std::endl;
        std::cout << "> L2 norm = " << norm_2(distancesVector) << std::endl;      
        std::cout << "> average value = " << sum(distancesVector)/distancesVector.size() << std::endl;      
        std::cout << "> points ignored: " << 100 * (double) outside_points_counter / points_counter << "%" << std::endl;      

    }

    // --------------------------------------------------------------------------
    void ComputeDistancesFromGaussPoints( boost::python::list rParameterResolution, int integration_degree, int max_iterations, double projection_tolerance,
                                          DoubleVector& results_vector )
    {
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        PatchVector& patch_vector = mrReconstructionDataBase.GetPatchVector();
        ModelPart& fe_model_part = mrReconstructionDataBase.GetFEModelPart();                       // INITIALIZATION
        results_vector.clear();
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
                  results_vector.push_back(distance);
                }
                else
                  outside_points_counter++;
                points_counter++;

            }
        }
    }

    // --------------------------------------------------------------------------
    void EvaluateDisplacementCoupling()
    {
      std::cout << "EvaluateDisplacementCoupling CALLED" << std::endl;

      double integral_d = 0;
      double length = 0;
      DoubleVector distances;

      BREPElementVector& brep_elements_vector = mrReconstructionDataBase.GetBREPElements();
      for(auto & brep_element_i : brep_elements_vector)
      {
        if(brep_element_i.HasCouplingCondition())
        {
          BREPGaussPointVector& coupling_gauss_points = brep_element_i.GetGaussPoints();
          for(auto & gauss_point_i : coupling_gauss_points)
          {
						// Read information from Gauss point
              unsigned int master_patch_id = gauss_point_i.GetMasterPatchId();
              Patch& master_patch = mrReconstructionDataBase.GetPatchFromPatchId( master_patch_id );
              unsigned int slave_patch_id = gauss_point_i.GetSlavePatchId();
              Patch& slave_patch = mrReconstructionDataBase.GetPatchFromPatchId( slave_patch_id );                
              double gp_i_weight = gauss_point_i.GetWeight();
              array_1d<double,2> location_on_master_patch = gauss_point_i.GetLocationOnMasterInParameterSpace();
              array_1d<double,2> location_on_slave_patch = gauss_point_i.GetLocationOnSlaveInParameterSpace();
              array_1d<double,2> tangent_on_master_patch = gauss_point_i.GetTangentOnMasterInParameterSpace();
              array_1d<double,2> tangent_on_slave_patch = gauss_point_i.GetTangentOnSlaveInParameterSpace();
              array_1d<int, 2> knotspans_on_master_patch = master_patch.ComputeSurfaceKnotSpans(location_on_master_patch);
              array_1d<int, 2> knotspans_on_slave_patch = slave_patch.ComputeSurfaceKnotSpans(location_on_slave_patch);

            // Compute Jacobian J1
              matrix<double> g_master = master_patch.ComputeBaseVectors(knotspans_on_master_patch, location_on_master_patch);
              Vector g1 = ZeroVector(3);
              g1(0) = g_master(0,0);
              g1(1) = g_master(1,0);
              g1(2) = g_master(2,0);
              Vector g2 = ZeroVector(3);
              g2(0) = g_master(0,1);
              g2(1) = g_master(1,1);
              g2(2) = g_master(2,1);
              double J1 = norm_2( g1* tangent_on_master_patch(0) + g2* tangent_on_master_patch(1) );

            // Compute distance between master and slave Gauss point
              Point<3> slave_point;
              Point<3> master_point;
              slave_patch.EvaluateSurfacePoint(location_on_slave_patch, slave_point);
              master_patch.EvaluateSurfacePoint(location_on_master_patch, master_point);
              double distance;
              EvaluateDistanceBetweenPoints(master_point, slave_point, distance);
              distances.push_back(distance);

            integral_d += distance * J1 * gp_i_weight;
            length += J1 * gp_i_weight;
          }
        }
      }
      
      // analyze the global quality of reconstruction
        Vector distancesVector = ZeroVector(distances.size());
        for (unsigned int i=0; i<distances.size(); ++i)
        {
          distancesVector(i) = distances[i];
        }
        std::cout << "\n> displacement coupling" << std::endl;
        std::cout << "> integral mean = " << integral_d/length << std::endl;      
        std::cout << "> Max value = " << norm_inf(distancesVector) << std::endl;
        std::cout << "> L2 norm = " << norm_2(distancesVector) << std::endl;      
        std::cout << "> average value = " << sum(distancesVector)/distancesVector.size() << std::endl;      
      
    }

    // --------------------------------------------------------------------------
    void EvaluateRotationCoupling()
    {
      std::cout << "EvaluateRotationCoupling CALLED" << std::endl;

      double integral_angle = 0;
      double length = 0;
      DoubleVector angles;

      BREPElementVector& brep_elements_vector = mrReconstructionDataBase.GetBREPElements();
      for(auto & brep_element_i : brep_elements_vector)
      {
        if(brep_element_i.HasCouplingCondition())
        {
          BREPGaussPointVector& coupling_gauss_points = brep_element_i.GetGaussPoints();
          for(auto & gauss_point_i : coupling_gauss_points)
          {
						// Read information from Gauss point
              unsigned int master_patch_id = gauss_point_i.GetMasterPatchId();
              Patch& master_patch = mrReconstructionDataBase.GetPatchFromPatchId( master_patch_id );
              unsigned int slave_patch_id = gauss_point_i.GetSlavePatchId();
              Patch& slave_patch = mrReconstructionDataBase.GetPatchFromPatchId( slave_patch_id );                
              double gp_i_weight = gauss_point_i.GetWeight();
              array_1d<double,2> location_on_master_patch = gauss_point_i.GetLocationOnMasterInParameterSpace();
              array_1d<double,2> location_on_slave_patch = gauss_point_i.GetLocationOnSlaveInParameterSpace();
              array_1d<double,2> tangent_on_master_patch = gauss_point_i.GetTangentOnMasterInParameterSpace();
              array_1d<double,2> tangent_on_slave_patch = gauss_point_i.GetTangentOnSlaveInParameterSpace();
              array_1d<int, 2> knotspans_on_master_patch = master_patch.ComputeSurfaceKnotSpans(location_on_master_patch);
              array_1d<int, 2> knotspans_on_slave_patch = slave_patch.ComputeSurfaceKnotSpans(location_on_slave_patch);

            // Compute Jacobian J1
              matrix<double> g_master = master_patch.ComputeBaseVectors(knotspans_on_master_patch, location_on_master_patch);
              Vector g1 = ZeroVector(3);
              g1(0) = g_master(0,0);
              g1(1) = g_master(1,0);
              g1(2) = g_master(2,0);
              Vector g2 = ZeroVector(3);
              g2(0) = g_master(0,1);
              g2(1) = g_master(1,1);
              g2(2) = g_master(2,1);
              double J1 = norm_2( g1* tangent_on_master_patch(0) + g2* tangent_on_master_patch(1) );

            // Compute normals on master and slave
              Vector g3_m = ZeroVector(3);
              g3_m(0) = g_master(0,2);
              g3_m(1) = g_master(1,2);
              g3_m(2) = g_master(2,2);
              matrix<double> g_slave = slave_patch.ComputeBaseVectors(knotspans_on_slave_patch, location_on_slave_patch);
              Vector g3_s = ZeroVector(3);
              g3_s(0) = g_slave(0,2);
              g3_s(1) = g_slave(1,2);
              g3_s(2) = g_slave(2,2);

            // Compute angle between normals on master and slave Gauss point
              double cosine = inner_prod(g3_m, g3_s);
              double angle_rad = acos(std::abs(cosine)); // abs() => ( 0 < angle_rad < PI/2 )
              double angle_deg = angle_rad * 180 / 3.1415926535;
              angles.push_back(angle_deg);

            // Sum contribution
              integral_angle += angle_deg * J1 * gp_i_weight;
              length += J1 * gp_i_weight;
          }
        }
      }
      
      // analyze the global quality of reconstruction
        Vector anglesVector = ZeroVector(angles.size());
        for (unsigned int i=0; i<angles.size(); ++i)
        {
          anglesVector(i) = angles[i];
        }
        std::cout << "\n> rotation coupling" << std::endl;
        std::cout << "> integral mean = " << integral_angle/length << std::endl;      
        std::cout << "> Max value = " << norm_inf(anglesVector) << std::endl;
        std::cout << "> L2 norm = " << norm_2(anglesVector) << std::endl;      
        std::cout << "> average value = " << sum(anglesVector)/anglesVector.size() << std::endl;      
      
    }

    // --------------------------------------------------------------------------
    void EvaluateDistanceBetweenNodes(NodeType::Pointer first, NodeType::Pointer second, double& rdistance)
    {
      Vector distance_vector = ZeroVector(3);
      distance_vector(0) = first->X() - second->X();
      distance_vector(1) = first->Y() - second->Y();
      distance_vector(2) = first->Z() - second->Z();
      rdistance = norm_2(distance_vector);
    }
    
    // --------------------------------------------------------------------------
    void EvaluateDistanceBetweenPoints(Point<3>& first, Point<3>& second, double& rdistance)
    {
      Vector distance_vector = ZeroVector(3);
      distance_vector(0) = first.X() - second.X();
      distance_vector(1) = first.Y() - second.Y();
      distance_vector(2) = first.Z() - second.Z();
      rdistance = norm_2(distance_vector);
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
    int points_counter=0;
    int outside_points_counter=0;
    /// Assignment operator.
    //      ReconstructionQualityEvaluationUtility& operator=(ReconstructionQualityEvaluationUtility const& rOther);

    /// Copy constructor.
    //      ReconstructionQualityEvaluationUtility(ReconstructionQualityEvaluationUtility const& rOther);

}; // Class ReconstructionQualityEvaluationUtility
} // namespace Kratos.

#endif // RECONSTRUCTION_QUALITY_EVALUATION_UTILITY_H
