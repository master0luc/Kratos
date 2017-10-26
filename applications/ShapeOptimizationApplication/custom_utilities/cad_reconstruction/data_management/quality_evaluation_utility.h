// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef QUALITY_EVALUATION_UTILITY_H
#define QUALITY_EVALUATION_UTILITY_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <string>
#include <cmath>
#include <iostream>
#include <numeric>

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
class QualityEvaluationUtility
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
    
    /// Pointer definition of QualityEvaluationUtility
    KRATOS_CLASS_POINTER_DEFINITION(QualityEvaluationUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    QualityEvaluationUtility( ReconstructionDataBase& reconstruction_data_base,
                              ReconstructionConditionContainer& condition_container):
     mrReconstructionDataBase( reconstruction_data_base ),
     mrReconstructionConditions( condition_container.GetReconstructionConditions() )
    {      
    }

    /// Destructor.
    virtual ~QualityEvaluationUtility()
    {
    }

    // --------------------------------------------------------------------------
    void EvaluateQuality(int max_iterations, double projection_tolerance, boost::python::list rParameterResolution, double projection_radius)
    {
      CADProjectionUtility FE2CADProjector( mrReconstructionDataBase.GetPatchVector(), max_iterations, projection_tolerance, "multiple_tree", projection_radius );
      FE2CADProjector.Initialize(rParameterResolution);

      std::vector<double> projection_distances_to_qk;
      std::vector<double> reconstruction_distances_to_qk;
      std::vector<double> reconstruction_distances_to_q;
      std::vector<unsigned int> qk_is_outside;
      
      for(auto & condition_i : mrReconstructionConditions)
      {
        // get data
            array_1d<double,3> fe_coordinates_undeformed;        condition_i->GetFECoordinatesInUndeformedConfiguration(fe_coordinates_undeformed);
            array_1d<double,3> fe_coordinates_deformed;          condition_i->GetFECoordinatesInDeformedConfiguration(fe_coordinates_deformed);
            array_1d<double,3> cad_coordinates_undeformed;       condition_i->GetCADCoordinatesInUndeformedConfiguration(cad_coordinates_undeformed);
            array_1d<double,3> cad_coordinates_deformed;         condition_i->GetCADCoordinatesInDeformedConfiguration(cad_coordinates_deformed);
            array_1d<double,3> nearest_cad_coordinates_deformed; FE2CADProjector.GetNearestCADCoordinates(fe_coordinates_deformed, nearest_cad_coordinates_deformed);
            bool qk_is_inside = condition_i->IsCADPointInside();
        
        // compute distances
          double projection_distance_to_qk =  EvaluateDistanceBetweenCoordinates(cad_coordinates_undeformed, fe_coordinates_undeformed);
          double reconstruction_distance_to_qk =  EvaluateDistanceBetweenCoordinates(cad_coordinates_deformed, fe_coordinates_deformed);
          double reconstruction_distance_to_q =  EvaluateDistanceBetweenCoordinates(nearest_cad_coordinates_deformed, fe_coordinates_deformed);

        // push_back results
          if(qk_is_inside)
          {
            projection_distances_to_qk.push_back(projection_distance_to_qk);
            reconstruction_distances_to_qk.push_back(reconstruction_distance_to_qk);
            qk_is_outside.push_back(0);
          }
          else
            qk_is_outside.push_back(1);          

          reconstruction_distances_to_q.push_back(reconstruction_distance_to_q);
      }


      // write statistics to console
      double average_dist = std::accumulate( projection_distances_to_qk.begin(), projection_distances_to_qk.end(), 0.0 ) / projection_distances_to_qk.size();
      double max = *std::max_element( projection_distances_to_qk.begin(), projection_distances_to_qk.end());
      int number_of_qk_outside = std::accumulate( qk_is_outside.begin(), qk_is_outside.end(), 0 );
      double percentage_of_qk_outside = (100.0 * number_of_qk_outside) / qk_is_outside.size();
      // double l2_norm = std::sqrt( std::inner_product(projection_distances_to_qk.begin(), projection_distances_to_qk.end(), projection_distances_to_qk.begin(), 0) );
      std::cout << "\n> Quality of projection:" << std::endl;
      std::cout << "\t average distance = " << std::fixed << std::scientific << std::setprecision(6) << average_dist << std::endl;
      std::cout << "\t max distance = " << std::fixed << std::setprecision(6) << max << std::endl;
      std::cout << "\t " << number_of_qk_outside << " out of " << qk_is_outside.size() << " points are outside of trimming boundaries (" << std::fixed << std::setprecision(2) << percentage_of_qk_outside << "%)" << std::endl;
      
      average_dist = std::accumulate( reconstruction_distances_to_qk.begin(), reconstruction_distances_to_qk.end(), 0.0 ) / reconstruction_distances_to_qk.size();
      max = *std::max_element( reconstruction_distances_to_qk.begin(), reconstruction_distances_to_qk.end());
      // l2_norm = std::sqrt( std::inner_product(reconstruction_distances_to_qk.begin(), reconstruction_distances_to_qk.end(), reconstruction_distances_to_qk.begin(), 0) );
      std::cout << "\n> Quality of reconstruction: | P - Q_k |" << std::endl;
      std::cout << "\t average distance = " << std::fixed << std::scientific << std::setprecision(6) << average_dist << std::endl;
      std::cout << "\t max distance = " << std::fixed << std::setprecision(6) << max << std::endl;
      
      average_dist = std::accumulate( reconstruction_distances_to_q.begin(), reconstruction_distances_to_q.end(), 0.0 ) / reconstruction_distances_to_q.size();
      max = *std::max_element( reconstruction_distances_to_q.begin(), reconstruction_distances_to_q.end());
      // l2_norm = std::sqrt( std::inner_product(reconstruction_distances_to_q.begin(), reconstruction_distances_to_q.end(), reconstruction_distances_to_q.begin(), 0) );
      std::cout << "\n> Quality of reconstruction: | P - Q |" << std::endl;
      std::cout << "\t average distance = " << std::fixed << std::scientific << std::setprecision(6) << average_dist << std::endl;
      std::cout << "\t max distance = " << std::fixed << std::setprecision(6) << max << std::endl;
      

    } 

    // --------------------------------------------------------------------------
    void EvaluateQualityAndOutputVtk(ModelPart& MdpaToEvaluateProjectionQuality, int max_iterations, double projection_tolerance, boost::python::list rParameterResolution, double projection_radius)
    {
      // add nodal solution step variables to the mdpa
        MdpaToEvaluateProjectionQuality.AddNodalSolutionStepVariable(PROJECTION_DISTANCE_TO_QK);
        MdpaToEvaluateProjectionQuality.AddNodalSolutionStepVariable(PATCH_ID_OF_QK);
        MdpaToEvaluateProjectionQuality.AddNodalSolutionStepVariable(INSIDE_TRIMMING_BOUNDARY_QK);
        MdpaToEvaluateProjectionQuality.AddNodalSolutionStepVariable(SHAPE_CHANGE_ABSOLUTE);
        MdpaToEvaluateProjectionQuality.AddNodalSolutionStepVariable(RECONSTRUCTION_DISTANCE_TO_QK);
        MdpaToEvaluateProjectionQuality.AddNodalSolutionStepVariable(RECONSTRUCTION_DISTANCE_TO_Q);

      CADProjectionUtility FE2CADProjector( mrReconstructionDataBase.GetPatchVector(), max_iterations, projection_tolerance, "multiple_tree", projection_radius );
      FE2CADProjector.Initialize(rParameterResolution);

      std::vector<double> projection_distances_to_qk;
      std::vector<double> reconstruction_distances_to_qk;
      std::vector<double> reconstruction_distances_to_q;
      std::vector<unsigned int> qk_is_outside;

      int counter=0;
      for(auto & condition_i : mrReconstructionConditions)
      {
        // get data
            array_1d<double,3> fe_coordinates_undeformed;        condition_i->GetFECoordinatesInUndeformedConfiguration(fe_coordinates_undeformed);
            array_1d<double,3> fe_coordinates_deformed;          condition_i->GetFECoordinatesInDeformedConfiguration(fe_coordinates_deformed);
            array_1d<double,3> cad_coordinates_undeformed;       condition_i->GetCADCoordinatesInUndeformedConfiguration(cad_coordinates_undeformed);
            array_1d<double,3> cad_coordinates_deformed;         condition_i->GetCADCoordinatesInDeformedConfiguration(cad_coordinates_deformed);
            array_1d<double,3> nearest_cad_coordinates_deformed; FE2CADProjector.GetNearestCADCoordinates(fe_coordinates_deformed, nearest_cad_coordinates_deformed);
            bool qk_is_inside = condition_i->IsCADPointInside();

        // compute distances
          double projection_distance_to_qk =  EvaluateDistanceBetweenCoordinates(cad_coordinates_undeformed, fe_coordinates_undeformed);
          double reconstruction_distance_to_qk =  EvaluateDistanceBetweenCoordinates(cad_coordinates_deformed, fe_coordinates_deformed);
          double reconstruction_distance_to_q =  EvaluateDistanceBetweenCoordinates(nearest_cad_coordinates_deformed, fe_coordinates_deformed);

        // push_back results
          if(qk_is_inside)
          {
            projection_distances_to_qk.push_back(projection_distance_to_qk);
            reconstruction_distances_to_qk.push_back(reconstruction_distance_to_qk);
            qk_is_outside.push_back(0);
          }
          else
            qk_is_outside.push_back(1);

          reconstruction_distances_to_q.push_back(reconstruction_distance_to_q);

        // create node at FE coordinates in undeformed configuration
          NodeType::Pointer new_node = Node <3>::Pointer(new Node<3>(counter, fe_coordinates_undeformed));
          new_node->SetSolutionStepVariablesList(&(MdpaToEvaluateProjectionQuality.GetNodalSolutionStepVariablesList()));
          MdpaToEvaluateProjectionQuality.AddNode(new_node);
          counter++;

        // Add nodal results
          new_node->FastGetSolutionStepValue(PROJECTION_DISTANCE_TO_QK) = cad_coordinates_undeformed - fe_coordinates_undeformed;
          new_node->FastGetSolutionStepValue(PATCH_ID_OF_QK) = condition_i->GetPatch().GetId();
          new_node->FastGetSolutionStepValue(INSIDE_TRIMMING_BOUNDARY_QK) = qk_is_inside;
          new_node->FastGetSolutionStepValue(SHAPE_CHANGE_ABSOLUTE) = fe_coordinates_deformed - fe_coordinates_undeformed;
          new_node->FastGetSolutionStepValue(RECONSTRUCTION_DISTANCE_TO_QK) = cad_coordinates_deformed - fe_coordinates_deformed;
          new_node->FastGetSolutionStepValue(RECONSTRUCTION_DISTANCE_TO_Q) = nearest_cad_coordinates_deformed - fe_coordinates_deformed;                 
      }

      // write statistics to console
      double average_dist = std::accumulate( projection_distances_to_qk.begin(), projection_distances_to_qk.end(), 0.0 ) / projection_distances_to_qk.size();
      double max = *std::max_element( projection_distances_to_qk.begin(), projection_distances_to_qk.end());
      int number_of_qk_outside = std::accumulate( qk_is_outside.begin(), qk_is_outside.end(), 0 );
      double percentage_of_qk_outside = (100.0 * number_of_qk_outside) / qk_is_outside.size();
      // double l2_norm = std::sqrt( std::inner_product(projection_distances_to_qk.begin(), projection_distances_to_qk.end(), projection_distances_to_qk.begin(), 0) );
      std::cout << "\n> Quality of projection:" << std::endl;
      std::cout << "\t average distance = " << std::fixed << std::scientific << std::setprecision(6) << average_dist << std::endl;
      std::cout << "\t max distance = " << std::fixed << std::setprecision(6) << max << std::endl;
      std::cout << "\t " << number_of_qk_outside << " out of " << qk_is_outside.size() << " points are outside of trimming boundaries (" << std::fixed << std::setprecision(2) << percentage_of_qk_outside << "%)" << std::endl;
      
      average_dist = std::accumulate( reconstruction_distances_to_qk.begin(), reconstruction_distances_to_qk.end(), 0.0 ) / reconstruction_distances_to_qk.size();
      max = *std::max_element( reconstruction_distances_to_qk.begin(), reconstruction_distances_to_qk.end());
      // l2_norm = std::sqrt( std::inner_product(reconstruction_distances_to_qk.begin(), reconstruction_distances_to_qk.end(), reconstruction_distances_to_qk.begin(), 0) );
      std::cout << "\n> Quality of reconstruction: | P - Q_k |" << std::endl;
      std::cout << "\t average distance = " << std::fixed << std::scientific << std::setprecision(6) << average_dist << std::endl;
      std::cout << "\t max distance = " << std::fixed << std::setprecision(6) << max << std::endl;
      
      average_dist = std::accumulate( reconstruction_distances_to_q.begin(), reconstruction_distances_to_q.end(), 0.0 ) / reconstruction_distances_to_q.size();
      max = *std::max_element( reconstruction_distances_to_q.begin(), reconstruction_distances_to_q.end());
      // l2_norm = std::sqrt( std::inner_product(reconstruction_distances_to_q.begin(), reconstruction_distances_to_q.end(), reconstruction_distances_to_q.begin(), 0) );
      std::cout << "\n> Quality of reconstruction: | P - Q |" << std::endl;
      std::cout << "\t average distance = " << std::fixed << std::scientific << std::setprecision(6) << average_dist << std::endl;
      std::cout << "\t max distance = " << std::fixed << std::setprecision(6) << max << std::endl;

    }    

    // --------------------------------------------------------------------------
    void EvaluateDisplacementCoupling()
    {
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
        double average_dist = std::accumulate(distances.begin(), distances.end(), 0.0) / distances.size();
        double max = *std::max_element(distances.begin(), distances.end());
        std::cout << "\n> Quality of reconstruction: G0 (displacement coupling)" << std::endl;
        std::cout << "\t average distance = " << std::fixed << std::scientific << std::setprecision(6) << average_dist << std::endl;
        std::cout << "\t max distance = " << std::fixed << std::setprecision(6) << max << std::endl;
        std::cout << "\t integral mean = " << std::fixed << std::setprecision(6) << integral_d/length << std::endl;
      
    }

    // --------------------------------------------------------------------------
    void EvaluateRotationCoupling()
    {
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
        double average_angle = std::accumulate(angles.begin(), angles.end(), 0.0) / angles.size();
        double max = *std::max_element(angles.begin(), angles.end());
        std::cout << "\n> Quality of reconstruction: G1 (rotation coupling)" << std::endl;
        std::cout << "\t average angle = " << std::fixed << std::scientific << std::setprecision(6) << average_angle << " deg" << std::endl;
        std::cout << "\t max angle = " << std::fixed << std::setprecision(6) << max << " deg" << std::endl;
        std::cout << "\t integral mean = " << std::fixed << std::setprecision(6) << integral_angle/length << " deg" << std::endl;
      
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

    // --------------------------------------------------------------------------
    double EvaluateDistanceBetweenCoordinates(array_1d<double,3>& first, array_1d<double,3>& second)
    {
      array_1d<double,3> d  = first - second;
      return std::sqrt( d[0]*d[0] + d[1]*d[1] + d[2]*d[2] );
    }

    // ==============================================================================

    /// Turn back information as a string.
    virtual std::string Info() const
    {
		return "QualityEvaluationUtility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
		rOStream << "QualityEvaluationUtility";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }

private:

    ReconstructionDataBase& mrReconstructionDataBase;
    std::vector<ReconstructionCondition::Pointer>& mrReconstructionConditions;
    int points_counter=0;
    int outside_points_counter=0;
    /// Assignment operator.
    //      QualityEvaluationUtility& operator=(QualityEvaluationUtility const& rOther);

    /// Copy constructor.
    //      QualityEvaluationUtility(QualityEvaluationUtility const& rOther);

}; // Class QualityEvaluationUtility
} // namespace Kratos.

#endif // QUALITY_EVALUATION_UTILITY_H
