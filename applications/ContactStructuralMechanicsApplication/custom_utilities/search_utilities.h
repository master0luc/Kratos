// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_SEARCH_UTILITIES)
#define KRATOS_SEARCH_UTILITIES

// System includes

// External includes

// Project includes
#include "contact_structural_mechanics_application_variables.h"
#include "includes/model_part.h"
#include "geometries/point.h"

/* Custom includes */
#include "includes/mortar_classes.h"
#include "custom_includes/point_item.h"

/* Custom utilities */
#include "utilities/exact_mortar_segmentation_utility.h"
#include "custom_utilities/contact_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    
///@}
///@name  Enum's
///@{
    
///@}
///@name  Functions
///@{
    
class SearchUtilities
{
public:
    ///@name Type Definitions
    ///@{
    
    typedef Point                                     PointType;
    typedef Node<3>                                       NodeType;
    typedef Geometry<NodeType>                        GeometryType;
    typedef Geometry<PointType>                  GeometryPointType;
    typedef ModelPart::NodesContainerType           NodesArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    
    // Type definition for integration methods
    typedef GeometryData::IntegrationMethod      IntegrationMethod;
    
    // Type definitions for the point item
    typedef PointItem                                PointItemType;
    typedef PointItemType::Pointer            PointItemTypePointer;
    typedef std::vector<PointItemTypePointer>      PointItemVector;
    
    // Auxiliar geometries
    typedef Line2D2<PointType>                            LineType;
    typedef Triangle3D3<PointType>                    TriangleType;
    
    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{
    
    ///@}
    ///@name Operations
    ///@{
    
    /**
     * This function checks if the geometry can be potentially in contact using exact integration
     * @param rVariables the kinematic variables
     * @param rThisMortarConditionMatrices The mortar operators
     * @param IntegrationUtility The exact integration utility
     * @param SlaveGeometry The geometry of the slave
     * @param MasterGeometry The geometry of the master
     * @param SlaveNormal The normals of the slave
     * @param MasterNormal The normals of the master
     * @param ActiveCheckLength The threshold distance to check the potential contact
     * @return condition_is_active: True if at least one node is active, false otherwise
     */
    template< const unsigned int TDim, const unsigned int TNumNodes, const bool TFill>
    static inline bool CheckExactIntegration(
        MortarKinematicVariables<TNumNodes> rVariables,
        MortarOperator<TNumNodes> rThisMortarConditionMatrices,   
        ExactMortarIntegrationUtility<TDim, TNumNodes> IntegrationUtility,
        GeometryType& SlaveGeometry,             // SLAVE
        GeometryType& MasterGeometry,            // MASTER
        const array_1d<double, 3>& SlaveNormal,  // SLAVE
        const array_1d<double, 3>& MasterNormal, // SLAVE
        const double ActiveCheckLength
        )
    {
        // Define the basic information
        bool condition_is_active = false;
        IntegrationMethod this_integration_method = GeometryData::GI_GAUSS_2;

        // Reading integration points
        std::vector<array_1d<PointType,TDim>> conditions_points_slave;
        const bool is_inside = IntegrationUtility.GetExactIntegration(SlaveGeometry, SlaveNormal, MasterGeometry, MasterNormal, conditions_points_slave);
        
        if (is_inside == true)
        {
            // Initialize general variables for the current master element
            rVariables.Initialize();
    
            // Initialize the mortar operators
            rThisMortarConditionMatrices.Initialize();
            
            for (unsigned int i_geom = 0; i_geom < conditions_points_slave.size(); i_geom++)
            {
                std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
                for (unsigned int i_node = 0; i_node < TDim; i_node++)
                {
                    PointType global_point;
                    SlaveGeometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                    points_array[i_node] = boost::make_shared<PointType>(global_point);
                }
                
                typename std::conditional<TDim == 2, LineType, TriangleType >::type decomp_geom( points_array );
                
                const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, SlaveGeometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);
                
                if (bad_shape == false)
                {
                    const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );
                    
                    // Integrating the mortar operators
                    for ( unsigned int point_number = 0; point_number < integration_points_slave.size(); point_number++ )
                    {
                        const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                        PointType local_point_parent;
                        PointType gp_global;
                        decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                        SlaveGeometry.PointLocalCoordinates(local_point_parent, gp_global);
                        
                        // Calculate the kinematic variables
                        const PointType& local_point = integration_points_slave[point_number].Coordinates();

                        /// SLAVE CONDITION ///
                        SlaveGeometry.ShapeFunctionsValues( rVariables.NSlave, local_point.Coordinates() );
                        rVariables.PhiLagrangeMultipliers = rVariables.NSlave;
                        rVariables.DetjSlave = 1.0; // NOTE: Using unitary area (we try to compute the gap, not the weighted gap) SlaveGeometry.DeterminantOfJacobian( local_point );
                        
                        /// MASTER CONDITION ///
                        PointType projected_gp_global;
                        const array_1d<double,3> gp_normal = MortarUtilities::GaussPointUnitNormal(rVariables.NSlave, SlaveGeometry);
                        
                        GeometryType::CoordinatesArrayType slave_gp_global;
                        SlaveGeometry.GlobalCoordinates( slave_gp_global, local_point );
                        MortarUtilities::FastProjectDirection( MasterGeometry, slave_gp_global, projected_gp_global, MasterNormal, -gp_normal ); // The opposite direction
                        
                        GeometryType::CoordinatesArrayType projected_gp_local;
                        
                        MasterGeometry.PointLocalCoordinates(projected_gp_local, projected_gp_global.Coordinates( ) ) ;
                        
                        // SHAPE FUNCTIONS 
                        MasterGeometry.ShapeFunctionsValues( rVariables.NMaster, projected_gp_local );    
                        
                        const double integration_weight = integration_points_slave[point_number].Weight(); // NOTE: Error in axisymmetric
                        
                        rThisMortarConditionMatrices.CalculateMortarOperators(rVariables, integration_weight);   
                    }

                    /* Setting the gap */

                    // Current coordinates 
                    const bounded_matrix<double, TNumNodes, TDim> x1 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(SlaveGeometry);
                    const bounded_matrix<double, TNumNodes, TDim> x2 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(MasterGeometry);
            
                    const bounded_matrix<double, TNumNodes, TDim> Dx1Mx2 = prod(rThisMortarConditionMatrices.DOperator, x1) - prod(rThisMortarConditionMatrices.MOperator, x2); 
                    
                    array_1d<double, TDim> aux_slave_normal;
                    for (unsigned int i_dim = 0; i_dim < TDim; i_dim++)
                    {
                        aux_slave_normal[i_dim] = SlaveNormal[i_dim];
                    }
                    
                    if (TFill == true)
                    {
                        for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
                        {
                            if (SlaveGeometry[i_node].Is(ACTIVE) == false)
                            {
                                const array_1d<double, 3> delta_disp = SlaveGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT, 0) - SlaveGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT, 1);
                                
                                // We check if the movement is in the direction of the normal
                                const bool moving_gap_direction = norm_2(delta_disp) > 0.0 ? (inner_prod(delta_disp, SlaveNormal) > 0 ) : true;
                                
                                const array_1d<double, TDim> aux_array = row(Dx1Mx2, i_node);
                                
                                const double nodal_gap = inner_prod(aux_array, - aux_slave_normal); 
                                
                                if ((nodal_gap < ActiveCheckLength) && (moving_gap_direction == true))
                                {                                    
                                    SlaveGeometry[i_node].Set(ACTIVE, true);
                                    condition_is_active = true;
                                }
                            }
                            else
                            {
                                condition_is_active = true;
                            }
                        }
                    }
                    else
                    {
                        for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
                        {
                            const array_1d<double, TDim> aux_array = row(Dx1Mx2, i_node);
                            
                            const double nodal_gap = inner_prod(aux_array, - aux_slave_normal); 
                            
                            if (nodal_gap < ActiveCheckLength)
                            {
                                return true;
                            }
                        }
                    }
                }
            }
        } 
        
        return condition_is_active;
    }
    
protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{       

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
};// class SearchUtilities

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
#endif /* KRATOS_SEARCH_UTILITIES defined */
