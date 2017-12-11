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

#if !defined(KRATOS_CONTACT_UTILITIES)
#define KRATOS_CONTACT_UTILITIES

#include "utilities/math_utils.h"
#include "contact_structural_mechanics_application_variables.h"
#include "includes/model_part.h"

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
    
class ContactUtilities
{
public:
    ///@name Type Definitions
    ///@{
    
    // General type definitions
    typedef Node<3>                                              NodeType;
    typedef Point                                               PointType;
    typedef PointType::CoordinatesArrayType          CoordinatesArrayType;
    typedef Geometry<NodeType>                               GeometryType;
    typedef ModelPart::NodesContainerType                  NodesArrayType;
    typedef ModelPart::ConditionsContainerType        ConditionsArrayType;
    
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
     * This function scales the points according to a factor (to increase the bounding box)
     * @param PointToScale The point to scale
     * @param Normal The normal of the point
     * @param LengthSearch The factor considered to "grow" the node
     */
    
    template<class TPointType>
    static inline void ScaleNode(
        TPointType& PointToScale,
        const array_1d<double, 3>& Normal,
        const double LengthSearch
        )
    {        
        PointToScale.Coordinates() = PointToScale.Coordinates() + Normal * LengthSearch;
    }
    
    /**
     * Calculates the distance between nodes
     * @param PointOrigin The first node
     * @param PointDestiny The second node
     */
    
    static inline double DistancePoints(
        const GeometryType::CoordinatesArrayType& PointOrigin,
        const GeometryType::CoordinatesArrayType& PointDestiny
        )
    {
        return std::sqrt((PointOrigin[0] - PointDestiny[0]) * (PointOrigin[0] - PointDestiny[0])
                       + (PointOrigin[1] - PointDestiny[1]) * (PointOrigin[1] - PointDestiny[1])
                       + (PointOrigin[2] - PointDestiny[2]) * (PointOrigin[2] - PointDestiny[2]));
    }
    
    /**
     * It computes the mean of the normal in the condition in all the nodes
     * @param rModelPart The model part to compute
     */
    
    static inline void ComputeNodesMeanNormalModelPart(ModelPart& rModelPart) 
    {
        // Tolerance
        const double tolerance = std::numeric_limits<double>::epsilon();
        
        NodesArrayType& nodes_array = rModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size()); 
        
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
        for(int i = 0; i < num_nodes; ++i) 
            noalias((nodes_array.begin() + i)->FastGetSolutionStepValue(NORMAL)) = ZeroVector(3);
        
        // Sum all the nodes normals
        ConditionsArrayType& conditions_array = rModelPart.Conditions();
        
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
        for(int i = 0; i < static_cast<int>(conditions_array.size()); ++i) 
        {
            auto it_cond = conditions_array.begin() + i;
            GeometryType& this_geometry = it_cond->GetGeometry();
            
            // Aux coordinates
            CoordinatesArrayType aux_coords;
            aux_coords = this_geometry.PointLocalCoordinates(aux_coords, this_geometry.Center());
            
            it_cond->SetValue(NORMAL, this_geometry.UnitNormal(aux_coords));
            
            const unsigned int number_nodes = this_geometry.PointsNumber();
            
            for (unsigned int i = 0; i < number_nodes; ++i)
            {
                auto& this_node = this_geometry[i];
                aux_coords = this_geometry.PointLocalCoordinates(aux_coords, this_node.Coordinates());
                const array_1d<double, 3>& normal = this_geometry.UnitNormal(aux_coords);
                auto& aux_normal = this_node.FastGetSolutionStepValue(NORMAL);
                for (unsigned int index = 0; index < 3; ++index)
                {
                #ifdef _OPENMP
                    #pragma omp atomic
                #endif
                    aux_normal[index] += normal[index];
                }
            }
        }

    #ifdef _OPENMP
        #pragma omp parallel for 
    #endif
        for(int i = 0; i < num_nodes; ++i) 
        {
            auto it_node = nodes_array.begin() + i;

            array_1d<double, 3>& normal = it_node->FastGetSolutionStepValue(NORMAL);
            const double norm_normal = norm_2(normal);
            if (norm_normal > tolerance) normal /= norm_normal;
            else KRATOS_ERROR << "WARNING:: ZERO NORM NORMAL IN NODE: " << it_node->Id() << std::endl;
        }
    }

    /**
     * It calculates the center updated in u_n+1/2
     * @param ThisGeometry The geometry to calculate
     * @param DeltaTime The increment of time considered
     * @return point: The center in u_n+1/2 (Newmark)
     */
    
    static inline array_1d<double, 3> GetHalfJumpCenter(
        GeometryType& ThisGeometry,
        const double DeltaTime
        )
    {
        PointType center = ThisGeometry.Center();
        
        // Initialize variables
        Vector N;
        GeometryType::CoordinatesArrayType local_point;
        
        // Get shape functions
        ThisGeometry.PointLocalCoordinates( local_point, center.Coordinates() );
        ThisGeometry.ShapeFunctionsValues( N, local_point );
        
        const Matrix new_delta_disp = 0.25 * DeltaTime * (GetVariableMatrix(ThisGeometry, VELOCITY, 0) + GetVariableMatrix(ThisGeometry, VELOCITY, 1)) + 0.125 * DeltaTime * DeltaTime * GetVariableMatrix(ThisGeometry, ACCELERATION, 1);
        
        const Vector new_delta_disp_center = prod(trans(new_delta_disp), N);
        
        for (unsigned int i = 0; i < new_delta_disp_center.size(); ++i)
            center.Coordinates()[i] += new_delta_disp_center[i];
        
        return center.Coordinates();
    }
    
         
    /** 
     * It calculates the matrix of a variable of a geometry 
     * @param Nodes The geometry to calculate 
     * @param rVarName The name of the variable to calculate 
     * @param Step The step where calculate 
     * @return var_matrix: The matrix containing the variables of the geometry 
     */ 
     
    static inline Matrix GetVariableMatrix( 
        const GeometryType& Nodes, 
        const Variable<array_1d<double,3> >& rVarName, 
        const unsigned int Step 
        ) 
    { 
        /* DEFINITIONS */         
        const std::size_t num_nodes = Nodes.size(); 
        const std::size_t dim = Nodes.WorkingSpaceDimension(); 
        Matrix var_matrix(num_nodes, dim); 
         
        for (unsigned int i_node = 0; i_node < num_nodes; i_node++) 
        { 
            const array_1d<double, 3> value = Nodes[i_node].FastGetSolutionStepValue(rVarName, Step); 
            for (unsigned int i_dof = 0; i_dof < dim; i_dof++) 
                var_matrix(i_node, i_dof) = value[i_dof]; 
        } 
         
        return var_matrix; 
    } 
    
private:
};// class ContactUtilities

///@name Explicit Specializations
///@{

}
#endif /* KRATOS_CONTACT_UTILITIES defined */
