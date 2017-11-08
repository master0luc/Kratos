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

#if !defined(KRATOS_DERIVATIVES_UTILITIES)
#define KRATOS_DERIVATIVES_UTILITIES

#include "contact_structural_mechanics_application_variables.h"

/* Includes */
#include "includes/model_part.h"
#include "includes/mortar_classes.h"

/* Utilities */
#include "utilities/mortar_utilities.h"
#include "utilities/math_utils.h"

/* Geometries */
#include "geometries/point.h"
#include "geometries/line_2d_2.h"
#include "geometries/triangle_3d_3.h"

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
    
template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
class DerivativesUtilities
{
public:
    ///@name Type Definitions
    ///@{
    
    // General type definitions
    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef Geometry<NodeType> GeometryType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Properties PropertiesType;
    
    typedef typename std::conditional<TNumNodes == 2, PointBelongsLine2D2N, typename std::conditional<TNumNodes == 3, PointBelongsTriangle3D3N, PointBelongsQuadrilateral3D4N>::type>::type BelongType;
    
    typedef PointBelong<TNumNodes>                                                 PointBelongType;
    
    typedef Geometry<PointBelongType>                                      GeometryPointBelongType;
    
    typedef array_1d<PointBelongType,TDim>                                      ConditionArrayType;
    
    typedef Line2D2<PointType>                                                            LineType;
    
    typedef Triangle3D3<PointType>                                                    TriangleType;
    
    typedef typename std::conditional<TDim == 2, LineType, TriangleType >::type  DecompositionType;
    
    typedef typename std::conditional<TFrictional == true, DerivativeDataFrictional<TDim, TNumNodes>, DerivativeData<TDim, TNumNodes> >::type DerivativeDataType;
    
    static constexpr unsigned int MatrixSize = TFrictional == true ? TDim * (TNumNodes + TNumNodes + TNumNodes) : TDim * (TNumNodes + TNumNodes) + TNumNodes;
    
    typedef MortarKinematicVariablesWithDerivatives<TDim, TNumNodes>              GeneralVariables;
    
    typedef DualLagrangeMultiplierOperatorsWithDerivatives<TDim, TNumNodes, TFrictional>    AeData;
    
    typedef MortarOperatorWithDerivatives<TDim, TNumNodes, TFrictional>    MortarConditionMatrices;
    
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
     * This method is used to compute the directional derivatives of the jacobian determinant
     */
    static inline void CalculateDeltaDetjSlave(
        GeneralVariables& rVariables,
        DerivativeDataType& rDerivativeData
        )
    {
        if (TDim == 2)
        {
            // Fill up the elements corresponding to the slave DOFs - the rest remains zero
            for ( unsigned int i_slave = 0, i = 0; i_slave < TNumNodes; ++i_slave, i += TDim )
            {
                rDerivativeData.DeltaDetjSlave[i    ] = rVariables.jSlave( 0, 0 ) * rVariables.DNDeSlave( i_slave, 0) / rVariables.DetjSlave;
                rDerivativeData.DeltaDetjSlave[i + 1] = rVariables.jSlave( 1, 0 ) * rVariables.DNDeSlave( i_slave, 0) / rVariables.DetjSlave;
            }
        }
        else
        {
            const array_1d<double,TNumNodes>& DNDxi  = column( rVariables.DNDeSlave, 0 );
            const array_1d<double,TNumNodes>& DNDeta = column( rVariables.DNDeSlave, 1 );
            
            const array_1d<double,TDim>& Jxi  = column( rVariables.jSlave, 0 );
            const array_1d<double,TDim>& Jeta = column( rVariables.jSlave, 1 );
            
            const array_1d<double,TDim>& normal = prod(trans(rDerivativeData.NormalMaster), rVariables.NSlave);
            
            bounded_matrix<double, TDim, TDim> DeltaJxixJeta;
            
            for ( unsigned int i_slave = 0, i = 0; i_slave < TNumNodes; ++i_slave, i += TDim )
            {
                DeltaJxixJeta(0,0) = 0.0;
                DeltaJxixJeta(0,1) =  Jeta(2) * DNDxi(i_slave) - Jxi(2) * DNDeta(i_slave); 
                DeltaJxixJeta(0,2) = -Jeta(1) * DNDxi(i_slave) + Jxi(1) * DNDeta(i_slave); 
                DeltaJxixJeta(1,0) = -Jeta(2) * DNDxi(i_slave) + Jxi(2) * DNDeta(i_slave); 
                DeltaJxixJeta(1,1) = 0.0;
                DeltaJxixJeta(1,2) =  Jeta(0) * DNDxi(i_slave) - Jxi(0) * DNDeta(i_slave);
                DeltaJxixJeta(2,0) =  Jeta(1) * DNDxi(i_slave) - Jxi(1) * DNDeta(i_slave); 
                DeltaJxixJeta(2,1) = -Jeta(0) * DNDxi(i_slave) + Jxi(0) * DNDeta(i_slave); 
                DeltaJxixJeta(2,2) = 0.0;
                
                rDerivativeData.DeltaDetjSlave[i    ] = inner_prod( normal, column( DeltaJxixJeta, 0 ) );
                rDerivativeData.DeltaDetjSlave[i + 1] = inner_prod( normal, column( DeltaJxixJeta, 1 ) );
                rDerivativeData.DeltaDetjSlave[i + 2] = inner_prod( normal, column( DeltaJxixJeta, 2 ) );
            }
        }
    }
    
    /**
     * This method is used to compute the local increment of the normal
     * NOTE: Not the mean, look in the contact utilities 
     */
    static inline bounded_matrix<double, TDim, TDim> LocalDeltaNormal(
        const GeometryType& CondGeometry,
        const unsigned int NodeIndex
        )
    {
        // Tolerance
        const double tolerance = std::numeric_limits<double>::epsilon();
            
        bounded_matrix<double, TDim, TDim> DeltaNeAdj;
        bounded_matrix<double, TDim, TDim> Ce;
        
        const bounded_matrix<double, TDim, TDim> I = IdentityMatrix(TDim, TDim);
        
        bounded_matrix<double, TDim, TDim> delta_normal = ZeroMatrix(TDim,TDim);
        
        // Normalized condition normal
        PointType auxiliar_center;
        auxiliar_center.Coordinates() = ZeroVector(3);
        const array_1d<double, 3>& Ne = CondGeometry.UnitNormal(auxiliar_center);
        bounded_matrix<double, TDim, TDim> NeoNe = subrange( outer_prod( Ne, Ne ), 0, TDim, 0, TDim );
        
        // Auxiliar value
        // The norm of a geometry's normal is its characteristic dimension - length for 2D and area for 3D 
        double NeNorm = (TDim == 2) ? CondGeometry.Length( ) : CondGeometry.Area( );
        
        if (TDim == 2)
        {                
            DeltaNeAdj( 0, 0 ) =  0.0;
            DeltaNeAdj( 0, 1 ) = -1.0;
            DeltaNeAdj( 1, 0 ) =  1.0;
            DeltaNeAdj( 1, 1 ) =  0.0;
            
            const double DNDej = (NodeIndex == 0) ? - 0.5 : 0.5;
            
            Ce = prod( I - NeoNe, DeltaNeAdj ) / NeNorm; // In 2D, DeltaNeAdj is node-independent => evaluated outside the nodes loop
            
            delta_normal = - 2.0 * Ce * DNDej; // NOTE: Check why - 2???!!!, it was the only wayto ensure the same value as the symbolic. You will need to repeat this in 3D            
        //         delta_normal = Ce * DNDej;     
        }
        else
        {
            MatrixType J = ZeroMatrix( 3, 2 ); // Jacobian [ 3D global x 2D local ]
            array_1d<double, 2> DNDej;
            array_1d<double, 3> LocalCoordsj;
            
            if( TNumNodes == 3 )    // linear triangle element
            {
                if( NodeIndex == 0 )
                {
                    LocalCoordsj[0] = 0.0;
                    LocalCoordsj[1] = 0.0;
                    DNDej[0] = - 1.0;
                    DNDej[1] = - 1.0;
                }
                else if( NodeIndex == 1 )
                {
                    LocalCoordsj[0] = 1.0;
                    LocalCoordsj[1] = 0.0;
                    DNDej[0] = 1.0;
                    DNDej[1] = 0.0;
                }
                else // NodeIndex == 2
                {
                    LocalCoordsj[0] = 0.0;
                    LocalCoordsj[1] = 1.0;
                    DNDej[0] = 0.0;
                    DNDej[1] = 1.0;
                }
            }
            else if( TNumNodes == 4 )    // linear quad element 
            {
                if( NodeIndex == 0 )
                {
                    LocalCoordsj[0] = - 1.0;
                    LocalCoordsj[1] = - 1.0;
                    DNDej[0] = - 0.5;
                    DNDej[1] = - 0.5;
                }
                else if( NodeIndex == 1 )
                {
                    LocalCoordsj[0] =   1.0;
                    LocalCoordsj[1] = - 1.0;
                    DNDej[0] =   0.5;
                    DNDej[1] = - 0.5;
                }
                else if( NodeIndex == 2 )
                {
                    LocalCoordsj[0] =  1.0;
                    LocalCoordsj[1] =  1.0;
                    DNDej[0] =  0.5;
                    DNDej[1] =  0.5;
                }
                else // NodeIndex == 3
                {
                    LocalCoordsj[0] = - 1.0;
                    LocalCoordsj[1] =   1.0;
                    DNDej[0] = - 0.5;
                    DNDej[1] =   0.5;
                }
            }
            
            CondGeometry.Jacobian( J, LocalCoordsj );
            
            DeltaNeAdj(0,0) = 0.0;
            DeltaNeAdj(0,1) = +J(2,1) * DNDej[0] - J(2,0) * DNDej[1]; 
            DeltaNeAdj(0,2) = -J(1,1) * DNDej[0] + J(1,0) * DNDej[1]; 
            DeltaNeAdj(1,0) = -J(2,1) * DNDej[0] + J(2,0) * DNDej[1]; 
            DeltaNeAdj(1,1) = 0.0;                   
            DeltaNeAdj(1,2) = +J(0,1) * DNDej[0] - J(0,0) * DNDej[1]; 
            DeltaNeAdj(2,0) = +J(1,1) * DNDej[0] - J(1,0) * DNDej[1]; 
            DeltaNeAdj(2,1) = -J(0,1) * DNDej[0] + J(0,0) * DNDej[1]; 
            DeltaNeAdj(2,2) = 0.0;
            
            Ce = prod( I - NeoNe, DeltaNeAdj ) / NeNorm;
            delta_normal = Ce;
        }
        
        NeNorm = norm_2( Ne );
        const double NeNorm3 = NeNorm * NeNorm * NeNorm;
        
        if ( NeNorm3 > tolerance )
        {
            const bounded_matrix<double, TDim, TDim> Cj = I / NeNorm - NeoNe / NeNorm3;
            delta_normal = prod( Cj, delta_normal );
        }
            
        return delta_normal; 
    }
    
    /**
     * Calculates the increment of the normal in the slave condition
     */
    
    static inline void CalculateDeltaNormalSlave(
        DerivativeDataType& rDerivativeData,
        const GeometryType& SlaveGeometry
        )
    {
        if (TDim == 2)
        {
            for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave)
            {
                const bounded_matrix<double, TDim, TDim>& delta_normal = SlaveGeometry[i_slave].GetValue(DELTA_NORMAL);
//                 const bounded_matrix<double, TDim, TDim> delta_normal = this->LocalDeltaNormal(SlaveGeometry, i_slave);
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof) 
                {
                    row(rDerivativeData.DeltaNormalSlave[i_slave * TDim + i_dof], i_slave) = trans(column(delta_normal, i_dof)); 
                }
            }
        }
    }
    
    /**
     * Calculates the increment of the normal and in the master condition
     */
    
    static inline void CalculateDeltaNormalMaster(
        DerivativeDataType& rDerivativeData,
        const GeometryType& MasterGeometry
        )
    {
        if (TDim == 2)
        {
            for ( unsigned int i_master = 0; i_master < 2; ++i_master )
            {
//                const bounded_matrix<double, 2, 2>& delta_normal = GetGeometry[i_master].GetValue(DELTA_NORMAL);
                const bounded_matrix<double, 2, 2> delta_normal = LocalDeltaNormal(MasterGeometry, i_master);
                for (unsigned i_dof = 0; i_dof < 2; ++i_dof) 
                {
                    row(rDerivativeData.DeltaNormalMaster[i_master * 2 + i_dof], i_master) = trans(column(delta_normal, i_dof));
                }
            }
        }
    }
    
    /**
     * This method is used to compute the directional derivatives of the cell vertex
     * @param rVariables: The kinematic variables
     * @param rDerivativeData: The derivatives container
     * @param TheseBelongs: The belongs list used in the derivatives
     * @param ConsiderNormalVariation: If consider the normal derivative
     * @param SlaveGeometry: The slave geometry
     * @param MasterGeometry: The master geometry
     * @param Normal: The normal vector
     */
    static inline void CalculateDeltaCellVertex(
        GeneralVariables& rVariables,
        DerivativeDataType& rDerivativeData,
        const array_1d<BelongType, TDim>& TheseBelongs,
        const bool ConsiderNormalVariation,
        GeometryType& SlaveGeometry,
        GeometryType& MasterGeometry,
        const array_1d<double, 3>& Normal
        )
    {
        // The Normal and delta Normal in the center of the element
        bounded_matrix<double, TDim, TDim> delta_normal = ZeroMatrix(TDim, TDim);
        
        const double aux_nodes_coeff = static_cast<double>(TNumNodes);
        
        const Point& slave_center = SlaveGeometry.Center().Coordinates();
        
//     #ifdef KRATOS_DEBUG
//         for (unsigned i_triangle = 0; i_triangle < 3; ++i_triangle) 
//         {
//             KRATOS_WATCH(static_cast<unsigned int>(TheseBelongs[i_triangle]));
//         }
//     #endif
        
        for (unsigned i_triangle = 0; i_triangle < 3; ++i_triangle) 
        {
            if (TheseBelongs[i_triangle] >= 2 * TNumNodes) // It belongs to an intersection
            {    
                // We compute the indexes
                unsigned int belong_index_slave_start, belong_index_slave_end, belong_index_master_start, belong_index_master_end;
                ConvertAuxHashIndex(static_cast<unsigned int>(TheseBelongs[i_triangle]), belong_index_slave_start, belong_index_slave_end, belong_index_master_start, belong_index_master_end);
                
                // The coordinates should be in the projected plane
                const array_1d<double, 3>& xs1 = MortarUtilities::FastProject(slave_center, SlaveGeometry[belong_index_slave_start], Normal).Coordinates(); // Start coordinates of the first segment
                const array_1d<double, 3>& xe1 = MortarUtilities::FastProject(slave_center, SlaveGeometry[belong_index_slave_end], Normal).Coordinates(); // End coordinates of the first segment
                const array_1d<double, 3>& xs2 = MortarUtilities::FastProject(slave_center, MasterGeometry[belong_index_master_start], Normal).Coordinates(); // Start coordinates of the second segment
                const array_1d<double, 3>& xe2 = MortarUtilities::FastProject(slave_center, MasterGeometry[belong_index_master_end], Normal).Coordinates(); // End coordinates of the second segment
                
//                 const array_1d<double, 3>& xs1 = SlaveGeometry[belong_index_slave_start].Coordinates();   // Start coordinates of the first segment
//                 const array_1d<double, 3>& xe1 = SlaveGeometry[belong_index_slave_end].Coordinates();     // End coordinates of the first segment
//                 const array_1d<double, 3>& xs2 = MasterGeometry[belong_index_master_start].Coordinates(); // Start coordinates of the second segment
//                 const array_1d<double, 3>& xe2 = MasterGeometry[belong_index_master_end].Coordinates();   // End coordinates of the second segment
                
//                 std::cout << belong_index_slave_start << "\t" << belong_index_slave_end << "\t" << belong_index_master_start << "\t" << belong_index_master_end << std::endl;
                
                array_1d<double, 3> aux_num, aux_denom;
                MathUtils<double>::CrossProduct(aux_num,   xs1 - xs2, xe2 - xs2);
                MathUtils<double>::CrossProduct(aux_denom, xe1 - xs1, xe2 - xs2);
                const double num   = inner_prod(aux_num,   Normal);
                const double denom = inner_prod(aux_denom, Normal);
                
                // We compute the first part
                const array_1d<double, 3>& aux_coords = (xe1 - xs1);
                array_1d<double, 3> aux_vertex_vector0, aux_vertex_vector1, aux_cross_product;
                
                // First term 
                aux_vertex_vector0 = ZeroVector(3);
                aux_vertex_vector1 = ZeroVector(3);
                
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof)
                {
                    aux_vertex_vector0 += LocalDeltaVertex(Normal, delta_normal, i_dof, belong_index_slave_start, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, 1.0);
                    
                    aux_vertex_vector1 += LocalDeltaVertex(Normal, delta_normal, i_dof, belong_index_slave_end, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, - 1.0);   
                }
                
                MathUtils<double>::CrossProduct(aux_cross_product, aux_vertex_vector0, xe2 - xs2);
                const double coeff1a = inner_prod(aux_cross_product, Normal);
                MathUtils<double>::CrossProduct(aux_cross_product, aux_vertex_vector1, xe2 - xs2);
                const double coeff1b = inner_prod(aux_cross_product, Normal);
                
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof)
                {
                    bounded_matrix<double, 3, 3>& local_delta_vertexa = rDerivativeData.DeltaCellVertex[belong_index_slave_start * TDim + i_dof];
                    bounded_matrix<double, 3, 3>& local_delta_vertexb = rDerivativeData.DeltaCellVertex[belong_index_slave_end * TDim + i_dof];
                    row(local_delta_vertexa, i_triangle) -= coeff1a/denom * aux_coords;
                    row(local_delta_vertexb, i_triangle) -= coeff1b/denom * aux_coords; 
                }
                
                // Second term
                aux_vertex_vector0 = ZeroVector(3);
                aux_vertex_vector1 = ZeroVector(3);
                
                if (ConsiderNormalVariation == true && belong_index_master_end < TNumNodes) delta_normal = LocalDeltaNormal(SlaveGeometry, belong_index_master_end) * (1.0/aux_nodes_coeff);
                else delta_normal = ZeroMatrix(3, 3);
                
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof)
                {
                    aux_vertex_vector0 += LocalDeltaVertex(Normal, delta_normal, i_dof, (belong_index_master_end + TNumNodes), ConsiderNormalVariation, SlaveGeometry, MasterGeometry, 1.0);
                }
                
                if (ConsiderNormalVariation == true && belong_index_master_start < TNumNodes) delta_normal = LocalDeltaNormal(SlaveGeometry, belong_index_master_start) * (1.0/aux_nodes_coeff);
                else delta_normal = ZeroMatrix(3, 3);
                
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof)
                {
                    aux_vertex_vector1 += LocalDeltaVertex(Normal, delta_normal, i_dof, (belong_index_master_start + TNumNodes), ConsiderNormalVariation, SlaveGeometry, MasterGeometry, - 1.0);   
                }

                MathUtils<double>::CrossProduct(aux_cross_product, xs1 - xs2, aux_vertex_vector0);
                const double coeff2a = inner_prod(aux_cross_product, Normal);
                MathUtils<double>::CrossProduct(aux_cross_product, xs1 - xs2, aux_vertex_vector1);
                const double coeff2b = inner_prod(aux_cross_product, Normal);
                
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof)
                {
                    bounded_matrix<double, 3, 3>& local_delta_vertexa = rDerivativeData.DeltaCellVertex[(belong_index_master_end + TNumNodes) * TDim + i_dof];
                    bounded_matrix<double, 3, 3>& local_delta_vertexb = rDerivativeData.DeltaCellVertex[(belong_index_master_start + TNumNodes) * TDim + i_dof];
                    row(local_delta_vertexa, i_triangle) -= coeff2a/denom * aux_coords;
                    row(local_delta_vertexb, i_triangle) -= coeff2b/denom * aux_coords; 
                }
                
                // Third term
                aux_vertex_vector0 = ZeroVector(3);
                aux_vertex_vector1 = ZeroVector(3);
                
                if (ConsiderNormalVariation == true && belong_index_slave_end < TNumNodes) delta_normal = LocalDeltaNormal(SlaveGeometry, belong_index_slave_end) * (1.0/aux_nodes_coeff);
                else delta_normal = ZeroMatrix(3, 3);
                
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof)
                {
                    aux_vertex_vector0 += LocalDeltaVertex(Normal, delta_normal, i_dof, belong_index_slave_end, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, - 1.0);
                }
                    
                if (ConsiderNormalVariation == true && belong_index_master_end < TNumNodes) delta_normal = LocalDeltaNormal(SlaveGeometry, belong_index_master_end) * (1.0/aux_nodes_coeff);
                else delta_normal = ZeroMatrix(3, 3);

                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof)
                {
                    aux_vertex_vector1 += LocalDeltaVertex(Normal, delta_normal, i_dof, belong_index_slave_start, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, 1.0);   
                }

                MathUtils<double>::CrossProduct(aux_cross_product, aux_vertex_vector0, xe2 - xs2);
                const double coeff3a = inner_prod(aux_cross_product, Normal);
                MathUtils<double>::CrossProduct(aux_cross_product, aux_vertex_vector1, xe2 - xs2);
                const double coeff3b = inner_prod(aux_cross_product, Normal);
                
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof)
                {
                    bounded_matrix<double, 3, 3>& local_delta_vertexa = rDerivativeData.DeltaCellVertex[belong_index_slave_end * TDim + i_dof];
                    bounded_matrix<double, 3, 3>& local_delta_vertexb = rDerivativeData.DeltaCellVertex[belong_index_slave_start * TDim + i_dof];
                    row(local_delta_vertexa, i_triangle) -= num * coeff3a/std::pow(denom, 2) * aux_coords; 
                    row(local_delta_vertexb, i_triangle) -= num * coeff3b/std::pow(denom, 2) * aux_coords; 
                }
                
                // Fourth term
                aux_vertex_vector0 = ZeroVector(3);
                aux_vertex_vector1 = ZeroVector(3);

                if (ConsiderNormalVariation == true && belong_index_master_end < TNumNodes) delta_normal = LocalDeltaNormal(SlaveGeometry, belong_index_master_end) * (1.0/aux_nodes_coeff);
                else delta_normal = ZeroMatrix(3, 3);
                
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof)
                {
                    aux_vertex_vector0 += LocalDeltaVertex(Normal, delta_normal, i_dof, (belong_index_master_end + TNumNodes), ConsiderNormalVariation, SlaveGeometry, MasterGeometry, - 1.0);
                }
                    
                if (ConsiderNormalVariation == true && belong_index_master_end < TNumNodes) delta_normal = LocalDeltaNormal(SlaveGeometry, belong_index_master_end) * (1.0/aux_nodes_coeff);
                else delta_normal = ZeroMatrix(3, 3);
                    
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof)
                {
                    aux_vertex_vector0 += LocalDeltaVertex(Normal, delta_normal, i_dof, (belong_index_master_start + TNumNodes), ConsiderNormalVariation, SlaveGeometry, MasterGeometry, 1.0);   
                }

                MathUtils<double>::CrossProduct(aux_cross_product, xe1 - xs1, aux_vertex_vector0);
                const double coeff4a = inner_prod(aux_cross_product, Normal);
                MathUtils<double>::CrossProduct(aux_cross_product, xe1 - xs1, aux_vertex_vector1);
                const double coeff4b = inner_prod(aux_cross_product, Normal);
                
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof)
                {
                    bounded_matrix<double, 3, 3>& local_delta_vertexa =  rDerivativeData.DeltaCellVertex[(belong_index_master_end + TNumNodes) * TDim + i_dof];
                    bounded_matrix<double, 3, 3>& local_delta_vertexb =  rDerivativeData.DeltaCellVertex[(belong_index_master_start + TNumNodes) * TDim + i_dof];
                    row(local_delta_vertexa, i_triangle) -= num * coeff4a/std::pow(denom, 2) * aux_coords; 
                    row(local_delta_vertexb, i_triangle) -= num * coeff4b/std::pow(denom, 2) * aux_coords; 
                }
                
                // Part corresponding to the delta Normal
                // We compute the first part
                if (ConsiderNormalVariation == true)
                {
                    for ( unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave)
                    {
                        delta_normal = LocalDeltaNormal(SlaveGeometry, i_slave) * (1.0/aux_nodes_coeff);
                        
                        for (unsigned i_dof = 0; i_dof < TDim; ++i_dof)
                        {
                            bounded_matrix<double, 3, 3>& local_delta_vertex =  rDerivativeData.DeltaCellVertex[i_slave * TDim + i_dof];
                            
                            row(local_delta_vertex, i_triangle) += inner_prod(aux_num,   trans(column(delta_normal, i_dof)))/std::pow(denom, 2) * aux_coords; 
                            
                            row(local_delta_vertex, i_triangle) -= inner_prod(aux_denom, trans(column(delta_normal, i_dof)))/std::pow(denom, 2) * aux_coords; 
                        }
                    }
                }
                
                // We compute the second part
                const double coeff0 = num/denom;
                
                if (ConsiderNormalVariation == true && belong_index_slave_end < TNumNodes) delta_normal = LocalDeltaNormal(SlaveGeometry, belong_index_slave_end) * (1.0/aux_nodes_coeff);
                else delta_normal = ZeroMatrix(3, 3);
                    
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof)
                {
                    // The term corresponding to xe1
                    bounded_matrix<double, 3, 3>& delta_cell_vertex_slave_end = rDerivativeData.DeltaCellVertex[belong_index_slave_end * TDim + i_dof];
                    
                    LocalDeltaVertex(delta_cell_vertex_slave_end, Normal, delta_normal, i_dof, i_triangle, belong_index_slave_end, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, - coeff0);
                }
                
                if (ConsiderNormalVariation == true && belong_index_slave_start < TNumNodes) delta_normal = LocalDeltaNormal(SlaveGeometry, belong_index_slave_start) * (1.0/aux_nodes_coeff);
                else delta_normal = ZeroMatrix(3, 3);

                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof)
                {
                    // The term corresponding to xs1
                    bounded_matrix<double, 3, 3>& delta_cell_vertex_slave_start = rDerivativeData.DeltaCellVertex[belong_index_slave_start * TDim + i_dof];
                    
                    LocalDeltaVertex(delta_cell_vertex_slave_start, Normal, delta_normal, i_dof, i_triangle, belong_index_slave_start, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, 1.0 + coeff0);
                }
            }
            else // It belongs to a master/slave node
            {
                const unsigned int belong_index = static_cast<unsigned int>(TheseBelongs[i_triangle]);
                
                if (ConsiderNormalVariation == true && belong_index < TNumNodes) delta_normal = LocalDeltaNormal(SlaveGeometry, belong_index) * (1.0/aux_nodes_coeff);
                else delta_normal = ZeroMatrix(3, 3);
                
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof)
                {
                    bounded_matrix<double, 3, 3>& local_delta_vertex = rDerivativeData.DeltaCellVertex[belong_index * TDim + i_dof];
                    
                    LocalDeltaVertex(local_delta_vertex, Normal, delta_normal, i_dof, i_triangle, belong_index, ConsiderNormalVariation, SlaveGeometry, MasterGeometry);
                }
            }
        }
    }

    /**
     * Calculates the increment of the shape functions
     * @param rVariables: The kinematic variables
     * @param rDerivativeData: The derivatives container
     * @param SlaveGeometry: The geometry of the slave side
     * @param MasterGeometry: The geometry of the master side
     * @param SlaveNormal: The normal of the slave side
     * @param MasterNormal: The normal of the master side
     * @param DecompGeom: The triangle used to decompose the geometry
     * @param LocalPointDecomp: The local coordinates in the decomposed geometry
     * @param LocalPointParent: The local coordinates in the slave geometry
     * @param ConsiderNormalVariation: If consider the normal derivative
     */
    
    static inline void CalculateDeltaN(
        const GeneralVariables& rVariables,
        DerivativeDataType& rDerivativeData,
        GeometryType& SlaveGeometry,
        GeometryType& MasterGeometry,
        const array_1d<double, 3> SlaveNormal,
        const array_1d<double, 3> MasterNormal,
        const DecompositionType& DecompGeom,
        const PointType& LocalPointDecomp,
        const PointType& LocalPointParent,
        const bool ConsiderNormalVariation = false
        )
    {
        /* Shape functions */
        const VectorType& N1 = rVariables.NSlave;
        const VectorType& N2 = rVariables.NMaster;
        
        /* Local gradients */
        const MatrixType& DNDe1 = rVariables.DNDeSlave;
        const MatrixType& DNDe2 = rVariables.DNDeMaster;
        
        // The Normal and delta Normal in the center of the element
        bounded_matrix<double, TDim, TDim> delta_normal = ZeroMatrix(TDim, TDim);
        
        /* Shape function decomposition */
        VectorType N_decomp;
        DecompGeom.ShapeFunctionsValues( N_decomp, LocalPointDecomp.Coordinates() );
        
        if (TDim == 2)
        {
            // TODO: Finish this!!!!
        }
        else
        {
            for ( unsigned int i_node = 0; i_node < 2 * TNumNodes; ++i_node)
            {
                for (unsigned i_dof = 0; i_dof < TDim; ++i_dof) 
                {
                    array_1d<double, 3> aux_RHS1 = ZeroVector(3);
                    
                    // The vertex cell contribution
                    const auto& local_delta_cell = rDerivativeData.DeltaCellVertex[i_node * TDim + i_dof]; 
                    for(std::size_t i_belong = 0; i_belong < 3; ++i_belong)
                    {
                        aux_RHS1 += N_decomp[i_belong] * row(local_delta_cell, i_belong);
                    }
                    array_1d<double, 3> aux_RHS2 = aux_RHS1;
                    
                    // Local contribution
                    const array_1d<double, 3>& aux_delta_node = LocalDeltaVertex( SlaveNormal,  delta_normal, i_dof, i_node, ConsiderNormalVariation, SlaveGeometry, MasterGeometry );
                    if (i_node < TNumNodes) aux_RHS1 -= N1[i_node] * aux_delta_node;
                    else aux_RHS2 -= N2[i_node - TNumNodes] * aux_delta_node;
                    
                    // We compute the delta coordinates 
                    array_1d<double, 2> aux_delta_coords1;
                    array_1d<double, 2> aux_delta_coords2;
                    DeltaPointLocalCoordinates(aux_delta_coords1, aux_RHS1, rVariables.DNDeSlave, SlaveGeometry);
                    DeltaPointLocalCoordinates(aux_delta_coords2, aux_RHS2, rVariables.DNDeMaster, MasterGeometry);
                    
                    // Now we can compute the delta shape functions
                    
                    auto& delta_n1 = rDerivativeData.DeltaN1[i_node * TDim + i_dof];
                    noalias(delta_n1) = (aux_delta_coords1[0] * column(DNDe1, 0) + aux_delta_coords1[1] * column(DNDe1, 1));
                    
                    auto& delta_n2 = rDerivativeData.DeltaN2[i_node * TDim + i_dof];
                    noalias(delta_n2) = (aux_delta_coords2[0] * column(DNDe2, 0) + aux_delta_coords2[1] * column(DNDe2, 1));
                }
            }
        }
    }
    
    /**
     * Calculates the increment of Phi (the dual shape function)
     * @param rVariables: The kinematic variables
     * @param rDerivativeData: The derivatives container
     */
    
    static inline void CalculateDeltaPhi(
        const GeneralVariables& rVariables,
        DerivativeDataType& rDerivativeData
        )
    {
        // Shape functions
        const VectorType& N1 = rVariables.NSlave;
        
        for (unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave)
        {
            for (unsigned int i_dim = 0; i_dim < TDim; ++i_dim)
            {
                const unsigned int i_dof = i_slave * TDim + i_dim;
                
                noalias(rDerivativeData.DeltaPhi[i_dof]) = prod(rDerivativeData.DeltaAe[i_dof], N1);;
            }
        }
    }
    
    /**
     * Returns a matrix with the increment of displacements, that can be used for compute the Jacobian reference (current) configuration
     * @param DeltaPosition: The matrix with the increment of displacements 
     * @param LocalCoordinates: The array containing the local coordinates of the exact integration segment
     */
    
    Matrix& CalculateDeltaPosition(
        Matrix& DeltaPosition,
        const GeometryType& ThisGeometry,
        const ConditionArrayType& LocalCoordinates
        )
    {
        KRATOS_TRY;

        DeltaPosition = ZeroMatrix(TDim, TDim);

        for ( unsigned int i_node = 0; i_node < TNumNodes; ++i_node )
        {
            const array_1d<double, 3 > delta_displacement = ThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT) - ThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT,1);
            
            for ( unsigned int j_node = 0; j_node < TDim; ++j_node )
            {
                Vector N;
                ThisGeometry.ShapeFunctionsValues( N, LocalCoordinates[j_node].Coordinates() );

                for ( unsigned int j_dim = 0; j_dim < TDim; ++j_dim )
                {
                    DeltaPosition(j_node, j_dim) += N[i_node] * delta_displacement[j_dim];
                }
            }
        }

        return DeltaPosition;

        KRATOS_CATCH( "" );
    }
    
    /**
     * Returns a matrix with the increment of displacements
     * @param DeltaPosition: The matrix with the increment of displacements 
     * @param ThisGeometry: The geometry considered 
     */
    
    static inline Matrix& CalculateDeltaPosition(
        Matrix& DeltaPosition,
        const GeometryType& ThisGeometry
        )
    {
        KRATOS_TRY;

        DeltaPosition = ZeroMatrix(TNumNodes, TDim);

        for ( unsigned int i_node = 0; i_node < TNumNodes; ++i_node )
        {
            const array_1d<double, 3 > delta_displacement = ThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT) - ThisGeometry[i_node].FastGetSolutionStepValue(DISPLACEMENT,1);
            
            for ( unsigned int i_dim = 0; i_dim < TDim; ++i_dim )
            {
                DeltaPosition(i_node, i_dim) += delta_displacement[i_dim];
            }
        }

        return DeltaPosition;

        KRATOS_CATCH( "" );
    }
    
    /**
     * Returns a vector with the increment of displacements
     */
    
    static inline void CalculateDeltaPosition(
        VectorType& DeltaPosition,
        const GeometryType& SlaveGeometry,
        const GeometryType& MasterGeometry,
        const unsigned int IndexNode
        )
    {
        KRATOS_TRY;

        if (IndexNode < TNumNodes)
        {
            DeltaPosition = SlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT) - SlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT,1);
        }
        else
        {
            const unsigned int index_master = IndexNode - TNumNodes;
            DeltaPosition = MasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT) - MasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT,1);
        }

        KRATOS_CATCH( "" );
    }
    
    /**
     * Returns a vector with the increment of displacements
     */
    
    static inline void CalculateDeltaPosition(
        VectorType& DeltaPosition,
        const GeometryType& SlaveGeometry,
        const GeometryType& MasterGeometry,
        const unsigned int IndexNode,
        const unsigned int iDoF
        )
    {
        KRATOS_TRY;

        DeltaPosition = ZeroVector(3);
        
        if (IndexNode < TNumNodes)
        {
            DeltaPosition[iDoF] = (SlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT) - SlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
        }
        else
        {
            const unsigned int index_master = IndexNode - TNumNodes;
            DeltaPosition[iDoF] = (MasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT) - MasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
        }

        KRATOS_CATCH( "" );
    }
    
    /**
     * Returns a double with the increment of displacements
     */
    
    static inline void CalculateDeltaPosition(
        double& DeltaPosition,
        const GeometryType& SlaveGeometry,
        const GeometryType& MasterGeometry,
        const unsigned int IndexNode,
        const unsigned int iDoF
        )
    {
        KRATOS_TRY;
    
        if (IndexNode < TNumNodes)
        {
            DeltaPosition = (SlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT) - SlaveGeometry[IndexNode].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
        }
        else
        {
            const unsigned int index_master = IndexNode - TNumNodes;
            DeltaPosition = (MasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT) - MasterGeometry[index_master].FastGetSolutionStepValue(DISPLACEMENT,1))[iDoF];
        }

        KRATOS_CATCH( "" );
    }
    
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
    
    /**
     * This method is used to compute the directional derivatives of the cell vertex (locally)
     * @param Normal: The normal of the slave surface
     * @param DeltaNormal: The derivative of the normal vector
     * @param N1: The shape function of the slave side
     * @param N2: The shape function of the master side
     * @param iDoF: The DoF computed index
     * @param iBelong: The belong (intersection, node, etc..) index
     * @param ConsiderNormalVariation: If the normal variation is considered
     * @param SlaveGeometry: The geometry of the slave side
     * @param MasterGeometry: The geometry of the master side
     * @param Coeff: The coefficient considered in proportion
     * @return The local vertex derivative
     */
    static inline array_1d<double, 3> LocalDeltaVertex(
        const array_1d<double, 3>& Normal,
        const bounded_matrix<double, TDim, TDim>& DeltaNormal,
        const unsigned int iDoF,
        const unsigned int iBelong,
        const bool ConsiderNormalVariation,
        const GeometryType& SlaveGeometry,
        const GeometryType& MasterGeometry,
        const double Coeff = 1.0
        )
    {
        // We create the auxiliar array
        array_1d<double, 3> aux_delta_vertex = ZeroVector(3);
        
        // This is the coefficient of the center contribution
        const double auxiliar_coeff = 1.0/static_cast<double>(TNumNodes);
        
        //  We initialize some values
        const array_1d<double, 3>& coords_center = SlaveGeometry.Center().Coordinates();
        const array_1d<double, 3>& coords_node = (iBelong < TNumNodes) ? SlaveGeometry[iBelong].Coordinates() : MasterGeometry[iBelong - TNumNodes].Coordinates();
    
        // The corresponding part to the nodal coordinates
        array_1d<double, 3> aux_der = ZeroVector(3);
        aux_der[iDoF] = 1.0;
        aux_delta_vertex += Coeff * aux_der;
        
        // The corresponding part to the normal
        const double coordsxdeltanormal = (ConsiderNormalVariation == true) ? inner_prod(coords_node - coords_center, column(DeltaNormal, iDoF)) : 0.0;
        
        const double factor_belong = (iBelong < TNumNodes) ? (1.0 - auxiliar_coeff) : 1.0; 
        const double deltacoordsxnormal =  factor_belong * Normal[iDoF];
        aux_delta_vertex += - Coeff * Normal * (deltacoordsxnormal + coordsxdeltanormal);
        
        // The corresponding part to delta normal
        const double coordsxnormal = - inner_prod(coords_node - coords_center, Normal);
        if (ConsiderNormalVariation == true) aux_delta_vertex += Coeff * coordsxnormal * trans(column(DeltaNormal, iDoF));
        
        return aux_delta_vertex;
    }
    
    /**
     * This method is used to compute the directional derivatives of the cell vertex (locally)
     * @param DeltaVertexMatrix: The whole delta vertex matrix
     * @param Normal: The normal of the slave surface
     * @param DeltaNormal: The derivative of the normal vector
     * @param N1: The shape function of the slave side
     * @param N2: The shape function of the master side
     * @param iDoF: The DoF computed index
     * @param iTriangle: The triangle point index
     * @param iBelong: The belong (intersection, node, etc..) index
     * @param ConsiderNormalVariation: If the normal variation is considered
     * @param SlaveGeometry: The geometry of the slave side
     * @param MasterGeometry: The geometry of the master side
     * @param Coeff: The coefficient considered in proportion
     */
    static inline void LocalDeltaVertex(
        bounded_matrix<double, 3, 3>& DeltaVertexMatrix,
        const array_1d<double, 3>& Normal,
        const bounded_matrix<double, TDim, TDim>& DeltaNormal,
        const unsigned int iDoF,
        const unsigned int iTriangle,
        const unsigned int iBelong,
        const bool ConsiderNormalVariation,
        const GeometryType& SlaveGeometry,
        const GeometryType& MasterGeometry,
        const double Coeff = 1.0
        )
    {
        noalias(row(DeltaVertexMatrix, iTriangle)) += LocalDeltaVertex( Normal,  DeltaNormal, iDoF, iBelong, ConsiderNormalVariation, SlaveGeometry, MasterGeometry, Coeff);
    }
    
    
    /**
     * This method computes the equivalent indexes to the auxiliar hash
     * @param AuxIndex: The auxiliar index to decompose
     * @param iBelongSlaveStart: The index of the first/slave segment and first node
     * @param iBelongSlaveEnd: The index of the first/slave segment and end node
     * @param iBelongMasterStart: The index of the second/master segment and first node
     * @param iBelongMasterEnd: The index of the second/master segment and end node
     */
    static inline void ConvertAuxHashIndex(
        const unsigned int AuxIndex,
        unsigned int& iBelongSlaveStart, 
        unsigned int& iBelongSlaveEnd, 
        unsigned int& iBelongMasterStart, 
        unsigned int& iBelongMasterEnd
        )
    {
        unsigned int index_to_decompose = AuxIndex - 2 * TNumNodes;
    
        iBelongMasterEnd = index_to_decompose/10000;
        index_to_decompose = std::fmod(index_to_decompose, 10000);
        iBelongMasterStart = index_to_decompose/1000;
        index_to_decompose = std::fmod(index_to_decompose, 1000);
        iBelongSlaveEnd = index_to_decompose/100;
        index_to_decompose = std::fmod(index_to_decompose, 100);
        iBelongSlaveStart = index_to_decompose/10;
    }
    
    /**
     * This method computes the increment of local coordinates
     * @param rResult: The solution obtained
     * @param DeltaPoint: The increment of position in the points
     * @param rVariables: The kinematic variables
     * @param ThisGeometry: The geometry considered
     */
    static inline void DeltaPointLocalCoordinates(
        array_1d<double, 2>& rResult,
        const array_1d<double, 3>& DeltaPoint,
        const MatrixType& rDNDe,
        const GeometryType& ThisGeometry
        )
    {
        bounded_matrix<double, 3, TNumNodes> X;
        for(unsigned int i = 0; i < TNumNodes; ++i)
        {
            X(0, i) = ThisGeometry[i].X();
            X(1, i) = ThisGeometry[i].Y();
            X(2, i) = ThisGeometry[i].Z();
        }
        
        const bounded_matrix<double, 3, 2> DN = prod(X, rDNDe);
        
        double det_j;
        const bounded_matrix<double, 2, 2> J = prod(trans(DN),DN);
        const bounded_matrix<double, 2, 2> invJ = MathUtils<double>::InvertMatrix<2>(J, det_j);
        
        const array_1d<double, 2> res = prod(trans(DN), DeltaPoint);
        noalias(rResult) = prod(invJ, res);

//         // NOTE: Testing, remove later
//         bounded_matrix<double, 3, 3> L;
//         L(0, 0) = DN(0, 0);
//         L(0, 1) = DN(0, 1);
//         L(1, 0) = DN(1, 0);
//         L(1, 1) = DN(1, 1);
//         L(2, 0) = DN(2, 0);
//         L(2, 1) = DN(2, 1);
//         L(0, 2) = 0.0; // NOTE: Use here the normal
//         L(1, 2) = 0.0;
//         L(2, 2) = 1.0;
//         double det_L;
//         const bounded_matrix<double, 3, 3> invL = MathUtils<double>::InvertMatrix<3>(L, det_L);
//         array_1d<double, 3> aux = prod(invL, DeltaPoint);
//         rResult[0] = aux[0];
//         rResult[1] = aux[1];
    }
    
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
};// class DerivativesUtilities

}
#endif /* KRATOS_DERIVATIVES_UTILITIES defined */
 
