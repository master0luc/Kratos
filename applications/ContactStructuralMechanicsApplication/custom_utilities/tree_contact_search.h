// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
// 

#if !defined(KRATOS_TREE_CONTACT_SEARCH_H_INCLUDED )
#define  KRATOS_TREE_CONTACT_SEARCH_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"

/* Custom includes*/
#include "custom_includes/point_item.h"
#include "custom_conditions/paired_condition.h"

/* Custom utilities*/
#include "custom_utilities/search_utilities.h"

/* Tree structures */
// #include "spatial_containers/bounding_volume_tree.h" // k-DOP
#include "spatial_containers/spatial_containers.h" // kd-tree 

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
    
    enum SearchTreeType {KdtreeInRadius = 0, KdtreeInBox = 1, Kdop = 2};
    
    enum CheckResult {Fail = 0, AlreadyInTheMap = 1, OK = 2};
    
///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
    
/** \brief TreeContactSearch
 * This utilitiy has as objective to create the contact conditions.
 * The conditions that can be created are Mortar conditions (or segment to segment) conditions: The created conditions will be between two segments
 * The utility employs the projection.h from MeshingApplication, which works internally using a kd-tree 
 */
template<unsigned int TDim>
class TreeContactSearch
{
public:
    ///@name Type Definitions
    ///@{
    
    // General type definitions
    typedef ModelPart::NodesContainerType                    NodesArrayType;
    typedef ModelPart::ConditionsContainerType          ConditionsArrayType;
    typedef Node<3>                                                NodeType;
    typedef Geometry<NodeType>                                 GeometryType;
    
    // Type definitions for the tree
    typedef PointItem                                             PointType;
    typedef PointType::Pointer                             PointTypePointer;
    typedef std::vector<PointTypePointer>                       PointVector;
    typedef PointVector::iterator                             PointIterator;
    typedef std::vector<double>                              DistanceVector;
    typedef DistanceVector::iterator                       DistanceIterator;
    
    // KDtree definitions
    typedef Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;

    /// Pointer definition of TreeContactSearch
    KRATOS_CLASS_POINTER_DEFINITION( TreeContactSearch );
      
    ///@}
    ///@name Life Cycle
    ///@{

    // Class Constructor
    
    /**
     * The constructor of the search utility uses the following inputs:
     * @param rMainModelPart The model part to be considered
     * @param ThisParameters The condiguration parameters, it includes:
     *                       - The allocation considered in the search
     *                       - The factor considered to check if active or not
     *                       - The integration order considered
     *                       - The size of the bucket
     *                       - The proportion increased of the Radius/Bounding-box volume for the search
     *                       - TypeSearch: 0 means search in radius, 1 means search in box // TODO: Add more types of bounding boxes, as kdops, look bounding_volume_tree.h
     * NOTE: Use an InterfacePreprocess object to create such a model part from a regular one:
     * InterfaceMapper = InterfacePreprocess()
     * InterfacePart = InterfaceMapper.GenerateInterfacePart(Complete_Model_Part)
     */
    
    TreeContactSearch( ModelPart& rMainModelPart, Parameters ThisParameters =  Parameters(R"({})") ); 
    
    virtual ~TreeContactSearch()= default;;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * This function initializes the ALM frictionless mortar conditions already created 
     */
    
    void InitializeMortarConditions();
    
    /**
     * This function clears the mortar conditions already created 
     */
    
    void TotalClearScalarMortarConditions();
    
    /**
     * This function clears the mortar conditions already created 
     */
    
    void TotalClearComponentsMortarConditions();
    
    /**
     * This function clears the ALM frictionless mortar conditions already created 
     */
    
    void TotalClearALMFrictionlessMortarConditions();
    
    /**
     * This function clears the mortar conditions already created (scalar version)
     */
    
    void PartialClearScalarMortarConditions();
    
    /**
     * This function clears the mortar conditions already created (components version)
     */
        
    void PartialClearComponentsMortarConditions();
    
    /**
     * This function clears the ALM frictionless mortar conditions already created 
     */
    
    void PartialClearALMFrictionlessMortarConditions();
      
    /**
     * This function creates a lists  points ready for the Mortar method
     */
    
    void CreatePointListMortar();

    /**
     * This function updates a lists  points ready for the Mortar method
     */
    
    void UpdatePointListMortar();

    /**
     * This function has as pourpose to find potential contact conditions and fill the mortar conditions with the necessary pointers
     */
    
    void UpdateMortarConditions();
    
    /**
     * This function has as pourpose to clean the existing pairs
     */
    
    void CleanMortarConditions();
    
    /**
     * It checks the current mortar conditions
     */
    
    void CheckMortarConditions();
    
    /**
     * It sets if the search is inverted
     */
    
    void InvertSearch();
    
    /**
     * This resets the contact operators
     */
        
    void ResetContactOperators();
    
    /**
     * This resets all the contact operators
     */
        
    void TotalResetContactOperators();
    
    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /************************************ GET INFO *************************************/
    /***********************************************************************************/
    
    virtual std::string Info() const
    {
        return "TreeContactSearch";
    }

    /************************************ PRINT INFO ***********************************/
    /***********************************************************************************/
    
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{
            
    /**
     * It check the conditions if they are correctly detected
     * @return ConditionPointers1: A vector containing the pointers to the conditions 
     * @param pCond1 The pointer to the condition in the destination model part
     * @param pCond2 The pointer to the condition in the destination model part  
     * @param InvertedSearch If the search is inverted
     */
    
    static inline CheckResult CheckCondition(
        IndexMap::Pointer IndexesMap,
        const Condition::Pointer pCond1,
        const Condition::Pointer pCond2,
        const bool InvertedSearch = false
        );
    
    /**
     * This method checks the potential pairing between two conditions/geometries
     */
    inline void CheckPotentialPairing(
        ModelPart& ComputingModelPart,
        std::size_t& ConditionId,
        Condition::Pointer pCondSlave,
        PointVector& PointsFound,
        const unsigned int NumberOfPointsFound,
        IndexMap::Pointer IndexesMap
        );
    
    /**
     * This method checks the potential pairing between two conditions/geometries
     */
    template<unsigned int TNumNodes>
    inline void AuxiliarCheckPotentialPairing(
        ModelPart& ComputingModelPart,
        std::size_t& ConditionId,
        Condition::Pointer pCondSlave,
        PointVector& PointsFound,
        const unsigned int NumberOfPointsFound,
        IndexMap::Pointer IndexesMap
        )
    {        
        // Some initial parameters
        const double active_check_factor = mrMainModelPart.GetProcessInfo()[ACTIVE_CHECK_FACTOR];
        const array_1d<double, 3>& contact_normal_origin = pCondSlave->GetValue(NORMAL);
        GeometryType& slave_geometry = pCondSlave->GetGeometry();
        const double active_check_length = slave_geometry.Length() * active_check_factor;
        
        // We update the base condition
        if (mCreateAuxiliarConditions && slave_geometry.GetGeometryType() != mGeometryType)
        {
            mConditionName = mThisParameters["mConditionName"].GetString(); 
            mConditionName.append("Condition"); 
            mConditionName.append(std::to_string(TDim));
            mConditionName.append("D"); 
            mConditionName.append(std::to_string(TNumNodes)); 
            mConditionName.append("N"); 
            mConditionName.append(mThisParameters["final_string"].GetString());
        }
        
        MortarKinematicVariables<TNumNodes> rVariables;
        MortarOperator<TNumNodes> rThisMortarConditionMatrices;
        ExactMortarIntegrationUtility<TDim, TNumNodes> integration_utility = ExactMortarIntegrationUtility<TDim, TNumNodes>(TDim);
        
        for(unsigned int i_pair = 0; i_pair < NumberOfPointsFound; ++i_pair)
        {   
            bool condition_is_active = false;
            
            Condition::Pointer p_cond_master = PointsFound[i_pair]->GetCondition(); // MASTER
            GeometryType& master_geometry = p_cond_master->GetGeometry();
            const array_1d<double, 3>& master_normal = p_cond_master->GetValue(NORMAL); 
                                
            const CheckResult condition_checked_right = CheckCondition(IndexesMap, pCondSlave, p_cond_master, mInvertedSearch);
            
            if (condition_checked_right == OK)
            {   
                condition_is_active = SearchUtilities::CheckExactIntegration<TDim, TNumNodes, true>(rVariables, rThisMortarConditionMatrices, integration_utility, slave_geometry, master_geometry, contact_normal_origin, master_normal, active_check_length);
                
                // If condition is active we add
                if (condition_is_active == true) AddPairing(ComputingModelPart, ConditionId, pCondSlave, p_cond_master, IndexesMap);
            }
            else if (condition_checked_right == AlreadyInTheMap)
            {
                condition_is_active = SearchUtilities::CheckExactIntegration<TDim, TNumNodes, false>(rVariables, rThisMortarConditionMatrices, integration_utility, slave_geometry, master_geometry, contact_normal_origin, master_normal, active_check_length);
                
                if (condition_is_active == false) 
                {
                    const std::size_t index_slave = p_cond_master->Id();
                    const std::size_t index_master = IndexesMap->GetAuxiliarIndex(index_slave);
                    IndexesMap->RemoveId(p_cond_master->Id());
                    if (mCreateAuxiliarConditions) ComputingModelPart.pGetCondition(index_master)->Set(TO_ERASE, true);
                }
            }
        }
    }
    
    /**
     * This method checks the potential pairing between two conditions/geometries
     */
    inline void AddPotentialPairing(
        ModelPart& ComputingModelPart,
        std::size_t& ConditionId,
        Condition::Pointer pCondSlave,
        PointVector& PointsFound,
        const unsigned int NumberOfPointsFound,
        IndexMap::Pointer IndexesMap
        );
    
    /**
     * This method add a new pair to the computing model part
     * @param ComputingModelPart The modelpart  used in the assemble of the system
     * @param ConditionId The ID of the new condition to be created
     * @param pCondSlave The pointer to the slave condition
     * @param pCondMaster The pointer to the master condition
     * @param IndexesMap The map of indexes considered
     */
    inline void AddPairing(
        ModelPart& ComputingModelPart,
        std::size_t& ConditionId,
        Condition::Pointer pCondSlave,
        Condition::Pointer pCondMaster,
        IndexMap::Pointer IndexesMap
        );
    
    /**
     * This method checks the current pairing between two conditions/geometries
     * @param pCondSlave The pointer to the slave condition
     * @param IndexesMap The map of indexes considered
     */
    inline void CheckCurrentPairing(
        Condition::Pointer pCondSlave,
        IndexMap::Pointer IndexesMap
        );
    
    /**
     * This method checks the current pairing between two conditions/geometries
     * @param pCondSlave The pointer to the slave condition
     * @param IndexesMap The map of indexes considered
     */
    template<unsigned int TNumNodes>
    inline void AuxiliarCheckCurrentPairing(
        Condition::Pointer pCondSlave,
        IndexMap::Pointer IndexesMap
        )
    {
        // Some initial parameters
        const double active_check_factor = mrMainModelPart.GetProcessInfo()[ACTIVE_CHECK_FACTOR];
        const array_1d<double, 3>& contact_normal_origin = pCondSlave->GetValue(NORMAL);
        const GeometryType& this_geometry = pCondSlave->GetGeometry();
        const double active_check_length = this_geometry.Length() * active_check_factor;
        
        MortarKinematicVariables<TNumNodes> rVariables;
        MortarOperator<TNumNodes> rThisMortarConditionMatrices;
        ExactMortarIntegrationUtility<TDim, TNumNodes> integration_utility = ExactMortarIntegrationUtility<TDim, TNumNodes>(TDim);
        
        for (auto it_pair = IndexesMap->begin(); it_pair != IndexesMap->end(); ++it_pair )
        {                                   
            Condition::Pointer p_cond_master = mrMainModelPart.pGetCondition(it_pair->first);
            const array_1d<double, 3>& master_normal = p_cond_master->GetValue(NORMAL); 

            const bool condition_is_active = SearchUtilities::CheckExactIntegration<TDim, TNumNodes, false>(rVariables, rThisMortarConditionMatrices, integration_utility, pCondSlave->GetGeometry(), p_cond_master->GetGeometry(), contact_normal_origin, master_normal, active_check_length);
            
            if (condition_is_active == false) 
            {
                mrMainModelPart.pGetCondition(it_pair->second)->Set(TO_ERASE, true);
                IndexesMap->RemoveId(p_cond_master->Id());
            }
        }
    }
    
    /**
     * This method checks the current pairing between two conditions/geometries
     * @param pCondSlave The pointer to the slave condition
     * @param IndexesMap The map of indexes considered
     */
    inline void CheckAllActivePairing(
        Condition::Pointer pCondSlave,
        IndexMap::Pointer IndexesMap
        );
    
    /**  
     * Calculates the minimal distance between one node and its center 
     * @return The radius of the geometry 
     */ 
    
    static inline double Radius(GeometryType& ThisGeometry);
    
    /**
     * This converts the framework string to an enum
     * @param str The string
     * @return SearchTreeType: The equivalent enum
     */
    
    SearchTreeType ConvertSearchTree(const std::string& str);
    
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
  
    ModelPart& mrMainModelPart;                      // The main model part
    GeometryData::KratosGeometryType mGeometryType;  // The current size of the geometry considered
    Parameters mThisParameters;                      // The configuration parameters
    bool mCheckGap;                                  // If the gap is checked during the search
    bool mInvertedSearch;                            // The search will be done inverting the way master and slave/master is assigned
    std::string mConditionName;                      // The name of the condition to be created
    bool mCreateAuxiliarConditions;                  // If the auxiliar conditions are created or not
    PointVector mPointListDestination;               // A list that contents the all the points (from nodes) from the modelpart 

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
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class TreeContactSearch

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

// /****************************** INPUT STREAM FUNCTION ******************************/
// /***********************************************************************************/
// 
// template<class TPointType, class TPointerType>
// inline std::istream& operator >> (std::istream& rIStream,
//                                   TreeContactSearch& rThis);
// 
// /***************************** OUTPUT STREAM FUNCTION ******************************/
// /***********************************************************************************/
// 
// template<class TPointType, class TPointerType>
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const TreeContactSearch& rThis)
// {
//     return rOStream;
// }

///@}

}  // namespace Kratos.

#endif // KRATOS_TREE_CONTACT_SEARCH_H_INCLUDED  defined 
