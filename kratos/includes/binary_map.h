//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_BINARY_MAP_H_INCLUDED)
#define KRATOS_BINARY_MAP_H_INCLUDED

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "define.h"

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
    
///@}
///@name Kratos Classes
///@{

/** @brief Custom Point container to be used by the mapper
*/
template<typename T1, typename T2> 
class BinaryMap : public std::unordered_map<T1, T2>
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::unordered_map<T1, T2> BaseType;
    
    /// Counted pointer of BinaryMap
    KRATOS_CLASS_POINTER_DEFINITION( BinaryMap );

    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Default constructors
    BinaryMap(){}

    /// Destructor
    virtual ~BinaryMap(){}
    
    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
    
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

    template <class Pair>
    struct Reverse
    {
        typedef typename Pair::first_type  second_type;
        typedef typename Pair::second_type first_type;
        second_type second;
        first_type first;
    };

    template <class Pair>
    Reverse<Pair> & mutate(Pair & p)
    {
        return reinterpret_cast<Reverse<Pair> &>(p);
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
}; // Class BinaryMap 
} // namespace Kratos 
