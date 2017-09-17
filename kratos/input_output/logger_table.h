//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

#if !defined(KRATOS_LOGGER_TABLE_H_INCLUDED )
#define  KRATOS_LOGGER_TABLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <map>
#include <chrono>

// External includes

// Project includes
#include "includes/kratos_export_api.h"
#include "includes/table_stream.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    
    typedef TableStream TableStreamType;
    
///@}
///@name  Enum's
///@{
    
///@}
///@name  Functions
///@{
    
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// LoggerTable class holdes message and the properties of the message in table format.
/** LoggerTable holds the origin of the message, severity, level and 
    the category of it.
    Most of the methods are defined in header to be inlined in order to
    increase the performance. 
*/
class KRATOS_API(KRATOS_CORE) LoggerTable
{
public:
    ///@name Type Definitions
    ///@{

    using TimePointType = std::chrono::steady_clock::time_point;

    ///@}
    ///@name Enums
    ///@{

    enum class Severity 
    {
        ERROR,
        WARNING,
        INFO,
        DETAIL,
        DEBUG,
        TRACE,
    };

    enum class Category 
    {
        STATUS,
        CRITICAL,
        STATISTICS,
        PROFILING
    };

    ///@}
    ///@name Life Cycle
    ///@{

    LoggerTable() : mMessage(), mLevel(1), mSeverity(Severity::INFO), mCategory(Category::STATUS) {}

    LoggerTable(std::string const& TheMessage) 
        : mMessage(TheMessage), mLevel(1), mSeverity(Severity::INFO), mCategory(Category::STATUS) {}

    LoggerTable(LoggerTable const& Other) 
        : mMessage(Other.mMessage), mLevel(Other.mLevel), mSeverity(Other.mSeverity), mCategory(Other.mCategory) {}

    /// Destructor.
    virtual ~LoggerTable() {}


    ///@}
    ///@name Operators
    ///@{

    LoggerTable& operator=(LoggerTable const& Other) 
    {
        mMessage = Other.mMessage;
        mLevel = Other.mLevel;
        mSeverity = Other.mSeverity;
        mCategory = Other.mCategory;

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{

    void SetMessage(std::string const& TheMessage) +
    {
        mMessage = TheMessage;
    }

    std::string const& GetMessage() const 
    {
        return mMessage;
    }

    void SetLevel(std::size_t TheLevel) 
    {
        mLevel = TheLevel;
    }

    std::size_t GetLevel() const 
    {
        return mLevel;
    }

    void SetSeverity(Severity const& TheSeverity) 
    {
        mSeverity = TheSeverity;
    }

    Severity GetSeverity() const 
    {
        return mSeverity;
    }

    void SetCategory(Category const& TheCategory) 
    {
        mCategory = TheCategory;
    }

    Category GetCategory() const 
    {
        return mCategory;
    }

    void SetTime() 
    {
        mTime = std::chrono::steady_clock::now();
    }

    TimePointType const& GetTime() const 
    {
        return mTime;
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;

    /// string stream function
    template<class StreamValueType>
    LoggerTable& operator << (StreamValueType const& rValue)
    {
        std::stringstream buffer;
        buffer << rValue;

        mMessage.append(buffer.str());

        return *this;
    }

    /// Manipulator stream function
    LoggerTable& operator << (std::ostream& (*pf)(std::ostream&));

    /// char stream function
    LoggerTable& operator << (const char * rString);

    /// Severity stream function
    LoggerTable& operator << (Severity const& TheSeverity);

    /// Category stream function
    LoggerTable& operator << (Category const& TheCategory);

    ///@}

private:
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    std::string mMessage;
    std::size_t mLevel;
    Severity mSeverity;
    Category mCategory;
    TimePointType mTime;
    TableStreamType mTable;

    ///@}
}; // Class LoggerTable

///@}

///@name Input and output
///@{

/// output stream function
std::ostream& operator << (std::ostream& rOStream,
    const LoggerTable& rThis);

///@}
///@name macros
///@{

///@}

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_LOGGER_TABLE_H_INCLUDED  defined
