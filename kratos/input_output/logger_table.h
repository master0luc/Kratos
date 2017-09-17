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

    typedef TableStream TableStreamType;
    
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
    
    // Class Constructor
    
    /**
     * The default constructor
     */
    
    LoggerTable(const bool UseBoldFont = true) 
        : mTable(&std::cout, "|", UseBoldFont), 
          mLevel(1), 
          mSeverity(Severity::INFO), 
          mCategory(Category::STATUS)
    {
        // TODO: Add something if necessary
    }

    /**
     * The constructor with table
     */
    
    LoggerTable(
        TableStream const& TheTable,
        const bool UseBoldFont = true
        ) 
        : mTable(TheTable), 
          mLevel(1), 
          mSeverity(Severity::INFO), 
          mCategory(Category::STATUS)
    {
        // TODO: Add something if necessary
    }

    /**
     * The copy constructor
     */
    LoggerTable(LoggerTable const& Other) 
        : mTable(Other.mTable), 
          mLevel(Other.mLevel), 
          mSeverity(Other.mSeverity), 
          mCategory(Other.mCategory)
    {
        // TODO: Add something if necessary
    }

    /// Destructor.
    virtual ~LoggerTable() {}

    ///@}
    ///@name Operators
    ///@{

    /**
     * The assigning operator 
     * @param Other: The logger to copy
     */
    LoggerTable& operator=(LoggerTable const& Other) 
    {
        mTable = Other.mTable;
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

    /**
     * This function sets the table
     * @param TheTable: The table to set
     */
    void GetTable(TableStream const& TheTable)
    {
        mTable = TheTable;
    }

    /**
     * This function recovers the table
     * @return mMessage: The table to set
     */
    TableStream const& GetMessage() const 
    {
        return mTable;
    }

    /**
     * This function sets the level
     * @param TheLevel: The level
     */
    void SetLevel(std::size_t TheLevel) 
    {
        mLevel = TheLevel;
    }

    /**
     * This function gets the level
     * @return mLevel: The level
     */
    std::size_t GetLevel() const 
    {
        return mLevel;
    }

    /**
     * This function sets the severity
     * @param TheSeverity: The severity
     */
    void SetSeverity(Severity const& TheSeverity) 
    {
        mSeverity = TheSeverity;
    }

    /**
     * This function gets the severity
     * @return mSeverity: The severity
     */
    Severity GetSeverity() const 
    {
        return mSeverity;
    }

    /**
     * This function sets the category
     * @param TheCategory: The category
     */
    void SetCategory(Category const& TheCategory) 
    {
        mCategory = TheCategory;
    }

    /**
     * This function gets the category
     * @return mCategory: The category
     */
    Category GetCategory() const 
    {
        return mCategory;
    }

    /**
     * This function sets the time
     */
    void SetTime() 
    {
        mTime = std::chrono::steady_clock::now();
    }

    /**
     * This function gets the time
     * @return mTime: The time
     */
    TimePointType const& GetTime() const 
    {
        return mTime;
    }

    /**
     * This function prints the header of the table
     */
    void PrintTableHeader()
    {
        mTable.PrintHeader();
    }
    
    /**
     * This function prints the footer of the table
     */
    void PrintTableFooter()
    {
        mTable.PrintFooter();
    }
    
    /**
     * It adds a column to the table
     * @param ThisName: The name of the variable
     * @param ThisSpaces: The number of spaces to consider
     */
    void AddColumnToTable(        
        std::string ThisName, 
        const unsigned int ThisSpaces = 10
        )
    {
        mTable.AddColumn(ThisName, ThisSpaces);
    }
    
    /**
     * This function sets if the table uses the bold UseBoldFont
     * @param UseBoldFont: If the bold font is used
     */
    void SetBoldTable(const bool UseBoldFont) 
    {
        mTable.SetBold(UseBoldFont);
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
    LoggerTable& operator << (StreamValueType const& rValue);

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

    TableStreamType mTable; // The table used as ouput
    std::size_t mLevel;     // The level of the message
    Severity mSeverity;     // The severity of the message
    Category mCategory;     // The category of the message
    TimePointType mTime;    // The time

    ///@}
}; // Class LoggerTable

///@}

///@name Input and output
///@{

/// output stream function
std::ostream& operator << (std::ostream& rOStream,
    const LoggerTable& rThis);

///@}
///@name Kratos Macros
///@{

///@}

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_LOGGER_TABLE_H_INCLUDED  defined
