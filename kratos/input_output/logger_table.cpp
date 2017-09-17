 
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

// System includes
#include <sstream>

// External includes 

// Project includes
#include "input_output/logger_table.h"

namespace Kratos
{
    std::string LoggerTable::Info() const
    {
        return "LoggerTable";
    }

    /// Print information about this object.
    void LoggerTable::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    void LoggerTable::PrintData(std::ostream& rOStream) const
    {
        rOStream << mTable;
    }

    /// char stream function
    LoggerTable& LoggerTable::operator << (const char * rString)
    {
        mTable << rString;

        return *this;
    }

    LoggerTable& LoggerTable::operator << (std::ostream& (*pf)(std::ostream&))
    {
        std::stringstream buffer;
        pf(buffer);

        mTable << buffer.str();

        return *this;
    }

    LoggerTable& LoggerTable::operator << (Severity const& TheSeverity)
    {
        mSeverity = TheSeverity;

        return *this;
    }

    LoggerTable& LoggerTable::operator << (Category const& TheCategory) 
    {
        mCategory = TheCategory;

        return *this;
    }


    /// output stream function
    std::ostream& operator << (std::ostream& rOStream,
        const LoggerTable& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

  
}  // namespace Kratos.
