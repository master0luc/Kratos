// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef RECONSTRUCTION_CONDITION_H
#define RECONSTRUCTION_CONDITION_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "../basic_nurbs_brep_handling/patch.h"

// ==============================================================================

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

/// Short class definition.
/** Detail class definition.

*/

class ReconstructionCondition
{
public:
    ///@name Type Definitions
    ///@{
        
    /// Pointer definition of ReconstructionCondition
    KRATOS_CLASS_POINTER_DEFINITION(ReconstructionCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ReconstructionCondition()
    {
    }

    /// Destructor.
    virtual ~ReconstructionCondition()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    virtual void FlagControlPointsRelevantForReconstruction() = 0;
    
    // --------------------------------------------------------------------------
    virtual void Initialize() = 0;    
    
    // --------------------------------------------------------------------------
    virtual void ComputeAndAddLHSContribution( CompressedMatrix& LHS ) = 0;

    // --------------------------------------------------------------------------
    virtual void ComputeAndAddRHSContribution( Vector& RHS ) = 0;

    // --------------------------------------------------------------------------
    virtual void EvaluateProjection( double &residual_distance, bool &is_outside_trimmed_area) = 0;

    // ==============================================================================

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ReconstructionCondition";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "ReconstructionCondition";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
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
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //      ReconstructionCondition& operator=(ReconstructionCondition const& rOther);

    /// Copy constructor.
    //      ReconstructionCondition(ReconstructionCondition const& rOther);

    ///@}

}; // Class ReconstructionCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // RECONSTRUCTION_CONDITION_H
