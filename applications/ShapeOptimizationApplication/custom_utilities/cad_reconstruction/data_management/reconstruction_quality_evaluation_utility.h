// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef RECONSTRUCTION_QUALITY_EVALUATION_UTILITY_H
#define RECONSTRUCTION_QUALITY_EVALUATION_UTILITY_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <string>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "reconstruction_data_base.h"

// ==============================================================================

namespace Kratos
{
class ReconstructionQualityEvaluationUtility
{
public:
    ///@name Type Definitions
    ///@{

    
    /// Pointer definition of ReconstructionQualityEvaluationUtility
    KRATOS_CLASS_POINTER_DEFINITION(ReconstructionQualityEvaluationUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ReconstructionQualityEvaluationUtility( ReconstructionDataBase& reconstruction_data_base):
     mrReconstructionDataBase( reconstruction_data_base )
    {      
    }

    /// Destructor.
    virtual ~ReconstructionQualityEvaluationUtility()
    {
    }

    // --------------------------------------------------------------------------
    void EvaluateGlobalQuality()
    {
      std::cout << "EvaluateGlobalQuality CALLED" << std::endl;
    }

    // --------------------------------------------------------------------------
    void EvaluateDisplacementCoupling()
    {
      std::cout << "EvaluateDisplacementCoupling CALLED" << std::endl;
    }

    // --------------------------------------------------------------------------
    void EvaluateRotationCoupling()
    {
      std::cout << "EvaluateRotationCoupling CALLED" << std::endl;
    }

    // ==============================================================================

    /// Turn back information as a string.
    virtual std::string Info() const
    {
		return "ReconstructionQualityEvaluationUtility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
		rOStream << "ReconstructionQualityEvaluationUtility";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }

private:

    ReconstructionDataBase& mrReconstructionDataBase;

    /// Assignment operator.
    //      ReconstructionQualityEvaluationUtility& operator=(ReconstructionQualityEvaluationUtility const& rOther);

    /// Copy constructor.
    //      ReconstructionQualityEvaluationUtility(ReconstructionQualityEvaluationUtility const& rOther);

}; // Class ReconstructionQualityEvaluationUtility
} // namespace Kratos.

#endif // RECONSTRUCTION_QUALITY_EVALUATION_UTILITY_H
