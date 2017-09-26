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

#if !defined(KRATOS_MORTAR_AND_CRITERIA_H)
#define  KRATOS_MORTAR_AND_CRITERIA_H

/* System includes */
#include <algorithm>    // std::sort

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/table_stream_utility.h"
#include "solving_strategies/convergencecriterias/and_criteria.h"
#if !defined(_WIN32)
    #include "utilities/color_utilities.h"
#endif
#include "utilities/svd_utils.h"
#if defined(INCLUDE_FEAST_CONTACT)
    #include "../ExternalSolversApplication/external_includes/feast_solver.h"
#endif
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

/** Short class definition.
Detail class definition.

\URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

\URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

\URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

\URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}

*/

template<class TSparseSpace,
         class TDenseSpace
         >
class MortarAndConvergenceCriteria : public And_Criteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /** Counted pointer of MortarAndConvergenceCriteria */

    KRATOS_CLASS_POINTER_DEFINITION(MortarAndConvergenceCriteria );

    typedef And_Criteria< TSparseSpace, TDenseSpace >            BaseType;

    typedef TSparseSpace                                  SparseSpaceType;

    typedef typename TSparseSpace::MatrixType            SparseMatrixType;

    typedef typename TSparseSpace::VectorType            SparseVectorType;

    typedef typename TDenseSpace::MatrixType              DenseMatrixType;

    typedef typename TDenseSpace::VectorType              DenseVectorType;
    
    typedef typename BaseType::TDataType                        TDataType;

    typedef typename BaseType::DofsArrayType                DofsArrayType;

    typedef typename BaseType::TSystemMatrixType        TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType        TSystemVectorType;

    typedef boost::shared_ptr<TableStreamUtility> TablePrinterPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /** 
     * Constructor.
     */
    MortarAndConvergenceCriteria(
        typename ConvergenceCriteria < TSparseSpace, TDenseSpace >::Pointer pFirstCriterion,
        typename ConvergenceCriteria < TSparseSpace, TDenseSpace >::Pointer pSecondCriterion,
        TablePrinterPointerType pTable = nullptr,
        const bool PrintingOutput = false,
        const bool ComputeConditionNumber = false
        )
        :And_Criteria< TSparseSpace, TDenseSpace >(pFirstCriterion, pSecondCriterion),
        mpTable(pTable),
        mPrintingOutput(PrintingOutput),
        mComputeConditionNumber(ComputeConditionNumber),
        mTableIsInitialized(false)
    {
    }

    /**
     * Copy constructor.
     */
    MortarAndConvergenceCriteria(MortarAndConvergenceCriteria const& rOther)
      :BaseType(rOther)
      ,mpTable(rOther.mpTable)
      ,mPrintingOutput(rOther.mPrintingOutput)
      ,mTableIsInitialized(rOther.mTableIsInitialized)
     {
         BaseType::mpFirstCriterion   =  rOther.mpFirstCriterion;
         BaseType::mpSecondCriterion  =  rOther.mpSecondCriterion;      
     }

    /** Destructor.
    */
    ~MortarAndConvergenceCriteria () override = default;

    ///@}
    ///@name Operators
    ///@{

    /**
     * Criteria that need to be called after getting the solution
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */

    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
        {
            if (mpTable != nullptr)
            {
                const unsigned int nl_iteration = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
                mpTable->AddToRow<unsigned int>(nl_iteration);
            }
        }
        
        bool criterion_result = BaseType::PostCriteria(rModelPart,rDofSet,A,Dx,b);
        
        if (mComputeConditionNumber == true)
        {
        #if defined(INCLUDE_FEAST_CONTACT)
            typedef FEASTSolver<TSparseSpace, TDenseSpace> FEASTSolverType;
            
            Parameters default_params(R"(
            {
                "solver_type": "FEAST",
                "print_feast_output": false,
                "perform_stochastic_estimate": true,
                "solve_eigenvalue_problem": true,
                "lambda_min": 0.0,
                "lambda_max": 1.0,
                "echo_level": 0,
                "number_of_eigenvalues": 0,
                "search_dimension": 10,
                "linear_solver_settings": {
                    "solver_type": "skyline_lu"
                }
            })");
            
            const std::size_t size = A.size1();
            
            const double normA = SparseSpaceType::TwoNorm(A);
            default_params["lambda_max"].SetDouble(normA);
            default_params["lambda_min"].SetDouble(-normA);
            default_params["number_of_eigenvalues"].SetInt(size * 2/3 - 1);
            default_params["search_dimension"].SetInt(3/2 * size + 1);
            SparseMatrixType copy_matrix = A;
            SparseMatrixType identity_matrix = IdentityMatrix(size, size);
            
            // Create the auxilary eigen system
            DenseMatrixType eigen_vectors;
            DenseVectorType eigen_values;
            
            // Create the FEAST solver
            FEASTSolverType FEASTSolver(boost::make_shared<Parameters>(default_params));
            
            // Solve the problem
            FEASTSolver.Solve(copy_matrix, identity_matrix, eigen_values, eigen_vectors);
            
            // Size of the eigen values vector
            const int dim_eigen_values = eigen_values.size();
            
            // We get the moduli of the eigen values
            #pragma omp parallel for 
            for (int i = 0; i < dim_eigen_values; i++)
            {
                eigen_values[i] = std::abs(eigen_values[i]);
            }
            
            // Now we sort the vector
            std::sort(eigen_values.begin(), eigen_values.end());
            
            // We compute the eigen value
            double condition_number = 0.0;
            if (dim_eigen_values > 0) condition_number = eigen_values[dim_eigen_values - 1]/eigen_values[0];
        #else
            const double condition_number = SVDUtils<double>::SVDConditionNumber(A);
        #endif
            
            if (mpTable != nullptr)
            {
                std::cout.precision(4);
                auto& Table = mpTable->GetTable();
                Table  << condition_number;
            }
            else
            {
                if (mPrintingOutput == false)
                {
                #if !defined(_WIN32)
                    std::cout << "\n" << BOLDFONT("CONDITION NUMBER:") << "\t " << std::scientific << condition_number << std::endl;
                #else
                    std::cout << "\n" << "CONDITION NUMBER:" << "\t" << std::scientific << condition_number << std::endl;
                #endif
                }
                else
                {
                    std::cout << "\n" << "CONDITION NUMBER:" << "\t" << std::scientific << condition_number << std::endl;
                }
            }
        }
        
        if (criterion_result == true && rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
        {
            if (mpTable != nullptr)
            {
                mpTable->PrintFooter();
            }
        }
        
        return criterion_result;
    }

    /**
     * This function initialize the convergence criteria
     * @param rModelPart: The model part of interest
     */ 
    
    void Initialize(ModelPart& rModelPart) override
    {
        if (mpTable != nullptr && mTableIsInitialized == false)
        {
            auto& table = mpTable->GetTable();
            table.AddColumn("ITER", 4);
            mTableIsInitialized = true;
        }
        
        BaseType::Initialize(rModelPart);
         
        if (mpTable != nullptr && mComputeConditionNumber == true)
        {
            auto& table = mpTable->GetTable();
            table.AddColumn("COND.NUM.", 10);
        }
    }

    /**
     * This function initializes the solution step
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
        {
            std::cout.precision(4);
            if (mPrintingOutput == false)
            {
            #if !defined(_WIN32)
                std::cout << "\n\n" << BOLDFONT("CONVERGENCE CHECK") << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tTIME: " << std::scientific << rModelPart.GetProcessInfo()[TIME] << "\tDELTA TIME: " << std::scientific << rModelPart.GetProcessInfo()[DELTA_TIME] << std::endl;
            #else
                std::cout << "\n\n" << "CONVERGENCE CHECK" << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tTIME: " << std::scientific << rModelPart.GetProcessInfo()[TIME] << "\tDELTA TIME: " << std::scientific << rModelPart.GetProcessInfo()[DELTA_TIME] << std::endl;
            #endif
            }
            else
            {
                std::cout << "\n\n" << "CONVERGENCE CHECK" << "\tSTEP: " << rModelPart.GetProcessInfo()[TIME_STEPS] << "\tTIME: " << std::scientific << rModelPart.GetProcessInfo()[TIME] << "\tDELTA TIME: " << std::scientific << rModelPart.GetProcessInfo()[DELTA_TIME] << std::endl;
            }
                
            if (mpTable != nullptr)
            {
                mpTable->PrintHeader();
            }
        }
        
        BaseType::InitializeSolutionStep(rModelPart,rDofSet,A,Dx,b);
    }

    /**
     * This function finalizes the solution step
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
        
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& A,
        const TSystemVectorType& Dx,
        const TSystemVectorType& b
        ) override
    {
        BaseType::FinalizeSolutionStep(rModelPart,rDofSet,A,Dx,b);
    }

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access 
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

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
    
    TablePrinterPointerType mpTable; // Pointer to the fancy table 
    bool mPrintingOutput;            // If the colors and bold are printed
    bool mComputeConditionNumber;    // If the condition number is computed
    bool mTableIsInitialized;        // If the table is already initialized
    
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

}; /* Class ClassName */

///@}

///@name Type Definitions */
///@{

///@}

}  /* namespace Kratos.*/

#endif /* KRATOS_MORTAR_AND_CRITERIA_H  defined */

