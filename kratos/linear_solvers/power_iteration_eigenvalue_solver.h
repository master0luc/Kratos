//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//  Collaborator:    Vicente Mataix Ferrandiz
//
//

#if !defined(KRATOS_POWER_ITERATION_EIGENVALUE_SOLVER_H_INCLUDED )
#define  KRATOS_POWER_ITERATION_EIGENVALUE_SOLVER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <numeric>
#include <vector>
#include <random>

// External includes

// Project includes
#include "processes/process.h"
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"

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
template<class TSparseSpaceType, class TDenseSpaceType, class TLinearSolverType,
         class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class PowerIterationEigenvalueSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PowerIterationEigenvalueSolver
    KRATOS_CLASS_POINTER_DEFINITION(PowerIterationEigenvalueSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    typedef std::size_t SizeType;

    typedef std::size_t IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PowerIterationEigenvalueSolver() {}

    PowerIterationEigenvalueSolver(
        double MaxTolerance,
        unsigned int MaxIterationNumber,
        unsigned int RequiredEigenvalueNumber,
        typename TLinearSolverType::Pointer pLinearSolver
    ): BaseType(MaxTolerance, MaxIterationNumber),   
       mRequiredEigenvalueNumber(RequiredEigenvalueNumber),
       mpLinearSolver(pLinearSolver)
    {

    }

    PowerIterationEigenvalueSolver(
        Parameters ThisParameters,
        typename TLinearSolverType::Pointer pLinearSolver
        ): mpLinearSolver(pLinearSolver)
    {
        Parameters DefaultParameters = Parameters(R"(
        {
            "solver_type"             : "PowerIterationEigenvalueSolver",
            "max_iteration"           : 500,
            "tolerance"               : 1e-9,
            "required_eigen_number"   : 1,
            "shifting_convergence"    : 0.25,
            "verbosity"               : 1
        })" );

        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);

        mRequiredEigenvalueNumber = ThisParameters["required_eigen_number"].GetInt();
        mEchoLevel = ThisParameters["verbosity"].GetInt();
        BaseType::SetTolerance( ThisParameters["tolerance"].GetDouble() );
        BaseType::SetMaxIterationsNumber( ThisParameters["max_iteration"].GetInt() );
    }

    /// Copy constructor.
    PowerIterationEigenvalueSolver(const PowerIterationEigenvalueSolver& Other) : BaseType(Other)
    {

    }


    /// Destructor.
    ~PowerIterationEigenvalueSolver() override
    {

    }


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    PowerIterationEigenvalueSolver& operator=(const PowerIterationEigenvalueSolver& Other)
    {
        BaseType::operator=(Other);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    static void RandomInitialize(
        const SparseMatrixType& K,
        DenseVectorType& R
        )
    {
        // We create a random vector as seed of the method
        std::random_device this_random_device;
        std::mt19937 generator(this_random_device());
        
        const SizeType size = K.size1();
        const double normK = TSparseSpaceType::TwoNorm(K);
        std::normal_distribution<> normal_distribution(normK, 0.25 * normK);
        
        for (SizeType i = 0; i < size; i++)
        {
            R[i] = normal_distribution(generator);
        }
        
//         for(SizeType i = 0 ; i < R.size() ; i++)
//         {
//             R[i] = 1.00; //rand();
//         }

//         R /= norm_2(R);
    }


    /**
     * The power iteration algorithm
     * @param K: The stiffness matrix
     * @param M: The mass matrix
     * @param Eigenvalues: The vector containing the eigen values
     * @param Eigenvectors: The matrix containing the eigen vectors
     */
    void Solve(
        SparseMatrixType& K,
        SparseMatrixType& M,
        DenseVectorType& Eigenvalues,
        DenseMatrixType& Eigenvectors
        ) override
    {

        using boost::numeric::ublas::trans;

        const SizeType size = K.size1();
        const SizeType max_iteration = BaseType::GetMaxIterationsNumber();
        const double tolerance = BaseType::GetTolerance();

        VectorType x = ZeroVector(size);
        VectorType y = ZeroVector(size);

        RandomInitialize(K, y);

        if(Eigenvalues.size() < 1)
        {
            Eigenvalues.resize(1, 0.0);
        }

        // Starting with first step
        double beta = 0.0;
        double ro = 0.0;
        double old_ro = Eigenvalues[0];

        if (mEchoLevel > 1)
        {
            std::cout << "Iteration  beta \t\t ro \t\t convergence norm" << std::endl;
        }

        for(SizeType i = 0 ; i < max_iteration ; i++)
        {
            // K*x = y
            mpLinearSolver->Solve(K, x, y);
            
            ro = inner_prod(y, x);
            
            // y = M*x
            noalias(y) = prod(M, x);
            beta = inner_prod(x, y);
            
            if(beta <= 0.0)
            {
                KRATOS_ERROR << "M is not Positive-definite. beta = " << beta << std::endl;
            }

            ro /= beta;
            beta = std::sqrt(beta);
            y /= beta;

            if(ro == 0.0)
            {
                KRATOS_ERROR << "Perpendicular eigenvector to M" << std::endl;
            }

            const double convergence_norm = std::abs((ro - old_ro) / ro);

            if (mEchoLevel > 1)
            {
                std::cout << "Iteration: " << i << " \t beta: " << beta << "\tro: " << ro << " \tConvergence norm: " << convergence_norm << std::endl;
            }
            
            if(convergence_norm < tolerance)
            {
                break;
            }

            old_ro = ro;
        }

        if (mEchoLevel > 0)
        {
            KRATOS_WATCH(ro);
            KRATOS_WATCH(y);
        }

        Eigenvalues[0] = ro;

        if((Eigenvectors.size1() < 1) || (Eigenvectors.size2() < size))
        {
            Eigenvectors.resize(1,size);
        }

        for(SizeType i = 0 ; i < size ; i++)
        {
            Eigenvectors(0,i) = y[i];
        }
    }
    
    /**
     * This method returns directly the first eigen value obtained
     * @param K: The stiffness matrix
     * @param M: The mass matrix
     * @return The first eigenvalue
     */
    double GetEigenValue(
        SparseMatrixType& K,
        SparseMatrixType& M
        )
    {
        DenseVectorType eigen_values;
        DenseMatrixType eigen_vectors;
        
        Solve(K, M, eigen_values, eigen_vectors);
        
        return eigen_values[0];
    }
    
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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Power iteration eigenvalue solver with " << BaseType::GetPreconditioner()->Info();
        return  buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
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

    unsigned int mRequiredEigenvalueNumber;

    unsigned int mEchoLevel;

    typename TLinearSolverType::Pointer mpLinearSolver;

    std::vector<DenseVectorType> mQVector;
    std::vector<DenseVectorType> mPVector;
    std::vector<DenseVectorType> mRVector;

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

}; // Class PowerIterationEigenvalueSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::istream& operator >> (std::istream& IStream,
                                  PowerIterationEigenvalueSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const PowerIterationEigenvalueSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_POWER_ITERATION_EIGENVALUE_SOLVER_H_INCLUDED defined 































