//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_POWER_METHOD_UTILS )
#define  KRATOS_POWER_METHOD_UTILS


/* System includes */
#include <random>

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "includes/ublas_interface.h"
#include "utilities/math_utils.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"

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

///Various mathematical utilities to compute SVD and the condition number of a matrix
/**
 * Defines several utility functions
 */
template<class TDataType>
class PowerMethodUtils
{
public:

    ///@name Type Definitions
    ///@{
    
    typedef Matrix MatrixType;

    typedef Vector VectorType;

    typedef std::size_t SizeType;
    
    typedef unsigned int IndexType;

    typedef UblasSpace<TDataType, CompressedMatrix, Vector> SparseSpaceType;
    
    typedef UblasSpace<TDataType, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType,  LocalSpaceType> LinearSolverType;

    typedef Reorderer<SparseSpaceType,  LocalSpaceType > ReordererType;

    typedef DirectSolver<SparseSpaceType,  LocalSpaceType, ReordererType > DirectSolverType;

    typedef SkylineLUFactorizationSolver<SparseSpaceType,  LocalSpaceType, ReordererType > SkylineLUFactorizationSolverType;

    ///@}
    ///@name Life Cycle
    ///@{

    /* Constructor */


    /** Destructor */

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    /**
     * This function computes using the power method the maximal eigenvalue
     * @param InputMatrix: The matrix to compute the eigenvalue
     * @param MaxIter: The max number of iterations
     * @return eigen: The maximal eigen value
     */
    
    static inline TDataType PowerMethod(
        const MatrixType& InputMatrix, 
        const SizeType MaxIter = 1000
        )
    {
        // We create a random vector as seed of the method
        std::random_device this_random_device;
        std::mt19937 generator(this_random_device());
        
        const SizeType size = InputMatrix.size1();
        const TDataType normA = SparseSpaceType::TwoNorm(InputMatrix);
        std::normal_distribution<> normal_distribution(normA, 0.25 * normA);
        
        VectorType eigen_vector(size);
        for (SizeType i = 0; i < size; i++)
        {
            eigen_vector[i] = normal_distribution(generator);
        }
        
        // We initialize values
        TDataType tolerance = 1.0e-10;
        TDataType eigen = 0.0;
        TDataType previous_eigen = 100;
        VectorType aux_vector = ZeroVector(size);
        
        SizeType iter = 0;
        while(std::abs(eigen-previous_eigen) > tolerance || SparseSpaceType::TwoNorm(eigen_vector - aux_vector) > tolerance || iter < MaxIter)
        {
            // We save the previous information
            previous_eigen = eigen;
            aux_vector = eigen_vector;
            
            // We update the vector 
            noalias(eigen_vector) = prod(InputMatrix, aux_vector);
            eigen = std::max_element(eigen_vector.begin(), eigen_vector.end());
            eigen_vector /= eigen;
            
            iter += 1;
        }
        
        return eigen;
    }
    
    /**
     * This function computes using the inverse power method the minimal eigenvalue
     * @param InputMatrix: The matrix to compute the eigenvalue
     * @param MaxIter: The max number of iterations
     * @return eigen: The maximal eigen value
     */
    
    static inline TDataType InvPowerMethod(
        const MatrixType& InputMatrix, 
        const SizeType MaxIter = 1000,
        LinearSolverType::Pointer pLinearSolver = nullptr
        )
    {
        if (pLinearSolver == nullptr)
        {
            pLinearSolver = boost::make_shared<SkylineLUFactorizationSolver()>;
        }
        
        // We create a random vector as seed of the method
        std::random_device this_random_device;
        std::mt19937 generator(this_random_device());
        
        const SizeType size = InputMatrix.size1();
        const TDataType normA = SparseSpaceType::TwoNorm(InputMatrix);
        std::normal_distribution<> normal_distribution(normA, 0.25 * normA);
        
        VectorType eigen_vector(size);
        for (SizeType i = 0; i < size; i++)
        {
            eigen_vector[i] = normal_distribution(generator);
        }
        
        // We initialize values
        TDataType tolerance = 1.0e-10;
        TDataType eigen = 0.0;
        TDataType previous_eigen = 100;
        VectorType aux_vector = ZeroVector(size);
        
        SizeType iter = 0;
        while(std::abs(eigen-previous_eigen) > tolerance || SparseSpaceType::TwoNorm(eigen_vector - aux_vector) > tolerance || iter < MaxIter)
        {
            // We save the previous information
            previous_eigen = eigen;
            aux_vector = eigen_vector;
            
            // We update the vector 
            pLinearSolver->Solve(InputMatrix, eigen_vector, aux_vector);
            eigen = 1.0/SparseSpaceType::TwoNorm(eigen_vector);
            eigen_vector *= eigen;
            
            iter += 1;
        }
        
        return eigen;
    }
    
    /**
     * This function computes using the inverse power method the minimal eigenvalue
     * @param InputMatrix: The matrix to compute the eigenvalue
     * @param MaxIter: The max number of iterations
     * @return eigen: The maximal eigen value
     */
    
    static inline TDataType ConditionNumber(
        const MatrixType& InputMatrix, 
        const SizeType MaxIter = 1000,
        LinearSolverType::Pointer pLinearSolver = nullptr
        )
    {
        TDataType condition_number;
        
        const TDataType max_lambda = PowerMethod(InputMatrix, MaxIter);
        const TDataType min_lambda = InvPowerMethod(InputMatrix, MaxIter, pLinearSolver);
        
        condition_number = std::abs(max_lambda)/std::abs(min_lambda); 
        
        return condition_number;
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

    ///@}
    ///@name Friends
    ///@{

private:
    
    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
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
    ///@name Private LifeCycle
    ///@{

    ///@}
    ///@name Unaccessible methods
    ///@{

    PowerMethodUtils(void);

    PowerMethodUtils(PowerMethodUtils& rSource);

}; /* Class PowerMethodUtils */

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/

#endif /* KRATOS_POWER_METHOD_UTILS  defined */

