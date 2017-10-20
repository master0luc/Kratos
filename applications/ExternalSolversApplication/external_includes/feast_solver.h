//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:		 BSD License
//   Kratos default license: kratos/license.txt
//
//   Project Name:        $ExternalSolversApplication   $
//   Last modified by:    $Author: michael.andre@tum.de $
//   Date:                $Date:             April 2017 $
//   Revision:            $Revision:                0.0 $
//
//

// System includes
#include <iostream>
#include <complex>
#include <vector>
#include <unordered_set>
#include <algorithm>

// External includes
#include <boost/smart_ptr.hpp>
extern "C" {
#include "feast.h"
}

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/skyline_lu_custom_scalar_solver.h"
#include "includes/ublas_interface.h"
#include "includes/ublas_complex_interface.h"
#include "spaces/ublas_space.h"

#if !defined(KRATOS_FEAST_SOLVER)
#define  KRATOS_FEAST_SOLVER

namespace Kratos {

///@name Kratos Classes
///@{

/// Adapter to FEAST eigenvalue problem solver.
template<class TSparseSpaceType, class TDenseSpaceType,
        class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class FEASTSolver: public LinearSolver<TSparseSpaceType, TDenseSpaceType,
        TReordererType> {

public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( FEASTSolver );

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType SparseVectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    typedef UblasSpace<std::complex<double>, ComplexCompressedMatrix, ComplexVector> ComplexSparseSpaceType;

    typedef UblasSpace<std::complex<double>, ComplexMatrix, ComplexVector> ComplexDenseSpaceType;

    typedef LinearSolver<ComplexSparseSpaceType, ComplexDenseSpaceType> ComplexLinearSolverType;

    typedef ComplexLinearSolverType::SparseMatrixType ComplexSparseMatrixType;

    typedef ComplexLinearSolverType::VectorType ComplexVectorType;
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for built-in linear solver.
    /**
     * Parameters let the user control the settings of the FEAST library.
     */
    FEASTSolver(Parameters::Pointer pParam) : mpParam(pParam)
    {
        Parameters default_params(R"(
        {
            "solver_type": "FEAST",
            "print_feast_output": false,
            "perform_stochastic_estimate": true,
            "solve_eigenvalue_problem": true,
            "lambda_min": 0.0,
            "lambda_max": 1.0,
            "echo_level": 1,
            "number_of_eigenvalues": 0,
            "search_dimension": 10,
            "linear_solver_settings": {
                "solver_type": "complex_skyline_lu_solver"
            }
        })");

        mpParam->RecursivelyValidateAndAssignDefaults(default_params);

        if (mpParam->GetValue("linear_solver_settings")["solver_type"].GetString() != "complex_skyline_lu_solver")
            KRATOS_ERROR << "built-in solver type must be used with this constructor" << std::endl;
            
        mpLinearSolver = boost::make_shared<SkylineLUCustomScalarSolver<ComplexSparseSpaceType, ComplexDenseSpaceType>>();
    }

    /// Constructor for externally provided linear solver.
    /**
     * Parameters let the user control the settings of the FEAST library.
     * Warning: For iterative solvers, very small tolerances (~1e-15)
     *          may be needed for FEAST to work properly. Common iterative 
     *          solvers normally don't perform efficiently with FEAST 
     *          (M. Galgon et al., Parallel Computing (49) 2015 153-163).
     */
    FEASTSolver(
        Parameters::Pointer pParam, 
        ComplexLinearSolverType::Pointer pLinearSolver
        ): mpParam(pParam)
    {
        Parameters default_params(R"(
        {
            "solver_type": "FEAST",
            "print_feast_output": false,
            "perform_stochastic_estimate": true,
            "solve_eigenvalue_problem": true,
            "lambda_min": 0.0,
            "lambda_max": 1.0,
            "echo_level": 1,
            "number_of_eigenvalues": 0,
            "search_dimension": 10,
            "linear_solver_settings": {}
        })");

        // Don't validate linear_solver_settings here
        mpParam->ValidateAndAssignDefaults(default_params);
        
        if (pLinearSolver != nullptr)
        {
            mpLinearSolver = pLinearSolver;
        }
        else
        {
            mpLinearSolver = boost::make_shared<SkylineLUSolver<ComplexSparseSpaceType, ComplexDenseSpaceType>>();
        }
        
    }

    /// Deleted copy constructor.
    FEASTSolver(const FEASTSolver& Other) = delete;

    /// Destructor.
    ~FEASTSolver() override {}

    ///@}
    ///@name Operators
    ///@{

    /// Deleted assignment operator.
    FEASTSolver& operator=(const FEASTSolver& Other) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// Solve the generalized eigenvalue problem.
    /**
     * K is a symmetric matrix. M is a symmetric positive-definite matrix.
     */
    void Solve(
            SparseMatrixType& K,
            SparseMatrixType& M,
            DenseVectorType& Eigenvalues,
            DenseMatrixType& Eigenvectors) override
    {
        const auto system_size = K.size1();

        Parameters& feast_settings = *mpParam;
        const double eigen_value_range_min = feast_settings["lambda_min"].GetDouble();
        const double eigen_value_range_max = feast_settings["lambda_max"].GetDouble();

        int search_dimension = feast_settings["search_dimension"].GetInt();
        int num_eigenvalues = feast_settings["number_of_eigenvalues"].GetInt();
        const int echo_level = feast_settings["echo_level"].GetInt();

        Eigenvalues.resize(search_dimension,false);
        Eigenvectors.resize(search_dimension,system_size,false);

        if (feast_settings["perform_stochastic_estimate"].GetBool())
        {
            // This estimates the number of eigenvalues in the interval [lambda_min, lambda_max]
            Calculate(M,K,eigen_value_range_min,eigen_value_range_max,search_dimension,
                    num_eigenvalues,Eigenvalues,Eigenvectors,true);

            if (echo_level > 0)
            {
                std::cout << "Estimated number of eigenvalues = " << num_eigenvalues << std::endl;
            }

            // Recommended estimate of search dimension from FEAST documentation
            search_dimension = num_eigenvalues + num_eigenvalues/2 + 1;
            feast_settings["search_dimension"].SetInt(search_dimension);
        }
        if (feast_settings["solve_eigenvalue_problem"].GetBool())
        {
            // This attempts to solve the generalized eigenvalue problem
            Calculate(M,K,eigen_value_range_min,eigen_value_range_max,search_dimension,
                    num_eigenvalues,Eigenvalues,Eigenvectors,false);

            Eigenvalues.resize(num_eigenvalues,true);
            Eigenvectors.resize(num_eigenvalues,system_size,true);
        }
        feast_settings["number_of_eigenvalues"].SetInt(num_eigenvalues);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FEAST solver.";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    Parameters::Pointer mpParam;

    ComplexLinearSolverType::Pointer mpLinearSolver;

    ///@}
    ///@name Private Operations
    ///@{

    /// Wrapper for FEAST library.
    void Calculate(
            SparseMatrixType& rMassMatrix,
            SparseMatrixType& rStiffnessMatrix,
            double EigenvalueRangeMin,
            double EigenvalueRangeMax,
            int SearchDimension,
            int& rNumEigenvalues,
            DenseVectorType& rEigenvalues,
            DenseMatrixType& rEigenvectors,
            bool PerformStochasticEstimate)
    {
        KRATOS_TRY

        int feast_params[64] = {};
        int num_iter, info, system_size;
        double eps_out;
        DenseVectorType residual(SearchDimension);
        std::vector<std::complex<double> > integration_nodes, integration_weights;
        system_size = static_cast<int>(rMassMatrix.size1());
        matrix<double,column_major> work(system_size,SearchDimension);
        matrix<std::complex<double>,column_major> zwork(system_size,SearchDimension);
        matrix<double,column_major> Aq(SearchDimension,SearchDimension);
        matrix<double,column_major> Bq(SearchDimension,SearchDimension);
        std::complex<double> Ze;
        ComplexSparseMatrixType Az;
        ComplexVectorType b(system_size);
        ComplexVectorType x(system_size);

        this->InitializeFEASTSystemMatrix(rMassMatrix, rStiffnessMatrix, Az);

        Parameters& feast_settings = *mpParam;

        // initialize FEAST eigenvalue solver (see FEAST documentation for details)
        feastinit(feast_params);
        if (feast_settings["print_feast_output"].GetBool())
            feast_params[0] = 1;
        feast_params[2] = 8; // stopping convergence criteria 10^-feast_params[2]
        feast_params[28] = 1;// not sure if this is needed
        if (PerformStochasticEstimate)
        {
            feast_params[1] = 4; // number of quadrature points (default: 8)
            feast_params[13] = 2;
        }

        integration_nodes.resize(feast_params[1]);
        integration_weights.resize(feast_params[1]);

        // get quadrature nodes and weights
        zfeast_contour(&EigenvalueRangeMin,
                &EigenvalueRangeMax,
                &feast_params[1],
                &feast_params[15],
                &feast_params[17],
                (double *)integration_nodes.data(),
                (double *)integration_weights.data());

        int ijob = -1;
        // solve the eigenvalue problem
        while (ijob != 0)
        {
            // FEAST's reverse communication interface
            dfeast_srcix(&ijob,&system_size,(double *)&Ze,(double *)work.data().begin(),
                    (double *)zwork.data().begin(),(double *)Aq.data().begin(),
                    (double *)Bq.data().begin(),feast_params,&eps_out,&num_iter,
                    &EigenvalueRangeMin,&EigenvalueRangeMax,&SearchDimension,
                    (double *)rEigenvalues.data().begin(),
                    (double *)rEigenvectors.data().begin(),
                    &rNumEigenvalues,(double *)residual.data().begin(),&info,
                    (double *)integration_nodes.data(),
                    (double *)integration_weights.data());

            switch (ijob)
            {
                case 10:
                {
                    // set up quadrature matrix (ZeM-K) and solver
                    this->CalculateFEASTSystemMatrix(Ze, rMassMatrix, rStiffnessMatrix, Az);
                    mpLinearSolver->Clear();
                    mpLinearSolver->Initialize(Az,x,b);
                    mpLinearSolver->InitializeSolutionStep(Az, x, b);
                } break;
                case 11:
                {
                    // solve the linear system for one quadrature point
                    for (int j=0; j < feast_params[22]; j++)
                    {
                        for (int i=0; i < system_size; i++)
                            b[i] = zwork(i,j);
                        mpLinearSolver->Solve(Az,x,b);
                        for (int i=0; i < system_size; i++)
                            zwork(i,j) = x[i];
                    }
                } break;
                case 30:
                {
                    // multiply Kx
                    for (int i=0; i < feast_params[24]; i++)
                    {
                        int k = feast_params[23]-1+i;
                        noalias(column(work,k)) = prod(rStiffnessMatrix,row(rEigenvectors,k));
                    }
                } break;
                case 40:
                {
                    // multiply Mx
                    for (int i=0; i < feast_params[24]; i++)
                    {
                        int k = feast_params[23]-1+i;
                        noalias(column(work,k)) = prod(rMassMatrix,row(rEigenvectors,k));
                    }
                }
            } // switch
        } // while

        KRATOS_CATCH("")
    }

    /**
     * Initialize CSR matrix structure for FEAST system matrix: C = z * B - A.
     */
    void InitializeFEASTSystemMatrix(
        const SparseMatrixType& B,
        const SparseMatrixType& A,
        ComplexSparseMatrixType& C
        )
    {
        C.resize(B.size1(), B.size2(), false);

        std::vector<std::unordered_set<std::size_t> > indices(C.size1());

        // indices for row begin / end
        C.index1_data()[0] = 0;
        for (std::size_t i = 0; i < C.size1(); ++i)
        {
            std::size_t row_begin, row_end;
            indices[i].reserve(40); // initialize C's indices

            row_begin = B.index1_data()[i];
            row_end = B.index1_data()[i + 1];
            indices[i].insert(B.index2_data().begin() + row_begin,
                    B.index2_data().begin() + row_end); // insert B's column indices for row i

            row_begin = A.index1_data()[i];
            row_end = A.index1_data()[i + 1];
            indices[i].insert(A.index2_data().begin() + row_begin,
                    A.index2_data().begin() + row_end); // insert A's column indices for row i

            // C.index1_data()[i+1] = number of non-zeros in rows <= i
            C.index1_data()[i + 1] = C.index1_data()[i] + indices[i].size();
        }

        // C.index1_data()[C.size1()] = number of non-zeros
        C.reserve(C.index1_data()[C.size1()]);

        // column indices
        std::size_t k = 0;
        for (std::size_t i = 0; i < C.size1(); ++i)
        {
            for (std::size_t j : indices[i])
                C.index2_data()[k++] = j; // fill C's column indices

            indices[i].clear();

            std::sort(C.index2_data().begin() + C.index1_data()[i],
                    C.index2_data().begin() + C.index1_data()[i + 1]);
        }

        C.set_filled(C.size1() + 1, C.index1_data()[C.size1()]);
    }

    /**
     * Calculate FEAST system matrix: C = z * B - A. Similar to FEAST's zdaddcsr subroutine.
     */
    void CalculateFEASTSystemMatrix(
        std::complex<double> z,
        SparseMatrixType& B,
        SparseMatrixType& A,
        ComplexSparseMatrixType& C
        )
    {
        std::size_t jb, ja;
        const std::size_t dimension = B.size1();

        std::size_t ptr = 0;
        for (std::size_t i = 0; i < dimension; ++i)
        {
            std::size_t b_ptr = B.index1_data()[i];
            std::size_t a_ptr = A.index1_data()[i];
            while (b_ptr < B.index1_data()[i + 1] || a_ptr < A.index1_data()[i + 1])
            {
                jb = (b_ptr < B.index1_data()[i + 1]) ?
                        B.index2_data()[b_ptr] : dimension;
                ja = (a_ptr < A.index1_data()[i + 1]) ?
                        A.index2_data()[a_ptr] : dimension;

                if (jb < ja)
                {
                    C.value_data()[ptr] = z * B(i, jb).ref();
                    b_ptr++;
                }
                else if (jb > ja)
                {
                    C.value_data()[ptr] = -A(i, ja).ref();
                    a_ptr++;
                }
                else
                { // jb == ja
                    C.value_data()[ptr] = z * B(i, jb).ref() - A(i, ja).ref();
                    b_ptr++;
                    a_ptr++;
                }
                ptr++;
            }
        }
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class FEASTSolver

///@}

///@name Input and output
///@{

/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >>(std::istream& rIStream,
        FEASTSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    return rIStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator <<(std::ostream& rOStream,
        const FEASTSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}// namespace Kratos.

#endif // KRATOS_FEAST_SOLVER  defined
