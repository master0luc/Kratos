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
//
#if !defined(KRATOS_CONTACT_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER )
#define  KRATOS_CONTACT_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER

/* System includes  */

/* External includes  */

/* Project includes  */
#include "includes/key_hash.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

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

Current class provides an implementation for standard builder and solving operations.

the RHS is constituted by the unbalanced loads (residual)

Degrees of freedom are reordered putting the restrained degrees of freedom at
the end of the system ordered in reverse order with respect to the DofSet.

Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
this information.

Calculation of the reactions involves a cost very similiar to the calculation of the total residual

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
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ContactResidualBasedBlockBuilderAndSolver
    : public ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(ContactResidualBasedBlockBuilderAndSolver);

    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BSBaseType;
    
    typedef ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType                                            TSchemeType;

    typedef typename BaseType::TDataType                                                TDataType;

    typedef typename BaseType::DofsArrayType                                        DofsArrayType;

    typedef typename BaseType::TSystemMatrixType                                TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType                                TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType                        LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType                        LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType                  TSystemMatrixPointerType;
    
    typedef typename BaseType::TSystemVectorPointerType                  TSystemVectorPointerType;

    typedef typename BaseType::NodesArrayType                                      NodesArrayType;
    
    typedef typename BaseType::ElementsArrayType                                ElementsArrayType;
    
    typedef typename BaseType::ConditionsArrayType                            ConditionsArrayType;

    typedef typename BaseType::ElementsContainerType                        ElementsContainerType;

    ///@}
    ///@name Life Cycle
    
    ///@{

    /** Constructor. */
    
    ContactResidualBasedBlockBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
    {
//         std::cout << "Using the standard builder and solver " << std::endl;
    }

    /** Destructor. */
    
    ~ContactResidualBasedBlockBuilderAndSolver() override
    {
    }

    ///@}
    ///@name Operators
    
    ///@{

    //**************************************************************************
    //**************************************************************************
    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& b) override
    {
        KRATOS_TRY
        
        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        //getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        //getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

        // Contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        // Vector containing the localization in the system of the different terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        const double start_build = OpenMPUtils::GetCurrentTime();

        #pragma omp parallel firstprivate(nelements,nconditions, LHS_Contribution, RHS_Contribution, EquationId )
        {
            # pragma omp for  schedule(guided, 512) nowait
            for (int k = 0; k < nelements; k++)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + k;

                // Detect if the element is active or not. If the user did not make any choice the element
                // Is active by default
                bool element_is_active = true;
                if (it->IsDefined(ACTIVE))
                    element_is_active = it->Is(ACTIVE);

                if (element_is_active)
                {
                    // Calculate elemental contribution
                    pScheme->CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    // Assemble the elemental contribution
#ifdef USE_LOCKS_IN_ASSEMBLY
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, BaseType::mlock_array);
#else
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
#endif
                    // Clean local elemental memory
                    pScheme->CleanMemory(*(it.base()));
                }

            }

            //#pragma omp parallel for firstprivate(nconditions, LHS_Contribution, RHS_Contribution, EquationId ) schedule(dynamic, 1024)
            #pragma omp for  schedule(guided, 512)
            for (int k = 0; k < nconditions; k++)
            {
                ModelPart::ConditionsContainerType::iterator it = cond_begin + k;

                // Detect if the element is active or not. If the user did not make any choice the element is active by default
                bool condition_is_active = true;
                if (it->IsDefined(ACTIVE))
                    condition_is_active = it->Is(ACTIVE);

                if (condition_is_active)
                {
                    // Calculate elemental contribution
                    pScheme->Condition_CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    // Assemble the elemental contribution
#ifdef USE_LOCKS_IN_ASSEMBLY
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, BaseType::mlock_array);
#else
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
#endif
                    // Clean local elemental memory
                    pScheme->CleanMemory(*(it.base()));
                }
            }
        }

        const double stop_build = OpenMPUtils::GetCurrentTime();
        if (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)
            std::cout << "build time: " << stop_build - start_build << std::endl;

        if (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)
        {
            std::cout << "Finished parallel building" << std::endl;
        }

        KRATOS_CATCH("")

    }

    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart
        ) override
    {
        KRATOS_TRY;

        if( this->GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0)
        {
            std::cout << "Setting up the dofs" << std::endl;
        }

        //Gets the array of elements from the modeler
        ElementsArrayType& pElements = rModelPart.Elements();
        const int nelements = static_cast<int>(pElements.size());

        Element::DofsVectorType ElementalDofList;

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        unsigned int nthreads = OpenMPUtils::GetNumThreads();

#ifdef USE_GOOGLE_HASH
        typedef google::dense_hash_set < Node<3>::DofType::Pointer, DoFIteratorHash>  set_type;
#else
        typedef std::unordered_set < Node<3>::DofType::Pointer, DoFIteratorHash>  set_type;
#endif
        std::vector<set_type> dofs_aux_list(nthreads);

        if( this->GetEchoLevel() > 2)
        {
            std::cout << "Number of threads" << nthreads << "\n" << std::endl;
        }
        for (int i = 0; i < static_cast<int>(nthreads); i++)
        {
#ifdef USE_GOOGLE_HASH
            dofs_aux_list[i].set_empty_key(Node<3>::DofType::Pointer());
#else
            dofs_aux_list[i].reserve(nelements);
#endif
        }
        if( this->GetEchoLevel() > 2)
        {
            std::cout << "Initializing element loop" << std::endl;
        }
        #pragma omp parallel firstprivate(nelements, ElementalDofList)
        {
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < nelements; i++)
            {
                typename ElementsArrayType::iterator it = pElements.begin() + i;
                const unsigned int this_thread_id = OpenMPUtils::ThisThread();

                // gets list of Dof involved on every element
                pScheme->GetElementalDofList(*(it.base()), ElementalDofList, CurrentProcessInfo);

                dofs_aux_list[this_thread_id].insert(ElementalDofList.begin(), ElementalDofList.end());
            }
            if( this->GetEchoLevel() > 2)
            {
                std::cout << "Initializing condition loop\n" << std::endl;
            }
            ConditionsArrayType& pConditions = rModelPart.Conditions();
            const int nconditions = static_cast<int>(pConditions.size());
            #pragma omp for  schedule(guided, 512)
            for (int i = 0; i < nconditions; i++)
            {
                typename ConditionsArrayType::iterator it = pConditions.begin() + i;
                const unsigned int this_thread_id = OpenMPUtils::ThisThread();

                // Gets list of Dof involved on every element
                pScheme->GetConditionDofList(*(it.base()), ElementalDofList, CurrentProcessInfo);
                dofs_aux_list[this_thread_id].insert(ElementalDofList.begin(), ElementalDofList.end());
            }
        }

        if( this->GetEchoLevel() > 2)
        {
            std::cout << "Initializing tree reduction\n" << std::endl;
        }
        // Here we do a reduction in a tree so to have everything on thread 0
        unsigned int old_max = nthreads;
        unsigned int new_max = ceil(0.5*static_cast<double>(old_max));
        while (new_max>=1 && new_max != old_max)
        {
            if( this->GetEchoLevel() > 2)
            {
                //just for debugging
                std::cout << "old_max" << old_max << " new_max:" << new_max << std::endl;
                for (int i = 0; i < static_cast<int>(new_max); i++)
                {
                    if (i + new_max < old_max)
                    {
                        std::cout << i << " - " << i+new_max << std::endl;
                    }
                }
                std::cout << "********************" << std::endl;
            }

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(new_max); i++)
            {
                if (i + new_max < old_max)
                {
                    dofs_aux_list[i].insert(dofs_aux_list[i+new_max].begin(), dofs_aux_list[i+new_max].end());
                    dofs_aux_list[i+new_max].clear();
                }
            }

            old_max = new_max;
            new_max = ceil(0.5*static_cast<double>(old_max));
        }

        if( this->GetEchoLevel() > 2)
        {
            std::cout << "Initializing ordered array filling\n" << std::endl;
        }

        DofsArrayType Doftemp;
        BaseType::mDofSet = DofsArrayType();

        Doftemp.reserve(dofs_aux_list[0].size());
        for (auto it= dofs_aux_list[0].begin(); it!= dofs_aux_list[0].end(); it++)
        {
            Doftemp.push_back( it->get() );
        }
        Doftemp.Sort();

        BaseType::mDofSet = Doftemp;

        //Throws an exception if there are no Degrees Of Freedom involved in the analysis
        KRATOS_ERROR_IF(BaseType::mDofSet.size() == 0) << "No degrees of freedom!" << std::endl;

        if( this->GetEchoLevel() > 2)
        {
            KRATOS_WATCH(BaseType::mDofSet.size())
        }
        BaseType::mDofSetIsInitialized = true;
        if( this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)
        {
            std::cout << "Finished setting up the dofs" << std::endl;
        }

        KRATOS_ERROR_IF(this->GetEchoLevel() > 2) << "Initializing lock array" << std::endl;

#ifdef _OPENMP
        if (BaseType::mlock_array.size() != 0)
        {
            for (int i = 0; i < static_cast<int>(BaseType::mlock_array.size()); i++)
            {
                omp_destroy_lock(&BaseType::mlock_array[i]);
            }
        }
        BaseType::mlock_array.resize(BaseType::mDofSet.size());

        for (int i = 0; i < static_cast<int>(BaseType::mlock_array.size()); i++)
        {
            omp_init_lock(&BaseType::mlock_array[i]);
        }
#endif
        if( this->GetEchoLevel() > 2)
        {
            std::cout << "End of setupdofset\n" << std::endl;
        }

        KRATOS_CATCH("");
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
    ///@name Protected Operators*/
    ///@{


    virtual void ConstructMatrixStructure(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixType& A,
        ElementsContainerType& rElements,
        ConditionsArrayType& rConditions,
        ProcessInfo& CurrentProcessInfo)
    {
        //filling with zero the matrix (creating the structure)
        Timer::Start("MatrixStructure");

        const std::size_t equation_size = BaseType::mEquationSystemSize;

#ifdef USE_GOOGLE_HASH
        std::vector<google::dense_hash_set<std::size_t> > indices(equation_size);
        const std::size_t empty_key = 2*equation_size + 10;
#else
        std::vector<std::unordered_set<std::size_t> > indices(equation_size);
#endif

        #pragma omp parallel for firstprivate(equation_size)
        for (int iii = 0; iii < static_cast<int>(equation_size); iii++)
        {
#ifdef USE_GOOGLE_HASH
        indices[iii].set_empty_key(empty_key);
#else
        indices[iii].reserve(40);
#endif
        }

        Element::EquationIdVectorType ids(3, 0);

        const int nelements = static_cast<int>(rElements.size());
        #pragma omp parallel for firstprivate(nelements, ids)
        for(int iii=0; iii<nelements; iii++)
        {
            typename ElementsContainerType::iterator i_element = rElements.begin() + iii;
            pScheme->EquationId( *(i_element.base()) , ids, CurrentProcessInfo);
            for (std::size_t i = 0; i < ids.size(); i++)
            {
#ifdef _OPENMP
                omp_set_lock(&BaseType::mlock_array[ids[i]]);
#endif
                auto& row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());

#ifdef _OPENMP
                omp_unset_lock(&BaseType::mlock_array[ids[i]]);
#endif
            }
        }

        const int nconditions = static_cast<int>(rConditions.size());
        #pragma omp parallel for firstprivate(nconditions, ids)
        for (int iii = 0; iii<nconditions; iii++)
        {
            typename ConditionsArrayType::iterator i_condition = rConditions.begin() + iii;
            pScheme->Condition_EquationId( *(i_condition.base()), ids, CurrentProcessInfo);
            for (std::size_t i = 0; i < ids.size(); i++)
            {
#ifdef _OPENMP
                omp_set_lock(&BaseType::mlock_array[ids[i]]);
#endif
                auto& row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());
#ifdef _OPENMP
                omp_unset_lock(&BaseType::mlock_array[ids[i]]);
#endif
            }
        }

        //count the row sizes
        unsigned int nnz = 0;
        for (unsigned int i = 0; i < indices.size(); i++)
            nnz += indices[i].size();

        A = compressed_matrix<double>(indices.size(), indices.size(), nnz);

        double* Avalues = A.value_data().begin();
        std::size_t* Arow_indices = A.index1_data().begin();
        std::size_t* Acol_indices = A.index2_data().begin();

        //filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
        Arow_indices[0] = 0;
        for (int i = 0; i < static_cast<int>(A.size1()); i++)
            Arow_indices[i+1] = Arow_indices[i] + indices[i].size();

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(A.size1()); i++)
        {
            const unsigned int row_begin = Arow_indices[i];
            const unsigned int row_end = Arow_indices[i+1];
            unsigned int k = row_begin;
            for (auto it = indices[i].begin(); it != indices[i].end(); it++)
            {
                Acol_indices[k] = *it;
                Avalues[k] = 0.0;
                k++;
            }

            indices[i].clear(); //deallocating the memory

            std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);
        }

        A.set_filled(indices.size()+1, nnz);

        Timer::Stop("MatrixStructure");
    }

    ///@}
    ///@name Protected Operations*/
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

    void BuildRHSNoDirichlet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        //Getting the Elements
        ElementsArrayType& pElements = rModelPart.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = rModelPart.Conditions();

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // Assemble all elements
        const int nelements = static_cast<int>(pElements.size());
        #pragma omp parallel firstprivate(nelements, RHS_Contribution, EquationId)
        {
            #pragma omp for schedule(guided, 512) nowait
            for(int i=0; i<nelements; i++)
            {
                typename ElementsArrayType::iterator it = pElements.begin() + i;
                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool element_is_active = true;
                if( it->IsDefined(ACTIVE) )
                    element_is_active = it->Is(ACTIVE);

                if(element_is_active)
                {
                    // Calculate elemental Right Hand Side Contribution
                    pScheme->Calculate_RHS_Contribution(*(it.base()), RHS_Contribution, EquationId, CurrentProcessInfo);

                    // Assemble the elemental contribution
                    BaseType::AssembleRHS(b, RHS_Contribution, EquationId);
                }
            }

            LHS_Contribution.resize(0, 0, false);
            RHS_Contribution.resize(0, false);

            // Assemble all conditions
            const int nconditions = static_cast<int>(ConditionsArray.size());

            #pragma omp for schedule(guided, 512)
            for (int i = 0; i<nconditions; i++)
            {
                auto it = ConditionsArray.begin() + i;
                // Detect if the element is active or not. If the user did not make any choice the element is active by default
                bool condition_is_active = true;
                if( it->IsDefined(ACTIVE) )
                    condition_is_active = it->Is(ACTIVE);

                if(condition_is_active)
                {
                    // Calculate elemental contribution
                    pScheme->Condition_Calculate_RHS_Contribution(*(it.base()), RHS_Contribution, EquationId, CurrentProcessInfo);

                    // Assemble the elemental contribution
                    BaseType::AssembleRHS(b, RHS_Contribution, EquationId);
                }
            }
        }

        KRATOS_CATCH("")
    }

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

}; /* Class ContactResidualBasedBlockBuilderAndSolver */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_CONTACT_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER  defined */
