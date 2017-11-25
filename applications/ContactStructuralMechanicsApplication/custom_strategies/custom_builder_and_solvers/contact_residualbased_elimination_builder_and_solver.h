//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
// //					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#if !defined(KRATOS_CONTACT_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER )
#define  KRATOS_CONTACT_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER

/* System includes */

/* External includes */

/* Project includes */
#include "includes/key_hash.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

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
class ContactResidualBasedEliminationBuilderAndSolver
    : public ResidualBasedEliminationBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(ContactResidualBasedEliminationBuilderAndSolver);

    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BSBaseType;
    
    typedef ResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType                                                  TSchemeType;

    typedef typename BaseType::TDataType                                                      TDataType;

    typedef typename BaseType::DofsArrayType                                              DofsArrayType;

    typedef typename BaseType::TSystemMatrixType                                      TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType                                      TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType                              LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType                              LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType                        TSystemMatrixPointerType;
    
    typedef typename BaseType::TSystemVectorPointerType                        TSystemVectorPointerType;

    typedef typename BaseType::NodesArrayType                                            NodesArrayType;
    
    typedef typename BaseType::ElementsArrayType                                      ElementsArrayType;
    
    typedef typename BaseType::ConditionsArrayType                                  ConditionsArrayType;

    typedef typename BaseType::ElementsContainerType                              ElementsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     */
    ContactResidualBasedEliminationBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : ResidualBasedEliminationBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
    {
//             std::cout << "using the standard builder and solver " << std::endl;
    }

    /** Destructor.
     */
    ~ContactResidualBasedEliminationBuilderAndSolver() override
    {
    }
    
    ///@}
    ///@name Operators
    ///@{

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

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        double start_build = OpenMPUtils::GetCurrentTime();

        #pragma omp parallel firstprivate(nelements, nconditions,  LHS_Contribution, RHS_Contribution, EquationId )
        {
            #pragma omp  for schedule(guided, 512) nowait
            for (int k = 0; k < nelements; k++)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + k;

                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool element_is_active = true;
                if (it->IsDefined(ACTIVE))
                    element_is_active = it->Is(ACTIVE);

                if (element_is_active)
                {
                    //calculate elemental contribution
                    pScheme->CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
#ifdef _OPENMP
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, mlock_array);
#else
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
#endif
                    // clean local elemental memory
                    pScheme->CleanMemory(*(it.base()));
                }
            }

            #pragma omp  for schedule(guided, 512)
            for (int k = 0; k < nconditions; k++)
            {
                    ModelPart::ConditionsContainerType::iterator it = cond_begin + k;

                    //detect if the element is active or not. If the user did not make any choice the element
                    //is active by default
                    bool condition_is_active = true;
                    if (it->IsDefined(ACTIVE))
                            condition_is_active = it->Is(ACTIVE);

                    if (condition_is_active)
                    {
                        // Calculate elemental contribution
                        pScheme->Condition_CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

#ifdef _OPENMP
                        BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, mlock_array);
#else
                        BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
#endif
                        // Clean local elemental memory
                        pScheme->CleanMemory(*(it.base()));
                    }
            }
        }

        double stop_build = OpenMPUtils::GetCurrentTime();
        if (this->GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0)
                std::cout << "build time: " << stop_build - start_build << std::endl;

        if (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)
        {
                KRATOS_WATCH("finished building");
        }

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    void BuildLHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A) override
    {
        KRATOS_TRY

        //getting the elements from the model
        ElementsArrayType& rElements = rModelPart.Elements();

        //getting the array of the conditions
        ConditionsArrayType& rConditions = rModelPart.Conditions();

        //resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        // assemble all elements
        for (typename ElementsArrayType::ptr_iterator it = rElements.ptr_begin(); it != rElements.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            BaseType::AssembleLHS(A, LHS_Contribution, EquationId);

            // clean local elemental memory
            pScheme->CleanMemory(*it);
        }

        LHS_Contribution.resize(0, 0, false);

        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = rConditions.ptr_begin(); it != rConditions.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            BaseType::AssembleLHS(A, LHS_Contribution, EquationId);
        }

        KRATOS_CATCH("")

    }

    //**************************************************************************
    //**************************************************************************

    void BuildLHS_CompleteOnFreeRows(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A) override
    {
        KRATOS_TRY

        //getting the elements from the model
        ElementsArrayType& rElements = rModelPart.Elements();

        //getting the array of the conditions
        ConditionsArrayType& rConditions = rModelPart.Conditions();

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        //resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // Assemble all elements
        for (typename ElementsArrayType::ptr_iterator it = rElements.ptr_begin(); it != rElements.ptr_end(); ++it)
        {
            // Calculate elemental contribution
            pScheme->Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, CurrentProcessInfo);

            // Assemble the elemental contribution
            AssembleLHS_CompleteOnFreeRows(A, LHS_Contribution, EquationId);

            // clean local elemental memory
            pScheme->CleanMemory(*it);
        }

        LHS_Contribution.resize(0, 0, false);
        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = rConditions.ptr_begin(); it != rConditions.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS_CompleteOnFreeRows(A, LHS_Contribution, EquationId);
        }


        KRATOS_CATCH("")

    }

    //**************************************************************************
    //**************************************************************************

    void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        //resetting to zero the vector of reactions

        if(BaseType::mCalculateReactionsFlag)
        {
            TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));
        }

        //Getting the Elements
        ElementsArrayType& pElements = rModelPart.Elements();

        //getting the array of the conditions
        ConditionsArrayType& pConditions = rModelPart.Conditions();

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        //contributions to the system
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements

        #pragma omp parallel firstprivate( RHS_Contribution, EquationId)
        {
            const int nelements = static_cast<int>(pElements.size());
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i<nelements; i++)
            {
                typename ElementsArrayType::iterator it = pElements.begin() + i;
                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool element_is_active = true;
                if (it->IsDefined(ACTIVE))
                        element_is_active = it->Is(ACTIVE);

                if (element_is_active)
                {
                        //calculate elemental Right Hand Side Contribution
                        pScheme->Calculate_RHS_Contribution(*(it.base()), RHS_Contribution, EquationId, CurrentProcessInfo);

                        //assemble the elemental contribution
                        AssembleRHS(b, RHS_Contribution, EquationId);
                }
            }

            // assemble all conditions
            const int nconditions = static_cast<int>(pConditions.size());
            #pragma omp  for schedule(guided, 512)
            for (int i = 0; i<nconditions; i++)
            {
                auto it = pConditions.begin() + i;
                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool condition_is_active = true;
                if (it->IsDefined(ACTIVE))
                        condition_is_active = it->Is(ACTIVE);

                if (condition_is_active)
                {

                        //calculate elemental contribution
                        pScheme->Condition_Calculate_RHS_Contribution(*(it.base()), RHS_Contribution, EquationId, CurrentProcessInfo);

                        //assemble the elemental contribution
                        AssembleRHS(b, RHS_Contribution, EquationId);
                }
            }
        }


        KRATOS_CATCH("")
    }

    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart
        ) override
    {
        KRATOS_TRY;

        if (this->GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0)
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

        for (int i = 0; i < static_cast<int>(nthreads); i++)
        {
#ifdef USE_GOOGLE_HASH
            dofs_aux_list[i].set_empty_key(Node<3>::DofType::Pointer());
#else
            dofs_aux_list[i].reserve(nelements);
#endif
        }

#pragma omp parallel for firstprivate(nelements, ElementalDofList)
        for (int i = 0; i < static_cast<int>(nelements); i++)
        {
            typename ElementsArrayType::iterator it = pElements.begin() + i;
            const unsigned int this_thread_id = OpenMPUtils::ThisThread();

            // gets list of Dof involved on every element
            pScheme->GetElementalDofList(*(it.base()), ElementalDofList, CurrentProcessInfo);

            dofs_aux_list[this_thread_id].insert(ElementalDofList.begin(), ElementalDofList.end());
        }

        ConditionsArrayType& pConditions = rModelPart.Conditions();
        const int nconditions = static_cast<int>(pConditions.size());
#pragma omp parallel for firstprivate(nconditions, ElementalDofList)
        for (int i = 0; i < nconditions; i++)
        {
            typename ConditionsArrayType::iterator it = pConditions.begin() + i;
            const unsigned int this_thread_id = OpenMPUtils::ThisThread();

            // gets list of Dof involved on every element
            pScheme->GetConditionDofList(*(it.base()), ElementalDofList, CurrentProcessInfo);
            dofs_aux_list[this_thread_id].insert(ElementalDofList.begin(), ElementalDofList.end());
        }

        //here we do a reduction in a tree so to have everything on thread 0
        unsigned int old_max = nthreads;
        unsigned int new_max = ceil(0.5*static_cast<double>(old_max));
        while (new_max >= 1 && new_max != old_max)
        {
#pragma omp parallel for
            for (int i = 0; i < static_cast<int>(new_max); i++)
            {
                if (i + new_max < old_max)
                {
                    dofs_aux_list[i].insert(dofs_aux_list[i + new_max].begin(), dofs_aux_list[i + new_max].end());
                    dofs_aux_list[i + new_max].clear();
                }
            }

            old_max = new_max;
            new_max = ceil(0.5*static_cast<double>(old_max));
        }

        DofsArrayType Doftemp;
        BaseType::mDofSet = DofsArrayType();

        Doftemp.reserve(dofs_aux_list[0].size());
        for (auto it = dofs_aux_list[0].begin(); it != dofs_aux_list[0].end(); it++)
        {
                Doftemp.push_back(it->get());
        }
        Doftemp.Sort();

        BaseType::mDofSet = Doftemp;

        //throws an execption if there are no Degrees of freedom involved in the analysis
        if (BaseType::mDofSet.size() == 0)
                KRATOS_THROW_ERROR(std::logic_error, "No degrees of freedom!", "");

        BaseType::mDofSetIsInitialized = true;
        if (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)
        {
            std::cout << "finished setting up the dofs" << std::endl;
        }

#ifdef _OPENMP
        if (mlock_array.size() != 0)
        {
                for (int i = 0; i < static_cast<int>(mlock_array.size()); i++)
                        omp_destroy_lock(&mlock_array[i]);
        }

        mlock_array.resize(BaseType::mDofSet.size());

        for (int i = 0; i < static_cast<int>(mlock_array.size()); i++)
                omp_init_lock(&mlock_array[i]);
#endif

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
    ///@name Protected Operators
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
            const std::size_t empty_key = 2 * equation_size + 10;
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
            for (int iii = 0; iii<nelements; iii++)
            {
                typename ElementsContainerType::iterator i_element = rElements.begin() + iii;
                pScheme->EquationId( *(i_element.base()), ids, CurrentProcessInfo);

                for (std::size_t i = 0; i < ids.size(); i++)
                {
                    if (ids[i] < BaseType::mEquationSystemSize)
                    {
#ifdef _OPENMP
                        omp_set_lock(&mlock_array[ids[i]]);
#endif
                        auto& row_indices = indices[ids[i]];
                        for (auto it = ids.begin(); it != ids.end(); it++)
                        {
                                if (*it < BaseType::mEquationSystemSize)
                                        row_indices.insert(*it);
                        }
#ifdef _OPENMP
                        omp_unset_lock(&mlock_array[ids[i]]);
#endif
                    }
                }
            }

            const int nconditions = static_cast<int>(rConditions.size());
#pragma omp parallel for firstprivate(nconditions, ids)
            for (int iii = 0; iii<nconditions; iii++)
            {
                typename ConditionsArrayType::iterator i_condition = rConditions.begin() + iii;
                pScheme->Condition_EquationId( *(i_condition.base()) , ids, CurrentProcessInfo);
                for (std::size_t i = 0; i < ids.size(); i++)
                {
                    if (ids[i] < BaseType::mEquationSystemSize)
                    {
#ifdef _OPENMP
                        omp_set_lock(&mlock_array[ids[i]]);
#endif
                        auto& row_indices = indices[ids[i]];
                        for (auto it = ids.begin(); it != ids.end(); it++)
                        {
                                if (*it < BaseType::mEquationSystemSize)
                                        row_indices.insert(*it);
                        }
#ifdef _OPENMP
                            omp_unset_lock(&mlock_array[ids[i]]);
#endif
                    }
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
                    Arow_indices[i + 1] = Arow_indices[i] + indices[i].size();

#pragma omp parallel for
            for (int i = 0; i < static_cast<int>(A.size1()); i++)
            {
                    const unsigned int row_begin = Arow_indices[i];
                    const unsigned int row_end = Arow_indices[i + 1];
                    unsigned int k = row_begin;
                    for (auto it = indices[i].begin(); it != indices[i].end(); it++)
                    {
                            Acol_indices[k] = *it;
                            Avalues[k] = 0.0;
                            k++;
                    }

                    std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);

            }

            A.set_filled(indices.size() + 1, nnz);

            Timer::Stop("MatrixStructure");
    }

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

#ifdef _OPENMP
    std::vector< omp_lock_t > mlock_array; // NOTE: In the block builder and solver the variable is protected (not private)
#endif
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void AssembleRHS(
        TSystemVectorType& b,
        const LocalSystemVectorType& RHS_Contribution,
        const Element::EquationIdVectorType& EquationId
        )
    {
        unsigned int local_size = RHS_Contribution.size();

        if (BaseType::mCalculateReactionsFlag == false)
        {
            for (unsigned int i_local = 0; i_local < local_size; i_local++)
            {
                const unsigned int i_global = EquationId[i_local];

                if (i_global < BaseType::mEquationSystemSize) //free dof
                {
                    // ASSEMBLING THE SYSTEM VECTOR
                    double& b_value = b[i_global];
                    const double& rhs_value = RHS_Contribution[i_local];

                    #pragma omp atomic
                    b_value += rhs_value;
                }
            }
        }
        else
        {
            TSystemVectorType& ReactionsVector = *BaseType::mpReactionsVector;
            for (unsigned int i_local = 0; i_local < local_size; i_local++)
            {
                const unsigned int i_global = EquationId[i_local];

                if (i_global < BaseType::mEquationSystemSize) //free dof
                {
                    // ASSEMBLING THE SYSTEM VECTOR
                    double& b_value = b[i_global];
                    const double& rhs_value = RHS_Contribution[i_local];

                    #pragma omp atomic
                    b_value += rhs_value;
                }
                else //fixed dof
                {
                    double& b_value = ReactionsVector[i_global - BaseType::mEquationSystemSize];
                    const double& rhs_value = RHS_Contribution[i_local];

                    #pragma omp atomic
                    b_value += rhs_value;
                }
            }
        }
    }
    
    void AssembleLHS_CompleteOnFreeRows( // NOTE: In the block builder and solver the method is protected (not private)
        TSystemMatrixType& A,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId
        )
    {
        unsigned int local_size = LHS_Contribution.size1();
        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];
            if (i_global < BaseType::mEquationSystemSize)
            {
                for (unsigned int j_local = 0; j_local < local_size; j_local++)
                {
                    int j_global = EquationId[j_local];
                    A(i_global, j_global) += LHS_Contribution(i_local, j_local);
                }
            }
        }
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

}; /* Class ContactResidualBasedEliminationBuilderAndSolver */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_CONTACT_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER  defined */
