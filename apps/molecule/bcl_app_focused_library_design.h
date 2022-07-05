//// (c) Copyright BCL @ Vanderbilt University 2014
//// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
//// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
//// (c)
//// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
//// (c)
//// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
//// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
//// (c)
//// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
//// (c)
//// (c) This file is part of the BCL software suite and is made available under the MIT license.
//// (c)
//
//#ifndef BCL_APP_FOCUSED_LIBRARY_DESIGN_H_
//#define BCL_APP_FOCUSED_LIBRARY_DESIGN_H_
//
//// include headers from the bcl - sorted alphabetically
//#include "app/bcl_app_apps.h"
//#include "chemistry/bcl_chemistry_configuration_set.h"
//#include "chemistry/bcl_chemistry_constitution_graph_converter.h"
//#include "chemistry/bcl_chemistry_constitution_set.h"
//#include "chemistry/bcl_chemistry_fragment_add_med_chem.h"
//#include "chemistry/bcl_chemistry_fragment_alchemy.h"
//#include "chemistry/bcl_chemistry_fragment_configuration_shared.h"
//#include "chemistry/bcl_chemistry_fragment_constitution_shared.h"
//#include "chemistry/bcl_chemistry_fragment_cyclize.h"
//#include "chemistry/bcl_chemistry_fragment_evolve_base.h"
//#include "chemistry/bcl_chemistry_fragment_extend_with_linker.h"
//#include "chemistry/bcl_chemistry_fragment_fluorinate.h"
//#include "chemistry/bcl_chemistry_fragment_grow.h"
//#include "chemistry/bcl_chemistry_fragment_halogenate.h"
//#include "chemistry/bcl_chemistry_fragment_mutate_mcm.h"
//#include "chemistry/bcl_chemistry_fragment_remove_atom.h"
//#include "chemistry/bcl_chemistry_fragment_remove_bond.h"
//#include "chemistry/bcl_chemistry_fragment_ring_swap.h"
//#include "chemistry/bcl_chemistry_fragment_split_interface.h"
//#include "chemistry/bcl_chemistry_fragment_track_mutable_atoms.h"
//#include "chemistry/bcl_chemistry_pick_atom_random.h"
//#include "chemistry/bcl_chemistry_pick_fragment_random.h"
//#include "chemistry/bcl_chemistry_rotamer_library_file.h"
//#include "chemistry/bcl_chemistry_sample_conformations.h"
//#include "chemistry/bcl_chemistry_score_function_generic.h"
//#include "command/bcl_command_app_default_flags.h"
//#include "command/bcl_command_command.h"
//#include "command/bcl_command_flag_dynamic.h"
//#include "command/bcl_command_flag_static.h"
//#include "command/bcl_command_parameter.h"
//#include "command/bcl_command_parameter_check_allowed.h"
//#include "command/bcl_command_parameter_check_file_existence.h"
//#include "command/bcl_command_parameter_check_ranged.h"
//#include "command/bcl_command_parameter_check_serializable.h"
//#include "descriptor/bcl_descriptor_cheminfo_properties.h"
//#include "descriptor/bcl_descriptor_combine.h"
//#include "io/bcl_io_file.h"
//#include "math/bcl_math_const_function.h"
//#include "math/bcl_math_mutate_combine.h"
//#include "math/bcl_math_mutate_decision_node.h"
//#include "math/bcl_math_mutate_repeat.h"
//#include "math/bcl_math_template_instantiations.h"
//#include "mc/bcl_mc_approximator.h"
//#include "mc/bcl_mc_temperature_accepted.h"
//#include "mc/bcl_mc_temperature_default.h"
//#include "mc/bcl_mc_temperature_interface.h"
//#include "model/bcl_model_retrieve_interface.h"
//#include "opti/bcl_opti_criterion_combine.h"
//#include "opti/bcl_opti_criterion_function.h"
//#include "opti/bcl_opti_criterion_number_iterations.h"
//#include "opti/bcl_opti_criterion_skipped_steps.h"
//#include "opti/bcl_opti_criterion_unimproved.h"
//#include "random/bcl_random_uniform_distribution.h"
//#include "sched/bcl_sched_scheduler_interface.h"
//#include "sched/bcl_sched_thunk_job.h"
//#include "sdf/bcl_sdf_fragment_factory.h"
//#include "sdf/bcl_sdf_mdl_handler.h"
//#include "util/bcl_util_format.h"
//#include "util/bcl_util_implementation.h"
//#include "util/bcl_util_sh_ptr.h"
//
//namespace bcl
//{
//  namespace app
//  {
//
//    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    //!
//    //! @class FocusedLibraryDesign
//    //! @brief Application for generating libraries for synthesis using QSAR models and a MCM structure generator
//    //!
//    //! @author brownbp1, mendenjl, loweew, geanesar
//    //! @date 05/09/2020
//    //!
//    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    class BCL_API FocusedLibraryDesign :
//        public InterfaceRelease
//    {
//
//    private:
//
//      //////////
//      // data //
//      //////////
//
//        //! flag to control number of molecules to be generated
//        util::ShPtr< command::FlagInterface> m_NumberMoleculesFlag;
//
//        //! flag to control the number of MC iterations in molecule optimization
//        util::ShPtr< command::FlagInterface> m_NumberIterationsFlag;
//
//        //! flag to control the number of maximum allowed consecutive unimproved MC iterations
//        util::ShPtr< command::FlagInterface> m_NumberUnimprovedFlag;
//
//        //! flag to control the number of maximum allowed skipped MC iterations
//        util::ShPtr< command::FlagInterface> m_NumberSkippedFlag;
//
//        //! flag to control the temperature for the Metropolis criterion
//        util::ShPtr< command::FlagInterface> m_MetropolisTemperatureFlag;
//
//        //! flag to control input base fragment
//        util::ShPtr< command::FlagInterface> m_StartFragmentFlag;
//
//        //! flag to control input mutable fragment within base fragment
//        util::ShPtr< command::FlagInterface> m_MutableFragmentFlag;
//
//        //! flag to control input mutable atoms within base fragment
//        util::ShPtr< command::FlagInterface> m_MutableAtomsFlag;
//
//        //! flag for defining output filename,
//        util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;
//
//        //! flag for defining input fragments
//        util::ShPtr< command::FlagInterface> m_GrowFragmentsFlag;
//
//        //! flag for an alternative score function to just the trained model
//        util::ShPtr< command::FlagInterface> m_PropertyScoringFunctionFlag;
//
//        //! flag for the druglikeness filter to use
//        util::ShPtr< command::FlagInterface> m_DrugLikenessTypeFlag;
//
//        //! flag to split molecules
//        util::ShPtr< command::FlagInterface> m_SplitImplementationFlag;
//
//        //! flag to do an internal MCM optimization
//        util::ShPtr< command::FlagInterface> m_SimulatedAnnealingFlag;
//
//        //! flag to maximize score istead of minimize
//        util::ShPtr< command::FlagInterface> m_LargerIsBetterFlag;
//
//        //! flag to save all molecules accepted or improved by main MCM
//        util::ShPtr< command::FlagInterface> m_SaveAllAcceptedImprovedFlag;
//
//        //! flag to use corina to generate starting conformer
//        util::ShPtr< command::FlagInterface> m_Corina;
//
//        //! flag to set 3D VDW score cutoff
//        util::ShPtr< command::FlagInterface> m_VDWClashCutoffFlag;
//
//        //! flag to enable pose-dependent scoring (default is ligand-based scoring)
//        util::ShPtr< command::FlagInterface> m_PoseDependentFlag;
//
//        //! flag controlling the maximum possible number of sequential mutates that can occur between MCM evaluation
//        util::ShPtr< command::FlagInterface> m_MaxSequentialMutatesFlag;
//
//        //! flags controling relative probabilities of different mutate objects
//        util::ShPtr< command::FlagInterface> m_RingSwapProbFlag;
//        util::ShPtr< command::FlagInterface> m_CyclizeProbFlag;
//        util::ShPtr< command::FlagInterface> m_AlchemyProbFlag;
//        util::ShPtr< command::FlagInterface> m_RemoveAtomProbFlag;
//        util::ShPtr< command::FlagInterface> m_RemoveBondProbFlag;
//        util::ShPtr< command::FlagInterface> m_AddMedChemProbFlag;
//        util::ShPtr< command::FlagInterface> m_FluorinateProbFlag;
//        util::ShPtr< command::FlagInterface> m_HalogenateProbFlag;
//        util::ShPtr< command::FlagInterface> m_ExtendWithLinkerProbFlag;
//
//      ///////////////////////////////////
//      // construction and destruction //
//      ///////////////////////////////////
//
//        //! default constructor
//        FocusedLibraryDesign();
//
//        // instantiate enumerator for PrepareSmallMoleculeEnsemble class
//        static const ApplicationType FocusedLibraryDesign_Instance;
//
//        // ThreadManager needs access to private nested classes
//        friend class ThreadManager;
//        friend class Worker;
//
//      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//      //!
//      //! @class ThreadManager
//      //! @brief manages threads for multithreaded structure generation
//      //!
//      //! @author mendenjl, geanesar, brownbp1
//      //! @date Nov 7, 2013
//      //!
//      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      class ThreadManager :
//          public util::ObjectInterface
//      {
//
//      private:
//
//          const size_t                                            m_NumberOfMoleculesRequested; // Number of molecules to build
//          size_t                                                  m_NumberOfMoleculesBuilt; // Number of molecules already built
//          const size_t                                            m_NumberMCIterations; // Number of iterations in the MC approximator
//          const size_t                                            m_NumberMCUnimproved; // Number of allowed consecutive unimproved MC iterations
//          const size_t                                            m_NumberMCSkipped; // Number of allowed skipped MC iterations
//          const float                                             m_MetropolisTemperature; // Tempterature during Metropolis criterion evaluation
//          const size_t                                            m_Threads; // Number of threads
//          chemistry::FragmentEnsemble                             m_Molecules; // The molecules which have been built and are ready for output
//          chemistry::ConstitutionSet                              m_UniqueConsts; // The unique molecules which have been built
//          chemistry::ConfigurationSet                             m_UniqueConfigs; // The unique molecules which have been built
//          io::OFStream                                            m_OutputStream; // Output file to write molecules to
//          const std::string                                       m_DrugLikenessType; // type of druglikeness filter to use for skipping MCM steps
//          const float                                             m_VDWScoreCutoff; // internal VDW score cutoff for 3D conformer (used to check for mols with reasonable substitutions)
//          const util::Implementation< chemistry::FragmentSplitInterface> m_SplitImplementation; // splitter to use when making fragments for internal MCM optimization
//          const std::string                                       m_PoseDependentMDLProperty; // enable pose-dependent scoring with the receptor indicated by this property
//          const std::string                                       m_PoseDependentResolveClashes; // resolve clashes between ligand and receptor
//          const size_t                                            m_MaxSequentialMutates;
//          const float                                             m_RingSwapProb;
//          const float                                             m_CyclizeProb;
//          const float                                             m_AlchemyProb;
//          const float                                             m_RemoveAtomProb;
//          const float                                             m_RemoveBondProb;
//          const float                                             m_AddMedChemProb;
//          const float                                             m_FluorinateProb;
//          const float                                             m_HalogenateProb;
//          const float                                             m_ExtendWithLinkerProb;
//          sched::Mutex                                            m_Mutex; // Lock for updating Workers
//
//          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//          //!
//          //! @class Worker
//          //! @brief runs the threads for Worker - builds molecules using metropolis monte-carlo routines
//          //!
//          //! @author brownbp1, mendenjl, geanesar
//          //! @date May 09, 2020
//          //!
//          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//          struct Worker
//          {
//            // Rotamer library to use - read in at Main()
//            util::ShPtr< chemistry::FragmentComplete>                                                m_StartFragment; // Base fragment to use
//            chemistry::FragmentEnsemble                                                              m_MutableFragment; // mutable fragment in base fragment
//            storage::Vector< size_t>                                                                 m_MutableAtomIndices; // mutable atoms in base fragment
//            descriptor::CheminfoProperty                                                             m_PropertyScorer; // Set objective function with property instead of model
//            size_t                                                                                   m_NumberMCIterations; // Number of MC iterations
//            size_t                                                                                   m_NumberMCUnimproved; // Number of allowed consecutive unimproved MC iterations
//            size_t                                                                                   m_NumberMCSkipped; // Number of allowed consecutive unimproved MC iterations
//            float                                                                                    m_MetropolisTemperature; // Temperature during Metropolis criterion evaluation
//            opti::Tracker< chemistry::FragmentComplete, double>                                      m_OptiGoal;
//            bool                                                                                     m_SaveAllAcceptedImproved;
//            std::string                                                                              m_ConformationComparer; // Conformation comparer
//            util::ShPtr< math::FunctionInterfaceSerializable< chemistry::FragmentComplete, double> > m_Score; // Objective function
//            util::ShPtr< math::MutateInterface< chemistry::FragmentComplete> >                       m_Mutate; // Grow molecules from scaffold
//            bool                                                                                     m_Corina; // enables corina conformers during cleaning
//            util::SiPtr< ThreadManager>                                                              m_ThreadManager; // Pointer to the thread manager, needed so Worker can be updated
//
//            // Builds and score the molecule
//            void RunThread();
//
//            // Print to screen the properties of our last accepted molecule
//            void ReportThreadMoleculeProgress();
//
//          }; // struct Worker
//
//      public:
//
//          //! @brief constructor
//          ThreadManager(
//            util::ShPtr< chemistry::FragmentComplete>              START_FRAGMENT, // Base fragment to use
//            chemistry::FragmentEnsemble                            MUTABLE_FRAGMENT, // mutable fragment in base fragment
//            storage::Vector< size_t>                               MUTABLE_ATOM_INDICES, // mutable atom indices in base fragment
//            util::ShPtr< chemistry::FragmentEnsemble>              FRAGMENT_POOL, // Fragments to add to base fragment
//            descriptor::CheminfoProperty                           PROPERTY_SCORER, // alternative scorer
//            bool                                                   INTERNAL_MCM_OPTI,
//            opti::Tracker< chemistry::FragmentComplete, double>    MCM_OPTI_GOAL,
//            bool                                                   SAVE_ALL_ACCEPTED_IMPROVED,
//            const size_t                                           &NUMBER_OF_MOLECULES, // Number to build
//            const size_t                                           &NUMBER_OF_ITERATIONS, // Number of MC iterations
//            const size_t                                           &NUMBER_UNIMPROVED_ITERATIONS, // Number of allowed consecutive unimproved MC iterations
//            const size_t                                           &NUMBER_SKIPPED_ITERATIONS, // Number of allowed consecutive unimproved MC iterations
//            const float                                            &METROPOLIS_TEMPERATURE, // Temperature during Metropolis criterion evaluation
//            const size_t                                           &NUMBER_THREADS, // Number of threads (from scheduler)
//            const std::string                                      &OUTPUT_FILENAME,
//            const std::string                                      &DRUG_LIKENESS_TYPE,
//            const float                                            &VDW_SCORE_CUTOFF,
//            const util::Implementation< chemistry::FragmentSplitInterface> &SPLIT_IMPLEMENTATION,
//            const size_t                                           &MAX_SEQUENTIAL_MUTATES,
//            const float                                            &RING_SWAP_PROB,
//            const float                                            &CYCLIZE_PROB,
//            const float                                            &ALCHEMY_PROB,
//            const float                                            &REMOVE_ATOM_PROB,
//            const float                                            &REMOVE_BOND_PROB,
//            const float                                            &ADD_MEDCHEM_PROB,
//            const float                                            &FLUORINATE_PROB,
//            const float                                            &HALOGENATE_PROB,
//            const float                                            &EXTEND_WITH_LINKER_PROB,
//            const std::string                                      &POSE_DEPENDENT_MDL_PROPERTY, // pose-dependent scoring
//            const std::string                                      &POSE_DEPENDENT_RESOLVE_CLASHES, // resolve clashes
//            const bool                                             &CORINA_CONFS // enables cornina conformers during cleaning
//          );
//
//          //! @brief clone function
//          ThreadManager *Clone() const;
//
//          //! @brief Get class identifier string
//          const std::string &GetClassIdentifier() const;
//
//          //! @brief
//          size_t GetNumberMoleculesBuilt();
//
//          //! @brief
//          size_t GetNumberMoleculesToBuild();
//
//          //! @brief
//          size_t GetNumberMCIterations();
//
//          //! @brief
//          size_t GetNumberMCUnimproved();
//
//          //! @brief
//          size_t GetNumberMCSkipped();
//
//          //! @brief Return FragmentEnsemble of the generated molecules
//          chemistry::FragmentEnsemble &GetMolecules()
//
//          //! @brief Tests to see if the worker should keep running
//          bool UpdateWorker( Worker &WORKER);
//
//          //! @brief Increase the number of molecules that have been built
//          void IncreaseMoleculeBuiltCount();
//
//          //! @brief Return true if this molecule is unique from the generated molecules at the constitutional level
//          bool CheckUniqueConstitution( const chemistry::FragmentComplete &MOLECULE);
//
//          //! @brief Return true if this molecule is unique from the generated molecules at the configurational level
//          bool CheckUniqueConfiguration( const chemistry::FragmentComplete &MOLECULE);
//
//          //! @brief Save a molecule to the growing ensemble
//          void AddMolecule( const chemistry::FragmentComplete &MOLECULE);
//
//      protected:
//
//          std::istream &Read( std::istream &INSTREAM);
//
//          std::ostream &Write( std::ostream &OUTSTREAM, const size_t INDENT) const;
//
//      }; // class ThreadManager
//
//    public:
//
//      //! @brief Clone function
//      //! @return pointer to new FoldProtein
//      FocusedLibraryDesign *Clone() const;
//
//    /////////////////
//    // data access //
//    /////////////////
//
//      //! @brief returns class name of the object behind a pointer or the current object
//      //! @return the class name
//      const std::string &GetClassIdentifier() const;
//
//      //! @brief returns readme information
//      //! @return string containing information about application
//      const std::string &GetReadMe() const;
//
//      //! @brief get a description for the app
//      //! @return a brief (no more than 3 line) description for the application
//      std::string GetDescription() const;
//
//    //////////////////////
//    // helper functions //
//    //////////////////////
//
//    ////////////////
//    // operations //
//    ////////////////
//
//      //! @brief initializes the command object for that executable
//      util::ShPtr< command::Command> InitializeCommand() const;
//
//      //! @brief the Main function
//      //! @return error code - 0 for success
//      int Main() const;
//
//    //////////////////////
//    // input and output //
//    //////////////////////
//
//    protected:
//
//      //! @brief read from std::istream
//      //! @param ISTREAM input stream
//      //! @return istream which was read from
//      std::istream &Read( std::istream &ISTREAM);
//
//      //! @brief write to std::ostream
//      //! @param OSTREAM output stream to write to
//      //! @param INDENT number of indentations
//      //! @return output stream which was written to
//      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;
//
//    }; // FocusedLibraryDesign
//
//  } // namespace app
//} // namespace bcl
//
//#endif BCL_APP_FOCUSED_LIBRARY_DESIGN_H_
