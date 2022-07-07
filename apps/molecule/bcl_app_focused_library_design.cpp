// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header for this class
//#include "molecule/bcl_app_focused_library_design.h"

// include headers from the bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "chemistry/bcl_chemistry_configuration_set.h"
#include "chemistry/bcl_chemistry_constitution_graph_converter.h"
#include "chemistry/bcl_chemistry_constitution_set.h"
#include "chemistry/bcl_chemistry_fragment_configuration_shared.h"
#include "chemistry/bcl_chemistry_fragment_constitution_shared.h"
#include "chemistry/bcl_chemistry_fragment_mutate_interface.h"
#include "chemistry/bcl_chemistry_score_function_generic.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "descriptor/bcl_descriptor_combine.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_const_function.h"
#include "math/bcl_math_mutate_combine.h"
#include "math/bcl_math_mutate_decision_node.h"
#include "math/bcl_math_mutate_repeat.h"
#include "math/bcl_math_template_instantiations.h"
#include "mc/bcl_mc_approximator.h"
#include "mc/bcl_mc_temperature_accepted.h"
#include "mc/bcl_mc_temperature_default.h"
#include "mc/bcl_mc_temperature_interface.h"
#include "model/bcl_model_retrieve_interface.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_function.h"
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "opti/bcl_opti_criterion_skipped_steps.h"
#include "opti/bcl_opti_criterion_unimproved.h"
#include "random/bcl_random_uniform_distribution.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_thunk_job.h"
#include "sdf/bcl_sdf_fragment_factory.h"
#include "sdf/bcl_sdf_mdl_handler.h"
#include "util/bcl_util_format.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FocusedLibraryDesign
    //! @brief Application for generating libraries for synthesis using QSAR models and a MCM structure generator
    //!
    //! @author brownbp1, mendenjl, loweew, geanesar
    //! @date 05/09/2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FocusedLibraryDesign :
      public InterfaceRelease
    {

    private:

      // typedefs for convenience
      typedef util::ShPtr< math::MutateDecisionNode< chemistry::FragmentComplete> > Mutates;
      typedef util::ShPtr< math::FunctionInterfaceSerializable< chemistry::FragmentComplete, double> > Score;
      typedef util::ShPtrVector< descriptor::CheminfoProperty> Properties;
      typedef opti::Tracker< chemistry::FragmentComplete, double> Tracker;
      typedef util::ShPtr< chemistry::FragmentComplete> FragmentComplete_p;

    //////////
    // data //
    //////////

      //! flag to control input base fragment
      util::ShPtr< command::FlagInterface> m_StartFragmentFlag;

      //! flag for defining output filename,
      util::ShPtr< command::FlagInterface> m_OutputFilenameFlag;

      //! flag to specify mutate implementations
      util::ShPtr< command::FlagInterface> m_MutateFlag;

      //! flag to specify mutate implementation probabilities
      util::ShPtr< command::FlagInterface> m_MutateProbabilityFlag;

      //! flag controlling the maximum possible number of sequential mutates that can occur between MCM evaluation
      util::ShPtr< command::FlagInterface> m_MaxSequentialMutatesFlag;

      //! flag to specify the MCM score function as a descriptor
      util::ShPtr< command::FlagInterface> m_PropertyScoringFunctionFlag;

      //! flag for descriptors to compute for the final ensemble
      util::ShPtr< command::FlagInterface> m_FinalMetricsFlag;

      //! flag to control number of molecules to be generated
      util::ShPtr< command::FlagInterface> m_NumberMoleculesFlag;

      //! flag to control the number of MC iterations in molecule optimization
      util::ShPtr< command::FlagInterface> m_NumberIterationsFlag;

      //! flag to control the number of maximum allowed consecutive unimproved MC iterations
      util::ShPtr< command::FlagInterface> m_NumberUnimprovedFlag;

      //! flag to control the number of maximum allowed skipped MC iterations
      util::ShPtr< command::FlagInterface> m_NumberSkippedFlag;

      //! flag to control the temperature for the Metropolis criterion
      util::ShPtr< command::FlagInterface> m_MetropolisTemperatureFlag;

      //! flag to maximize score istead of minimize
      util::ShPtr< command::FlagInterface> m_LargerIsBetterFlag;

      //! flag to save all molecules accepted or improved by main MCM
      util::ShPtr< command::FlagInterface> m_SaveAllAcceptedImprovedFlag;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      FocusedLibraryDesign();

      // ThreadManager needs access to private nested classes
      friend class ThreadManager;
      friend class Worker;

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class ThreadManager
      //! @brief manages threads for multithreaded structure generation
      //!
      //! @author mendenjl, geanesar, brownbp1
      //! @date Nov 7, 2013
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class ThreadManager :
        public util::ObjectInterface
      {

      private:

      ///////////
      // Data //
      ///////////

          const size_t                  m_NumberOfMoleculesRequested; //< Number of molecules to build
          size_t                        m_NumberOfMoleculesBuilt;     //< Number of molecules already built
          const size_t                  m_Threads;                    //< Number of threads
          chemistry::FragmentEnsemble   m_Molecules;                  //< The molecules which have been built and are ready for output
          chemistry::ConstitutionSet    m_UniqueConsts;               //< The unique molecules which have been built
          chemistry::ConfigurationSet   m_UniqueConfigs;              //< The unique molecules which have been built
          io::OFStream                  m_OutputStream;               //< Output file to write molecules to
          sched::Mutex                  m_Mutex;                      //< Lock for updating Workers

          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //!
          //! @class Worker
          //! @brief runs the threads for Worker - builds molecules using metropolis monte-carlo routines
          //!
          //! @author brownbp1, mendenjl, geanesar
          //! @date May 09, 2020
          //!
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          struct Worker
          {
            util::SiPtr< ThreadManager>    m_ThreadManager;           //< Pointer to the thread manager, needed so Worker can be updated
            FragmentComplete_p             m_StartFragment;           //< Base fragment to use
            descriptor::CheminfoProperty   m_PropertyScorer;          //< Set objective function with property instead of model
            Mutates                        m_Mutate;                  //< Grow molecules from scaffold
            Score                          m_Score;                   //< Objective function
            Properties                     m_FinalMetrics;            //< Metrics to compute on final ensemble
            Tracker                        m_OptiGoal;                //< Tracker tracking up or down
            size_t                         m_NumberMCIterations;      //< Number of MC iterations
            size_t                         m_NumberMCUnimproved;      //< Number of allowed consecutive unimproved MC iterations
            size_t                         m_NumberMCSkipped;         //< Number of allowed consecutive unimproved MC iterations
            float                          m_MetropolisTemperature;   //< Temperature during Metropolis criterion evaluation
            bool                           m_SaveAllAcceptedImproved; //< Collect every molecule that is accepted or improved
            descriptor::CheminfoProperty   m_BondDruglikeness;        //< MoleculeTotalDruglikeBondEnergy

           // Builds and score the molecule
            void RunThread()
            {
              util::ShPtr< storage::Pair< chemistry::FragmentComplete, double> > last_accepted;
              do
              {
                // create the temperature control
                util::ShPtr< mc::TemperatureInterface> sp_temperature( new mc::TemperatureDefault( float( m_MetropolisTemperature)));

                // create the metropolis
                mc::Metropolis< double> metropolis( sp_temperature, true);

                // create the termination criterion
                opti::CriterionCombine< chemistry::FragmentComplete, double> criterion_combine;

                // insert termination criteria that depends on the total number of MC iterations
                opti::CriterionNumberIterations< chemistry::FragmentComplete, double> maximum_number_iterations( m_NumberMCIterations);
                criterion_combine.InsertCriteria( maximum_number_iterations);

                // insert termination criteria that depends on the total number of unimproved MC iterations
                opti::CriterionUnimproved< chemistry::FragmentComplete, double> maximum_number_unimproved_iterations( m_NumberMCUnimproved);
                criterion_combine.InsertCriteria( maximum_number_unimproved_iterations);

                // insert termination criteria that depends on the total number of skipped MC iterations
                opti::CriterionSkippedSteps< chemistry::FragmentComplete, double> maximum_number_skipped_iterations( m_NumberMCSkipped);
                criterion_combine.InsertCriteria( maximum_number_skipped_iterations);

                // Set up MC method
                mc::Approximator< chemistry::FragmentComplete, double> approximator
                (
                  *m_Score,
                  *m_Mutate,
                  metropolis,
                  criterion_combine,
                  *m_StartFragment,
                  m_OptiGoal
                );

                // assume we start with druglike molecule
                double druglike_mol_activity( ( *m_Score)( approximator.GetTracker().GetCurrent()->First()));

                // tell me about the scaffold
                BCL_MessageStd( "Scaffold properties");
                ReportThreadMoleculeProgress( approximator.GetTracker().GetCurrent()->First());
                BCL_MessageStd( "FLD_Score: " + util::Format()( druglike_mol_activity));

                // run the approximator
                BCL_MessageStd( "MCM BEGIN");
                while( approximator.CanContinue() && approximator.ShouldContinue() && m_ThreadManager->GetNumberMoleculesBuilt() + 1 <= m_ThreadManager->GetNumberMoleculesToBuild())
                {
                  // do next step in approximation and get new molecule
                  approximator.Next();
                  const util::ShPtr< storage::Pair< chemistry::FragmentComplete, double> > &current_mol( approximator.GetTracker().GetCurrent());

                  // Check for undruglike properties of the current molecule
                  if( approximator.GetTracker().GetStatusOfLastStep() == opti::e_Accepted || approximator.GetTracker().GetStatusOfLastStep() == opti::e_Improved)
                  {
                    // save the new molecule
                    last_accepted = current_mol;

                    // tell me about the new mol
                    BCL_MessageStd( "Molecule tracker updated at iteration: " + util::Format()( approximator.GetTracker().GetIteration()));
//                    ReportThreadMoleculeProgress( last_accepted->First());
                    BCL_MessageStd( "FLD_Score: " + util::Format()( last_accepted->Second()));
                  }
                  else
                  {
                    BCL_MessageStd( "FLD_Score: " + util::Format()( current_mol->Second()));
                    BCL_MessageStd( "MCM Rejected");
                  }
                  // save every accepted/improved step of MCM
                  // TODO: add this to approximator to avoid having to manually do the MCM
                  if( last_accepted.IsDefined() && m_SaveAllAcceptedImproved)
                  {
                    m_ThreadManager->m_Mutex.Lock();
                    if( m_ThreadManager->GetNumberMoleculesBuilt() + 1 <= m_ThreadManager->GetNumberMoleculesToBuild())
                    {
//                      // get best molecule and best score
//                      BCL_MessageStd
//                      (
//                        "Chemical perturbations from most recent to least recent: "
//                      );
//                      approximator.GetMu

                      // grab the last accepted molecule
                      chemistry::FragmentComplete best_mol( last_accepted->First());

                      // save the final MCM molecule
                      if( m_ThreadManager->CheckUniqueConfiguration( best_mol))
                      {

                        // assign final FLD score to last accepted molecule
                        linal::Vector< double> best_score( 1, last_accepted->Second());
                        best_mol.StoreProperty( "FLD_Score", best_score);

                        // compute the final metrics and save on molecule
                        for( size_t i( 0), sz( m_FinalMetrics.GetSize()); i < sz; ++i)
                        {
                          const std::string &property_name( ( *m_FinalMetrics( i))->GetAlias());
                          const linal::Vector< float> &property_value( ( *m_FinalMetrics( i))->SumOverObject( best_mol));
                          best_mol.StoreProperty( property_name, property_value);
                        }

                        // add final molecule to ensemble
                        m_ThreadManager->AddMolecule( best_mol);
                        m_ThreadManager->IncreaseMoleculeBuiltCount();
                      }
                    }
                    m_ThreadManager->m_Mutex.Unlock();
                  }
                }
                BCL_MessageStd( "MCM END");

                // save molecules to output, lock it down with a mutex
                if( last_accepted.IsDefined())
                {
                  m_ThreadManager->m_Mutex.Lock();
                  if( m_ThreadManager->GetNumberMoleculesBuilt() + 1 <= m_ThreadManager->GetNumberMoleculesToBuild())
                  {
                    // print approximator endhook
                    // since the PrintEndHook() method is protected and we're sort of
                    // hacking the approximator here, I just copied this.
                    // TODO: add the conditional output accepted/improved molecules to approximator function so that
                    // I can get rid of this hack
                    BCL_MessageStd( "MC Minimization ended");
                    const storage::Vector< size_t> &counts( approximator.GetTracker().GetCounts());
                    const size_t &tot_nr_steps( approximator.GetTracker().GetIteration());
                    const size_t &nr_improved( counts( opti::e_Improved));
                    const size_t &nr_accepted( counts( opti::e_Accepted));
                    const size_t &nr_rejected( counts( opti::e_Rejected));
                    const size_t &nr_skipped( counts( opti::e_Skipped));
                    util::Format format_a, format_b;
                    format_a.W( 5);
                    format_b.W( 5).FFP( 2);
                    BCL_MessageStd( "#MC steps: " + util::Format()( tot_nr_steps));
                    BCL_MessageStd
                    (
                      "#MC steps improved:\t" + format_a( nr_improved) + "\t%" + format_b( 100.0 * nr_improved / tot_nr_steps)
                    );
                    BCL_MessageStd
                    (
                      "#MC steps accepted:\t" + format_a( nr_accepted) + "\t%" + format_b( 100.0 * nr_accepted / tot_nr_steps)
                    );
                    BCL_MessageStd
                    (
                      "#MC steps rejected:\t" + format_a( nr_rejected) + "\t%" + format_b( 100.0 * nr_rejected / tot_nr_steps)
                    );
                    BCL_MessageStd
                    (
                      "#MC steps skipped:\t" + format_a( nr_skipped) + "\t%" + format_b( 100.0 * nr_skipped / tot_nr_steps)
                    );

                    // get best molecule and best score
                    chemistry::FragmentComplete best_mol( last_accepted->First());
                    linal::Vector< double> best_score( 1, last_accepted->Second());
                    best_mol.StoreProperty( "FLD_Score", best_score);

                    // save the final MCM molecule
                    if( m_ThreadManager->CheckUniqueConfiguration( best_mol))
                    {
                      m_ThreadManager->AddMolecule( best_mol);
                      m_ThreadManager->IncreaseMoleculeBuiltCount();
                    }
                  }
                  m_ThreadManager->m_Mutex.Unlock();
                }
              } while( m_ThreadManager->UpdateWorker( *this));
            } // RunThread()

            void ReportThreadMoleculeProgress( const chemistry::FragmentComplete &MOLECULE)
            {
              BCL_MessageStd( "MolWeight: " + util::Format()( descriptor::GetCheminfoProperties().calc_MolWeight->SumOverObject( MOLECULE)( 0)));
              BCL_MessageStd( "# of HBondAcceptors + HBondDonors: " +
                util::Format()( descriptor::GetCheminfoProperties().calc_HbondAcceptor->SumOverObject( MOLECULE)( 0)
                    + descriptor::GetCheminfoProperties().calc_HbondDonor->SumOverObject( MOLECULE)( 0)));
              BCL_MessageStd( "# of NRotBonds: " + util::Format()( descriptor::GetCheminfoProperties().calc_NRotBond->SumOverObject( MOLECULE)( 0)));
              BCL_MessageStd( "LogP: " + util::Format()( descriptor::GetCheminfoProperties().calc_XLogP->SumOverObject( MOLECULE)( 0)));
              BCL_MessageStd( "Bond energy and atom propensity score: " + util::Format()( m_BondDruglikeness->SumOverObject( MOLECULE)( 3)));
              BCL_MessageStd( "# of F: " + util::Format()( descriptor::GetCheminfoProperties().calc_IsF->SumOverObject( MOLECULE)( 0)));
              BCL_MessageStd( "# of Cl: " + util::Format()(descriptor::GetCheminfoProperties().calc_IsCl->SumOverObject( MOLECULE)( 0)));
              BCL_MessageStd( "# of Br: " + util::Format()(descriptor::GetCheminfoProperties().calc_IsBr->SumOverObject( MOLECULE)( 0)));
              BCL_MessageStd( "# of I: " + util::Format()( descriptor::GetCheminfoProperties().calc_IsI->SumOverObject( MOLECULE)( 0)));
              BCL_MessageStd( "# of Halogens: " + util::Format()( descriptor::GetCheminfoProperties().calc_IsHalogen->SumOverObject( MOLECULE)( 0)));
              BCL_MessageStd( "Complexity : " + util::Format()( descriptor::GetCheminfoProperties().calc_MolComplexity->SumOverObject( MOLECULE)( 0)));
            }

          }; // struct Worker

          // Tests to see if the worker should keep running
          bool UpdateWorker( ThreadManager::Worker &WORKER)
          {
            // Lock structure during the modification
            m_Mutex.Lock();

            if( m_NumberOfMoleculesBuilt >= m_NumberOfMoleculesRequested)
            {
              m_Mutex.Unlock();
              return false;
            }

            m_Mutex.Unlock();
            return true;
          } // UpdateWorker()

      public:

          //! param
          ThreadManager(
            const size_t                        NUMBER_THREADS,               //< Number of threads (from scheduler)
            const FragmentComplete_p            &START_FRAGMENT,              //< Base fragment to use
            const std::string                   &OUTPUT_FILENAME,             //< Output filename
            const Mutates                       &MUTATES,                     //< Mutates with probabilities to apply to chemical scaffold
            const descriptor::CheminfoProperty  &PROPERTY_SCORER,             //< Alternative scorer
            const Properties                    &FINAL_METRICS,               //< Metrics to compute on final ensemble
            const Tracker                       &MCM_OPTI_GOAL,               //< Whether to optimize up or down
            const size_t                        NUMBER_OF_MOLECULES,          //< Number to build
            const size_t                        NUMBER_OF_ITERATIONS,         //< Number of MC iterations
            const size_t                        NUMBER_UNIMPROVED_ITERATIONS, //< Number of allowed consecutive unimproved MC iterations
            const size_t                        NUMBER_SKIPPED_ITERATIONS,    //< Number of allowed consecutive unimproved MC iterations
            const float                         METROPOLIS_TEMPERATURE,       //< Temperature during Metropolis criterion evaluation
            bool                                SAVE_ALL_ACCEPTED_IMPROVED,   //< Save every accepted/improved molecule
            const size_t                        MAX_SEQUENTIAL_MUTATES        //< Maximum number of times to apply the mutate prior to scoring
          ) :
            m_NumberOfMoleculesRequested( NUMBER_OF_MOLECULES),
            m_NumberOfMoleculesBuilt( 0),
            m_Threads( std::min( NUMBER_THREADS, NUMBER_OF_MOLECULES))
          {
            // prepare output filestream
            io::File::MustOpenOFStream( m_OutputStream, OUTPUT_FILENAME);

            // set up sequential mutate to perform 1 to N mutates in a row prior to scoring (does not bypass druglikeness filtering)
//            util::ShPtr< math::MutateInterface< chemistry::FragmentComplete> > mutate_repeater
//            (
//              new math::MutateRepeat< chemistry::FragmentComplete>
//              (
//                MUTATES,
//                1,
//                MAX_SEQUENTIAL_MUTATES
//              )
//            );

            // Set up workers
            std::vector< Worker> workers( m_Threads);
            for(
                std::vector< Worker>::iterator itr( workers.begin()), end( workers.end());
                itr != end;
                ++itr
            )
            {
              Worker &worker_ref( *itr);
              worker_ref.m_ThreadManager                = this;
              worker_ref.m_StartFragment                = START_FRAGMENT.HardCopy();
              worker_ref.m_StartFragment->GetCacheMap() = util::ShPtr< descriptor::CacheMap>( new descriptor::CacheMap);
//              worker_ref.m_Mutate                       = mutate_repeater;
              worker_ref.m_Mutate                       = MUTATES;
              worker_ref.m_PropertyScorer               = PROPERTY_SCORER.HardCopy();
              worker_ref.m_Score                        = Score( new chemistry::ScoreFunctionGeneric( worker_ref.m_PropertyScorer) );
              worker_ref.m_FinalMetrics                 = FINAL_METRICS;
              worker_ref.m_OptiGoal                     = MCM_OPTI_GOAL;
              worker_ref.m_NumberMCIterations           = NUMBER_OF_ITERATIONS;
              worker_ref.m_NumberMCUnimproved           = NUMBER_UNIMPROVED_ITERATIONS;
              worker_ref.m_NumberMCSkipped              = NUMBER_SKIPPED_ITERATIONS;
              worker_ref.m_MetropolisTemperature        = METROPOLIS_TEMPERATURE;
              worker_ref.m_SaveAllAcceptedImproved      = SAVE_ALL_ACCEPTED_IMPROVED;
              worker_ref.m_BondDruglikeness             = descriptor::CheminfoProperty( "MoleculeTotalBondEnergy");
            }

            // Allocate space for jobs
            util::ShPtrVector< sched::JobInterface> jobs;
            jobs.AllocateMemory( m_Threads);

            const size_t group_id( 1);
            for( size_t proc_number( 0); proc_number < m_Threads; ++proc_number)
            {
              Worker &worker_ref( workers[ proc_number]);
              jobs.PushBack
              (
                util::ShPtr< sched::JobInterface>
                (
                  new sched::ThunkJob< Worker, void>
                  (
                    group_id,
                    worker_ref,
                    &Worker::RunThread,
                    sched::JobInterface::e_READY,
                    NULL
                  )
                )
              );

              // Submit the jobs to the scheduler
              sched::GetScheduler().RunJob( jobs.LastElement());
            }

            // join all the jobs
            for( size_t proc_number( 0); proc_number < m_Threads; ++proc_number)
            {
              sched::GetScheduler().Join( jobs( proc_number));
            }

            // Close output
            io::File::CloseClearFStream( m_OutputStream);
          }; // ThreadManager()

          // Increase the number of molecules that have been built
          void IncreaseMoleculeBuiltCount()
          {
            BCL_MessageVrb( "Number of molecules built: " + util::Format()( m_NumberOfMoleculesBuilt + 1));
            m_NumberOfMoleculesBuilt++;
          }

          // Return the number of molecules built
          size_t GetNumberMoleculesBuilt()
          {
            return m_NumberOfMoleculesBuilt;
          }

          // Return the number of molecules requested be built
          size_t GetNumberMoleculesToBuild()
          {
            return m_NumberOfMoleculesRequested;
          }

          // Return true if the molecule is unique among built molecules
          bool CheckUniqueConstitution( const chemistry::FragmentComplete &MOLECULE)
          {
            bool unique( m_UniqueConsts.Insert( chemistry::FragmentConstitutionShared( MOLECULE)).second);
            BCL_MessageStd( "Number in ConstitutionSet: " + util::Format()( m_UniqueConsts.GetSize()));
            BCL_MessageVrb( "Unique? : " + util::Format()( unique));
            return unique;
          }

          // Return true if the molecule is unique among built molecules
          bool CheckUniqueConfiguration( const chemistry::FragmentComplete &MOLECULE)
          {
            bool unique( m_UniqueConfigs.Insert( chemistry::FragmentConfigurationShared( MOLECULE)).second);
            BCL_MessageStd( "Number in ConfigurationSet: " + util::Format()( m_UniqueConfigs.GetSize()));
            BCL_MessageVrb( "Unique? : " + util::Format()( unique));
            return unique;
          }

          // Add molecule to final ensemble
          void AddMolecule( const chemistry::FragmentComplete &MOLECULE)
          {
            m_Molecules.PushBack( MOLECULE);
            MOLECULE.WriteMDL( m_OutputStream);
          }

          // Return FragmentEnsemble of the generated molecules
          chemistry::FragmentEnsemble &GetMolecules()
          {
            return m_Molecules;
          }

          //! @brief clone function
          ThreadManager *Clone() const
          {
            BCL_Exit( "ThreadManager cannot be cloned.", -1);
            return NULL;
          }

          //! @brief Get class identifier string
          const std::string &GetClassIdentifier() const
          {
            return GetStaticClassName( *this);
          }

      protected:

          std::istream &Read( std::istream &INSTREAM)
          {
            return INSTREAM;
          }

          std::ostream &Write( std::ostream &OUTSTREAM, const size_t INDENT) const
          {
            return OUTSTREAM;
          }

      }; // class ThreadManager

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      FocusedLibraryDesign *Clone() const
      {
        return new FocusedLibraryDesign( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const
      {
        static std::string s_read_me =
          "FocusedLibraryDesign generates a distribution of new molecules by applying alchemical and medicinal chemistry-like transformations"
          " to a starting scaffold or molecule. A quantitative structure-activity relationship (QSAR) model scores generated molecules "
          " and a Monte Carlo - Metropolis (MCM) algorithm samples the distribution of generated molecules based on QSAR score. Optionally,"
          " an internal MCM-simulated annealing (SA) optimization trial can be performed for each individual transformation, localized to a randomly"
          " selected or manually specified substructure(s). After every transformation, the new molecule is evaluated for drug-likeness according to"
          " a user-specified composite metric. If the molecule is deemed non-druglike, then the molecule is excluded from further analysis and the MCM"
          " move is repeated. When the MCM termination criteria are met, the new molecule is output into an SDF file. Optionally, all generated molecules"
          " that are accepted and/or improved by the MCM engine can be output to the SDF file instead of just the final molecule (note that this will"
          " reduce the QSAR score distribution of the final library). The program terminates when the requested number of unique molecules have been"
          " generated.";
        return s_read_me;
      }

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const
      {
        return "Generates distributions of molecules utilizing chemical perturbations ('mutations') and a property-based score metric, such as a QSAR model.";
      }

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const
      {
        // initialize command to be returned
        util::ShPtr< command::Command> sp_cmd( new command::Command());

      ////////////////////
      // common options //
      ////////////////////

        // add command line options to add/remove hydrogens
        sdf::AddMoleculeIOPrefFlags( *sp_cmd);

        // add AtomMdlLine to molecule
        sdf::MdlHandler::AddAtomMdlLineFlag( *sp_cmd);

      /////////////////////
      // special options //
      /////////////////////

        // insert all the flags and params
        sp_cmd->AddFlag( m_StartFragmentFlag);
        sp_cmd->AddFlag( m_OutputFilenameFlag);
        sp_cmd->AddFlag( m_MutateFlag);
        sp_cmd->AddFlag( m_MutateProbabilityFlag);
        sp_cmd->AddFlag( m_MaxSequentialMutatesFlag);
        sp_cmd->AddFlag( m_PropertyScoringFunctionFlag);
        sp_cmd->AddFlag( m_FinalMetricsFlag);
        sp_cmd->AddFlag( m_NumberMoleculesFlag);
        sp_cmd->AddFlag( m_NumberIterationsFlag);
        sp_cmd->AddFlag( m_NumberUnimprovedFlag);
        sp_cmd->AddFlag( m_NumberSkippedFlag);
        sp_cmd->AddFlag( m_MetropolisTemperatureFlag);
        sp_cmd->AddFlag( m_LargerIsBetterFlag);
        sp_cmd->AddFlag( m_SaveAllAcceptedImprovedFlag);

      ///////////////////
      // default flags //
      ///////////////////

        // default flags
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        // return assembled Command object
        return sp_cmd;
      }

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const
      {

        // setup the starting fragment
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_StartFragmentFlag->GetFirstParameter()->GetValue());
        util::ShPtr< chemistry::FragmentComplete> sp_startfragment
        (
          new chemistry::FragmentComplete( sdf::FragmentFactory::MakeFragment( input))
        );
        io::File::CloseClearFStream( input);

        // setup the scorer
        descriptor::CheminfoProperty property_scorer
        (
          m_PropertyScoringFunctionFlag->GetFlag() ?
              m_PropertyScoringFunctionFlag->GetFirstParameter()->GetValue() :
              "Constant(0.0)"
        );

        // setup the mutates
        auto mutate_input( m_MutateFlag->GetStringList());
        auto mutate_probs( m_MutateProbabilityFlag->GetNumericalList< float>());
        if( mutate_probs.GetSize() != mutate_input.GetSize())
        {
          mutate_probs = storage::Vector< float>( mutate_input.GetSize(), 1.0);
        }
        Mutates mutates( new math::MutateDecisionNode< chemistry::FragmentComplete>() );
        for
        (
            size_t mutate_i( 0), mutate_sz( mutate_input.GetSize());
            mutate_i < mutate_sz;
            ++mutate_i
        )
        {
          util::Implementation< chemistry::FragmentMutateInterface> implementation( mutate_input( mutate_i));
          mutates->AddMutate( *implementation, mutate_probs( mutate_i));
        }

        // setup the final metrics
        auto final_metrics_input( m_FinalMetricsFlag->GetStringList());
        Properties final_metrics;
        for
        (
            size_t prop_i( 0), prop_sz( final_metrics_input.GetSize());
            prop_i < prop_sz;
            ++prop_i
        )
        {
          util::ShPtr< descriptor::CheminfoProperty> property( new descriptor::CheminfoProperty( final_metrics_input( prop_i)));
          final_metrics.PushBack( property);
        }

        // set MCM optimization goal
        opti::Tracker< chemistry::FragmentComplete, double> mcm_opti_goal( opti::e_SmallerIsBetter);
        if( m_LargerIsBetterFlag->GetFlag())
        {
          mcm_opti_goal = opti::e_LargerIsBetter;
        }

        // set save options for intermediate molecules
        bool save_all_accepted_improved( false);
        if( m_SaveAllAcceptedImprovedFlag->GetFlag())
        {
          save_all_accepted_improved = true;
        }

      /////////////////////////
      // parse the arguments //
      /////////////////////////

      /////////////////////////////
      // Prepare rotamer library //
      /////////////////////////////

        // Start track time
        util::Stopwatch threadmanager_timer( "Molecule Building", util::Time( 1, 0), util::Message::e_Standard, true, false);
        threadmanager_timer.Start();

        // Build the molecules using metropolis monte-carlo
        ThreadManager thread_manager
        (
          sched::GetNumberCPUs(),
          sp_startfragment,
          m_OutputFilenameFlag->GetFirstParameter()->GetValue(),
          mutates,
          property_scorer,
          final_metrics,
          mcm_opti_goal,
          m_NumberMoleculesFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
          m_NumberIterationsFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
          m_NumberUnimprovedFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
          m_NumberSkippedFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
          m_MetropolisTemperatureFlag->GetFirstParameter()->GetNumericalValue< float>(),
          save_all_accepted_improved,
          m_MaxSequentialMutatesFlag->GetFirstParameter()->GetNumericalValue< size_t>()
        );

        // End track time
        threadmanager_timer.Stop();
        return 0;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      // instantiate enumerator for PrepareSmallMoleculeEnsemble class
      static const ApplicationType FocusedLibraryDesign_Instance;

    }; // FocusedLibraryDesign

    //! @brief standard constructor
    FocusedLibraryDesign::FocusedLibraryDesign() :
      m_StartFragmentFlag
      (
        new command::FlagStatic
        (
          "start_fragment", "filename for input starting fragment",
          command::Parameter
          (
            "fragment_filename", "filename for input sdf of molecules", ""
          )
        )
      ),
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output_filename", "flag selecting the output file name",
          command::Parameter
          (
            "output_filename_param", "filename for output sdf of molecules"
          )
        )
      ),
      m_MutateFlag
      (
        new command::FlagDynamic
        (
          "mutates",
          "methods with which to mutate molecules; "
          "multiple mutates may be provided; the mutate performed at each iteration will be chosen at random "
          "based on the weighted probability of each mutation.",
          command::Parameter
          (
            "mutate",
            "",
            command::ParameterCheckSerializable
            (
              util::Implementation< chemistry::FragmentMutateInterface>()
            )
          )
        )
      ),
      m_MutateProbabilityFlag
      (
        new command::FlagDynamic
        (
          "mutate_probs",
          "relative probability of performing a mutate; the number of values passed here must "
          "be equal to the number of mutates passed via the 'mutates' flag, otherwise "
          "each mutate will be initialized with an equal probability.",
          command::Parameter
          (
            "mutate",
            "",
            command::ParameterCheckRanged< float>( 0.0, std::numeric_limits< float>::max())
          )
        )
      ),
      m_MaxSequentialMutatesFlag
      (
        new command::FlagStatic
        (
          "max_sequential_mutates", "flag for the maximum number of mutates that can occur between MCM evaluations; "
              "a number between 1 and <max_sequential_mutates> will be randomly selected with uniform probability",
          command::Parameter
          (
            "max_sequential_mutates", "perform multiple mutates in a row prior to druglikeness filtering and MCM evaluation",
            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()), "1"
          )
        )
      ),
      m_PropertyScoringFunctionFlag
      (
        new command::FlagDynamic
        (
          "score_function",
          "the scoring function to use",
          command::Parameter
          (
            "function",
            "the scoring function implementation to use",
            command::ParameterCheckSerializable
            (
//              chemistry::ScoreFunctionGeneric()
              descriptor::CheminfoProperty()
            )
          )
        )
      ),
      m_FinalMetricsFlag
      (
        new command::FlagDynamic
        (
          "final_metrics",
          "descriptors to compute on the final molecules",
          command::Parameter
          (
            "descriptors",
            "",
            command::ParameterCheckSerializable
            (
              descriptor::CheminfoProperty()
            )
          )
        )
      ),
      m_NumberMoleculesFlag
      (
        new command::FlagStatic
        (
          "number_molecules", "flag for number of molecules to generate",
          command::Parameter
          (
            "number_molecules", "total number of molecules",
            command::ParameterCheckRanged< int>( 1, std::numeric_limits< int>::max()), "10"
          )
        )
      ),
      m_NumberIterationsFlag
      (
        new command::FlagStatic
        (
          "number_iterations", "flag for number of MC iterations",
          command::Parameter
          (
            "number_iterations", "maximum number of MC iterations",
            command::ParameterCheckRanged< int>( 1, std::numeric_limits< int>::max()), "100"
          )
        )
      ),
      m_NumberUnimprovedFlag
      (
        new command::FlagStatic
        (
          "number_unimproved", "flag for number of maximum allowed consecutive unimproved MC iterations",
          command::Parameter
          (
            "number_unimproved", "maximum number of allowed consecutive unimproved MC iterations",
            command::ParameterCheckRanged< int>( 1, std::numeric_limits< int>::max()), "100"
          )
        )
      ),
      m_NumberSkippedFlag
      (
        new command::FlagStatic
        (
          "number_skipped", "flag for number of maximum allowed skipped MC iterations",
          command::Parameter
          (
            "number_skipped", "maximum number of allowed skipped MC iterations",
            command::ParameterCheckRanged< int>( 1, std::numeric_limits< int>::max()), "100"
          )
        )
      ),
      m_MetropolisTemperatureFlag
      (
        new command::FlagStatic
        (
          "temperature", "flag for the temperature used in the Metropolis criterion;"
          " units match the units of the score function",
          command::Parameter
          (
            "temperature", "temperature for MCM evaluation",
            command::ParameterCheckRanged< float>( 1.0, std::numeric_limits< float>::max()), "1.0"
          )
        )
      ),
      m_LargerIsBetterFlag
      (
        new command::FlagStatic
        (
          "larger_score_better",
          "sets all MCM objects to increase the score during optimization; "
          "default behavior attempts to decrease the score during optimization"
        )
      ),
      m_SaveAllAcceptedImprovedFlag
      (
        new command::FlagStatic
        (
          "save_all_accepted_improved",
          "save all molecules that are accepted or improved according to the main (outer) MCM; "
          "overall distribution of molecule scores will be skewed worse, but intermediate "
          "structures will be available for analysis"
        )
      )
    {
    }

    const ApplicationType FocusedLibraryDesign::FocusedLibraryDesign_Instance
    (
      GetAppGroups().AddAppToGroup( new FocusedLibraryDesign(), GetAppGroups().e_ChemInfo)
    );

  } // namespace app
} // namespace bcl
