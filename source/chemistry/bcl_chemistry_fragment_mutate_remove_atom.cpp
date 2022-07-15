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

// include header of this class
#include "chemistry/bcl_chemistry_fragment_mutate_remove_atom.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_split_largest_component.h"
#include "chemistry/bcl_chemistry_fragment_track_mutable_atoms.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "command/bcl_command_command_state.h"
#include "find/bcl_find_collector_interface.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "random/bcl_random_uniform_distribution.h"

// external includes - sorted alphabetically
#include <iostream>
#include <iterator>
#include <vector>

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentMutateRemoveAtom::s_Instance
    (
      util::Enumerated< FragmentMutateInterface>::AddInstance( new FragmentMutateRemoveAtom())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentMutateRemoveAtom::FragmentMutateRemoveAtom()
    {
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief druglikeness constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    FragmentMutateRemoveAtom::FragmentMutateRemoveAtom
    (
      const descriptor::CheminfoProperty &DRUG_LIKENESS_TYPE,
      const bool &CORINA_CONFS
    )
    {
      m_DrugLikenessType = DRUG_LIKENESS_TYPE;
      m_Corina = CORINA_CONFS;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief full constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
    //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
    //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
    FragmentMutateRemoveAtom::FragmentMutateRemoveAtom
    (
      const descriptor::CheminfoProperty &DRUG_LIKENESS_TYPE,
      const FragmentComplete &SCAFFOLD_FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
      const bool &CORINA_CONFS
    )
    {
      m_DrugLikenessType = DRUG_LIKENESS_TYPE;
      m_ScaffoldFragment = SCAFFOLD_FRAGMENT;
      m_MutableFragments = MUTABLE_FRAGMENTS;
      m_MutableAtomIndices = MUTABLE_ATOM_INDICES;
      m_Corina = CORINA_CONFS;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief local mutate pose-sensitive constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
    //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
    //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
    //! @param MDL property label containing path to protein binding pocket PDB file
    //! @param PROPERTY_SCORER property that will be used to score interactions with protein pocket
    //! @param RESOLVE_CLASHES if true, resolve clashes with specified protein pocket after mutatation
    //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
    FragmentMutateRemoveAtom::FragmentMutateRemoveAtom
    (
      const descriptor::CheminfoProperty &DRUG_LIKENESS_TYPE,
      const FragmentComplete &SCAFFOLD_FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
      const std::string &MDL,
      const descriptor::CheminfoProperty &PROPERTY_SCORER,
      const bool &RESOLVE_CLASHES,
      const storage::Vector< float> &BFACTORS,
      const bool &CORINA_CONFS
    )
    {
      m_DrugLikenessType = DRUG_LIKENESS_TYPE;
      m_ScaffoldFragment = SCAFFOLD_FRAGMENT;
      m_MutableFragments = MUTABLE_FRAGMENTS;
      m_MutableAtomIndices = MUTABLE_ATOM_INDICES;
      m_MDL = MDL;
      m_PropertyScorer = PROPERTY_SCORER;
      m_ResolveClashes = RESOLVE_CLASHES;
      m_BFactors = BFACTORS;
      m_Corina = CORINA_CONFS;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief local clash resolver constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
    //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
    //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
    //! @param MDL property label containing path to protein binding pocket PDB file
    //! @param PROPERTY_SCORER property that will be used to score interactions with protein pocket
    //! @param RESOLVE_CLASHES if true, resolve clashes with specified protein pocket after mutatation
    //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
    FragmentMutateRemoveAtom::FragmentMutateRemoveAtom
    (
      const descriptor::CheminfoProperty &DRUG_LIKENESS_TYPE,
      const FragmentComplete &SCAFFOLD_FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
      const std::string &MDL,
      const bool &RESOLVE_CLASHES,
      const storage::Vector< float> &BFACTORS,
      const bool &CORINA_CONFS
    )
    {
    }

    //! @brief clone constructor
    FragmentMutateRemoveAtom *FragmentMutateRemoveAtom::Clone() const
    {
      return new FragmentMutateRemoveAtom( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentMutateRemoveAtom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentMutateRemoveAtom::GetAlias() const
    {
      static const std::string s_name( "RemoveAtom");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an fragment and generating a new fragment by growing on a valence
    //! @param FRAGMENT small molecule of interest
    //! @return MutateResult with Constitution after the mutate
    math::MutateResult< FragmentComplete> FragmentMutateRemoveAtom::operator()( const FragmentComplete &FRAGMENT) const
    {
      BCL_MessageStd( "RemoveAtom!");

      // these two options are incompatible if enabled simultaneously
      if( ( m_RestrictToBondedHeavy == m_RestrictToBondedH) == true)
      {
        BCL_MessageStd( "Cannot simultaneously restrict atom selection to heavy atoms and hydrogen atoms; returning NULL");
        return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
      }

      // pick an atom to remove
      util::SiPtr< const AtomConformationalInterface> picked_atom;
      if( m_MutableAtomIndices.GetSize() || m_MutableElements.GetSize() || m_MutableFragments.GetSize())
      {
        picked_atom = this->PickAtom( FRAGMENT, false);
      }
      else
      {
        picked_atom = this->PickAtom( FRAGMENT, true);
      }

      // if atom is hydrogen atom, grab the atom to which it is connected
      if( picked_atom->GetElementType() == GetElementTypes().e_Hydrogen && m_RestrictToBondedHeavy)
      {
        if( !picked_atom->GetBonds().GetSize())
        {
          BCL_MessageStd
          (
            "Hydrogen atom selected and restrict to bonded heavy atoms is enabled. "
            "However, no bonded heavy atom is found. Returning NULL."
          );
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }
        picked_atom = util::SiPtr<const AtomConformationalInterface>( picked_atom->GetBonds().Begin()->GetTargetAtom());
      }

      // if restricted to hydrogen atoms get the new picked atom
      else if( m_RestrictToBondedH && picked_atom->GetElementType() != GetElementTypes().e_Hydrogen)
      {
        // loop over bonds and find the hydrogen atom indices
        storage::Vector< size_t> h_indices;
        for
        (
            auto bond_itr( picked_atom->GetBonds().Begin()),
            bond_itr_end( picked_atom->GetBonds().End());
            bond_itr != bond_itr_end;
            ++bond_itr
        )
        {
          if( bond_itr->GetTargetAtom().GetElementType() == GetElementTypes().e_Hydrogen)
          {
            h_indices.PushBack( FRAGMENT.GetAtomVector().GetAtomIndex( bond_itr->GetTargetAtom()));
          }
        }

        // require some hydrogen atoms that can be mutated
        if( !h_indices.GetSize())
        {
          BCL_MessageStd
          (
            "Heavy atom selected and restrict to bonded hydrogen atom is enabled. "
            "However, no bonded hydrogen atom is found. Returning NULL."
          );
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }
        // set new picked atom
        else if( h_indices.GetSize() > size_t( 1))
        {
          h_indices.Shuffle();
        }
        picked_atom = util::SiPtr< const AtomConformationalInterface>( FRAGMENT.GetAtomVector()( h_indices( 0)));
      }

      // removal atom index
      size_t picked_atom_index( FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom));

      // create an atom vector and an index vector from input molecule
      AtomVector< AtomComplete> atoms( FRAGMENT.GetAtomVector());
      storage::Vector< size_t> atom_indices( storage::CreateIndexVector( FRAGMENT.GetSize()));

      // Remove a random atom from the atom vector
      atom_indices.RemoveElements( picked_atom_index, size_t( 1));
      atoms.Reorder( atom_indices);

      // for cleaning and optimizing the new molecule conformer
      FragmentMapConformer cleaner
      (
        m_DrugLikenessType,
        m_MDL,
        FRAGMENT.GetMDLProperty( m_MDL),
        m_PropertyScorer,
        m_ResolveClashes,
        m_BFactors,
        m_Corina,
        storage::Vector< size_t>(),
        m_ChooseBestAlignedConf,
        m_FixGeometry,
        m_ExtendAdjacentAtoms,
        m_ExtendRingAtoms
      );

      // standardize and return
      HydrogensHandler::Remove( atoms);
      return math::MutateResult< FragmentComplete>
      (
        cleaner.Clean
        (
          atoms,
          m_ScaffoldFragment.GetSize() ? m_ScaffoldFragment : FRAGMENT,
          m_DrugLikenessType,
          m_SkipNeutralization,
          m_SkipSaturateH,
          m_SkipSplit
        ),
        *this
      );

      // TODO: consider adding a step to have a certain probability of closing the gap created by this atom removal

      // failed all tries; return null
      return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
    }

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // helper functions //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    io::Serializer FragmentMutateRemoveAtom::GetSerializer() const
    {
      io::Serializer parameters( FragmentMutateInterface::GetSerializer());
      parameters.SetClassDescription
      (
        "Removes an atom from a molecule. If the atom removal splits the molecule into "
        "two or more distinct non-bonded components, save only the largest component."
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool FragmentMutateRemoveAtom::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      // static initialization check
      if( command::CommandState::IsInStaticInitialization())
      {
        return true;
      }

      // call RISH function of the base class
      if( !FragmentMutateInterface::ReadInitializerSuccessHook( LABEL, ERROR_STREAM))
      {
        return false;
      }

      // done
      return true;
    }

  } // namespace chemistry
} // namespace bcl
