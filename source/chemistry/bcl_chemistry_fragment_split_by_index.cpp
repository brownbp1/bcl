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
#include "chemistry/bcl_chemistry_fragment_split_by_index.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_fragment_split_isolate.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_binary_function_stl_wrapper.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentSplitByIndex::s_Instance
    (
      util::Enumerated< FragmentSplitInterface>::AddInstance( new FragmentSplitByIndex)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor, sets default steps to 4
    FragmentSplitByIndex::FragmentSplitByIndex
    (
      const storage::Vector< size_t> &ATOM_INDICES,
      const bool REMOVE_BONDED_H,
      const bool INVERT,
      const bool BREAK
    ) :
      m_AtomIndices( ATOM_INDICES),
      m_RemoveBondedH( REMOVE_BONDED_H),
      m_Invert( INVERT),
      m_Break( BREAK)
    {
    }

    //! virtual copy constructor
    FragmentSplitByIndex *FragmentSplitByIndex::Clone() const
    {
      return new FragmentSplitByIndex( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &FragmentSplitByIndex::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentSplitByIndex::GetAlias() const
    {
      static const std::string s_name( "Index");
      return s_name;
    }

    //! @brief Get a description for what this class does (used when writing help)
    //! @return a description for what this class does (used when writing help)
    const std::string &FragmentSplitByIndex::GetClassDescription() const
    {
      return GetStaticClassName( *this);
    }

    //! @return the minimum size of fragments
    const size_t FragmentSplitByIndex::GetMinSize() const
    {
      return 0;
    }

    //! @return the indices of all hydrogen atoms bonded to target atoms
    storage::Vector< size_t> FragmentSplitByIndex::GetBondedHydrogenAtoms
    (
      const storage::Vector< size_t> &ATOM_INDICES,
      const AtomVector< AtomComplete> &MOLECULE_ATOMS
    ) const
    {
      // initialize output
      storage::Vector< size_t> h_atoms;

      // find each target atom's bonded hydrogens
      for
      (
          auto atom_itr( ATOM_INDICES.Begin()), atom_itr_end( ATOM_INDICES.End());
          atom_itr != atom_itr_end;
          ++atom_itr
      )
      {
        // skip if selected atom is a hydrogen;
        // yes, this does prevent someone from specifying removal of a hydrogen atoms
        // and its bonded hydrogen if the entire fragment is H2; probably not the desired use-case
        if( MOLECULE_ATOMS( *atom_itr).GetElementType() == GetElementTypes().e_Hydrogen)
        {
          continue;
        }

        // if heavy atom, go over bonds
        for
        (
            auto bond_itr( MOLECULE_ATOMS( *atom_itr).GetBonds().Begin()),
            bond_itr_end( MOLECULE_ATOMS( *atom_itr).GetBonds().End());
            bond_itr != bond_itr_end;
            ++bond_itr
        )
        {
          if( bond_itr->GetTargetAtom().GetElementType() == GetElementTypes().e_Hydrogen)
          {
            h_atoms.PushBack( MOLECULE_ATOMS.GetAtomIndex( bond_itr->GetTargetAtom()));
          }
        }
      }

      // return hydrogen atom indices
      return h_atoms;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief returns an ensemble of fragments of a molecule
    //! @param CONFORMATION molecule of interest
    //! @return an ensemble of common substructures relative to those in a file
    //! TODO: Implement this
    FragmentEnsemble FragmentSplitByIndex::operator()( const ConformationInterface &CONFORMATION) const
    {
      // we will want to construct an atom vector to build our return ensemble
      storage::Vector< sdf::AtomInfo> atominfo( CONFORMATION.GetAtomInfo());
      storage::Vector< sdf::BondInfo> bondinfo( CONFORMATION.GetBondInfo());
      AtomVector< AtomComplete> atoms( atominfo, bondinfo);

      // invert atom indices if desired
      storage::Vector< size_t> keep_indices;
      if( m_Invert)
      {
        // we need to add the bonded hydrogens after assigning the atom indices
        if( m_RemoveBondedH)
        {
          auto h_atoms( GetBondedHydrogenAtoms( m_AtomIndices, atoms));
          storage::Set< size_t> keep_indices_set( m_AtomIndices.Begin(), m_AtomIndices.End());
          keep_indices_set.InsertElements( h_atoms.Begin(), h_atoms.End());
          keep_indices = storage::Vector< size_t>( keep_indices_set.Begin(), keep_indices_set.End());
        }
        else
        {
          keep_indices = m_AtomIndices;
        }
      }
      else
      {
        // if we are not inverting but want to remove bonded H, then we must do it to
        // the m_Atomindices vector before we identify keep indices
        if( m_RemoveBondedH)
        {
          auto h_atoms( GetBondedHydrogenAtoms( m_AtomIndices, atoms));
          storage::Set< size_t> keep_indices_set( m_AtomIndices.Begin(), m_AtomIndices.End());
          keep_indices_set.InsertElements( h_atoms.Begin(), h_atoms.End());
          m_AtomIndices = storage::Vector< size_t>( keep_indices_set.Begin(), keep_indices_set.End());
        }

        // each atom in original conformation
        for( size_t i( 0); i < CONFORMATION.GetSize(); ++i)
        {
          // if not found in removal indices
          if( m_AtomIndices.Find( i) >= m_AtomIndices.GetSize())
          {
            keep_indices.PushBack( i);
          }
        }
      }

      // remove atoms except those we designated to save
      atoms.Reorder( keep_indices);

      // build a new molecule
      FragmentComplete mol( atoms, CONFORMATION.GetName());

      // break if desired
      if( m_Break)
      {
        // accept fragments as small as 1 atom
        FragmentSplitIsolate isolater( 1);

        // return the isolated ensemble
        return isolater( mol);
      }

      return FragmentEnsemble( storage::List< FragmentComplete>( 1, mol));
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool FragmentSplitByIndex::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // read in atom indices
      if( m_AtomIndicesString.size())
      {
        m_AtomIndices.Reset();
        m_AtomIndices = util::SplitStringToNumerical< size_t>( m_AtomIndicesString);
        return true;
      }

      BCL_MessageStd( "No atom indices provided! Exiting without splitting molecules.");
      return false;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer FragmentSplitByIndex::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "splits molecules into fragments by removing specified indices"
      );
      parameters.AddInitializer
      (
        "atom_indices",
        "the 0-indexed atom indices to remove from the input molecules",
        io::Serialization::GetAgent( &m_AtomIndicesString),
        ""
      );
      parameters.AddInitializer
      (
        "include_bonded_h",
        "in addition to removing the specified 'atom_indices', also remove their bonded hydrogen atoms",
        io::Serialization::GetAgent( &m_RemoveBondedH),
        "true"
      );
      parameters.AddInitializer
      (
        "invert",
        "invert the atom index selection prior to atom removal",
        io::Serialization::GetAgent( &m_Invert),
        "false"
      );
      parameters.AddInitializer
      (
        "break",
        "if removing an atom separates a molecule into isolated components, return those isolated "
        "components as separate fragments; by default, a fragment complex is returned in such cases.",
        io::Serialization::GetAgent( &m_Break),
        "false"
      );

      return parameters;
    }

  } // namespace chemistry
} // namespace bcl
