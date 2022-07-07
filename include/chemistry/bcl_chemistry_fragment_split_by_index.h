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

#ifndef BCL_CHEMISTRY_FRAGMENT_SPLIT_BY_INDEX_H_
#define BCL_CHEMISTRY_FRAGMENT_SPLIT_BY_INDEX_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_ensemble.h"
#include "bcl_chemistry_fragment_split_interface.h"
#include "io/bcl_io_serializer.h"
#include "util/bcl_util_enumerated.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentSplitByIndex
    //! @brief Splits molecules into the largest common substructure they possess relative to molecules of an input file
    //!
    //! @see @link example_chemistry_fragment_split_largest_common_substructure.cpp @endlink
    //! @author ben
    //! @date Jul 07, 2022
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentSplitByIndex :
      public FragmentSplitInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the atom indices to remove from the input molecules
      std::string m_AtomIndicesString;
      storage::Vector< size_t> m_AtomIndices;

      //! invert the atom index selection prior to removal
      bool m_Invert;

      //! break a fragment complex into isolated fragments
      bool m_Break;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! virtual copy constructor
      FragmentSplitByIndex *Clone() const;

      //! @brief constructor
      FragmentSplitByIndex
      (
        const storage::Vector< size_t> &ATOM_INDICES = storage::Vector< size_t>(),
        const bool INVERT = false,
        const bool BREAK = false
      );

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief Get a description for what this class does (used when writing help)
      //! @return a description for what this class does (used when writing help)
      const std::string &GetClassDescription() const;

      //! @return the minimum size of fragments
      const size_t GetMinSize() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns an ensemble of fragments of a molecule
      //! @param CONFORMATION molecule of interest
      //! @return an ensemble of common substructures relative to those in a file
      FragmentEnsemble operator()( const ConformationInterface &CONFORMATION) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief reads in molecules from a given file if it is necessary
      void ReadFile() const;

    protected:

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    };

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_SPLIT_BY_INDEX_H_
