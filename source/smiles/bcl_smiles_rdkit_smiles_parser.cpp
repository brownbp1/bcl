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
#include <io/bcl_io_file.h>
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "smiles/bcl_smiles_rdkit_smiles_parser.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_rdkit_mol_utils.h"
#include "command/bcl_command_command.h"
#include "sdf/bcl_sdf_mdl_entry_types.h"
#include "sdf/bcl_sdf_mdl_header.h"
#include "sdf/bcl_sdf_mdl_line_types.h"
#include "sdf/bcl_sdf_mdl_property.h"

// external includes - sorted alphabetically
#include "ForceField/ForceField.h"
#include "GraphMol/DistGeomHelpers/Embedder.h"
#include "GraphMol/FileParsers/FileParsers.h"
#include "GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h"
#include "GraphMol/ForceFieldHelpers/MMFF/Builder.h"
#include "GraphMol/ForceFieldHelpers/UFF/AtomTyper.h"
#include "GraphMol/ForceFieldHelpers/UFF/Builder.h"
#include "GraphMol/GraphMol.h"
#include "GraphMol/ROMol.h"
#include "GraphMol/SmilesParse/SmartsWrite.h"
#include "GraphMol/SmilesParse/SmilesParse.h"
#include "GraphMol/SmilesParse/SmilesWrite.h"

// stl includes
#include <fstream>
#include <iostream>

namespace bcl
{
  namespace smiles
  {

    //! @brief standard constructor
    RdkitSmilesParser::RdkitSmilesParser() :
      m_Header( storage::Vector< std::string>()),
      m_FileData( storage::Vector< storage::Map< std::string, std::string>>()),
          m_Molecules( storage::Vector< chemistry::FragmentComplete>()),
          m_FileType( e_SMILES)
        {
        }

    //! @brief constructor from an input stream
    RdkitSmilesParser::RdkitSmilesParser( std::istream &ISTREAM) :
      m_Header( storage::Vector< std::string>()),
      m_FileData( storage::Vector< storage::Map< std::string, std::string>>()),
          m_Molecules( storage::Vector< chemistry::FragmentComplete>()),
          m_FileType( e_SMILES)
        {
          ReadFile( ISTREAM);
    }

    //! @brief constructor from a pre-read set of lines
    RdkitSmilesParser::RdkitSmilesParser( const storage::List< std::string> &LINES) :
      m_Header( storage::Vector< std::string>()),
      m_FileData( storage::Vector< storage::Map< std::string, std::string>>()),
          m_Molecules( storage::Vector< chemistry::FragmentComplete>()),
          m_FileType( e_SMILES)
        {
          ReadFile( LINES.Begin(), LINES.End());
    }

    //! @brief constructor from an input stream and file type specification
    RdkitSmilesParser::RdkitSmilesParser( std::istream &ISTREAM, const FileType &TYPE) :
      m_Header( storage::Vector< std::string>()),
      m_FileData( storage::Vector< storage::Map< std::string, std::string>>()),
          m_Molecules( storage::Vector< chemistry::FragmentComplete>()),
          m_FileType( TYPE)
        {
          ReadFile( ISTREAM);
    }

    //! @brief constructor from a pre-read set of lines and file type specification
    RdkitSmilesParser::RdkitSmilesParser( const storage::List< std::string> &LINES, const FileType &TYPE) :
      m_Header( storage::Vector< std::string>()),
      m_FileData( storage::Vector< storage::Map< std::string, std::string>>()),
          m_Molecules( storage::Vector< chemistry::FragmentComplete>()),
          m_FileType( TYPE)
        {
          ReadFile( LINES.Begin(), LINES.End());
    }

//    //! @brief write constructor
//    RdkitSmilesParser::RdkitSmilesParser
//    (
//      std::ostream &OSTREAM,
//      const storage::Vector< chemistry::FragmentComplete> &MOLECULES,
//      const FileType &TYPE
//    ) :
//    m_Header( storage::Vector< std::string>()),
//    m_FileData( storage::Vector< storage::Map< std::string, std::string>>()),
//    m_Molecules( MOLECULES),
//    m_FileType( TYPE)
//    {
//      WriteFile( OSTREAM);
//    }

    //! @brief virtual copy constructor
    RdkitSmilesParser *RdkitSmilesParser::Clone() const
    {
      return new RdkitSmilesParser( *this);
    }

    //! @brief return the file type
    //! returns the file type as a string
    const std::string &RdkitSmilesParser::GetFileTypeAsString( const FileType &ENUM)
    {
      static const std::string s_Names[ size_t( s_TotalFileTypes) + 1] =
      {
          "SMILES",
          "SMARTS",
          GetStaticClassName< FileType>()
      };
      return s_Names[ ENUM];
    }

    //! @brief Read the SMILES/SMARTS file using iterators to strings
    //! @param LINE_BEGIN a line that represents a header/counts line
    //! @param LINE_END one-past-end of possible lines
    //! @param READ_TO_DATA_ONLY if true, do not store molecules in memory, but populate
    //! the file data member; useful if you are reading many molecules but working on
    //! one at a time.
    //! @return an iterator one line past the end of the file
    storage::List< std::string>::const_iterator RdkitSmilesParser::ReadFile
    (
      const storage::List< std::string>::const_iterator &LINE_BEGIN,
      const storage::List< std::string>::const_iterator &LINE_END,
      const bool READ_TO_DATA_ONLY
    )
    {
      // clear data
      ResetAll();

      // sanity check on file entries
      storage::List< std::string>::const_iterator itr( LINE_BEGIN);
      if( LINE_BEGIN == LINE_END)
      {
        BCL_MessageStd( "RdkitSmilesParser: nothing to read");
        return LINE_BEGIN;
      }

      // first row is header; every other row is a molecule entry
      m_Header = util::SplitString( ( *itr), ",\n\t\r\0");

      // check if the first column of the header indicates that a SMILES or SMARTS string is expected in subsequent rows
      if( !ContainsValidHeader())
      {
        return LINE_BEGIN;
      }
      ++itr;

      // parse molecules and associated properties line-by-line
      size_t line_index( 0);
      for( ; itr != LINE_END; ++itr, ++line_index)
      {
        // skip blank lines
        if( !ContainsNonspaceCharacters( *itr))
        {
          BCL_MessageStd( "Empty line on row " + util::Format()( line_index) + " of input file - skipping");
          continue;
        }

        // skip empty SMILES/SMARTS columns
        const storage::Vector< std::string> line_str( util::SplitString( ( *itr), ",\n\t\r\0"));
        const std::string &mol_str( line_str( 0));
        if( mol_str.empty())
        {
          continue;
        }

        // add properties to our molecule
        size_t col_index( 0);
        for
        (
            storage::Vector< std::string>::const_iterator
            header_itr( m_Header.Begin()), header_itr_end( m_Header.End());
            header_itr != header_itr_end;
            ++header_itr, ++col_index
        )
        {
          // includes the SMILES or SMARTS string as well as the properties
          m_FileData.PushBack( storage::Map< std::string, std::string>());
          m_FileData.LastElement()[ *header_itr] = line_str( col_index);
        }

        // generate our molecule from the SMILES/SMARTS string
        if( !READ_TO_DATA_ONLY)
        {
          chemistry::FragmentComplete mol
          (
            m_FileType == e_SMILES ?
                ConvertSMILESToMOL( mol_str) :
                ConvertSMARTSToMOL( mol_str)
          );

          size_t col_index( 0);
          for
          (
              storage::Vector< std::string>::const_iterator
              header_itr( m_Header.Begin()), header_itr_end( m_Header.End());
              header_itr != header_itr_end;
              ++header_itr, ++col_index
          )
          {
            mol.GetStoredPropertiesNonConst().SetMDLProperty( *header_itr, line_str( col_index));
          }
          m_Molecules.PushBack( mol);
        }
      }
      return itr;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @param READ_TO_DATA_ONLY if true, do not store molecules in memory, but populate
    //! the file data member; useful if you are reading many molecules but working on
    //! one at a time.
    //! @return istream which was read from
    std::istream &RdkitSmilesParser::ReadFile( std::istream &ISTREAM, const bool READ_TO_DATA_ONLY)
    {
      // buffer for the stream
      storage::List< std::string> lines;

      // check that istream is good
      if( !ISTREAM.good() || ISTREAM.eof())
      {
        BCL_MessageCrt( "Error reading SMILES: passed bad istream");
        return ISTREAM;
      }

      while( !ISTREAM.eof())
      {
        lines.PushBack( std::string());
        std::getline( ISTREAM, lines.LastElement());
      }
      ReadFile( lines.Begin(), lines.End(), READ_TO_DATA_ONLY);
      return ISTREAM;
    }

    //! @brief write to std::ostream in mdl format
    //! @param OSTREAM the stream to write to
    //! @param FORCE_WRITE_ATOM_TYPES whether to force writing BCL atom types (and other info)
    std::ostream &RdkitSmilesParser::WriteFile
    (
      std::ostream &OSTREAM
    ) const
    {
      WriteHeaderLine( OSTREAM);
      for( auto mol_itr( m_Molecules.Begin()), mol_itr_end( m_Molecules.End()); mol_itr != mol_itr_end; ++mol_itr)
      {
        WriteNonHeaderLine( OSTREAM, *mol_itr, true);
      }
      if( !ContainsValidHeader())
      {
        OSTREAM.clear();
        return OSTREAM;
      }
      return OSTREAM;
    }

    //! @brief write header lines to std::ostream to top the output SMILES or SMARTS file
    //! @param OSTREAM the stream to write to
    std::ostream &RdkitSmilesParser::WriteHeaderLine( std::ostream &OSTREAM, const std::string &DELIMITER) const
    {
      storage::Vector< std::string>::const_iterator header_itr( m_Header.Begin()), header_itr_end( m_Header.End());
      OSTREAM << *header_itr;
      ++header_itr;
      for
      (
          ;
          header_itr != header_itr_end;
          ++header_itr
      )
      {
        OSTREAM << DELIMITER;
        OSTREAM << *header_itr;
      }
      OSTREAM << std::endl;
      return OSTREAM;
    }

    //! @brief write to std::ostream to fill molecule and property data into the output SMILES or SMARTS file
    //! @param OSTREAM the stream to write to
    //! @param MOL molecule data to write to file; only properties in m_Headers will be written
    std::ostream &RdkitSmilesParser::WriteNonHeaderLine
    (
      std::ostream &OSTREAM,
      const chemistry::FragmentComplete &MOL,
      const bool INCLUDE_PROPERTIES,
      const std::string &DELIMITER
    ) const
    {
      // write SMILES or SMARTS string from molecule
      std::string mol_str
      (
        m_FileType == e_SMILES ?
          ConvertMolToSMILES( MOL) :
          ConvertMolToSMARTS( MOL)
      );
      OSTREAM << mol_str;

      // add properties to subsequent columns
      if( INCLUDE_PROPERTIES)
      {
        storage::Vector< std::string>::const_iterator header_itr( m_Header.Begin()), header_itr_end( m_Header.End());
        ++header_itr; // first column should always be SMILES/SMARTS string
        for
        (
            ;
            header_itr != header_itr_end;
            ++header_itr
        )
        {
          OSTREAM << DELIMITER;
          OSTREAM << MOL.GetStoredProperties().GetMDLProperty( *header_itr);
        }
        OSTREAM << std::endl;
      }
      else
      {
        OSTREAM << std::endl;
      }
      return OSTREAM;
    }

    //! @brief build a molecule from the file data member
    //! @param MOL_INDEX the index of the file data that will be filled
    //! @returns molecule filled from data SMILES/SMARTS and any properties
    chemistry::FragmentComplete RdkitSmilesParser::FillMolFromData( const size_t MOL_INDEX) const
    {
      // get the molecule string id
      if( !m_Header.GetSize() || !ContainsValidHeader())
      {
        BCL_MessageStd( "Invalid header on input file. Returning null.");
        return chemistry::FragmentComplete();
      }
      const std::string &id_str( m_Header( 0));

      // get the appropriate line from the file data
      if( m_FileData.GetSize() <= MOL_INDEX)
      {
        BCL_MessageStd( "Invalid molecule index specified. Returning null.");
        return chemistry::FragmentComplete();
      }
      const std::string &mol_str( m_FileData( MOL_INDEX).Find( id_str)->second);

      // no empty molecule strings
      if( mol_str.empty())
      {
        BCL_MessageStd( "Empty molecule. Returning null.");
        return chemistry::FragmentComplete();
      }

      // build molecule from SMILES or SMARTS string
      chemistry::FragmentComplete mol
      (
        m_FileType == e_SMILES ?
            ConvertSMILESToMOL( mol_str) :
            ConvertSMARTSToMOL( mol_str)
      );

      // add properties to molecule
      for
      (
          storage::Vector< std::string>::const_iterator
          header_itr( m_Header.Begin()), header_itr_end( m_Header.End());
          header_itr != header_itr_end;
          ++header_itr
      )
      {
        mol.GetStoredPropertiesNonConst().SetMDLProperty( *header_itr, m_FileData( MOL_INDEX).Find( *header_itr)->second);
      }
      return mol;
    }

    //! @brief write to std::ostream directly from a molecule without object construction
    //! @param OSTREAM the stream to write to
    //! @param MOL the molecule whose data we will write
    //! @returns a stream with the new molecule SMILES; note that this write function
    //! does not write properties to file because no header line is produced
    std::ostream &RdkitSmilesParser::WriteSMILESFromMol
    (
      std::ostream &OSTREAM,
      const chemistry::FragmentComplete &MOL
    )
    {
      const std::string mol_str( ConvertMolToSMILES( MOL));
      OSTREAM << mol_str;
      return OSTREAM;
    }

    //! @brief write to std::ostream directly from a molecule without object construction
    //! @param OSTREAM the stream to write to
    //! @param NAME the molecule name
    //! @param ATOM_INFO the atominfo describing the molecule
    //! @param BOND_INFO the bondingo describing the molecule
    //! @returns a stream with the new molecule SMILES; note that this write function
    //! does not write properties to file because no header line is produced
    std::ostream &RdkitSmilesParser::WriteSMILESFromMolInfo
    (
      std::ostream &OSTREAM,
      const std::string &NAME,
      const storage::Vector< sdf::AtomInfo> &ATOM_INFO,
      const storage::Vector< sdf::BondInfo> &BOND_INFO
    )
    {
      const chemistry::AtomVector< chemistry::AtomComplete> atom_v( ATOM_INFO, BOND_INFO);
      const chemistry::FragmentComplete mol( atom_v, NAME);
      const std::string mol_str( ConvertMolToSMILES( mol));
      OSTREAM << mol_str;
      return OSTREAM;
    }

    //! @brief write to std::ostream directly from a molecule without object construction
    //! @param OSTREAM the stream to write to
    //! @param MOL the molecule whose data we will write
    //! @returns a stream with the new molecule SMARTS; note that this write function
    //! does not write properties to file because no header line is produced
    std::ostream &RdkitSmilesParser::WriteSMARTSFromMol
    (
      std::ostream &OSTREAM,
      const chemistry::FragmentComplete &MOL
    )
    {
      const std::string mol_str( ConvertMolToSMARTS( MOL));
      OSTREAM << mol_str;
      return OSTREAM;
    }

    //! @brief write to std::ostream directly from a molecule without object construction
    //! @param OSTREAM the stream to write to
    //! @param NAME the molecule name
    //! @param ATOM_INFO the atominfo describing the molecule
    //! @param BOND_INFO the bondingo describing the molecule
    //! @returns a stream with the new molecule SMARTS; note that this write function
    //! does not write properties to file because no header line is produced
    std::ostream &RdkitSmilesParser::WriteSMARTSFromMolInfo
    (
      std::ostream &OSTREAM,
      const std::string &NAME,
      const storage::Vector< sdf::AtomInfo> &ATOM_INFO,
      const storage::Vector< sdf::BondInfo> &BOND_INFO
    )
    {
      const chemistry::AtomVector< chemistry::AtomComplete> atom_v( ATOM_INFO, BOND_INFO);
      const chemistry::FragmentComplete mol( atom_v, NAME);
      const std::string mol_str( ConvertMolToSMARTS( mol));
      OSTREAM << mol_str;
      return OSTREAM;
    }

    //! @brief convert input molecule into SMILES string
    //! @param MOL molecule to be converted to SMILES string
    std::string RdkitSmilesParser::ConvertMolToSMILES( const chemistry::FragmentComplete &MOL)
    {
      // generate our molecule from the SMILES/SMARTS string
      auto rdkit_mol( chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( MOL, chemistry::GetElementTypes().e_Undefined, true, true));
      return ::RDKit::MolToSmiles( *rdkit_mol);
    }

    //! @brief convert input molecule into SMARTS string
    //! @param MOL molecule to be converted to SMARTS string
    std::string RdkitSmilesParser::ConvertMolToSMARTS( const chemistry::FragmentComplete &MOL)
    {
      auto rdkit_mol( chemistry::RdkitMolUtils::FragmentCompleteToRDKitRWMol( MOL, chemistry::GetElementTypes().e_Undefined, true, true));
      return ::RDKit::MolToSmarts( *rdkit_mol);
    }

    //! @brief convert SMILES string to molecule
    //! @param MOL molecule to be converted to SMILES string
    chemistry::FragmentComplete RdkitSmilesParser::ConvertSMILESToMOL
    (
      const std::string &SMILES,
      const bool IS_SMARTS,
      const bool ADD_H,
      const bool GEN3D,
      const size_t GEOOPT_ITERATIONS,
      const std::string &GEOOPT_FF,
      const float GEOOPT_NONBONDED_THRESH,
      const bool GEOOPT_IGNORE_INTERFRAG_INTERACTIONS
    )
    {
      // generate our molecule from the SMILES/SMARTS string
      ::RDKit::RWMol *rdkit_mol
      (
        IS_SMARTS ?
          ::RDKit::SmartsToMol( SMILES) :
          ::RDKit::SmilesToMol( SMILES)
      );

      // check validity of SMILES/SMARTS syntax
      ::RDKit::MolOps::sanitizeMol( *rdkit_mol);
      if( rdkit_mol == nullptr)
      {
        BCL_MessageStd( "RDKit could not sanitize the molecule upon conversion from SMILES. Returning null...");
        return chemistry::FragmentComplete();
      }

      // add explicit hydrogen atoms
      if( ADD_H)
      {
        ::RDKit::MolOps::addHs( *rdkit_mol);
      }

      // generate a 3D conformer
      if( GEN3D)
      {
        ::RDKit::DGeomHelpers::EmbedMolecule( *rdkit_mol);
        std::ofstream ofs( "foo.mol");
        ofs << ::RDKit::MolToMolBlock( *rdkit_mol);
      }

      // geometry optimize the structure with a molecular mechanics force field
      if( GEOOPT_ITERATIONS)
      {
        // initialize force field
        ::ForceFields::ForceField *ff;

        // setup force field
        if( GEOOPT_FF == "MMFF94" || GEOOPT_FF == "MMFF94s")
        {
          ::RDKit::MMFF::MMFFMolProperties mmff_mol_properties( *rdkit_mol, GEOOPT_FF == "MMFF94" ? "MMFF94" : "MMFF94s");
          if( !mmff_mol_properties.isValid())
          {
            BCL_MessageStd( "Invalid MMFF molecule properties. Skipping geometry optimization.");
            return *( chemistry::RdkitMolUtils::RDKitRWMolToFragmentComplete( *rdkit_mol, SMILES));
          }
          ff = ::RDKit::MMFF::constructForceField( *rdkit_mol, GEOOPT_NONBONDED_THRESH, -1, GEOOPT_IGNORE_INTERFRAG_INTERACTIONS);
        }
        else if( GEOOPT_FF == "UFF")
        {
          ff = ::RDKit::UFF::constructForceField( *rdkit_mol, GEOOPT_NONBONDED_THRESH, -1, GEOOPT_IGNORE_INTERFRAG_INTERACTIONS);
        }
        else
        {
          BCL_MessageStd( "[WARNING] Invalid force field specification - not performing geometry optimization!");
          return *( chemistry::RdkitMolUtils::RDKitRWMolToFragmentComplete( *rdkit_mol, SMILES));
        }
        BCL_MessageVrb( "Running short geometry optimization on the chemical fragment converted from SMILES/SMARTS...");
        ff->initialize();
        ff->minimize( GEOOPT_ITERATIONS);
      }
      return *( chemistry::RdkitMolUtils::RDKitRWMolToFragmentComplete( *rdkit_mol, SMILES));
    }

    //! @brief convert SMARTS string to molecule
    //! @param MOL molecule to be converted to SMARTS string
    chemistry::FragmentComplete RdkitSmilesParser::ConvertSMARTSToMOL
    (
      const std::string &SMILES,
      const bool ADD_H,
      const bool GEN3D,
      const size_t GEOOPT_ITERATIONS,
      const std::string &GEOOPT_FF,
      const float GEOOPT_NONBONDED_THRESH,
      const bool GEOOPT_IGNORE_INTERFRAG_INTERACTIONS
    )
    {
      return ConvertSMILESToMOL
          (
            SMILES,
            true,
            ADD_H,
            GEN3D,
            GEOOPT_ITERATIONS,
            GEOOPT_FF,
            GEOOPT_NONBONDED_THRESH,
            GEOOPT_IGNORE_INTERFRAG_INTERACTIONS
          );
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief test whether a line contains only spaces or is otherwise empty
    //! @param STRING the string to test
    bool RdkitSmilesParser::ContainsNonspaceCharacters( const std::string &STRING)
    {
      for( std::string::const_iterator itr( STRING.begin()), itr_end( STRING.end()); itr != itr_end; ++itr)
      {
        if( !isspace( *itr))
        {
          return true;
        }
      }
      return false;
    }

    //! @brief check validity of first column header name; this column MUST
    //! indicate the output file type (SMILES or SMARTS)
    //! @returns true if valid; false otherwise
    bool RdkitSmilesParser::ContainsValidHeader() const
    {
      const std::string &mol_str( m_Header( 0));

      // first column of header labels must indicate the molecule string type
      if
      (
          mol_str == "SMILES" ||
          mol_str == "smiles" ||
          mol_str == "Smiles" ||
          mol_str == "smi" ||
          mol_str == "SMARTS" ||
          mol_str == "smarts" ||
          mol_str == "Smarts" ||
          mol_str == "sma"
      )
      {
        return true;
      }
      // alternatively, you can get by with no actual header so long as there is only one column;
      // this assumes that
      else if( m_Header.GetSize() == 1)
      {
        return true;
      }
      BCL_MessageCrt
      (
        "Invalid SMILES/SMARTS header name! "
        "The first column of the the first row (header row) of the SMILES/SMARTS file must contain one of the following: "
        "'SMILES', 'smiles', 'Smiles', 'smi', 'SMARTS', 'smarts', 'Smarts', or 'sma'. "
        "Here, the following value was found: '" + util::Format()( mol_str) + "'. Consider checking your delimiters "
        "to verify that if your file contains multiple columns that they are comma- or tab-separated."
      );
      return false;
    }

    //! @brief read RdkitSmilesParser object from std::istream
    //! @param ISTREAM istream that contains RdkitSmilesParser object
    //! @return istream after RdkitSmilesParser object was extracted
    std::istream &RdkitSmilesParser::Read( std::istream &ISTREAM)
    {
      return ReadFile( ISTREAM);
    }

    //! @brief write RdkitSmilesParser into std::ostream
    //! @param OSTREAM ostream that gets RdkitSmilesParser object
    //! @param INDENT number of indentations
    //! @return ostream after RdkitSmilesParser object was inserted
    std::ostream &RdkitSmilesParser::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return
      return OSTREAM;
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> RdkitSmilesParser::s_Instance
    (
      GetObjectInstances().AddInstance( new RdkitSmilesParser())
    );

  } // namespace smiles
} // namespace bcl

