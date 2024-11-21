#ifndef MOLEREADER_H
#define MOLEREADER_H

#include "headers/ICMole/global.h"


namespace ICMole
{



  class MoleReader
  {
       /*!< For debbugging purpose only*/
      static unsigned int  debug;


  private:
    /*!< Path of the file to read */
    std::string    fileName;
    /*!< Stream for parsing file  */
    std::ifstream  finput;
    /*!< File format (PDB/MOL2/SDF) */
    unsigned int   format;
    /*!< Set to true if a file is open , false otherwise */
    bool           isOpen;
    /*!< Set to false at the beginning of a new molecule reading.
     * Set to true if the read went well */
    bool           moleCheck;
    /*!< Boolean telling if we are at the end of the file*/
    bool           EOFile;
    /*!< Boolean telling if the file is already loaded as string vector or not*/
    bool           lineinfile;
    /*!< For PDB only, set the processing of the PDB File */
    bool           wAnalyse;
    /*!< For all format, split the molecular file into protein/ligand/cofactor ...*/
    bool          wSplit;

    /*!< List of lines that represent the file */
    std::vector<std::string>            lineFile;
    /*!< Iterator over lines */
    std::vector<std::string>::iterator  Itlines;
    /*!< Makes the relationship between residu names and their corresponding
     *   object*/
    std::map<std::string,Residu*> mapRes;
    /*!< Makes the link between Calpha and residues */
    std::map<unsigned int,Residu*> mapResRoot;
    bool open();
    void close();
    void getLine(std::string &ligne) ;
    void readMOL2Substruct(const unsigned int NbrSubst, Molecule& mole) throw(MoleExcept);
    void readMOL2Atom(const unsigned int NbrAt,Molecule &mole) throw(MoleExcept);
    void readMOL2Bond(const unsigned int NbrBd, Molecule& mole) throw(MoleExcept);
    void tokenize( const std::string& str,std:: map<std::string,std::string>& tokens);

public:
    MoleReader();

    MoleReader(const std::vector<std::string> &linefile, const unsigned int& fformat=FileFormat::UNDEFINED) throw(MoleExcept);

    MoleReader( const std::string& filename, const unsigned int& fformat=FileFormat::UNDEFINED) throw(MoleExcept);
    
    const bool& isEOF() const {return EOFile;}

    void loadNewFile( const std::string& filename, const unsigned int& fformat=FileFormat::UNDEFINED, const bool &with_split=false) throw(MoleExcept);

    void loadNextMolecule(Molecule &molecule, const unsigned int& typeMole=MoleType::UNDEFINED);

    void loadInComplex(Complex &cp, const unsigned int& TypeMol);

    void loadPDB(Molecule& molecule,const unsigned int& TypeMol) throw(MoleExcept);

    size_t getNumMolecules();

    void loadNextChm(Molecule &molecule, const unsigned int& typeMole=MoleType::UNDEFINED);

    std::vector<std::string> split(std::string str, char delimiter);

  };
}

#endif // MOLEREADER_H
