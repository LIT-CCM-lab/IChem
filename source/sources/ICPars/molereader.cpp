#include <iostream>
#include "headers/ICPars/molereader.h"
#include "headers/ICMole/complex.h"
#include "headers/ICMole/molecule.h"
using namespace ICMole;
using namespace std;

unsigned int  MoleReader::debug=0;




/*! \fn MoleReader::MoleReader()
  * \brief Standard constructor
  *
  */
MoleReader::MoleReader()
{
    fileName="";
    EOFile=false;
    isOpen=false;
    lineinfile=false;
    moleCheck=false;
    wAnalyse=false;
    wSplit=false;
    format=MoleType::UNDEFINED;
}




/**
 * @brief MoleReader::MoleReader
 * @param filename : Path of the molecular file
 * @param fformat : file format (see FileFormat namespace)
 * @throw 2010101 : No file given
 * @throw 2010102 : Unrecognized file format
 * @throw 2010103 : Unable to open file
 */
MoleReader::MoleReader(const std::string&  filename,
                       const unsigned int& fformat)
throw(MoleExcept)
{
    if (filename.length()== 0)
        throw MoleExcept(2010101,
                         "MoleReader::CONSTRUCTOR",
                         "No file given");
    if (fformat>=FileFormat::UNDEFINED){
        const string extension = filename.substr(filename.length()-4,4);
        if (extension == "mol2" ||
                extension == "MOL2") format = FileFormat::MOL2;
        else if (extension == ".pdb" ||
                 extension == ".PDB") format = FileFormat::PDB;
        else if (extension == ".sdf" ||
                 extension == ".SDF") format = FileFormat::SDF;
        else throw MoleExcept(2010102,
                              "MoleReader::CONSTRUCTOR",
                              "Unrecognized file format : "+extension);
    }
    else format = fformat;
    EOFile=false;
    isOpen=false;
    lineinfile=false;
    wAnalyse=false;
    wSplit=false;
    fileName=filename;
    if (!open()) throw MoleExcept(2010103,
                                  "MoleReader::CONSTRUCTOR",
                                  "Unable to open "+filename);
}


/**
 * @brief MoleReader::MoleReader
 * @param linefile : list of line of parse
 * @param fformat : Molecular type
 * @throw 2010201 : No line given
 */
MoleReader::MoleReader(const vector<string> &linefile,
                       const unsigned int& fformat)
throw(MoleExcept)
{
    if (linefile.size() ==0)
        throw MoleExcept(2010201,
                         "MoleReader::CONSTRUCTOR",
                         "No line given");
    format = fformat;
    wSplit=false;
    EOFile = false;
    isOpen = true;
    lineinfile=true;
    wAnalyse=false;
    lineFile = linefile;
    Itlines = lineFile.begin();
}



/*! \fn bool MoleReader::open()
  * \return True if the opening went well. False otherwise
  *
  * Close an already opened file and open the new one.
  *
  */
bool MoleReader::open()
{
    if (finput) close();
    if (debug) { cout << "DEBUG -- Opening file " << fileName<<endl;}
    finput.open(fileName.c_str(),ios::binary);

    if (!finput.is_open()){ return false;  }
    EOFile = false;
    isOpen=true;
    return true;
}




/*! \fn void MoleReader::Close()
  *
  * Close an already opened file.
  */
void MoleReader::close()
{
    if (isOpen && !lineinfile) finput.close();
}




/**
 * @brief MoleReader::loadNewFile
 * @param filename
 * @param fformat
 * @param with_split
 * @throw 2010301 : No file given
 * @throw 2010302 : Unrecognized file format
 * @throw 2010303 : Unable to open file
 */
void MoleReader::loadNewFile( const std::string& filename,
                              const unsigned int& fformat,
                              const bool &with_split)
throw(MoleExcept)
{
    if (filename.length()== 0)
        throw MoleExcept(2010301,
                         "MoleReader::loadNewFile",
                         "No file given");

    if (fformat>=FileFormat::UNDEFINED){
        const string extension = filename.substr(filename.length()-4,4);
        if (extension == "mol2" ||
                extension == "MOL2") format = FileFormat::MOL2;
        else if (extension == ".pdb" ||
                 extension == ".PDB") format = FileFormat::PDB;
        else if (extension == ".sdf" ||
                 extension == ".SDF") format = FileFormat::SDF;
        else if (extension == ".chm" ||
                 extension == ".CHM") format = FileFormat::CHM;
        else throw MoleExcept(2010302,
                              "MoleReader::loadNewFile",
                              "Unrecognized file format "+extension
                              +" in "+filename);}
    else format = fformat;
    if (isOpen) 
        finput.close();
    EOFile = false;
    isOpen=false;
    fileName = filename;
    wSplit=with_split;
    if (!open()) throw MoleExcept(2010303,
                                  "MoleReader::loadNewFile",
                                  "Unable to open "+filename);

}


void MoleReader::getLine(std::string &ligne)  {
    ligne= "";
    if (lineinfile)
    {
        if (Itlines == lineFile.end()
                || (Itlines+1) == lineFile.end())
        { EOFile=true;return;  }
        else {Itlines++; ligne= *Itlines;}
        return;
    }
    else
    {
        if (std::getline( finput, ligne ))
        {
            const size_t len=ligne.length();
            for (size_t i=0;i<len;i++)
                if (ligne[i]==char(13)) ligne.replace(i,1,"");
            return;
        }

        if (finput.eof())  EOFile=true;
        ligne="";
        return;
    }
}

size_t MoleReader::getNumMolecules()
{
    size_t nMole=0;
    std::string ligne;
    const size_t former_fpos= finput.tellg();
    finput.clear();
    finput.seekg(0,ios_base::beg);
    while (!finput.eof())
    {
        getLine(ligne);
        if (ligne.find("@<TRIPOS>MOLECULE") !=string::npos)
            nMole++;
    }
    finput.clear();
    finput.seekg(former_fpos,ios_base::beg);
    EOFile=false;
    return nMole;
}


void MoleReader::loadInComplex(Complex &cp, const unsigned int& TypeMol)
{
    while(!EOFile)
    {
        Molecule *mole = new Molecule(100,100,TypeMol);
        if (format == FileFormat::MOL2)
        {

            loadNextMolecule(*mole,TypeMol);
            mole->checkMOL2();
            mole->ringPerception();

        }
        else if  (format == FileFormat::CHM){
            loadNextChm(*mole,TypeMol);
        }else
        {
            loadPDB(*mole,TypeMol);
        }

        if (mole != (Molecule*)NULL)
        {
            if (mole->numAtom() > 0) cp.addMolecule(mole);
            else delete mole;


        }
    }
}
