#ifndef MOLEWRITER_H
#define MOLEWRITER_H
#include "headers/ICMole/global.h"

namespace ICMole
{

class MoleWriter
{

private:
    /**
     * @brief Name of the output file
     */
    std::string   fname;

    /**
     * @brief  Object used to save molecule
     */
    std::ofstream ofstr;

    /**
     * @brief File Format
     */
    short FFormat;
     /*!< \brief Open file */
    /**
     * @brief Open the given file set as fname
     * @return TRUE when the file has been opened, false otherwise
     */
    bool openFile();
    /**
     * @brief Output Hydrogen in the file. Default TRUE
     */
    bool withH;
    /**
     * @brief Use rotpos coordinates instead of fixpos. Default: FALSE
     */
    bool withRotPos;

    /**
     * @brief Include dummy atom describing aromatic center in the file.
     * Default: FALSE
     */
    bool withCycles;

    /**
     * @brief Include atom sets in the file. Default FALSE
     */
    bool withSets;

    /**
     * @brief Output only used atoms. Default FALSE
     */
    bool onlyUsed;


    void putDataInFile();

public:
    MoleWriter();
    MoleWriter(const std::string& name, const short& FFormat=FileFormat::MOL2);
    ~MoleWriter();
    void writeMOL2( Molecule *const  mole, const bool outints = false);
    void printMOL2( Molecule *const  mole);
    void writePharmaLig( Molecule *const mole);
    void writeMOL2ss(Molecule *const  mole, int const mod);
    void writeMOL2(const Complex &complex, const std::string &name);

    inline void setWithH     (const bool& b){withH=b;}
    inline void setWithRotPos(const bool& b){withRotPos=b;}
    inline void setWithCycles(const bool& b){withCycles=b;}
    inline void setWithSets  (const bool& b){withSets=b;}
    inline void setOnlyUsed  (const bool& b){onlyUsed=b;}
    bool newFile(const std::string& name, const unsigned int& FFormat=FileFormat::MOL2);

};

}
#endif // MOLEWRITER_H
