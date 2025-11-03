#ifndef PDBCONVERT_H
#define PDBCONVERT_H
#include "headers/ICMole/global.h"

namespace ICMole
{

class PDBConvert
{
     typedef std::vector<unsigned int> listLine;

    std::vector<std::string> filelines;

    unsigned int header_line;
     int master;

    /**
     * @brief All file pointer of REMARK Lines are stored
     */
    listLine remark;

    /**
     * @brief Keep file positions of every ATOM, HETATM, TER blocks
     */
    listLine atoms;


    listLine connect;
    listLine seqres;

    listLine coords;
    /*
    listLine dbref;
    listLine compnd;
    listLine sheet;
    listLine het;
    listLine helix;
    listLine site;
    listLine ssbond;
    listLine links;
    listLine cispep;
    listLine source;
    listLine title;
    listLine status;*/


    std::string PDBName;
    std::fstream filestream;


public:
    static bool for_scPDB;

    void loadPDB(const std::string& Path) ;
    void harmonizeSizeLine(const size_t& length=80);
    void changeMSEtoMET();
    void changeCSEtoCYS();
    void moveHETATMtoend();
    void selAltAtm();
    void renumAtms();
    void toComplex(Complex&cp) const;
    void setPDBName(const std::string &name) {PDBName=name;}

    PDBConvert(const std::string& name="");
};

}
#endif // PDBCONVERT_H
