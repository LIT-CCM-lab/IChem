#ifndef PHARMWRITER_H
#define PHARMWRITER_H


#include "headers/ICMole/global.h"
#include "headers/ICMole/pharmprop.h"


namespace ICMole
{

class PharmWriter
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



    void putDataInFile();

public:
    PharmWriter();
    PharmWriter(const std::string& name, const short& FFormat=FileFormat::CHM);
    ~PharmWriter();
    void writeCHMss(const PharmList &pharmacophore, std::string toto, bool boolWeigth=false);
    void writeCHM(const PharmList &pharmacophore, const char *namefile, bool boolWeigth=false);
    void writePML(const PharmList &pharmacophore, const char *namefile, bool boolWeigth=false);

    bool newFile(const std::string& name, const unsigned int& FFormat=FileFormat::CHM);

    void writeHBD(int nb, PharmProp prop, const char *namefile);
    void writeHBA(int nb, PharmProp prop, const char *namefile);
    void writePos(int nb, PharmProp prop, const char *namefile);
    void writeNeg(int nb, PharmProp prop, const char *namefile);
    void writeAro(int nb, PharmProp prop, const char *namefile);
    void writehyd(int nb, PharmProp prop, const char *namefile);
    void writemet(int nb, PharmProp prop, const char *namefile);

    };

}

#endif // PHARMWRITER_H
