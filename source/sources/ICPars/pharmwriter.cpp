#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <ctime>

#include "headers/ICPars/pharmwriter.h"
#include "headers/ICMole/molecule.h"
#include "headers/ICMole/pharmprop.h"
#include <sstream>


using namespace std;
using namespace ICMole;

PharmWriter::PharmWriter():
    fname(""),
    FFormat(FileFormat::CHM)
{}


PharmWriter::PharmWriter(const string &name,
                         const short &FFormat):
    fname(name),
    FFormat(FFormat)
{

    if (fname != "")  openFile();
}

/**
 * @brief PharmWriter::newFile
 * @param newfname: path of the file
 * @param FForma : File format
 * @return True if the file has been open. False otherwise
 *
 * Close a potentially already opened file, and open the new one defined by fname
 */
bool PharmWriter::newFile(const string& newfname, const unsigned int& FForma)
{

    if (ofstr.is_open())ofstr.close();
    fname=newfname;FFormat=FForma;
    return openFile();
}


string IntToString (int a)
{
    ostringstream temp;
    temp<<a;
    return temp.str();
}
string DoubleToString (double a)
{
    ostringstream temp;
    temp<<a;
    return temp.str();
}

void comptable (const PharmList &pharmacophore, int nbprop[7])
{
    for (ItcPharm ItPha = pharmacophore.begin(); ItPha != pharmacophore.end();++ItPha ){
        PharmProp &Pharm=**ItPha;
        if (Pharm.prop()         == PharType::HYDROPHOBIC){
            nbprop[1]++;nbprop[0]++;
        }else if (Pharm.prop()   == PharType::HBACCEPTOR){
            nbprop[2]++;nbprop[0]++;
        }else if (Pharm.prop()   == PharType::HBDONOR){
            nbprop[3]++;nbprop[0]++;
        }else if (Pharm.prop()   == PharType::AROMATIC){
            nbprop[4]++;nbprop[0]++;
        } else   if (Pharm.prop()== PharType::IONNEG){
            nbprop[5]++;nbprop[0]++;
        } else   if (Pharm.prop()== PharType::IONPOS){
            nbprop[6]++;nbprop[0]++;
        } else   if (Pharm.prop()== PharType::METAL){
            nbprop[7]++;nbprop[0]++;
        } else if (Pharm.prop()  == PharType::EXCLUSION){
            nbprop[8]++;
        }
    }
}
bool PharmWriter::openFile()
{

    if (fname.length() ==0)return false;
    ofstr.open(fname.c_str(),ios::out|ios::binary);
    if (!ofstr.is_open())
    {    cout <<"Error Writing to "<< fname <<endl;return false;}
    return true;
}

PharmWriter::~PharmWriter() { if (ofstr.is_open()) ofstr.close();}

void PharmWriter::writeAro(int nb, PharmProp prop, const char * namefile){
    ofstream ofs;
    ofs.open(namefile,ios::out|ios::app);

    //ARO
    ofs <<"  AND Ring6A"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 37.5 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM R6A11."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM R6A21."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM R6A31."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM R6A41."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM R6A51."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM R6A61."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND R6A561."<<nb<<" ( R6A51."<<nb<<" R6A61."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND R6A451."<<nb<<" ( R6A41."<<nb<<" R6A51."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND R6A341."<<nb<<" ( R6A31."<<nb<<" R6A41."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND R6A231."<<nb<<" ( R6A21."<<nb<<" R6A31."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND R6A121."<<nb<<" ( R6A11."<<nb<<" R6A21."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND R6A611."<<nb<<" ( R6A61."<<nb<<" R6A11."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"    OBJECTS :"<< endl <<"      POINT normal-R6A1."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 300.0 0.0 0.0 ) ;"<< endl <<"      POINT centroid-R6A1."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      VECTOR vector-R6A1."<<nb<<" HEAD ( normal-R6A1."<<nb<<" ) TAIL ( centroid-R6A1."<<nb<<" ) ;"<< endl <<"      PLANE plane-R6A1."<<nb<<" ORIGIN ( centroid-R6A1."<<nb<<" ) NORMAL ( normal-R6A1."<<nb<<" )"<< endl <<"    LEAST_SQ_PLANE ( R6A11."<<nb<<" R6A21."<<nb<<" R6A31."<<nb<<" R6A41."<<nb<<" R6A51."<<nb<<" R6A61."<<nb<<" ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND Ring5A"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 42.9 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM R5A11."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM R5A21."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM R5A31."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM R5A41."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM R5A51."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND R5A451."<<nb<<" ( R5A41."<<nb<<" R5A51."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND R5A341."<<nb<<" ( R5A31."<<nb<<" R5A41."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND R5A231."<<nb<<" ( R5A21."<<nb<<" R5A31."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND R5A121."<<nb<<" ( R5A11."<<nb<<" R5A21."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND R5A511."<<nb<<" ( R5A51."<<nb<<" R5A11."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"    OBJECTS :"<< endl <<"      POINT normal-R5A1."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 300.0 0.0 0.0 ) ;"<< endl <<"      POINT centroid-R5A1."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      VECTOR vector-R5A1."<<nb<<" HEAD ( normal-R5A1."<<nb<<" ) TAIL ( centroid-R5A1."<<nb<<" ) ;"<< endl <<"      PLANE plane-R5A1."<<nb<<" ORIGIN ( centroid-R5A1."<<nb<<" ) NORMAL ( normal-R5A1."<<nb<<" )"<< endl <<"    LEAST_SQ_PLANE ( R5A11."<<nb<<" R5A21."<<nb<<" R5A31."<<nb<<" R5A41."<<nb<<" R5A51."<<nb<<" ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND Ring6B"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 37.5 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM R6B11."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM R6B21."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM R6B31."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM R6B41."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM R6B51."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM R6B61."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND R6B561."<<nb<<" ( R6B51."<<nb<<" R6B61."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND R6B451."<<nb<<" ( R6B41."<<nb<<" R6B51."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND R6B341."<<nb<<" ( R6B31."<<nb<<" R6B41."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND R6B231."<<nb<<" ( R6B21."<<nb<<" R6B31."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND R6B121."<<nb<<" ( R6B11."<<nb<<" R6B21."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND R6B611."<<nb<<" ( R6B61."<<nb<<" R6B11."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"    OBJECTS :"<< endl <<"      POINT normal-R6B1."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 300.0 0.0 0.0 ) ;"<< endl <<"      POINT centroid-R6B1."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      VECTOR vector-R6B1."<<nb<<" HEAD ( normal-R6B1."<<nb<<" ) TAIL ( centroid-R6B1."<<nb<<" ) ;"<< endl <<"      PLANE plane-R6B1."<<nb<<" ORIGIN ( centroid-R6B1."<<nb<<" ) NORMAL ( normal-R6B1."<<nb<<" )"<< endl <<"    LEAST_SQ_PLANE ( R6B61."<<nb<<" R6B51."<<nb<<" R6B41."<<nb<<" R6B31."<<nb<<" R6B21."<<nb<<" R6B11."<<nb<<" ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND Ring5B"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 42.9 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM R5B11."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM R5B21."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM R5B31."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM R5B41."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM R5B51."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND R5B451."<<nb<<" ( R5B41."<<nb<<" R5B51."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND R5B341."<<nb<<" ( R5B31."<<nb<<" R5B41."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND R5B231."<<nb<<" ( R5B21."<<nb<<" R5B31."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND R5B121."<<nb<<" ( R5B11."<<nb<<" R5B21."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND R5B511."<<nb<<" ( R5B51."<<nb<<" R5B11."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"    OBJECTS :"<< endl <<"      POINT normal-R5B1."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 300.0 0.0 0.0 ) ;"<< endl <<"      POINT centroid-R5B1."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      VECTOR vector-R5B1."<<nb<<" HEAD ( normal-R5B1."<<nb<<" ) TAIL ( centroid-R5B1."<<nb<<" ) ;"<< endl <<"      PLANE plane-R5B1."<<nb<<" ORIGIN ( centroid-R5B1."<<nb<<" ) NORMAL ( normal-R5B1."<<nb<<" )"<< endl <<"    LEAST_SQ_PLANE ( R5B51."<<nb<<" R5B41."<<nb<<" R5B31."<<nb<<" R5B21."<<nb<<" R5B11."<<nb<<" ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl ;

    ofs << setprecision(1)<< fixed<<"  OR RING_AROMATIC"<<nb<<" ( Ring6A"<<nb<<" Ring5A"<<nb<<" Ring6B"<<nb<<" Ring5B"<<nb<<" )"<< endl
        <<"  {"<< endl
       <<"    ATTRIBUTES :"<< endl
      <<"      COORD2D ( 0 0 ) ;"<< endl
     <<"      COORD3D ( "<<prop.get_1A().x*100<<" "<<prop.get_1A().y*100<<" "<<prop.get_1A().z*100<<" ) ;"<< endl
    <<"      PROP ( SOURCE \"BIOCAD_FUNCTION\" ) ;"<< endl
    <<"    OBJECTS :"<< endl
    <<"      POINT centroid"<<nb<< endl
    <<"    COORD2D ( 0 0 )"<< endl
    <<"    COORD3D ( "<<prop.getBase().x*100<<" "<<prop.getBase().y*100<<" "<<prop.getBase().z*100<<" )"<< endl
    <<"    MAP ( Ring6A"<<nb<<"/centroid-R6A1."<<nb<<" Ring5A"<<nb<<"/centroid-R5A1."<<nb<<" Ring6B"<<nb<<"/centroid-R6B1."<<nb<<" Ring5B"<<nb<<"/centroid-R5B1."<<nb<<" )"<< endl
    <<"    PROP ( INDEX_CODE \"20\" ) ;"<< endl
    <<"      POINT normal"<<nb<< endl
    <<"    COORD2D ( 0 0 )"<< endl
    <<"    COORD3D ( "<<prop.get_3A().x*100<<" "<<prop.get_3A().y*100<<" "<<prop.get_3A().z*100<<" )"<< endl
    <<"    MAP ( Ring6A"<<nb<<"/normal-R6A1."<<nb<<" Ring5A"<<nb<<"/normal-R5A1."<<nb<<" Ring6B"<<nb<<"/normal-R6B1."<<nb<<" Ring5B"<<nb<<"/normal-R5B1."<<nb<<" )"<< endl
    <<"    PROP ( INDEX_CODE \"21\" ) ;"<< endl
    <<"      VECTOR vector"<<nb<<" HEAD ( normal"<<nb<<" ) TAIL ( centroid"<<nb<<" )"<< endl
    <<"    MAP ( Ring6A"<<nb<<"/vector-R6A1."<<nb<<" Ring5A"<<nb<<"/vector-R5A1."<<nb<<" Ring6B"<<nb<<"/vector-R6B1."<<nb<<" Ring5B"<<nb<<"/vector-R5B1."<<nb<<" ) ;"<< endl
    <<"      PLANE plane"<<nb<<" ORIGIN ( centroid"<<nb<<" ) NORMAL ( normal"<<nb<<" )"<< endl
    <<"    MAP ( Ring6A"<<nb<<"/plane-R6A1."<<nb<<" Ring5A"<<nb<<"/plane-R5A1."<<nb<<" Ring6B"<<nb<<"/plane-R6B1."<<nb<<" Ring5B"<<nb<<"/plane-R5B1."<<nb<<" ) ;"<< endl
    <<"  }"<< endl;
    ofs.close();
}
#include <ctime>
void PharmWriter::writeHBD(int nb,PharmProp prop,const char * namefile){
    ofstream ofs;
    ofs.open(namefile,ios::out|ios::app);
    // HBB
    ofs <<"  AND Hydroxyl"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 100.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM Heavy#O#51."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #H#61."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 4-51."<<nb<<" ( Heavy#O#51."<<nb<<" #H#61."<<nb<<" ) ;"<< endl <<"    OBJECTS :"<< endl <<"      POINT DonorProj11."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 300.0 0.0 0.0 ) ;"<< endl <<"      VECTOR DonorVector11."<<nb<<" HEAD ( DonorProj11."<<nb<<" ) TAIL ( Heavy#O#51."<<nb<<" ) HBDONOR ;"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND Exclude_IonizedHydroxyl"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #CPS...#81."<<nb<<""<< endl <<"    ELEMENTS ( C P S )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM O#101."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 2 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #O#111."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #H#121."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 7-91."<<nb<<" ( #CPS...#81."<<nb<<" O#101."<<nb<<" ) DOUBLE ;"<< endl <<"      BOND 7-101."<<nb<<" ( #CPS...#81."<<nb<<" #O#111."<<nb<<" ) ;"<< endl <<"      BOND 10-111."<<nb<<" ( #O#111."<<nb<<" #H#121."<<nb<<" ) ;"<< endl <<"    OBJECTS :"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  EXCLUDE Hydroxyl_withExclusions"<<nb<<" ( Hydroxyl"<<nb<<" Exclude_IonizedHydroxyl"<<nb<<" )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    PRESENT ( Hydroxyl"<<nb<<" ) ;"<< endl <<"    ABSENT ( Exclude_IonizedHydroxyl"<<nb<<" ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"      MAP Hydroxyl"<<nb<<"/Heavy#O#51."<<nb<<""<< endl <<"    ( Exclude_IonizedHydroxyl"<<nb<<"/#O#111."<<nb<<" ) ;"<< endl <<"      MAP Hydroxyl"<<nb<<"/#H#61."<<nb<<""<< endl <<"    ( Exclude_IonizedHydroxyl"<<nb<<"/#H#121."<<nb<<" ) ;"<< endl <<"      MAP Hydroxyl"<<nb<<"/4-51."<<nb<<""<< endl <<"    ( Exclude_IonizedHydroxyl"<<nb<<"/10-111."<<nb<<" ) ;"<< endl <<"  }"<< endl <<"  AND ThiolAcetylene"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 100.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM Heavy#CS#141."<<nb<<""<< endl <<"    ELEMENTS ( C S )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 1 )"<< endl <<"    COORDINATION ( 2 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #H#151."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 13-141."<<nb<<" ( Heavy#CS#141."<<nb<<" #H#151."<<nb<<" ) ;"<< endl <<"    OBJECTS :"<< endl <<"      POINT DonorProj21."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 300.0 0.0 0.0 ) ;"<< endl <<"      VECTOR DonorVector21."<<nb<<" HEAD ( DonorProj21."<<nb<<" ) TAIL ( Heavy#CS#141."<<nb<<" ) HBDONOR ;"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND amine"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 100.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #H#191."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM Heavy#N#201."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    HCOUNT ( 1 2 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 19-181."<<nb<<" ( Heavy#N#201."<<nb<<" #H#191."<<nb<<" ) ;"<< endl <<"    OBJECTS :"<< endl <<"      POINT DonorProj31."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 300.0 0.0 0.0 ) ;"<< endl <<"      VECTOR DonorVector31."<<nb<<" HEAD ( DonorProj31."<<nb<<" ) TAIL ( Heavy#N#201."<<nb<<" ) HBDONOR ;"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND Exclude_tetrazolyl"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #H#281."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #C#231."<<nb<<""<< endl <<"    ELEMENTS ( C )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#241."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#251."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#261."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#271."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    HCOUNT ( 1 2 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 26-271."<<nb<<" ( #N#271."<<nb<<" #H#281."<<nb<<" ) ;"<< endl <<"      BOND 22-231."<<nb<<" ( #C#231."<<nb<<" #N#241."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 23-241."<<nb<<" ( #N#241."<<nb<<" #N#251."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 24-251."<<nb<<" ( #N#251."<<nb<<" #N#261."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 25-261."<<nb<<" ( #N#261."<<nb<<" #N#271."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 22-261."<<nb<<" ( #C#231."<<nb<<" #N#271."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"    OBJECTS :"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND Exclude_cf3Sulfonamide"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #F#301."<<nb<<""<< endl <<"    ELEMENTS ( F )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #F#31."<<nb<<"1."<<nb<<""<< endl <<"    ELEMENTS ( F )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #C#321."<<nb<<""<< endl <<"    ELEMENTS ( C )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #H#341."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#351."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    HCOUNT ( 1 2 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #O#361."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #O#371."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #F#381."<<nb<<""<< endl <<"    ELEMENTS ( F )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #S#391."<<nb<<""<< endl <<"    ELEMENTS ( S )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 33-341."<<nb<<" ( #H#341."<<nb<<" #N#351."<<nb<<" ) ;"<< endl <<"      BOND 38-341."<<nb<<" ( #S#391."<<nb<<" #N#351."<<nb<<" ) EXOCYCLIC ;"<< endl <<"      BOND 38-351."<<nb<<" ( #S#391."<<nb<<" #O#361."<<nb<<" ) DOUBLE EXOCYCLIC ;"<< endl <<"      BOND 38-361."<<nb<<" ( #S#391."<<nb<<" #O#371."<<nb<<" ) DOUBLE EXOCYCLIC ;"<< endl <<"      BOND 31-301."<<nb<<" ( #C#321."<<nb<<" #F#31."<<nb<<"1."<<nb<<" ) EXOCYCLIC ;"<< endl <<"      BOND 31-381."<<nb<<" ( #C#321."<<nb<<" #S#391."<<nb<<" ) EXOCYCLIC ;"<< endl <<"      BOND 31-371."<<nb<<" ( #C#321."<<nb<<" #F#381."<<nb<<" ) EXOCYCLIC ;"<< endl <<"      BOND 31-291."<<nb<<" ( #C#321."<<nb<<" #F#301."<<nb<<" ) EXOCYCLIC ;"<< endl <<"    OBJECTS :"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  EXCLUDE good_amine_withExclusions"<<nb<<" ( amine"<<nb<<" Exclude_tetrazolyl"<<nb<<" Exclude_cf3Sulfonamide"<<nb<<" )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    PRESENT ( amine"<<nb<<" ) ;"<< endl <<"    ABSENT ( Exclude_tetrazolyl"<<nb<<" Exclude_cf3Sulfonamide"<<nb<<" ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"      MAP amine"<<nb<<"/#H#191."<<nb<<""<< endl <<"    ( Exclude_tetrazolyl"<<nb<<"/#H#281."<<nb<<" Exclude_cf3Sulfonamide"<<nb<<"/#H#341."<<nb<<" ) ;"<< endl <<"      MAP amine"<<nb<<"/Heavy#N#201."<<nb<<""<< endl <<"    ( Exclude_tetrazolyl"<<nb<<"/#N#271."<<nb<<" Exclude_cf3Sulfonamide"<<nb<<"/#N#351."<<nb<<" ) ;"<< endl <<"      MAP amine"<<nb<<"/19-181."<<nb<<""<< endl <<"    ( Exclude_tetrazolyl"<<nb<<"/26-271."<<nb<<" Exclude_cf3Sulfonamide"<<nb<<"/33-341."<<nb<<" ) ;"<< endl <<"  }"<< endl;


    ofs << setprecision(1)<< fixed<<  "  OR HB_DONOR"<<nb<<" ( Hydroxyl_withExclusions"<<nb<<" ThiolAcetylene"<<nb<<" good_amine_withExclusions"<<nb<<" )"<< endl
        <<"  {"<< endl
       <<"    ATTRIBUTES :"<< endl
      <<"      COORD2D ( 0 0 ) ;"<< endl
     <<"      COORD3D ( "<<prop.get_1A().x*100<<" "<<prop.get_1A().y*100<<" "<<prop.get_1A().z*100<<" ) ;"<< endl
    <<"      PROP ( SOURCE \"BIOCAD_FUNCTION\" ) ;"<< endl
    <<"    OBJECTS :"<< endl
    <<"      ATOM HB_DONOR"<<nb<<"_TAIL"<< endl
    <<"    COORD2D ( 0 0 )"<< endl
    <<"    COORD3D ( "<<prop.getBase().x*100<<" "<<prop.getBase().y*100<<" "<<prop.getBase().z*100<<" )"<< endl
    <<"    MAP ( Hydroxyl_withExclusions"<<nb<<"/Hydroxyl"<<nb<<"/Heavy#O#51."<<nb<<" ThiolAcetylene"<<nb<<"/Heavy#CS#141."<<nb<<" good_amine_withExclusions"<<nb<<"/amine"<<nb<<"/Heavy#N#201."<<nb<<" )"<< endl
    <<"    PROP ( INDEX_CODE \"12\" ) ;"<< endl
    <<"      POINT HB_DONOR"<<nb<<"_HEAD"<< endl
    <<"    COORD2D ( 0 0 )"<< endl
    <<"    COORD3D ( "<<prop.get_3A().x*100<<" "<<prop.get_3A().y*100<<" "<<prop.get_3A().z*100<<" )"<< endl
    <<"    MAP ( Hydroxyl_withExclusions"<<nb<<"/Hydroxyl"<<nb<<"/DonorProj11."<<nb<<" ThiolAcetylene"<<nb<<"/DonorProj21."<<nb<<" good_amine_withExclusions"<<nb<<"/amine"<<nb<<"/DonorProj31."<<nb<<" )"<< endl
    <<"    PROP ( INDEX_CODE \"13\" ) ;"<< endl
    <<"      VECTOR DonorVector"<<nb<<" HEAD ( HB_DONOR"<<nb<<"_HEAD ) TAIL ( HB_DONOR"<<nb<<"_TAIL )"<< endl
    <<"    MAP ( Hydroxyl_withExclusions"<<nb<<"/Hydroxyl"<<nb<<"/DonorVector11."<<nb<<" ThiolAcetylene"<<nb<<"/DonorVector21."<<nb<<" good_amine_withExclusions"<<nb<<"/amine"<<nb<<"/DonorVector31."<<nb<<" ) ;"<< endl
    <<"  }"<<endl;
    ofs.close();
}

void PharmWriter::writeHBA(int nb,PharmProp prop,const char * namefile){
    ofstream ofs;
    ofs.open(namefile,ios::out|ios::app);
    ofs <<"  AND atomWithLP"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 150.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM Heavy#NOS#31."<<nb<<""<< endl <<"    ELEMENTS ( N O S )"<< endl <<"    CHARGE ( -7 -6 -5 -4 -3 -2 -1 0 )"<< endl <<"    VALENCE ( 2 3 )"<< endl <<"    LPCOUNT ( 1 2 3 )"<< endl <<"    COORDINATION ( 1 2 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    OBJECTS :"<< endl <<"      POINT AcceptorProj11."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 300.0 0.0 0.0 ) ;"<< endl <<"      VECTOR AcceptorVector11."<<nb<<" HEAD ( AcceptorProj11."<<nb<<" ) TAIL ( Heavy#NOS#31."<<nb<<" ) HBACCEPTOR ;"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND goodN"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 100.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM Heavy#N#31."<<nb<< endl <<"    ELEMENTS ( N )"<< endl <<"    LPCOUNT ( 1 )"<< endl <<"    COORDINATION ( 3 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #NonC#31."<<nb<<""<< endl <<"    ELEMENTS ( ANY_NONCARBON )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND N-NonC1."<<nb<<" ( Heavy#N#31."<<nb<<" #NonC#31."<<nb<<" ) SINGLE DOUBLE AROMATIC TRIPLE ;"<< endl <<"    OBJECTS :"<< endl <<"      POINT AcceptorProj21."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 300.0 0.0 0.0 ) ;"<< endl <<"      VECTOR AcceptorVector21."<<nb<<" HEAD ( AcceptorProj21."<<nb<<" ) TAIL ( Heavy#N#31."<<nb<<" ) HBACCEPTOR ;"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl ;

    ofs<< setprecision(1)<< fixed <<"  OR HB_ACCEPTOR"<<nb<<" ( atomWithLP"<<nb<<" goodN"<<nb<<" )"<< endl
       <<"  {"<< endl
      <<"    ATTRIBUTES :"<< endl
     <<"      COORD2D ( 0 0 ) ;"<< endl
    <<"      COORD3D ( "<<prop.get_1A().x*100<<" "<<prop.get_1A().y*100<<" "<<prop.get_1A().z*100<<" ) ;"<< endl
    <<"      PROP ( SOURCE \"BIOCAD_FUNCTION\" ) ;"<< endl
    <<"    OBJECTS :"<< endl
    <<"      ATOM HB_ACCEPTOR"<<nb<<"_TAIL"<< endl
    <<"    COORD2D ( 0 0 )"<< endl
    <<"    COORD3D ( "<<prop.getBase().x*100<<" "<<prop.getBase().y*100<<" "<<prop.getBase().z*100<<" )"<< endl
    <<"    MAP ( atomWithLP"<<nb<<"/Heavy#NOS#31."<<nb<<" goodN"<<nb<<"/Heavy#N#31."<<nb<<" )"<< endl
    <<"    PROP ( INDEX_CODE \"14\" ) ;"<< endl
    <<"      POINT HB_ACCEPTOR"<<nb<<"_HEAD"<< endl
    <<"    COORD2D ( 0 0 )"<< endl
    <<"    COORD3D ( "<<prop.get_3A().x*100<<" "<<prop.get_3A().y*100<<" "<<prop.get_3A().z*100<<" )"<< endl
    <<"    MAP ( atomWithLP"<<nb<<"/AcceptorProj11."<<nb<<" goodN"<<nb<<"/AcceptorProj21."<<nb<<" )"<< endl
    <<"    PROP ( INDEX_CODE \"15\" ) ;"<< endl
    <<"      VECTOR AcceptorVector"<<nb<<" HEAD ( HB_ACCEPTOR"<<nb<<"_HEAD ) TAIL ( HB_ACCEPTOR"<<nb<<"_TAIL )"<< endl
    <<"    MAP ( atomWithLP"<<nb<<"/AcceptorVector11."<<nb<<" goodN"<<nb<<"/AcceptorVector21."<<nb<<" ) ;"<< endl
    <<"  }"<< endl;
    ofs.close();

}

void PharmWriter::writePos(int nb,PharmProp prop,const char * namefile){
    ofstream ofs;
    ofs.open(namefile,ios::out|ios::app);
    //POS
    ofs <<"  AND primary_amine"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #N#31."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 2 )"<< endl <<"    LPCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #C#41."<<nb<<""<< endl <<"    ELEMENTS ( C )"<< endl <<"    ALIPHATIC"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORDINATION ( 4 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #H#51."<<nb<< endl <<"    ELEMENTS ( H )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #H#62."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 2-31."<<nb<<" ( #N#31."<<nb<<" #C#41."<<nb<<" ) ;"<< endl <<"      BOND 2-41."<<nb<<" ( #N#31."<<nb<<" #H#51."<<nb<<" ) ;"<< endl <<"      BOND 2-51."<<nb<<" ( #N#31."<<nb<<" #H#62."<<nb<<" ) ;"<< endl <<"    OBJECTS :"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND secondary_amine"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #N#81."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 1 )"<< endl <<"    LPCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #C#91."<<nb<<""<< endl <<"    ELEMENTS ( C )"<< endl <<"    ALIPHATIC"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORDINATION ( 4 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #C#101."<<nb<<""<< endl <<"    ELEMENTS ( C )"<< endl <<"    ALIPHATIC"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORDINATION ( 4 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #H#111."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 7-81."<<nb<<" ( #N#81."<<nb<<" #C#91."<<nb<<" ) ;"<< endl <<"      BOND 7-92."<<nb<<" ( #N#81."<<nb<<" #C#101."<<nb<<" ) ;"<< endl <<"      BOND 7-102."<<nb<<" ( #N#81."<<nb<<" #H#111."<<nb<<" ) ;"<< endl <<"    OBJECTS :"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND tertiary_amine"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #N#131."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #C#141."<<nb<<""<< endl <<"    ELEMENTS ( C )"<< endl <<"    ALIPHATIC"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORDINATION ( 4 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #C#151."<<nb<<""<< endl <<"    ELEMENTS ( C )"<< endl <<"    ALIPHATIC"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORDINATION ( 4 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #C#161."<<nb<<""<< endl <<"    ELEMENTS ( C )"<< endl <<"    ALIPHATIC"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORDINATION ( 4 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 12-131."<<nb<<" ( #N#131."<<nb<<" #C#141."<<nb<<" ) ;"<< endl <<"      BOND 12-141."<<nb<<" ( #N#131."<<nb<<" #C#151."<<nb<<" ) ;"<< endl <<"      BOND 12-151."<<nb<<" ( #N#131."<<nb<<" #C#161."<<nb<<" ) ;"<< endl <<"    OBJECTS :"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND amidine"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #N#191."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 0 1 )"<< endl <<"    LPCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #C#201."<<nb<<""<< endl <<"    ELEMENTS ( C )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#211."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 18-191."<<nb<<" ( #N#191."<<nb<<" #C#201."<<nb<<" ) DOUBLE ;"<< endl <<"      BOND 19-201."<<nb<<" ( #C#201."<<nb<<" #N#211."<<nb<<" ) ;"<< endl <<"    OBJECTS :"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND guanidino"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #H#241."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#252."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 1 2 )"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #C#261."<<nb<<""<< endl <<"    ELEMENTS ( C )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#272."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 0 1 )"<< endl <<"    LPCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#281."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 1 2 )"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #H#311."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 24-231."<<nb<<" ( #N#252."<<nb<<" #H#241."<<nb<<" ) ;"<< endl <<"      BOND 24-252."<<nb<<" ( #N#252."<<nb<<" #C#261."<<nb<<" ) ;"<< endl <<"      BOND 25-262."<<nb<<" ( #C#261."<<nb<<" #N#272."<<nb<<" ) DOUBLE ;"<< endl <<"      BOND 25-271."<<nb<<" ( #C#261."<<nb<<" #N#281."<<nb<<" ) ;"<< endl <<"      BOND 27-301."<<nb<<" ( #N#281."<<nb<<" #H#311."<<nb<<" ) ;"<< endl <<"    OBJECTS :"<< endl <<"      POINT pointP-11."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 )"<< endl <<"    CENTROID ( #N#252."<<nb<<" #N#272."<<nb<<" #N#281."<<nb<<" ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND amidine-H"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #H#361."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#381."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 0 1 )"<< endl <<"    LPCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #C#391."<<nb<<""<< endl <<"    ELEMENTS ( C )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#401."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 1 2 )"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 39-351."<<nb<<" ( #N#401."<<nb<<" #H#361."<<nb<<" ) ;"<< endl <<"      BOND 37-381."<<nb<<" ( #N#381."<<nb<<" #C#391."<<nb<<" ) DOUBLE ;"<< endl <<"      BOND 38-391."<<nb<<" ( #C#391."<<nb<<" #N#401."<<nb<<" ) ;"<< endl <<"    OBJECTS :"<< endl <<"      POINT pointP-21."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 )"<< endl <<"    CENTROID ( #N#381."<<nb<<" #N#401."<<nb<<" ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND Exclude_guanidino"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #N#431."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 1 2 )"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #C#441."<<nb<<""<< endl <<"    ELEMENTS ( C )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#451."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 0 1 )"<< endl <<"    LPCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#461."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    ALIPHATIC"<< endl <<"    HCOUNT ( 1 2 )"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #H#471."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #H#501."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND P42-431."<<nb<<" ( #N#431."<<nb<<" #C#441."<<nb<<" ) ;"<< endl <<"      BOND P43-441."<<nb<<" ( #C#441."<<nb<<" #N#451."<<nb<<" ) DOUBLE ;"<< endl <<"      BOND P43-451."<<nb<<" ( #C#441."<<nb<<" #N#461."<<nb<<" ) ;"<< endl <<"      BOND 42-461."<<nb<<" ( #N#431."<<nb<<" #H#471."<<nb<<" ) ;"<< endl <<"      BOND 45-491."<<nb<<" ( #N#461."<<nb<<" #H#501."<<nb<<" ) ;"<< endl <<"    OBJECTS :"<< endl <<"      POINT point-2ex1."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 )"<< endl <<"    CENTROID ( #N#431."<<nb<<" #N#451."<<nb<<" ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  EXCLUDE amidineH_withExclusions"<<nb<<" ( amidine-H"<<nb<<" Exclude_guanidino"<<nb<<" )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    PRESENT ( amidine-H"<<nb<<" ) ;"<< endl <<"    ABSENT ( Exclude_guanidino"<<nb<<" ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"      MAP amidine-H"<<nb<<"/#H#361."<<nb<<""<< endl <<"    ( Exclude_guanidino"<<nb<<"/#H#501."<<nb<<" ) ;"<< endl <<"      MAP amidine-H"<<nb<<"/#N#381."<<nb<<""<< endl <<"    ( Exclude_guanidino"<<nb<<"/#N#451."<<nb<<" ) ;"<< endl <<"      MAP amidine-H"<<nb<<"/#C#391."<<nb<<""<< endl <<"    ( Exclude_guanidino"<<nb<<"/#C#441."<<nb<<" ) ;"<< endl <<"      MAP amidine-H"<<nb<<"/#N#401."<<nb<<""<< endl <<"    ( Exclude_guanidino"<<nb<<"/#N#461."<<nb<<" ) ;"<< endl <<"      MAP amidine-H"<<nb<<"/39-351."<<nb<<""<< endl <<"    ( Exclude_guanidino"<<nb<<"/45-491."<<nb<<" ) ;"<< endl <<"      MAP amidine-H"<<nb<<"/37-381."<<nb<<""<< endl <<"    ( Exclude_guanidino"<<nb<<"/P43-441."<<nb<<" ) ;"<< endl <<"      MAP amidine-H"<<nb<<"/38-391."<<nb<<""<< endl <<"    ( Exclude_guanidino"<<nb<<"/P43-451."<<nb<<" ) ;"<< endl <<"  }"<< endl <<"  AND pos"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #Any_AtomP#31."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    POSITIVE"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    OBJECTS :"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND Exclude_Pos"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #Any_AtomP#51."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    POSITIVE"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #Any_AtomP#61."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    NEGATIVE"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 4-52."<<nb<<" ( #Any_AtomP#51."<<nb<<" #Any_AtomP#61."<<nb<<" ) SINGLE DOUBLE AROMATIC TRIPLE ;"<< endl <<"    OBJECTS :"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  EXCLUDE positive_withExclusions"<<nb<<" ( pos"<<nb<<" Exclude_Pos"<<nb<<" )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      PROP ( SOURCE \"BIOCAD_FUNCTION\" ) ;"<< endl <<"    PRESENT ( pos"<<nb<<" ) ;"<< endl <<"    ABSENT ( Exclude_Pos"<<nb<<" ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"      MAP pos"<<nb<<"/#Any_AtomP#31."<<nb<<""<< endl <<"    ( Exclude_Pos"<<nb<<"/#Any_AtomP#51."<<nb<<" ) ;"<< endl <<"  }"<< endl  ;

    ofs<< setprecision(1) << fixed<<"  OR Positive"<<nb<<" ( primary_amine"<<nb<<" secondary_amine"<<nb<<" tertiary_amine"<<nb<<" amidine"<<nb<<" guanidino"<<nb<<" amidineH_withExclusions"<<nb<<" positive_withExclusions"<<nb<<" )"<< endl
       <<"  {"<< endl
      <<"    ATTRIBUTES :"<< endl
     <<"      COORD2D ( 0 0 ) ;"<< endl
    <<"      COORD3D ( "<<prop.getBase().x*100<<" "<<prop.getBase().y*100<<" "<<prop.getBase().z*100<<" ) ;"<< endl
    <<"      PROP ( SOURCE \"BIOCAD_FUNCTION\" ) ;"<< endl
    <<"    OBJECTS :"<< endl
    <<"      POINT PosIon"<<nb<< endl
    <<"    COORD2D ( 0 0 )"<< endl
    <<"    COORD3D ( "<<prop.getBase().x*100<<" "<<prop.getBase().y*100<<" "<<prop.getBase().z*100<<" )"<< endl
    <<"    MAP ( primary_amine"<<nb<<"/#N#31."<<nb<<" secondary_amine"<<nb<<"/#N#81."<<nb<<" tertiary_amine"<<nb<<"/#N#131."<<nb<<" amidine"<<nb<<"/#N#191."<<nb<<" guanidino"<<nb<<"/pointP-11."<<nb<<" amidineH_withExclusions"<<nb<<"/amidine-H"<<nb<<"/pointP-21."<<nb<<" positive_withExclusions"<<nb<<"/pos"<<nb<<"/#Any_AtomP#31."<<nb<<" )"<< endl
    <<"    PROP ( INDEX_CODE \"18\" ) ;"<< endl
    <<"  }"<< endl;
    ofs.close();
}

void PharmWriter::writeNeg(int nb,PharmProp prop,const char * namefile){
    ofstream ofs;
    ofs.open(namefile,ios::out|ios::app);
    //NEG
    ofs <<"  AND trifluorosulfonamide"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #F#31."<<nb<<""<< endl <<"    ELEMENTS ( F )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #F#41."<<nb<<""<< endl <<"    ELEMENTS ( F )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #C#51."<<nb<<""<< endl <<"    ELEMENTS ( C )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #H#63."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#71."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #O#81."<<nb<<"."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #O#91."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #F#101."<<nb<<""<< endl <<"    ELEMENTS ( F )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #S#111."<<nb<<""<< endl <<"    ELEMENTS ( S )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 5-61."<<nb<<" ( #H#63."<<nb<<" #N#71."<<nb<<" ) ;"<< endl <<"      BOND 10-61."<<nb<<" ( #S#111."<<nb<<" #N#71."<<nb<<" ) EXOCYCLIC ;"<< endl <<"      BOND 10-71."<<nb<<" ( #S#111."<<nb<<" #O#81."<<nb<<"."<<nb<<" ) DOUBLE EXOCYCLIC ;"<< endl <<"      BOND 10-81."<<nb<<" ( #S#111."<<nb<<" #O#91."<<nb<<" ) DOUBLE EXOCYCLIC ;"<< endl <<"      BOND 4-31."<<nb<<" ( #C#51."<<nb<<" #F#41."<<nb<<" ) EXOCYCLIC ;"<< endl <<"      BOND 4-101."<<nb<<" ( #C#51."<<nb<<" #S#111."<<nb<<" ) EXOCYCLIC ;"<< endl <<"      BOND 4-91."<<nb<<" ( #C#51."<<nb<<" #F#101."<<nb<<" ) EXOCYCLIC ;"<< endl <<"      BOND 4-21."<<nb<<" ( #C#51."<<nb<<" #F#31."<<nb<<" ) EXOCYCLIC ;"<< endl <<"    OBJECTS :"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND sulfonic_acid"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #O#131."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #H#141."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #O#151."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #S#161."<<nb<<""<< endl <<"    ELEMENTS ( S )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #O#171."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 15-161."<<nb<<" ( #S#161."<<nb<<" #O#171."<<nb<<" ) DOUBLE EXOCYCLIC ;"<< endl <<"      BOND 14-151."<<nb<<" ( #O#151."<<nb<<" #S#161."<<nb<<" ) EXOCYCLIC ;"<< endl <<"      BOND 14-131."<<nb<<" ( #O#151."<<nb<<" #H#141."<<nb<<" ) ;"<< endl <<"      BOND 15-121."<<nb<<" ( #S#161."<<nb<<" #O#131."<<nb<<" ) DOUBLE EXOCYCLIC ;"<< endl <<"    OBJECTS :"<< endl <<"      POINT pointN-11."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 )"<< endl <<"    CENTROID ( #O#151."<<nb<<" #O#171."<<nb<<" #O#131."<<nb<<" ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND phosphoryl"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #H#201."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #O#211."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #H#221."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #O#231."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #P#241."<<nb<<""<< endl <<"    ELEMENTS ( P )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #O#251."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 24-191."<<nb<<" ( #O#251."<<nb<<" #H#201."<<nb<<" ) ;"<< endl <<"      BOND 23-243."<<nb<<" ( #P#241."<<nb<<" #O#251."<<nb<<" ) EXOCYCLIC ;"<< endl <<"      BOND 22-233."<<nb<<" ( #O#231."<<nb<<" #P#241."<<nb<<" ) EXOCYCLIC ;"<< endl <<"      BOND 22-211."<<nb<<" ( #O#231."<<nb<<" #H#221."<<nb<<" ) ;"<< endl <<"      BOND 23-201."<<nb<<" ( #P#241."<<nb<<" #O#211."<<nb<<" ) DOUBLE EXOCYCLIC ;"<< endl <<"    OBJECTS :"<< endl <<"      POINT pointN-21."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 )"<< endl <<"    CENTROID ( #O#231."<<nb<<" #O#211."<<nb<<" #O#251."<<nb<<" ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND acidic_OH-1"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #O#291."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #H#301."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #O#311."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #S#321."<<nb<<""<< endl <<"    ELEMENTS ( C P S )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 30-311 ( #O#311."<<nb<<" #S#321."<<nb<<" ) EXOCYCLIC ;"<< endl <<"      BOND 30-291 ( #O#311."<<nb<<" #H#301."<<nb<<" ) ;"<< endl <<"      BOND 31-281 ( #S#321."<<nb<<" #O#291."<<nb<<" ) DOUBLE EXOCYCLIC ;"<< endl <<"    OBJECTS :"<< endl <<"      POINT pnt-11."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 )"<< endl <<"    CENTROID ( #O#291."<<nb<<" #O#311."<<nb<<" ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND Exclude_sulfonic_acid"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #O#341."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #H#351."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #O#363."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #S#371."<<nb<<""<< endl <<"    ELEMENTS ( S )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #O#381."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 36-371."<<nb<<" ( #S#371."<<nb<<" #O#381."<<nb<<" ) DOUBLE EXOCYCLIC ;"<< endl <<"      BOND 35-361."<<nb<<" ( #O#363."<<nb<<" #S#371."<<nb<<" ) EXOCYCLIC ;"<< endl <<"      BOND 35-341."<<nb<<" ( #O#363."<<nb<<" #H#351."<<nb<<" ) ;"<< endl <<"      BOND 36-331."<<nb<<" ( #S#371."<<nb<<" #O#341."<<nb<<" ) DOUBLE EXOCYCLIC ;"<< endl <<"    OBJECTS :"<< endl <<"      POINT pnt-21."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 )"<< endl <<"    CENTROID ( #O#381."<<nb<<" #O#341."<<nb<<" ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND Exclude_phosphoryl"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #H#401."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #O#411."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #H#421."<<nb<<""<< endl <<"    ELEMENTS ( H )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #O#431."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #P#441."<<nb<<""<< endl <<"    ELEMENTS ( P )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #O#451."<<nb<<""<< endl <<"    ELEMENTS ( O )"<< endl <<"    ALIPHATIC"<< endl <<"    EXOCYCLIC"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 40-391."<<nb<<" ( #O#411."<<nb<<" #H#401."<<nb<<" ) ;"<< endl <<"      BOND N43-441."<<nb<<" ( #P#441."<<nb<<" #O#451."<<nb<<" ) DOUBLE EXOCYCLIC ;"<< endl <<"      BOND N42-431."<<nb<<" ( #O#431."<<nb<<" #P#441."<<nb<<" ) EXOCYCLIC ;"<< endl <<"      BOND 42-411."<<nb<<" ( #O#431."<<nb<<" #H#421."<<nb<<" ) ;"<< endl <<"      BOND 43-401."<<nb<<" ( #P#441."<<nb<<" #O#411."<<nb<<" ) EXOCYCLIC ;"<< endl <<"    OBJECTS :"<< endl <<"      POINT pnt-31."<<nb<<""<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 )"<< endl <<"    CENTROID ( #O#451."<<nb<<" #O#431."<<nb<<" ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  EXCLUDE acidic_OH_withExclusions"<<nb<<" ( acidic_OH-1"<<nb<<" Exclude_sulfonic_acid"<<nb<<" Exclude_phosphoryl"<<nb<<" )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    PRESENT ( acidic_OH-1"<<nb<<" ) ;"<< endl <<"    ABSENT ( Exclude_sulfonic_acid1 Exclude_phosphoryl"<<nb<<" ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"      MAP acidic_OH-1"<<nb<<"/#O#291."<<nb<<""<< endl <<"    ( Exclude_sulfonic_acid"<<nb<<"/#O#341."<<nb<<" Exclude_phosphoryl"<<nb<<"/#O#451."<<nb<<" ) ;"<< endl <<"      MAP acidic_OH-1"<<nb<<"/#H#301."<<nb<<""<< endl <<"    ( Exclude_sulfonic_acid"<<nb<<"/#H#351."<<nb<<" Exclude_phosphoryl"<<nb<<"/#H#421."<<nb<<" ) ;"<< endl <<"      MAP acidic_OH-1"<<nb<<"/#O#311."<<nb<<""<< endl <<"    ( Exclude_sulfonic_acid"<<nb<<"/#O#363."<<nb<<" Exclude_phosphoryl"<<nb<<"/#O#431."<<nb<<" ) ;"<< endl <<"      MAP acidic_OH-1"<<nb<<"/#S#321."<<nb<<""<< endl <<"    ( Exclude_sulfonic_acid"<<nb<<"/#S#371."<<nb<<" Exclude_phosphoryl"<<nb<<"/#P#441."<<nb<<" ) ;"<< endl <<"      MAP acidic_OH-1"<<nb<<"/30-311"<< endl <<"    ( Exclude_sulfonic_acid"<<nb<<"/35-361."<<nb<<" Exclude_phosphoryl"<<nb<<"/N42-431."<<nb<<" ) ;"<< endl <<"      MAP acidic_OH-1"<<nb<<"/30-291"<< endl <<"    ( Exclude_sulfonic_acid"<<nb<<"/35-341."<<nb<<" Exclude_phosphoryl"<<nb<<"/42-411."<<nb<<" ) ;"<< endl <<"      MAP acidic_OH-1"<<nb<<"/31-281"<< endl <<"    ( Exclude_sulfonic_acid"<<nb<<"/36-331."<<nb<<" Exclude_phosphoryl"<<nb<<"/N43-441."<<nb<<" ) ;"<< endl <<"      MAP acidic_OH-1"<<nb<<"/pnt-11."<<nb<<""<< endl <<"    ( Exclude_sulfonic_acid"<<nb<<"/pnt-21."<<nb<<" Exclude_phosphoryl"<<nb<<"/pnt-31."<<nb<<" ) ;"<< endl <<"  }"<< endl <<"  AND tetrazole_N-1"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #C#511."<<nb<<""<< endl <<"    ELEMENTS ( C )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#521."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 1 )"<< endl <<"    LPCOUNT ( 0 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#531."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 1 )"<< endl <<"    LPCOUNT ( 0 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#541."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 1 )"<< endl <<"    LPCOUNT ( 0 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#551."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 1 )"<< endl <<"    LPCOUNT ( 0 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 50-511."<<nb<<" ( #C#511."<<nb<<" #N#521."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 51-521."<<nb<<" ( #N#521."<<nb<<" #N#531."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 52-531."<<nb<<" ( #N#531."<<nb<<" #N#541."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 53-541."<<nb<<" ( #N#541."<<nb<<" #N#551."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 50-541."<<nb<<" ( #C#511."<<nb<<" #N#551."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"    OBJECTS :"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND Exclude_redundant_mappingN1"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #C#611."<<nb<<""<< endl <<"    ELEMENTS ( C )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#621."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#631."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#641."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#651."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 60-611."<<nb<<" ( #C#611."<<nb<<" #N#621."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 61-621."<<nb<<" ( #N#621."<<nb<<" #N#631."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 62-631."<<nb<<" ( #N#631."<<nb<<" #N#641."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 63-641."<<nb<<" ( #N#641."<<nb<<" #N#651."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 60-641."<<nb<<" ( #C#611."<<nb<<" #N#651."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"    OBJECTS :"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  EXCLUDE tetrazole_N1_withExclusions"<<nb<<" ( tetrazole_N-1"<<nb<<" Exclude_redundant_mappingN1"<<nb<<" )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    PRESENT ( tetrazole_N-11 ) ;"<< endl <<"    ABSENT ( Exclude_redundant_mappingN1"<<nb<<" ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"      MAP tetrazole_N-1"<<nb<<"/#C#511."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN1"<<nb<<"/#C#611."<<nb<<" ) ;"<< endl <<"      MAP tetrazole_N-1"<<nb<<"/#N#521."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN1"<<nb<<"/#N#621."<<nb<<" ) ;"<< endl <<"      MAP tetrazole_N-1"<<nb<<"/#N#531."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN1"<<nb<<"/#N#631."<<nb<<" ) ;"<< endl <<"      MAP tetrazole_N-1"<<nb<<"/#N#541."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN1"<<nb<<"/#N#641."<<nb<<" ) ;"<< endl <<"      MAP tetrazole_N-1"<<nb<<"/#N#551."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN1"<<nb<<"/#N#651."<<nb<<" ) ;"<< endl <<"      MAP tetrazole_N-1"<<nb<<"/50-511."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN1"<<nb<<"/60-611."<<nb<<" ) ;"<< endl <<"      MAP tetrazole_N-1"<<nb<<"/51-521."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN1"<<nb<<"/61-621."<<nb<<" ) ;"<< endl <<"      MAP tetrazole_N-1"<<nb<<"/52-531."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN1"<<nb<<"/62-631."<<nb<<" ) ;"<< endl <<"      MAP tetrazole_N-1"<<nb<<"/53-541."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN1"<<nb<<"/63-641."<<nb<<" ) ;"<< endl <<"      MAP tetrazole_N-1"<<nb<<"/50-541."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN1"<<nb<<"/60-641."<<nb<<" ) ;"<< endl <<"  }"<< endl <<"  AND tetrazole_N-2"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #C#721."<<nb<<""<< endl <<"    ELEMENTS ( C )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#731."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 1 )"<< endl <<"    LPCOUNT ( 0 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#741."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 1 )"<< endl <<"    LPCOUNT ( 0 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#751."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 1 )"<< endl <<"    LPCOUNT ( 0 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#761."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 1 )"<< endl <<"    LPCOUNT ( 0 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 71-721."<<nb<<" ( #C#721."<<nb<<" #N#731."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 72-731."<<nb<<" ( #N#731."<<nb<<" #N#741."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 73-741."<<nb<<" ( #N#741."<<nb<<" #N#751."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 74-751."<<nb<<" ( #N#751."<<nb<<" #N#761."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 71-751."<<nb<<" ( #C#721."<<nb<<" #N#761."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"    OBJECTS :"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND Exclude_redundant_mappingN2"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #C#821."<<nb<<""<< endl <<"    ELEMENTS ( C )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#831."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#841."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#851."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #N#861."<<nb<<""<< endl <<"    ELEMENTS ( N )"<< endl <<"    AROMATIC"<< endl <<"    ENDOCYCLIC"<< endl <<"    BRIDGEHEAD"<< endl <<"    HCOUNT ( 0 )"<< endl <<"    LPCOUNT ( 0 1 )"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 81-821."<<nb<<" ( #C#821."<<nb<<" #N#831."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 82-831."<<nb<<" ( #N#831."<<nb<<" #N#841."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 83-841."<<nb<<" ( #N#841."<<nb<<" #N#851."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 84-851."<<nb<<" ( #N#851."<<nb<<" #N#861."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      BOND 81-851."<<nb<<" ( #C#821."<<nb<<" #N#861."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"    OBJECTS :"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  EXCLUDE tetrazole_N2_withExclusions"<<nb<<" ( tetrazole_N-2"<<nb<<" Exclude_redundant_mappingN2"<<nb<<" )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    PRESENT ( tetrazole_N-21 ) ;"<< endl <<"    ABSENT ( Exclude_redundant_mappingN21 ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"      MAP tetrazole_N-2"<<nb<<"/#C#721."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN2"<<nb<<"/#C#821."<<nb<<" ) ;"<< endl <<"      MAP tetrazole_N-2"<<nb<<"/#N#731."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN2"<<nb<<"/#N#831."<<nb<<" ) ;"<< endl <<"      MAP tetrazole_N-2"<<nb<<"/#N#741."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN2"<<nb<<"/#N#841."<<nb<<" ) ;"<< endl <<"      MAP tetrazole_N-2"<<nb<<"/#N#751."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN2"<<nb<<"/#N#851."<<nb<<" ) ;"<< endl <<"      MAP tetrazole_N-2"<<nb<<"/#N#761."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN2"<<nb<<"/#N#861."<<nb<<" ) ;"<< endl <<"      MAP tetrazole_N-2"<<nb<<"/71-721."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN2"<<nb<<"/81-821."<<nb<<" ) ;"<< endl <<"      MAP tetrazole_N-2"<<nb<<"/72-731."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN2"<<nb<<"/82-831."<<nb<<" ) ;"<< endl <<"      MAP tetrazole_N-2"<<nb<<"/73-741."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN2"<<nb<<"/83-841."<<nb<<" ) ;"<< endl <<"      MAP tetrazole_N-2"<<nb<<"/74-751."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN2"<<nb<<"/84-851."<<nb<<" ) ;"<< endl <<"      MAP tetrazole_N-2"<<nb<<"/71-751."<<nb<<""<< endl <<"    ( Exclude_redundant_mappingN2"<<nb<<"/81-851."<<nb<<" ) ;"<< endl <<"  }"<< endl <<"  AND neg"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #Any_AtomN#31."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    NEGATIVE"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    OBJECTS :"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  AND Exclude_Neg"<<nb<<" ( )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"    TOPOLOGY :"<< endl <<"      ATOM #Any_AtomN#51."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    POSITIVE"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      ATOM #Any_AtomN#61."<<nb<<""<< endl <<"    ELEMENTS ( ANY )"<< endl <<"    NEGATIVE"<< endl <<"    COORD2D ( 0 0 )"<< endl <<"    COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      BOND 4-53 ( #Any_AtomN#51."<<nb<<" #Any_AtomN#61."<<nb<<" ) SINGLE DOUBLE AROMATIC TRIPLE ;"<< endl <<"    OBJECTS :"<< endl <<"    CONSTRAINTS :"<< endl <<"  }"<< endl <<"  EXCLUDE negative_withExclusions"<<nb<<" ( neg"<<nb<<" Exclude_Neg"<<nb<<" )"<< endl <<"  {"<< endl <<"    ATTRIBUTES :"<< endl <<"      COORD2D ( 0 0 ) ;"<< endl <<"      COORD3D ( 0.0 0.0 0.0 ) ;"<< endl <<"      PROP ( SOURCE \"BIOCAD_FUNCTION\" ) ;"<< endl <<"    PRESENT ( neg"<<nb<<" ) ;"<< endl <<"    ABSENT ( Exclude_Neg"<<nb<<" ) ;"<< endl <<"    CONSTRAINTS :"<< endl <<"      MAP neg"<<nb<<"/#Any_AtomN#31."<<nb<<""<< endl <<"    ( Exclude_Neg"<<nb<<"/#Any_AtomN#61."<<nb<<" ) ;"<< endl <<"  }"<< endl <<""<< endl ;

    ofs << setprecision(1) << fixed <<"  OR Negative"<<nb<<" ( trifluorosulfonamide"<<nb<<" sulfonic_acid"<<nb<<" phosphoryl"<<nb<<" acidic_OH_withExclusions"<<nb<<" tetrazole_N1_withExclusions"<<nb<<" tetrazole_N2_withExclusions"<<nb<<" negative_withExclusions"<<nb<<" )"<< endl
        <<"  {"<< endl
       <<"    ATTRIBUTES :"<< endl
      <<"      COORD2D ( 0 0 ) ;"<< endl
     <<"      COORD3D ( "<<prop.getBase().x*100<<" "<<prop.getBase().y*100<<" "<<prop.getBase().z*100<<" ) ;"<< endl
    <<"      PROP ( SOURCE \"BIOCAD_FUNCTION\" ) ;"<< endl
    <<"    OBJECTS :"<< endl
    <<"      POINT NegIon"<<nb<< endl
    <<"    COORD2D ( 0 0 )"<< endl
    <<"      COORD3D ( "<<prop.getBase().x*100<<" "<<prop.getBase().y*100<<" "<<prop.getBase().z*100<<" )"<< endl
    <<"    MAP ( trifluorosulfonamide"<<nb<<"/#N#71."<<nb<<" sulfonic_acid"<<nb<<"/pointN-11."<<nb<<" phosphoryl"<<nb<<"/pointN-21."<<nb<<" acidic_OH_withExclusions"<<nb<<"/acidic_OH-1"<<nb<<"/pnt-11."<<nb<<" tetrazole_N1_withExclusions"<<nb<<"/tetrazole_N-1"<<nb<<"/#N#551."<<nb<<" tetrazole_N2_withExclusions"<<nb<<"/tetrazole_N-2"<<nb<<"/#N#751."<<nb<<" negative_withExclusions"<<nb<<"/neg"<<nb<<"/#Any_AtomN#31."<<nb<<" )"<< endl
    <<"    PROP ( INDEX_CODE \"19\" ) ;"<< endl
    <<"  }"<< endl;
    ofs.close();
}


void PharmWriter::writehyd(int nb,PharmProp prop,const char * namefile){
    ofstream ofs;
    ofs.open(namefile,ios::out|ios::app);

    ofs<< setprecision(1) << fixed<<"  AND HYDROPHOBIC"<<nb<<" ( )"<< endl
       <<"  {"<< endl
      <<"    ATTRIBUTES :"<< endl
     <<"      COORD2D ( 0 0 ) ;"<< endl
    <<"      COORD3D ( "<<prop.getBase().x*100<<" "<<prop.getBase().y*100<<" "<<prop.getBase().z*100<<" ) ;"<< endl
    <<"      PROP ( SOURCE \"BIOCAD_FUNCTION\" ) ;"<< endl
    <<"    TOPOLOGY :"<< endl
    <<"    OBJECTS :"<< endl
    <<"      POINT Hydrophobic"<<nb<< endl
    <<"    COORD2D ( 0 0 )"<< endl
    <<"    COORD3D ( "<<prop.getBase().x*100<<" "<<prop.getBase().y*100<<" "<<prop.getBase().z*100<<" )"<< endl
    <<"    HYDROPHOBE"<< endl
    <<"    PROP ( INDEX_CODE \"11\" ) ;"<< endl
    <<"    CONSTRAINTS :"<< endl
    <<"  }"<< endl;
    ofs.close();

}
void PharmWriter::writemet(int nb,PharmProp prop,const char * namefile){
    ofstream ofs;
    ofs.open(namefile,ios::out|ios::app);
    //NEG
    ofs <<"  AND sulfur"<<nb<<" ( )"<< endl <<"    {"<< endl <<"      ATTRIBUTES :"<< endl <<"        COORD2D ( 0 0 ) ;"<< endl <<"        COORD3D ( -364.8 294.7 43.5 ) ;"<< endl <<"        PROP ( SOURCE \"BIOCAD_FUNCTION\" ) ;"<< endl <<"      TOPOLOGY :"<< endl <<"        ATOM #S#3."<<nb<<""<< endl <<"      ELEMENTS ( S )"<< endl <<"      ALIPHATIC"<< endl <<"      EXOCYCLIC"<< endl <<"      LPCOUNT ( 2 )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( -364.8 294.7 43.5 ) ;"<< endl <<"      OBJECTS :"<< endl <<"      CONSTRAINTS :"<< endl <<"    }"<< endl <<"    AND CONO"<<nb<<" ( )"<< endl <<"    {"<< endl <<"      ATTRIBUTES :"<< endl <<"        COORD2D ( 0 0 ) ;"<< endl <<"        COORD3D ( -304.1 -93.5 -5.2 ) ;"<< endl <<"        PROP ( SOURCE \"UNASSIGNED\" ) ;"<< endl <<"      TOPOLOGY :"<< endl <<"        ATOM #O#6."<<nb<<""<< endl <<"      ELEMENTS ( O )"<< endl <<"      ALIPHATIC"<< endl <<"      EXOCYCLIC"<< endl <<"      HCOUNT ( 1 )"<< endl <<"      LPCOUNT ( 2 )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( -349.9 -219.2 21.3 ) ;"<< endl <<"        ATOM #H#7."<<nb<<""<< endl <<"      ELEMENTS ( H )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( -295.9 -292.2 -9.9 ) ;"<< endl <<"        ATOM #H#8."<<nb<<""<< endl <<"      ELEMENTS ( H )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( -204.4 -101.5 -62.9 ) ;"<< endl <<"        ATOM #N#9."<<nb<<""<< endl <<"      ELEMENTS ( N )"<< endl <<"      ALIPHATIC"<< endl <<"      EXOCYCLIC"<< endl <<"      HCOUNT ( 1 )"<< endl <<"      LPCOUNT ( 0 )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( -291.0 -101.4 -12.8 ) ;"<< endl <<"        ATOM #O#10."<<nb<<""<< endl <<"      ELEMENTS ( O )"<< endl <<"      ALIPHATIC"<< endl <<"      EXOCYCLIC"<< endl <<"      HCOUNT ( 0 )"<< endl <<"      LPCOUNT ( 2 )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( -301.3 132.4 -6.9 ) ;"<< endl <<"        ATOM #C#11."<<nb<<""<< endl <<"      ELEMENTS ( C )"<< endl <<"      ALIPHATIC"<< endl <<"      EXOCYCLIC"<< endl <<"      LPCOUNT ( 0 )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( -354.6 25.9 24.0 ) ;"<< endl <<"        BOND 7-8."<<nb<<" ( #H#8."<<nb<<" #N#9."<<nb<<" ) SINGLE ;"<< endl <<"        BOND 10-8."<<nb<<" ( #C#11."<<nb<<" #N#9."<<nb<<" ) SINGLE EXOCYCLIC ;"<< endl <<"        BOND 10-9."<<nb<<" ( #C#11."<<nb<<" #O#10."<<nb<<" ) DOUBLE EXOCYCLIC ;"<< endl <<"        BOND 5-6."<<nb<<" ( #O#6."<<nb<<" #H#7."<<nb<<" ) SINGLE ;"<< endl <<"        BOND 5-8."<<nb<<" ( #O#6."<<nb<<" #N#9."<<nb<<" ) SINGLE EXOCYCLIC ;"<< endl <<"      OBJECTS :"<< endl <<"        POINT point-1"<<nb<<""<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( -331.8 -98.2 10.8 )"<< endl <<"      CENTROID ( #C#11."<<nb<<" #N#9."<<nb<<" #O#6."<<nb<<" ) ;"<< endl <<"      CONSTRAINTS :"<< endl <<"    }"<< endl <<"    AND carboxyl"<<nb<<" ( )"<< endl <<"    {"<< endl <<"      ATTRIBUTES :"<< endl <<"        COORD2D ( 0 0 ) ;"<< endl <<"        COORD3D ( 42.1 245.8 -5.5 ) ;"<< endl <<"        PROP ( SOURCE \"UNASSIGNED\" ) ;"<< endl <<"      TOPOLOGY :"<< endl <<"        ATOM #H#14."<<nb<<""<< endl <<"      ELEMENTS ( H )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 27.7 84.8 2.8 ) ;"<< endl <<"        ATOM #O#15."<<nb<<""<< endl <<"      ELEMENTS ( O )"<< endl <<"      ALIPHATIC"<< endl <<"      EXOCYCLIC"<< endl <<"      HCOUNT ( 1 )"<< endl <<"      LPCOUNT ( 1 )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 69.2 167.9 -21.2 ) ;"<< endl <<"        ATOM #O#16."<<nb<<""<< endl <<"      ELEMENTS ( O )"<< endl <<"      ALIPHATIC"<< endl <<"      EXOCYCLIC"<< endl <<"      HCOUNT ( 0 )"<< endl <<"      LPCOUNT ( 2 )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 60.6 398.3 -16.2 ) ;"<< endl <<"        ATOM #C#17."<<nb<<""<< endl <<"      ELEMENTS ( C )"<< endl <<"      ALIPHATIC"<< endl <<"      EXOCYCLIC"<< endl <<"      LPCOUNT ( 0 )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 7.3 291.8 14.6 ) ;"<< endl <<"        BOND 13-14."<<nb<<" ( #H#14."<<nb<<" #O#15."<<nb<<" ) SINGLE ;"<< endl <<"        BOND 16-14."<<nb<<" ( #C#17."<<nb<<" #O#15."<<nb<<" ) SINGLE EXOCYCLIC ;"<< endl <<"        BOND 16-15."<<nb<<" ( #C#17."<<nb<<" #O#16."<<nb<<" ) DOUBLE EXOCYCLIC ;"<< endl <<"      OBJECTS :"<< endl <<"        POINT point-2"<<nb<<""<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 45.7 286.0 -7.6 )"<< endl <<"      CENTROID ( #C#17."<<nb<<" #O#15."<<nb<<" #O#16."<<nb<<" ) ;"<< endl <<"      CONSTRAINTS :"<< endl <<"    }"<< endl <<"    AND NNC(OS)"<<nb<<" ( )"<< endl <<"    {"<< endl <<"      ATTRIBUTES :"<< endl <<"        COORD2D ( 0 0 ) ;"<< endl <<"        COORD3D ( 77.7 -231.9 -8.0 ) ;"<< endl <<"        PROP ( SOURCE \"UNASSIGNED\" ) ;"<< endl <<"      TOPOLOGY :"<< endl <<"        ATOM #N#20."<<nb<<""<< endl <<"      ELEMENTS ( N )"<< endl <<"      ALIPHATIC"<< endl <<"      EXOCYCLIC"<< endl <<"      LPCOUNT ( 1 )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 51.5 -399.9 7.1 ) ;"<< endl <<"        ATOM #N#21."<<nb<<""<< endl <<"      ELEMENTS ( N )"<< endl <<"      ALIPHATIC"<< endl <<"      EXOCYCLIC"<< endl <<"      HCOUNT ( 0 )"<< endl <<"      LPCOUNT ( 1 )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 112.1 -278.7 -27.8 ) ;"<< endl <<"        ATOM #O."<<nb<<",S#22."<<nb<<""<< endl <<"      ELEMENTS ( O S )"<< endl <<"      ALIPHATIC"<< endl <<"      EXOCYCLIC"<< endl <<"      HCOUNT ( 0 )"<< endl <<"      LPCOUNT ( 2 )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 103.5 -48.3 -22.9 ) ;"<< endl <<"        ATOM #C#23."<<nb<<""<< endl <<"      ELEMENTS ( C )"<< endl <<"      ALIPHATIC"<< endl <<"      EXOCYCLIC"<< endl <<"      LPCOUNT ( 0 )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 50.3 -154.8 7.9 ) ;"<< endl <<"        BOND 22-20."<<nb<<" ( #C#23."<<nb<<" #N#21."<<nb<<" ) SINGLE EXOCYCLIC ;"<< endl <<"        BOND 22-21."<<nb<<" ( #C#23."<<nb<<" #O."<<nb<<",S#22."<<nb<<" ) DOUBLE EXOCYCLIC ;"<< endl <<"        BOND 19-20."<<nb<<" ( #N#20."<<nb<<" #N#21."<<nb<<" ) DOUBLE EXOCYCLIC ;"<< endl <<"      OBJECTS :"<< endl <<"        POINT point-3"<<nb<<""<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 71.3 -277.8 -4.3 )"<< endl <<"      CENTROID ( #C#23."<<nb<<" #N#21."<<nb<<" #N#20."<<nb<<" ) ;"<< endl <<"      CONSTRAINTS :"<< endl <<"    }"<< endl <<"    AND phospho-group"<<nb<<" ( )"<< endl <<"    {"<< endl <<"      ATTRIBUTES :"<< endl <<"        COORD2D ( 0 0 ) ;"<< endl <<"        COORD3D ( 506.6 236.9 -18.2 ) ;"<< endl <<"        PROP ( SOURCE \"UNASSIGNED\" ) ;"<< endl <<"      TOPOLOGY :"<< endl <<"        ATOM #C."<<nb<<",N."<<nb<<",O#25."<<nb<<""<< endl <<"      ELEMENTS ( C N O )"<< endl <<"      ALIPHATIC"<< endl <<"      HCOUNT ( 0 )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 513.3 62.4 -20.1 ) ;"<< endl <<"        ATOM #O#26."<<nb<<""<< endl <<"      ELEMENTS ( O )"<< endl <<"      ALIPHATIC"<< endl <<"      EXOCYCLIC"<< endl <<"      LPCOUNT ( 2 )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 513.4 297.0 145.6 ) ;"<< endl <<"        ATOM #O#27."<<nb<<""<< endl <<"      ELEMENTS ( O )"<< endl <<"      ALIPHATIC"<< endl <<"      EXOCYCLIC"<< endl <<"      HCOUNT ( 0 )"<< endl <<"      LPCOUNT ( 2 )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 633.9 292.7 -97.3 ) ;"<< endl <<"        ATOM #P#28."<<nb<<""<< endl <<"      ELEMENTS ( P )"<< endl <<"      ALIPHATIC"<< endl <<"      EXOCYCLIC"<< endl <<"      HCOUNT ( 0 )"<< endl <<"      LPCOUNT ( 0 )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 509.7 238.3 -20.2 ) ;"<< endl <<"        ATOM #C."<<nb<<",N."<<nb<<",O#29."<<nb<<""<< endl <<"      ELEMENTS ( C N O )"<< endl <<"      ALIPHATIC"<< endl <<"      HCOUNT ( 0 )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 362.6 294.1 -99.0 ) ;"<< endl <<"        BOND 27-24."<<nb<<" ( #P#28."<<nb<<" #C."<<nb<<",N."<<nb<<",O#25."<<nb<<" ) SINGLE EXOCYCLIC ;"<< endl <<"        BOND 27-25."<<nb<<" ( #P#28."<<nb<<" #O#26."<<nb<<" ) SINGLE EXOCYCLIC ;"<< endl <<"        BOND 27-26."<<nb<<" ( #P#28."<<nb<<" #O#27."<<nb<<" ) DOUBLE EXOCYCLIC ;"<< endl <<"        BOND 28-27."<<nb<<" ( #C."<<nb<<",N."<<nb<<",O#29."<<nb<<" #P#28."<<nb<<" ) SINGLE EXOCYCLIC ;"<< endl <<"      OBJECTS :"<< endl <<"      CONSTRAINTS :"<< endl <<"    }"<< endl <<"    AND NNarom"<<nb<<" ( )"<< endl <<"    {"<< endl <<"      ATTRIBUTES :"<< endl <<"        COORD2D ( 0 0 ) ;"<< endl <<"        COORD3D ( 595.4 -235.6 14.5 ) ;"<< endl <<"        PROP ( SOURCE \"UNASSIGNED\" ) ;"<< endl <<"      TOPOLOGY :"<< endl <<"        ATOM #H#32."<<nb<<""<< endl <<"      ELEMENTS ( H )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 518.2 -67.4 -2.7 ) ;"<< endl <<"        ATOM #N#33."<<nb<<""<< endl <<"      ELEMENTS ( N )"<< endl <<"      AROMATIC"<< endl <<"      ENDOCYCLIC"<< endl <<"      BRIDGEHEAD"<< endl <<"      HCOUNT ( 1 )"<< endl <<"      LPCOUNT ( 0 )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 560.3 -157.6 6.0 ) ;"<< endl <<"        ATOM #C."<<nb<<",N#34."<<nb<<""<< endl <<"      ELEMENTS ( C N )"<< endl <<"      AROMATIC"<< endl <<"      ENDOCYCLIC"<< endl <<"      BRIDGEHEAD"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 678.3 -196.0 -47.9 ) ;"<< endl <<"        ATOM #C."<<nb<<",N#35."<<nb<<""<< endl <<"      ELEMENTS ( C N )"<< endl <<"      AROMATIC"<< endl <<"      ENDOCYCLIC"<< endl <<"      BRIDGEHEAD"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 699.3 -329.4 -11.1 ) ;"<< endl <<"        ATOM #N#36."<<nb<<""<< endl <<"      ELEMENTS ( N )"<< endl <<"      AROMATIC"<< endl <<"      ENDOCYCLIC"<< endl <<"      BRIDGEHEAD"<< endl <<"      HCOUNT ( 0 )"<< endl <<"      LPCOUNT ( 1 )"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 594.3 -371.1 64.8 ) ;"<< endl <<"        ATOM #C."<<nb<<",N#37."<<nb<<""<< endl <<"      ELEMENTS ( C N )"<< endl <<"      AROMATIC"<< endl <<"      ENDOCYCLIC"<< endl <<"      BRIDGEHEAD"<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 509.0 -264.1 74.8 ) ;"<< endl <<"        BOND 36-32."<<nb<<" ( #C."<<nb<<",N#37."<<nb<<" #N#33."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"        BOND 31-32."<<nb<<" ( #H#32."<<nb<<" #N#33."<<nb<<" ) SINGLE ;"<< endl <<"        BOND 33-32."<<nb<<" ( #C."<<nb<<",N#34."<<nb<<" #N#33."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"        BOND 34-33."<<nb<<" ( #C."<<nb<<",N#35."<<nb<<" #C."<<nb<<",N#34."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"        BOND 35-34."<<nb<<" ( #N#36."<<nb<<" #C."<<nb<<",N#35."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"        BOND 36-35."<<nb<<" ( #C."<<nb<<",N#37."<<nb<<" #N#36."<<nb<<" ) AROMATIC ENDOCYCLIC ;"<< endl <<"      OBJECTS :"<< endl <<"        POINT point-4"<<nb<<""<< endl <<"      COORD2D ( 0 0 )"<< endl <<"      COORD3D ( 608.2 -263.6 17.3 )"<< endl <<"      CENTROID ( #N#33."<<nb<<" #C."<<nb<<",N#34."<<nb<<" #C."<<nb<<",N#35."<<nb<<" #N#36."<<nb<<" #C."<<nb<<",N#37."<<nb<<" ) ;"<< endl <<"      CONSTRAINTS :"<< endl <<"    } "<< endl <<""<< endl;
    ofs << setprecision(1) << fixed << "  OR Zinc-Binder"<<nb<<" ( sulfur"<<nb<<" CONO"<<nb<<" carboxyl"<<nb<<" NNC(OS)"<<nb<<" phospho-group"<<nb<<" NNarom"<<nb<<" )"<<endl
        <<"  {" << endl
       <<"    ATTRIBUTES :"<<endl
      <<"      COORD2D ( 0 0 ) ;"<<endl
     <<"      COORD3D ( -701.6 272.3 0.0 ) ;"<<endl
    <<"    OBJECTS :"<<endl
    <<"      POINT point.1"<<nb<<endl
    <<"    COORD2D ( 0 0 )"<<endl
    <<"    COORD3D ( -701.6 272.3 0.0 )"<<endl
    <<"    MAP ( sulfur"<<nb<<"/#S#3."<<nb<<" CONO"<<nb<<"/point-1"<<nb<<" carboxyl"<<nb<<"/point-2"<<nb<<" NNC(OS)"<<nb<<"/point-3"<<nb<<" phospho-group"<<nb<<"/#P#28."<<nb<<" NNarom"<<nb<<"/point-4"<<nb<<" )"<<endl
    <<"    PROP ( INDEX_CODE \"22\" ) ;"<<endl
    <<"  }" << endl;

}


void PharmWriter::writeCHMss(const PharmList &pharmacophore,string toto,bool boolWeigth)
{



    double maxdist=0;
    double dist;
    // Feature importante, depart des pharmacophores
    Coords Eatm[6];
    for (int i =0;i<4;i++){
        Eatm[i]= Coords(0,0,0);
    }
    int flag=0;
    double div=1;

    int nbfeat =0;
    // Lectre du pharmacophore généré
    // Calcul des plus grande distances
    for (ItcPharm ItPha = pharmacophore.begin(); ItPha != pharmacophore.end();++ItPha ){
        PharmProp &Pharm=**ItPha;
        if ( !Pharm.prop()==PharType::EXCLUSION){
            nbfeat++;
            continue;}
        for (ItcPharm ItPha2 = ItPha+1; ItPha2 != pharmacophore.end();++ItPha2 ){
            PharmProp &Pharm2=**ItPha2;
            //            if ( !Pharm2.prop()==PharType::EXCLUSION){continue;}
            dist = Pharm.getBase().calcDist(Pharm2.getBase());
            //            cout <<" Distance  " << dist<<endl;
            if ( dist > maxdist){
                maxdist = dist;
                if (Eatm[0] != Coords(0,0,0)){
                    Coords temp1 = Eatm[0];
                    Coords temp2 = Eatm[1];
                    Eatm[3].setCoords(temp1);
                    Eatm[4].setCoords(temp2);
                }
                Eatm[0].setCoords(Pharm.getBase());
                Eatm[1].setCoords(Pharm2.getBase());
                //                    cout << "maxdist : "<<maxdist<< "  " << Eatm[0].toString() << " 2 different ? " << Eatm[3].toString() << endl;
                //                    cout << "maxdist : "<<maxdist<< "  " << Eatm[1].toString() << " 2 different ? " << Eatm[4].toString() << endl;
            }
        }

    }

    // Vérification du nombre de feature présente
    // Flag corespond au nombre de pharmacophore finalement créés
    //    cout <<nbfeat << endl;
    if (nbfeat<20){
        flag=0;
        div=1.0;
    }else if (nbfeat<30){
        flag=2;
        div=0.6;

    }else if (nbfeat<40){
        flag = 3;
        div = 0.5;
    }else {
        flag = 6;
        div = 0.5;
    }
    double div_temp = 0;
    Eatm[2].setCoords((Eatm[0].x+Eatm[1].x)/2,(Eatm[0].y+Eatm[1].y)/2,(Eatm[0].z+Eatm[1].z)/2);
    Eatm[5].setCoords((Eatm[3].x+Eatm[4].x)/2,(Eatm[3].y+Eatm[4].y)/2,(Eatm[3].z+Eatm[4].z)/2);
    int posi=0;
    for (posi=0 ; posi < flag; posi++){
        if (posi == 2 || posi == 5) {
            div_temp = div/2;
        }else {
            div_temp = div;
        }
        int nbHBD=1, nbHBA=1, nbPOS=1, nbNEG=1,nbARO=1,nbHYD=1, nbMET=1;
        int nbExclusion=1;
        double W[8]={1,1,1,1,1,1,1,1};
        int nbprop[8]={0,0,0,0,0,0,0,0};
//#ifdef ICHEM_DEBUG
//        cout << name.substr(0,4) << "|PHARMA|NBCAV \t" << size_cavs.size() << endl;
//#endif
        ofstream ofs;
        string totox = toto + IntToString(posi) +".chm";
        const char *namefile = totox.c_str();
        ofs.open(namefile,ios::out);
        namefile = namefile;
        time_t t = time(NULL);
        char mbstr[25];
        strftime(mbstr, sizeof(mbstr), "%a %b %d %T %Y", localtime(&t));
        string n(namefile, 0, 9);
        //    strftime(format, 24, "%a %b %d %T %Y", &timeInfo);

        ofs << "!BioCAD CHM format, all rights reserved"<<endl << "!" <<mbstr <<endl << "!VERSION 1.1" << endl << "QUERY "<< n << endl << "{"<<endl;
        ofs.close();
        string name=" ",weight="",blob="",exclu="",donor="",acceptor="",pos="",neg="",aro="",met="";





        for (ItcPharm ItPha = pharmacophore.begin(); ItPha != pharmacophore.end();++ItPha ){
            PharmProp &Pharm=**ItPha;

            if (Eatm[posi].calcDist(Pharm.getBase()) > (maxdist*div_temp) && Pharm.prop() != PharType::EXCLUSION){continue;}

            char buff[100];

            if (Pharm.prop()== PharType::HYDROPHOBIC){
                writehyd(nbHYD,Pharm,namefile);

                string nb=IntToString(nbHYD);
                sprintf(buff, "%.1f", Pharm.getBase().x*100);
                string CooX=buff;
                sprintf(buff, "%.1f", Pharm.getBase().y*100);
                string CooY=buff;
                sprintf(buff, "%.1f", Pharm.getBase().z*100);
                string CooZ=buff;
                sprintf(buff, "%.4f", W[1]);
                string Ws=buff;


                name =name +"HYDROPHOBIC"+nb+" ";
                weight=weight+"    WEIGHT ( HYDROPHOBIC"+nb+" "+Ws+" ) ;\n";
                string tolS =IntToString(160+Pharm.get_tol());
                blob=blob+"      BLOB HYDROPHOBIC"+nb+"-BLOB ( HYDROPHOBIC"+nb+"/Hydrophobic"+nb+" )\n    COORD2D ( 0 0 )\n    COORD3D ( "+CooX+" "+CooY+" "+CooZ+" )\n    TOLERANCE ( "+  tolS +".0000 ) ;\n";


                nbHYD++;
            } else if (Pharm.prop() == PharType::HBACCEPTOR){
                writeHBA(nbHBA,Pharm,namefile);

                string nb=IntToString(nbHBA);
                sprintf(buff, "%.1f", Pharm.getBase().x*100);
                string CooX=buff;
                sprintf(buff, "%.1f", Pharm.getBase().y*100);
                string CooY=buff;
                sprintf(buff, "%.1f", Pharm.getBase().z*100);
                string CooZ=buff;
                sprintf(buff, "%.1f", Pharm.get_3A().x*100);
                string Coo3X=buff;
                sprintf(buff, "%.1f", Pharm.get_3A().y*100);
                string Coo3Y=buff;
                sprintf(buff, "%.1f", Pharm.get_3A().z*100);
                string Coo3Z=buff;
                sprintf(buff, "%.4f", W[2]);
                string Ws=buff;

                name =name +"HB_ACCEPTOR"+nb+" ";
                weight=weight+"    WEIGHT ( HB_ACCEPTOR"+nb+" "+Ws+" ) ;\n";
                acceptor=acceptor+"      BLOB HB_ACCEPTOR"+nb+"_TAIL_BLOB ( HB_ACCEPTOR"+nb+"/HB_ACCEPTOR"+nb+"_TAIL )\n    COORD2D ( 0 0 )\n    COORD3D ( "+CooX+" "+CooY+" "+CooZ+" )\n    TOLERANCE ( "+IntToString(160+Pharm.get_tol())+".0000 ) ;\n";
                acceptor=acceptor+"      BLOB HB_ACCEPTOR"+nb+"_HEAD_BLOB ( HB_ACCEPTOR"+nb+"/HB_ACCEPTOR"+nb+"_HEAD )\n    COORD2D ( 0 0 )\n    COORD3D ( "+Coo3X+" "+Coo3Y+" "+Coo3Z+" )\n    TOLERANCE ( 220.0000 ) ;\n";

                nbHBA++;
            }else if (Pharm.prop()== PharType::HBDONOR){
                writeHBD(nbHBD,Pharm,namefile);

                string nb=IntToString(nbHBD);
                sprintf(buff, "%.1f", Pharm.getBase().x*100);
                string CooX=buff;
                sprintf(buff, "%.1f", Pharm.getBase().y*100);
                string CooY=buff;
                sprintf(buff, "%.1f", Pharm.getBase().z*100);
                string CooZ=buff;
                char buff[100];
                sprintf(buff, "%.1f", Pharm.get_3A().x*100);
                string Coo3X=buff;
                sprintf(buff, "%.1f", Pharm.get_3A().y*100);
                string Coo3Y=buff;
                sprintf(buff, "%.1f", Pharm.get_3A().z*100);
                string Coo3Z=buff;
                sprintf(buff, "%.4f", W[3]);
                string Ws=buff;

                name =name +"HB_DONOR"+nb+" ";
                weight=weight+"    WEIGHT ( HB_DONOR"+nb+" "+Ws+" ) ;\n";
                donor=donor+"      BLOB HB_DONOR"+nb+"_TAIL_BLOB ( HB_DONOR"+nb+"/HB_DONOR"+nb+"_TAIL )\n    COORD2D ( 0 0 )\n    COORD3D ( "+CooX+" "+CooY+" "+CooZ+" )\n    TOLERANCE ( "+IntToString(160+Pharm.get_tol())+".0000 ) ;\n";
                donor=donor+"      BLOB HB_DONOR"+nb+"_HEAD_BLOB ( HB_DONOR"+nb+"/HB_DONOR"+nb+"_HEAD )\n    COORD2D ( 0 0 )\n    COORD3D ( "+Coo3X+" "+Coo3Y+" "+Coo3Z+" )\n    TOLERANCE ( 220.0000 ) ;\n";

                nbHBD++;
            } else if (Pharm.prop()== PharType::AROMATIC){
                writeAro(nbARO,Pharm,namefile);

                string nb=IntToString(nbARO);
                sprintf(buff, "%.1f", Pharm.getBase().x*100);
                string CooX=buff;
                sprintf(buff, "%.1f", Pharm.getBase().y*100);
                string CooY=buff;
                sprintf(buff, "%.1f", Pharm.getBase().z*100);
                string CooZ=buff;
                char buff[100];
                sprintf(buff, "%.1f", Pharm.get_3A().x*100);
                string Coo3X=buff;
                sprintf(buff, "%.1f", Pharm.get_3A().y*100);
                string Coo3Y=buff;
                sprintf(buff, "%.1f", Pharm.get_3A().z*100);
                string Coo3Z=buff;
                sprintf(buff, "%.4f", W[4]);
                string Ws=buff;

                name =name +"RING_AROMATIC"+nb+" ";
                weight=weight+"    WEIGHT ( RING_AROMATIC"+nb+" "+Ws+" ) ;\n";
                aro=aro+"      BLOB RING_AROMATIQUE"+nb+"_HEAL_BLOB ( RING_AROMATIC"+nb+"/centroid"+nb+" )\n    COORD2D ( 0 0 )\n    COORD3D ( "+CooX+" "+CooY+" "+CooZ+" )\n    TOLERANCE ( "+IntToString(160+Pharm.get_tol())+".0000 ) ;\n";
                aro=aro+"      BLOB RING_AROMATIQUE"+nb+"_TAIL_BLOB ( RING_AROMATIC"+nb+"/normal"+nb+" )\n    COORD2D ( 0 0 )\n    COORD3D ( "+Coo3X+" "+Coo3Y+" "+Coo3Z+" )\n    TOLERANCE ( 220.0000 ) ;\n";

                nbARO++;
            } else   if (Pharm.prop()== PharType::IONNEG){
                writeNeg(nbNEG,Pharm,namefile);

                string nb=IntToString(nbNEG);
                sprintf(buff, "%.1f", Pharm.getBase().x*100);
                string CooX=buff;
                sprintf(buff, "%.1f", Pharm.getBase().y*100);
                string CooY=buff;
                sprintf(buff, "%.1f", Pharm.getBase().z*100);
                string CooZ=buff;
                sprintf(buff, "%.4f", W[5]);
                string Ws=buff;

                name =name +"NEG_IONISABLE"+nb+" ";
                weight=weight+"    WEIGHT ( NEG_IONISABLE"+nb+" "+Ws+" ) ;\n";
                neg=neg+"      BLOB NEG_IONISABLE"+nb+"_BLOB ( NEG_IONISABLE"+nb+"/NegIon"+nb+" )\n    COORD2D ( 0 0 )\n    COORD3D ( "+CooX+" "+CooY+" "+CooZ+" )\n    TOLERANCE ( "+IntToString(160+Pharm.get_tol())+".0000 ) ;\n";



                nbNEG++;

            } else   if (Pharm.prop()== PharType::IONPOS){
                writePos(nbPOS,Pharm,namefile);
                string nb=IntToString(nbPOS);
                sprintf(buff, "%.1f", Pharm.getBase().x*100);
                string CooX=buff;
                sprintf(buff, "%.1f", Pharm.getBase().y*100);
                string CooY=buff;
                sprintf(buff, "%.1f", Pharm.getBase().z*100);
                string CooZ=buff;
                sprintf(buff, "%.4f", W[6]);
                string Ws=buff;

                name =name +"POS_IONISABLE"+nb+" ";
                weight=weight+"    WEIGHT ( POS_IONISABLE"+nb+" "+Ws+" ) ;\n";
                pos=pos+"      BLOB POS_IONISABLE"+nb+"_BLOB ( POS_IONISABLE    "+nb+"/PosIon"+nb+" )\n    COORD2D ( 0 0 )\n    COORD3D ( "+CooX+" "+CooY+" "+CooZ+" )\n    TOLERANCE ( "+IntToString(160+Pharm.get_tol())+".0000 ) ;\n";

                nbPOS++;
            } else if (Pharm.prop() == PharType::EXCLUSION){

                string nb=IntToString(nbExclusion);
                sprintf(buff, "%.1f", Pharm.getBase().x*100);
                string CooX=buff;
                sprintf(buff, "%.1f", Pharm.getBase().y*100);
                string CooY=buff;
                sprintf(buff, "%.1f", Pharm.getBase().z*100);
                string CooZ=buff;
                exclu=exclu+"      EXCLUDED_VOLUME ExcludedVolume"+nb+"\n    COORD2D ( 0 0 )\n    COORD3D ( "+CooX+" "+CooY+" "+CooZ+" )\n    TOLERANCE ( "+IntToString(Pharm.get_tol())+".0000 ) ;\n";
                nbExclusion++;
            } else if ( Pharm.prop() == PharType::METAL){
                writemet(nbMET,Pharm,namefile);

                string nb=IntToString(nbMET);
                sprintf(buff, "%.1f", Pharm.getBase().x*100);
                string CooX=buff;
                sprintf(buff, "%.1f", Pharm.getBase().y*100);
                string CooY=buff;
                sprintf(buff, "%.1f", Pharm.getBase().z*100);
                string CooZ=buff;
                sprintf(buff, "%.4f", W[1]);
                string Ws=buff;


                name =name +"Zinc-Binder"+nb+" ";
                weight=weight+"    WEIGHT ( Zinc-Binder"+nb+" "+Ws+" ) ;\n";
                string tolS =IntToString(160+Pharm.get_tol());
                met=met+"      BLOB Metal"+nb+" ( Zinc-Binder"+nb+"/point.1"+nb+" )\n    COORD2D ( 0 0 )\n    COORD3D ( "+CooX+" "+CooY+" "+CooZ+" )\n    TOLERANCE ( "+  tolS +".0000 ) ;\n";

                nbMET++;
            }

        }

        ofs.open(namefile,ios::out|ios::app);
        ofs << "  FIT fit ("<<name<<")"<<endl <<"  {"<<endl
            <<"    ATTRIBUTES :"<<endl
           <<"      COORD2D ( 0 0 ) ;" <<endl << "      COORD3D ( 0.0 0.0 0.0 ) ;"<<endl<<"      ROOT ;"<<endl;
        ofs << weight;
        ofs << "    MINIMUM_FIT ( 0.0000 ) ;"<<endl<<"    CONSTRAINTS :"<<endl;
        ofs << acceptor;
        ofs << donor;
        ofs << pos;
        ofs << neg;
        ofs << aro;
        ofs << blob;
        ofs << met;
        ofs <<exclu<<"  }\n}\n";

        ofs.close();

    }

    //        int n=0, nbfeature=0;
    //        for (n=300 ; n< 800; n++ ){
    //            int nbHBD=1, nbHBA=1, nbPOS=1, nbNEG=1,nbARO=1,nbHYD=1, nbMET=1;
    //            int nbExclusion=1;
    //            double W[8]={0,0,0,0,0,0,0,0};
    //            int nbprop[8]={0,0,0,0,0,0,0,0};
    //            if (boolWeigth){
    //                comptable(pharmacophore,nbprop);
    //                for (int i = 0 ; i < 8 ; i++ ){
    //        //            cout << nbprop[i] << " / " << 100.0*nbprop[i]/nbprop[0]<<"%" <<endl;
    //                    if (100.0*nbprop[i]/nbprop[0] < 10) {
    //                        W[i]=2.0;
    //                    }else if (100.0*nbprop[i]/nbprop[0] > 30) {
    //                       W[i]=1.0;
    //                    }else {
    //                        W[i]=1.0+(30-100.0*nbprop[i]/nbprop[0])/20;
    //                    }
    //                }
    //            }

    //            nbfeature = n/100;
    //    //        cout << /*nbfeature << */endl;
    //            for (int feature = 0 ; feature < nbfeature; feature++){

    //                int toto = rand()%pharmacophore.size();
    //    //            cout << toto << endl;

    //            }
    //        }
}

void PharmWriter::writePML(const PharmList &pharmacophore,const char *namefile,bool boolWeigth)
{

    int nbHBD=1, nbHBA=1, nbPOS=1, nbNEG=1,nbARO=1,nbHYD=1, nbMET=1;
    int nbExclusion=1;
    int nbtot=1 ;
    double W[8]={1,1,1,1,1,1,1,1};
    int nbprop[8]={0,0,0,0,0,0,0,0};
    string tol3S = "2.200";
    string pmlname[20] = {"","","HBA","HBD","","PI","","NI","","","","H","","AR","","","","ZNB","","exclusion" };


//#ifdef ICHEM_DEBUG
//    cout << name.substr(0,4) << "|PHARMA|NBCAV \t" << size_cavs.size() << endl;
//#endif

    if (boolWeigth){
        comptable(pharmacophore,nbprop);
        for (int i = 0 ; i < 8 ; i++ ){
            //            cout << nbprop[i] << " / " << 100.0*nbprop[i]/nbprop[0]<<"%" <<endl;
            if (100.0*nbprop[i]/nbprop[0] < 10) {
                W[i]=1.0;
            }else if (100.0*nbprop[i]/nbprop[0] > 30) {
                W[i]=0.01;
            }else {
                W[i]=(30-100.0*nbprop[i]/nbprop[0])/20;
            }
        }
    }
    ofstream ofs;
    ofs.open(namefile,ios::out);

    time_t t = time(NULL);

    ofs << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl <<"<MolecularEnvironment version =\"1.4\" id=\"moleculeEnvironment0\" name=\"" << namefile << "\">"<<endl;
    ofs << "<pharmacophore name=\"" << namefile << ".1\" id=\"pharma0\" pharmacophoreType=\"LIGAND_SCOUT\">" << endl;


    for (ItcPharm ItPha = pharmacophore.begin(); ItPha != pharmacophore.end();++ItPha ){
        PharmProp &Pharm=**ItPha;
        char buff[100];

        string nb=IntToString(nbtot);
        sprintf(buff, "%.1f", Pharm.getBase().x);
        string CooX=buff;
        sprintf(buff, "%.1f", Pharm.getBase().y);
        string CooY=buff;
        sprintf(buff, "%.1f", Pharm.getBase().z);
        string CooZ=buff;
        sprintf(buff, "%.4f", W[1]);
        string Ws=buff;
        string tolS =DoubleToString(1.5+Pharm.get_tol());
        string tolSs2 =DoubleToString((1.94+Pharm.get_tol())/2);
        sprintf(buff, "%.1f", Pharm.get_tol()*0.01);
        string tolSEx = buff;

        if (Pharm.prop() == PharType::HBACCEPTOR || Pharm.prop() == PharType::HBDONOR){
            char buff[100];
            sprintf(buff, "%.1f", Pharm.get_3A().x);
            string base1X=buff;
            sprintf(buff, "%.1f", Pharm.get_3A().y);
            string base1Y=buff;
            sprintf(buff, "%.1f", Pharm.get_3A().z);
            string base1Z=buff;
            ofs << "<vector name=\""<< pmlname[Pharm.prop()]  << "\" featureId=\""<< nb << "\" pointsToLigand=\"false\" hasSyntheticProjectedPoint=\"false\" optional=\"true\" disabled=\"false\" weight=\"" << Ws << "\" id=\"feature" << nb << "\">" << endl;
            ofs << "<origin x3=\"" << CooX << "\" y3=\"" << CooY << "\" z3=\"" << CooZ << "\" tolerance=\"" << tolSs2 << "\" />" << endl;
            ofs << "<target x3=\"" << base1X << "\" y3=\"" << base1Y << "\" z3=\"" << base1Z << "\" tolerance=\"" << tolS << "\" />" << endl;
            ofs << "</vector>" << endl;

        } else if (Pharm.prop() == PharType::AROMATIC) {
            char buff[100];
            sprintf(buff, "%.1f", Pharm.get_3A().x);
            string base1X=buff;
            sprintf(buff, "%.1f", Pharm.get_3A().y);
            string base1Y=buff;
            sprintf(buff, "%.1f", Pharm.get_3A().z);
            string base1Z=buff;

            ofs << "<plane name=\""<< pmlname[Pharm.prop()]  << "\" featureId=\""<< nb << "\" optional=\"true\" disabled=\"false\" weight=\"" << Ws << "\" id=\"feature" << nb << "\">" << endl;
            ofs << "<position x3=\"" << CooX << "\" y3=\"" << CooY << "\" z3=\"" << CooZ << "\" tolerance=\"" << tolS << "\" />" << endl;
            ofs << "<normal x3=\"" << base1X << "\" y3=\"" << base1Y << "\" z3=\"" << base1Z << "\" tolerance=\"" << tolS << "\" />" << endl;
            ofs << "</plane>" <<endl;
        }
        else if (Pharm.prop() == PharType::EXCLUSION){
            ofs << "<volume name=\""<< pmlname[Pharm.prop()]  << "\" featureId=\""<< nb << "\" optional=\"true\" disabled=\"false\" weight=\"" << Ws << "\" id=\"feature" << nb << "\">" << endl;
            ofs << "<position x3=\"" << CooX << "\" y3=\"" << CooY << "\" z3=\"" << CooZ << "\" tolerance=\""<< tolSEx <<"\" />" << endl;
            ofs << "</volume>" << endl;
        }
        else {
            ofs << "<point name=\""<< pmlname[Pharm.prop()]  << "\" featureId=\""<< nb << "\" optional=\"true\" disabled=\"false\" weight=\"" << Ws << "\" id=\"feature" << nb << "\">" << endl;
            ofs << "<position x3=\"" << CooX << "\" y3=\"" << CooY << "\" z3=\"" << CooZ << "\" tolerance=\"" << tolS << "\" />" << endl;
            ofs << "</point>" << endl;
        }
        nbtot++;
    }

    ofs << "</pharmacophore>";
    ofs << "</MolecularEnvironment>"<< endl;
    ofs.close();


}
void PharmWriter::writeCHM(const PharmList &pharmacophore,const char *namefile,bool boolWeigth)
{
    int nbHBD=1, nbHBA=1, nbPOS=1, nbNEG=1,nbARO=1,nbHYD=1, nbMET=1;
    int nbExclusion=1;
    double W[8]={1,1,1,1,1,1,1,1};
    int nbprop[8]={0,0,0,0,0,0,0,0};
//#ifdef ICHEM_DEBUG
//    cout << name.substr(0,4) << "|PHARMA|NBCAV \t" << size_cavs.size() << endl;
//#endif

    if (boolWeigth){
        comptable(pharmacophore,nbprop);
        for (int i = 0 ; i < 8 ; i++ ){
            //            cout << nbprop[i] << " / " << 100.0*nbprop[i]/nbprop[0]<<"%" <<endl;
            if (100.0*nbprop[i]/nbprop[0] < 10) {
                W[i]=2.0;
            }else if (100.0*nbprop[i]/nbprop[0] > 30) {
                W[i]=1.0;
            }else {
                W[i]=1.0+(30-100.0*nbprop[i]/nbprop[0])/20;
            }
        }
    }
    ofstream ofs;
    ofs.open(namefile,ios::out);

    time_t t = time(NULL);
    char mbstr[25];
    strftime(mbstr, sizeof(mbstr), "%a %b %d %T %Y", localtime(&t));
    string n(namefile, 0, 9);
    //    strftime(format, 24, "%a %b %d %T %Y", &timeInfo);

    ofs << "!BioCAD CHM format, all rights reserved"<<endl << "!" <<mbstr <<endl << "!VERSION 1.1" << endl << "QUERY "<< n << endl << "{"<<endl;
    ofs.close();
    string name=" ",weight="",blob="",exclu="",donor="",acceptor="",pos="",neg="",aro="",met="";
    for (ItcPharm ItPha = pharmacophore.begin(); ItPha != pharmacophore.end();++ItPha ){
        PharmProp &Pharm=**ItPha;
        char buff[100];

        if (Pharm.prop() == PharType::HYDROPHOBIC){
            writehyd(nbHYD,Pharm,namefile);

            string nb=IntToString(nbHYD);
            sprintf(buff, "%.1f", Pharm.getBase().x*100);
            string CooX=buff;
            sprintf(buff, "%.1f", Pharm.getBase().y*100);
            string CooY=buff;
            sprintf(buff, "%.1f", Pharm.getBase().z*100);
            string CooZ=buff;
            sprintf(buff, "%.4f", W[1]);
            string Ws=buff;


            name =name +"HYDROPHOBIC"+nb+" ";
            weight=weight+"    WEIGHT ( HYDROPHOBIC"+nb+" "+Ws+" ) ;\n";
            string tolS =IntToString(160+Pharm.get_tol());
            blob=blob+"      BLOB HYDROPHOBIC"+nb+"-BLOB ( HYDROPHOBIC"+nb+"/Hydrophobic"+nb+" )\n    COORD2D ( 0 0 )\n    COORD3D ( "+CooX+" "+CooY+" "+CooZ+" )\n    TOLERANCE ( "+  tolS +".0000 ) ;\n";


            nbHYD++;
        } else if (Pharm.prop() == PharType::HBACCEPTOR){
            writeHBA(nbHBA,Pharm,namefile);

            string nb=IntToString(nbHBA);
            sprintf(buff, "%.1f", Pharm.getBase().x*100);
            string CooX=buff;
            sprintf(buff, "%.1f", Pharm.getBase().y*100);
            string CooY=buff;
            sprintf(buff, "%.1f", Pharm.getBase().z*100);
            string CooZ=buff;
            sprintf(buff, "%.1f", Pharm.get_3A().x*100);
            string Coo3X=buff;
            sprintf(buff, "%.1f", Pharm.get_3A().y*100);
            string Coo3Y=buff;
            sprintf(buff, "%.1f", Pharm.get_3A().z*100);
            string Coo3Z=buff;
            sprintf(buff, "%.4f", W[2]);
            string Ws=buff;

            name =name +"HB_ACCEPTOR"+nb+" ";
            weight=weight+"    WEIGHT ( HB_ACCEPTOR"+nb+" "+Ws+" ) ;\n";
            acceptor=acceptor+"      BLOB HB_ACCEPTOR"+nb+"_TAIL_BLOB ( HB_ACCEPTOR"+nb+"/HB_ACCEPTOR"+nb+"_TAIL )\n    COORD2D ( 0 0 )\n    COORD3D ( "+CooX+" "+CooY+" "+CooZ+" )\n    TOLERANCE ( "+IntToString(160+Pharm.get_tol())+".0000 ) ;\n";
            acceptor=acceptor+"      BLOB HB_ACCEPTOR"+nb+"_HEAD_BLOB ( HB_ACCEPTOR"+nb+"/HB_ACCEPTOR"+nb+"_HEAD )\n    COORD2D ( 0 0 )\n    COORD3D ( "+Coo3X+" "+Coo3Y+" "+Coo3Z+" )\n    TOLERANCE ( 220.0000 ) ;\n";

            nbHBA++;
        }else if (Pharm.prop() == PharType::HBDONOR){
            writeHBD(nbHBD,Pharm,namefile);

            string nb=IntToString(nbHBD);
            sprintf(buff, "%.1f", Pharm.getBase().x*100);
            string CooX=buff;
            sprintf(buff, "%.1f", Pharm.getBase().y*100);
            string CooY=buff;
            sprintf(buff, "%.1f", Pharm.getBase().z*100);
            string CooZ=buff;
            char buff[100];
            sprintf(buff, "%.1f", Pharm.get_3A().x*100);
            string Coo3X=buff;
            sprintf(buff, "%.1f", Pharm.get_3A().y*100);
            string Coo3Y=buff;
            sprintf(buff, "%.1f", Pharm.get_3A().z*100);
            string Coo3Z=buff;
            sprintf(buff, "%.4f", W[3]);
            string Ws=buff;

            name =name +"HB_DONOR"+nb+" ";
            weight=weight+"    WEIGHT ( HB_DONOR"+nb+" "+Ws+" ) ;\n";
            donor=donor+"      BLOB HB_DONOR"+nb+"_TAIL_BLOB ( HB_DONOR"+nb+"/HB_DONOR"+nb+"_TAIL )\n    COORD2D ( 0 0 )\n    COORD3D ( "+CooX+" "+CooY+" "+CooZ+" )\n    TOLERANCE ( "+IntToString(160+Pharm.get_tol())+".0000 ) ;\n";
            donor=donor+"      BLOB HB_DONOR"+nb+"_HEAD_BLOB ( HB_DONOR"+nb+"/HB_DONOR"+nb+"_HEAD )\n    COORD2D ( 0 0 )\n    COORD3D ( "+Coo3X+" "+Coo3Y+" "+Coo3Z+" )\n    TOLERANCE ( 220.0000 ) ;\n";

            nbHBD++;
        } else if (Pharm.prop() == PharType::AROMATIC){
            writeAro(nbARO,Pharm,namefile);

            string nb=IntToString(nbARO);
            sprintf(buff, "%.1f", Pharm.getBase().x*100);
            string CooX=buff;
            sprintf(buff, "%.1f", Pharm.getBase().y*100);
            string CooY=buff;
            sprintf(buff, "%.1f", Pharm.getBase().z*100);
            string CooZ=buff;
            char buff[100];
            sprintf(buff, "%.1f", Pharm.get_3A().x*100);
            string Coo3X=buff;
            sprintf(buff, "%.1f", Pharm.get_3A().y*100);
            string Coo3Y=buff;
            sprintf(buff, "%.1f", Pharm.get_3A().z*100);
            string Coo3Z=buff;
            sprintf(buff, "%.4f", W[4]);
            string Ws=buff;

            name =name +"RING_AROMATIC"+nb+" ";
            weight=weight+"    WEIGHT ( RING_AROMATIC"+nb+" "+Ws+" ) ;\n";
            aro=aro+"      BLOB RING_AROMATIQUE"+nb+"_HEAL_BLOB ( RING_AROMATIC"+nb+"/centroid"+nb+" )\n    COORD2D ( 0 0 )\n    COORD3D ( "+CooX+" "+CooY+" "+CooZ+" )\n    TOLERANCE ( "+IntToString(160+Pharm.get_tol())+".0000 ) ;\n";
            aro=aro+"      BLOB RING_AROMATIQUE"+nb+"_TAIL_BLOB ( RING_AROMATIC"+nb+"/normal"+nb+" )\n    COORD2D ( 0 0 )\n    COORD3D ( "+Coo3X+" "+Coo3Y+" "+Coo3Z+" )\n    TOLERANCE ( 220.0000 ) ;\n";

            nbARO++;
        } else   if (Pharm.prop() == PharType::IONNEG){
            writeNeg(nbNEG,Pharm,namefile);

            string nb=IntToString(nbNEG);
            sprintf(buff, "%.1f", Pharm.getBase().x*100);
            string CooX=buff;
            sprintf(buff, "%.1f", Pharm.getBase().y*100);
            string CooY=buff;
            sprintf(buff, "%.1f", Pharm.getBase().z*100);
            string CooZ=buff;
            sprintf(buff, "%.4f", W[5]);
            string Ws=buff;

            name =name +"NEG_IONISABLE"+nb+" ";
            weight=weight+"    WEIGHT ( NEG_IONISABLE"+nb+" "+Ws+" ) ;\n";
            neg=neg+"      BLOB NEG_IONISABLE"+nb+"_BLOB ( NEG_IONISABLE"+nb+"/NegIon"+nb+" )\n    COORD2D ( 0 0 )\n    COORD3D ( "+CooX+" "+CooY+" "+CooZ+" )\n    TOLERANCE ( "+IntToString(160+Pharm.get_tol())+".0000 ) ;\n";



            nbNEG++;

        } else   if (Pharm.prop() == PharType::IONPOS){
            writePos(nbPOS,Pharm,namefile);
            string nb=IntToString(nbPOS);
            sprintf(buff, "%.1f", Pharm.getBase().x*100);
            string CooX=buff;
            sprintf(buff, "%.1f", Pharm.getBase().y*100);
            string CooY=buff;
            sprintf(buff, "%.1f", Pharm.getBase().z*100);
            string CooZ=buff;
            sprintf(buff, "%.4f", W[6]);
            string Ws=buff;

            name =name +"POS_IONISABLE"+nb+" ";
            weight=weight+"    WEIGHT ( POS_IONISABLE"+nb+" "+Ws+" ) ;\n";
            pos=pos+"      BLOB POS_IONISABLE"+nb+"_BLOB ( POS_IONISABLE    "+nb+"/PosIon"+nb+" )\n    COORD2D ( 0 0 )\n    COORD3D ( "+CooX+" "+CooY+" "+CooZ+" )\n    TOLERANCE ( "+IntToString(160+Pharm.get_tol())+".0000 ) ;\n";

            nbPOS++;
        } else if (Pharm.prop() == PharType::EXCLUSION){

            string nb=IntToString(nbExclusion);
            sprintf(buff, "%.1f", Pharm.getBase().x*100);
            string CooX=buff;
            sprintf(buff, "%.1f", Pharm.getBase().y*100);
            string CooY=buff;
            sprintf(buff, "%.1f", Pharm.getBase().z*100);
            string CooZ=buff;
            exclu=exclu+"      EXCLUDED_VOLUME ExcludedVolume"+nb+"\n    COORD2D ( 0 0 )\n    COORD3D ( "+CooX+" "+CooY+" "+CooZ+" )\n    TOLERANCE ( "+IntToString(Pharm.get_tol())+".0000 ) ;\n";
            nbExclusion++;
        } else if ( Pharm.prop() == PharType::METAL){
            writemet(nbMET,Pharm,namefile);

            string nb=IntToString(nbMET);
            sprintf(buff, "%.1f", Pharm.getBase().x*100);
            string CooX=buff;
            sprintf(buff, "%.1f", Pharm.getBase().y*100);
            string CooY=buff;
            sprintf(buff, "%.1f", Pharm.getBase().z*100);
            string CooZ=buff;
            sprintf(buff, "%.4f", W[1]);
            string Ws=buff;


            name =name +"Zinc-Binder"+nb+" ";
            weight=weight+"    WEIGHT ( Zinc-Binder"+nb+" "+Ws+" ) ;\n";
            string tolS =IntToString(160+Pharm.get_tol());
            met=met+"      BLOB Metal"+nb+" ( Zinc-Binder"+nb+"/point.1"+nb+" )\n    COORD2D ( 0 0 )\n    COORD3D ( "+CooX+" "+CooY+" "+CooZ+" )\n    TOLERANCE ( "+  tolS +".0000 ) ;\n";

            nbMET++;
        }

    }

    ofs.open(namefile,ios::out|ios::app);
    ofs << "  FIT fit ("<<name<<")"<<endl <<"  {"<<endl
        <<"    ATTRIBUTES :"<<endl
       <<"      COORD2D ( 0 0 ) ;" <<endl << "      COORD3D ( 0.0 0.0 0.0 ) ;"<<endl<<"      ROOT ;"<<endl;
    ofs << weight;
    ofs << "    MINIMUM_FIT ( 0.0000 ) ;"<<endl<<"    CONSTRAINTS :"<<endl;
    ofs << acceptor;
    ofs << donor;
    ofs << pos;
    ofs << neg;
    ofs << aro;
    ofs << blob;
    ofs << met;
    ofs <<exclu<<"  }\n}\n";

    ofs.close();
}

//    Calcul des vecteurs :
//    Coords vectorD = atomlist[itPHA]->fixpos - boxliste[itPHA]->fixpos;
//    double norme = sqrt(pow(vectorD.x,2)+pow(vectorD.y,2)+pow(vectorD.z,2));
//    Coords Vector1 = Coords(vectorD.x/norme,vectorD.y/norme,vectorD.z/norme);
//    Coords pos1 = boxliste[itPHA]->fixpos+Vector1;
//    Coords pos2 = pos1+Vector1;
//    Coords pos3 = pos2+Vector1;
