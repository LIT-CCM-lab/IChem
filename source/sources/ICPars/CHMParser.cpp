#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include "headers/ICPars/molereader.h"
#include "headers/ICMole/molecule.h"
#include "headers/ICMole/atom.h"
//#define ICHEM_DEBUG Debug pas refait, merci Jeremy
//# ICHEM_DEBUG2
using namespace std;
using namespace ICMole;

vector<string> MoleReader::split(string str, char delimiter) {
    vector<string> internal;
    stringstream ss(str); // Turn the string into a stream.
    string tok;

    while(getline(ss, tok, delimiter)) {
        internal.push_back(tok);
    }

    return internal;
}

void MoleReader::loadNextChm(Molecule &molecule,
                             const unsigned int& typeMole)
{

#ifdef ICHEM_DEBUG
    if (debug)
        cout << "##################################################"<<endl
             << "##################################################"<<endl
             << "#################### NEW  CHM ####################"<<endl
             << "##################################################"<<endl
             << "##################################################"<<endl;
#endif

    molecule.setMoleType(typeMole);

    // CHECKS :
    if (!isOpen || !finput) throw MoleExcept(2010401,
                                             "MoleReader::loeadNextChm",
                                             "No file opened");

    if (finput.eof()) {EOFile=true;return;}


    string ligne, name; // Each line of file goes here

    unsigned int NbrAt   =100000,  // Number of atoms
            NbrFeat =0,       // Number of features

            ////////////////////////////
            ////// READING HEADER //////
            ////////////////////////////
            /// Objective: knowing how many feature we will have to
            /// handle
            size_start=finput.tellg();

    unsigned int atm_resId  =0, atm_Id     =0;
    string       atm_resName="", atm_Name   ="";

    Residu &resSER = molecule.addResidu("X",1,"SER");
    Residu &resALA = molecule.addResidu("X",2,"ALA");
    Residu &resASP = molecule.addResidu("X",3,"ASP");
    Residu &resLYS = molecule.addResidu("X",4,"LYS");
    Residu &resPHE = molecule.addResidu("X",5,"PHE");
    Residu &resGLY = molecule.addResidu("X",6,"GLY");
    while ( !EOFile)
    {
        getLine(ligne);

        if (EOFile) break;
        if (ligne.find("FIT fit", 0)!=string::npos)
        {
            ligne.substr(10);
            string toto = ligne;
            int chara = 0, space =0;
            for (chara;chara < toto.size(); chara++){
                if ( ligne[chara] == ' '  ){
                    space++;
                }
            }
//            cout << "Lecture chm" << endl;
//            cout << space -5 << " features"  << endl;
//            cout << ligne << endl;
            NbrAt = space -5;


            molecule.Atoms.reserve(NbrAt);
            molecule.name="pharmacophore";
        }
        if (ligne.find("BLOB",0)!=string::npos)
        {

            //            if (bestatm == (Atom*)NULL || bestdist >= 4)                        {molen.addAtom("H",bx.fixpos,"DU","H"   , &resCUB).setFormalCharge((double)projs[bx.getId()]); }
            //            else if (bestatm->props.isCation())                                 {molen.addAtom("O",bx.fixpos,"OD1","O.3", &resASP).setFormalCharge((double)projs[bx.getId()]);}
            //            else if (bestatm->props.isAnion())                                  {molen.addAtom("N",bx.fixpos,"NZ","N.2" , &resLYS).setFormalCharge((double)projs[bx.getId()]);}
            //            else if (bestatm->props.isDonor()&& bestatm->props.isAcceptor())    {molen.addAtom("O",bx.fixpos,"OG","O.2" , &resSER).setFormalCharge((double)projs[bx.getId()]);}
            //            else if (bestatm->props.isAcceptor())                               {molen.addAtom("N",bx.fixpos,"N","N.2"  , &resALA).setFormalCharge((double)projs[bx.getId()]);}
            //            else if (bestatm->props.isDonor())                                  {molen.addAtom("O",bx.fixpos,"O","O.2"  , &resALA).setFormalCharge((double)projs[bx.getId()]);}
            //            else if (bestatm->props.isAromatic())                               {molen.addAtom("C",bx.fixpos,"CZ","C.ar", &resPHE).setFormalCharge((double)projs[bx.getId()]);}
            //            else if (bestatm->props.isHydrophobic())                            {molen.addAtom("C",bx.fixpos,"CA","C.3" , &resGLY).setFormalCharge((double)projs[bx.getId()]);}


            if (ligne.find("HB_A",0)!=string::npos) {
                if (ligne.find("heavy",0)!=string::npos)
                {
                    getLine(ligne);
                    getLine(ligne);
//                    cout << ligne << endl;
                    Coords pos;
                    vector<string> sep = split(ligne, ' ');
                    if (atof(sep[2].c_str()) != 0){
                        pos.setCoords(atof(sep[2].c_str())/100,atof(sep[3].c_str())/100,atof(sep[4].c_str())/100);
//                        cout << "Acceptor (" <<atof(sep[2].c_str())/100<< " , "<< atof(sep[3].c_str())/100<<" , "<< atof(sep[4].c_str())/100 << ")"<< endl;
                    } else {
                        pos.setCoords(atof(sep[6].c_str())/100,atof(sep[7].c_str())/100,atof(sep[8].c_str())/100);
//                        cout << "Acceptor (" <<atof(sep[6].c_str())/100<< " , "<< atof(sep[7].c_str())/100<<" , "<< atof(sep[8].c_str())/100 << ")"<< endl;
                    }
                    molecule.addAtom("O",pos,"O","O.2"  , &resALA).setFormalCharge(2);

                }
                if (ligne.find("Proj",0)!=string::npos)
                {
                    getLine(ligne);
                    getLine(ligne);
//                    cout << ligne << endl;
                    Coords pos;
                    vector<string> sep = split(ligne, ' ');
                    if (atof(sep[2].c_str()) != 0){
                        pos.setCoords(atof(sep[2].c_str())/100,atof(sep[3].c_str())/100,atof(sep[4].c_str())/100);
//                        cout << "Acceptor (" <<atof(sep[2].c_str())/100<< " , "<< atof(sep[3].c_str())/100<<" , "<< atof(sep[4].c_str())/100 << ")"<< endl;
                    } else {
                        pos.setCoords(atof(sep[6].c_str())/100,atof(sep[7].c_str())/100,atof(sep[8].c_str())/100);
//                        cout << "Acceptor (" <<atof(sep[6].c_str())/100<< " , "<< atof(sep[7].c_str())/100<<" , "<< atof(sep[8].c_str())/100 << ")"<< endl;
                    }
                    molecule.addAtom("P",pos,"P","P"  , &resALA).setFormalCharge(2);

                }
            }else if (ligne.find("HB_D",0)!=string::npos) {
                if (ligne.find("heavy",0)!=string::npos)
                {
                    getLine(ligne);
                    getLine(ligne);
//                    cout << ligne << endl;

                    Coords pos;
                    vector<string> sep = split(ligne, ' ');

                    if (atof(sep[2].c_str()) != 0){
                        pos.setCoords(atof(sep[2].c_str())/100,atof(sep[3].c_str())/100,atof(sep[4].c_str())/100);
                    } else {
                        pos.setCoords(atof(sep[6].c_str())/100,atof(sep[7].c_str())/100,atof(sep[8].c_str())/100);
                    }

                    molecule.addAtom("N",pos,"N","N.2"  , &resALA).setFormalCharge(3);
//                    cout << "Donnor ("<<pos.x<< " , "<< pos.y <<" , "<< pos.z<< ")"<< endl;
                }
                if (ligne.find("Proj",0)!=string::npos)
                {
                    getLine(ligne);
                    getLine(ligne);
//                    cout << ligne << endl;
                    Coords pos;
                    vector<string> sep = split(ligne, ' ');
                    if (atof(sep[2].c_str()) != 0){
                        pos.setCoords(atof(sep[2].c_str())/100,atof(sep[3].c_str())/100,atof(sep[4].c_str())/100);
//                        cout << "Acceptor (" <<atof(sep[2].c_str())/100<< " , "<< atof(sep[3].c_str())/100<<" , "<< atof(sep[4].c_str())/100 << ")"<< endl;
                    } else {
                        pos.setCoords(atof(sep[6].c_str())/100,atof(sep[7].c_str())/100,atof(sep[8].c_str())/100);
//                        cout << "Acceptor (" <<atof(sep[6].c_str())/100<< " , "<< atof(sep[7].c_str())/100<<" , "<< atof(sep[8].c_str())/100 << ")"<< endl;
                    }
                    molecule.addAtom("P",pos,"P","P"  , &resALA).setFormalCharge(2);

                }
            }else if (ligne.find("POS_",0)!=string::npos) {
                if (ligne.find("PosI",0)!=string::npos)
                {
                    getLine(ligne);
                    getLine(ligne);
//                    cout << ligne << endl;

                    Coords pos;
                    vector<string> sep = split(ligne, ' ');

                    if (atof(sep[2].c_str()) != 0){
                        pos.setCoords(atof(sep[2].c_str())/100,atof(sep[3].c_str())/100,atof(sep[4].c_str())/100);
                    } else {
                        pos.setCoords(atof(sep[6].c_str())/100,atof(sep[7].c_str())/100,atof(sep[8].c_str())/100);
                    }                    molecule.addAtom("N",pos,"NZ","N.2"  , &resLYS).setFormalCharge(5);
//                    cout << "Positif ("<<pos.x<< " , "<< pos.y <<" , "<< pos.z<< ")"<< endl;
                }
            }else if (ligne.find("NEG_",0)!=string::npos) {
                if (ligne.find("NegI",0)!=string::npos)
                {
                    getLine(ligne);
                    getLine(ligne);
//                    cout << ligne << endl;
                    Coords pos;
                    vector<string> sep = split(ligne, ' ');

                    if (atof(sep[2].c_str()) != 0){
                        pos.setCoords(atof(sep[2].c_str())/100,atof(sep[3].c_str())/100,atof(sep[4].c_str())/100);
                    } else {
                        pos.setCoords(atof(sep[6].c_str())/100,atof(sep[7].c_str())/100,atof(sep[8].c_str())/100);
                    }                    molecule.addAtom("O",pos,"OD1","O.3"  , &resASP).setFormalCharge(7);
//                    cout << "Negatif (" <<pos.x<< " , "<< pos.y <<" , "<< pos.z<< ")"<< endl;
                }
            }else if (ligne.find("HYDR",0)!=string::npos) {
                if (ligne.find("Hydro",0)!=string::npos)
                {
                    getLine(ligne);
                    getLine(ligne);
//                    cout << ligne << endl;
                    Coords pos;
                    vector<string> sep = split(ligne, ' ');
                    if (atof(sep[2].c_str()) != 0){
                        pos.setCoords(atof(sep[2].c_str())/100,atof(sep[3].c_str())/100,atof(sep[4].c_str())/100);
                    } else {
                        pos.setCoords(atof(sep[6].c_str())/100,atof(sep[7].c_str())/100,atof(sep[8].c_str())/100);
                    }                    molecule.addAtom("C",pos,"CA","C.3"  , &resGLY).setFormalCharge(11);
//                    cout << "Hydrophobe ("<<pos.x<< " , "<< pos.y <<" , "<< pos.z<< ")"<< endl;
                }
            }else if (ligne.find("RING",0)!=string::npos){
                if (ligne.find("centroid",0)!=string::npos)
                {
                    getLine(ligne);
                    getLine(ligne);
//                    cout << ligne << endl;
                    Coords pos;
                    vector<string> sep = split(ligne, ' ');
                    if (atof(sep[2].c_str()) != 0){
                        pos.setCoords(atof(sep[2].c_str())/100,atof(sep[3].c_str())/100,atof(sep[4].c_str())/100);
                    } else {
                        pos.setCoords(atof(sep[6].c_str())/100,atof(sep[7].c_str())/100,atof(sep[8].c_str())/100);
                    }                    molecule.addAtom("C",pos,"CZ","C.ar"  , &resPHE).setFormalCharge(13);
//                    cout << "Aromatique ("<<pos.x<< " , "<< pos.y <<" , "<< pos.z<< ")"<< endl;
                }
                if (ligne.find("normal",0)!=string::npos)
                {
                    getLine(ligne);
                    getLine(ligne);
//                    cout << ligne << endl;
                    Coords pos;
                    vector<string> sep = split(ligne, ' ');
                    if (atof(sep[2].c_str()) != 0){
                        pos.setCoords(atof(sep[2].c_str())/100,atof(sep[3].c_str())/100,atof(sep[4].c_str())/100);
//                        cout << "Acceptor (" <<atof(sep[2].c_str())/100<< " , "<< atof(sep[3].c_str())/100<<" , "<< atof(sep[4].c_str())/100 << ")"<< endl;
                    } else {
                        pos.setCoords(atof(sep[6].c_str())/100,atof(sep[7].c_str())/100,atof(sep[8].c_str())/100);
//                        cout << "Acceptor (" <<atof(sep[6].c_str())/100<< " , "<< atof(sep[7].c_str())/100<<" , "<< atof(sep[8].c_str())/100 << ")"<< endl;
                    }
                    molecule.addAtom("P",pos,"P","P"  , &resPHE).setFormalCharge(2);

                }
            }
        }
    }

}




