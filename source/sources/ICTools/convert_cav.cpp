
#include <stdio.h>
#include "headers/ICTools/switch.h"
//#include "headers/ICCalcs/volsite.h"
#include "headers/ICPars/pharmwriter.h"

using namespace std;
using namespace ICMole;
void IChemSwitch::helpconvertCav() const
{
    cout << "LIG2CAV - ligand to cavity property"                                  <<endl
         << "    cavLig input output backup                                               "<<endl
         << endl
         << "###########################################################################"<<endl
         <<endl;

}

void IChemSwitch::convertCav() const throw(ICMole::MoleExcept)
{
    MoleReader mread;
    Complex cplx;
    bool convertpharma = false;

    const unsigned int InputSize = (unsigned int)Input_Values.size();
    if ( InputSize == 0 ) { helpconvertCav(); return;}

    if ( InputSize == 1) {
        convertpharma = true;
    }
    else if ( InputSize < 3 ) { helpconvertCav();throw MoleExcept(9010501,"IChem::runlig_cav","Not enough parameters");}


    if (convertpharma){
        const std::string RLig = Input_Values.at(0);
        mread.loadNewFile(RLig);
        mread.loadInComplex(cplx,MoleType::LIGAND);
        if (cplx.getNumMolecule(MoleType::LIGAND)==0) { helpconvertCav(); throw MoleExcept(9010503,"IChem::runlig_cav","No pharmacophore found in "+ RLig);}
        MoleWriter lipha;

//        for (ItCAtom itAtm= listAtoms.begin(); itAtm != listAtoms.end(); ++itAtm)
//        {
//            Atom &atm = **itAtm;


//            PharmProp *newPharm = new PharmProp(PharType::AROMATIC,boxliste[itPHA]->fixpos,pos1,pos3,0);
//            MonPharma.push_back(newPharm);
//            PharmProp* newPharm = new PharmProp(PharType::METAL,boxliste[itPHA]->fixpos,pos1,pos3,0);
//            MonPharma.push_back(newPharm);
//            PharmProp* newPharm = new PharmProp(PharType::HBDONOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
//            MonPharma.push_back(newPharm);
//            PharmProp* newPharm2 = new PharmProp(PharType::HBACCEPTOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
//            MonPharma.push_back(newPharm2);
//            PharmProp *newPharm = new PharmProp(PharType::HYDROPHOBIC,boxliste[itPHA]->fixpos,pos1,pos3,0);
//            MonPharma.push_back(newPharm);

//        }
//        // Ecris le pharmacophore en chm
//        ossname.str("");
//        if (name.length() != 0) ossname << name << "_Pharmacophore";
//        else ossname << "Pharmacophore";
//        ossname<<".chm";
//        PharmWriter pw(ossname.str(),FileFormat::CHM);
//        string name_str = ossname.str();
//        pw.writeCHM(pharmaco,name_str.c_str(),Bweight); //% Recupéré Bweightcelui de volsit
//        ossname.str("");
//        ossname << "Pharmacophore.pml";
//        string name_pml = ossname.str();
//        pw.writePML(pharmaco,name_pml.c_str(),Bweight);
//        ossname.str("");ossname << "Pharmacophore";
//        //            ossname << "2";   ossname<<".chm";
//        name_str = ossname.str();
//        pw.writeCHMss(pharmaco,name_str,Bweight); //% Recupéré Bweightcelui de volsite


    }
    else {

    const std::string RLig = Input_Values.at(0);
    const std::string OLig = Input_Values.at(1);
    const std::string backup = Input_Values.at(2);
    mread.loadNewFile(RLig);
    mread.loadInComplex(cplx,MoleType::LIGAND);
    if (cplx.getNumMolecule(MoleType::LIGAND)==0) { helpconvertCav(); throw MoleExcept(9010503,"IChem::runlig_cav","No ligand found in "+ RLig);}
    cout << "Ligand file : " << RLig << endl;

    string name="";
    ostringstream ossnamelipha;
    if (name.length() != 0) ossnamelipha << name << "_mod";
    else ossnamelipha << OLig;
    ossnamelipha << ".mol2";




    MoleReader iread;
    iread.loadNewFile((RLig));
    const size_t nLigand=iread.getNumMolecules();
    cout << "num  mol " <<  nLigand << endl;
    MoleWriter test;
    test.newFile(ossnamelipha.str(),FileFormat::MOL2);
    while (!iread.isEOF())
    {
        Molecule ligand;
        iread.loadNextMolecule(ligand,MoleType::LIGAND);
        ligand.checkMOL2();
        ligand.ringPerception();
        test.writePharmaLig(&ligand);
    }


    MoleWriter lipha;

//    lipha.newFile(ossnamelipha.str(),FileFormat::MOL2); lipha.writePharmaLig(cplx.getMole(MoleType::LIGAND));
    lipha.newFile(backup,FileFormat::MOL2); lipha.writeMOL2(cplx.getMole(MoleType::LIGAND));

    }

}


