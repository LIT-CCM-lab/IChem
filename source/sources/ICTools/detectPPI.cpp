#include "headers/ICTools/switch.h"
#include "headers/ICCalcs/ppi.h"
#include "headers/ICCalcs/pdbconvert.h"
#include "headers/ICMole/box.h"
using namespace std;
using namespace ICMole;
//using namespace ICInter;

void IChemSwitch::helpdetectPPI() const
{
    cout << "DetectPPI - Interface detection and characterization"                                                <<endl
         << "    detectppi PDB_name name_file                                                    (1)"<<endl
         << endl
         << "   [General options]"<<endl
         << "      -c  C    : The name of all chains that might be in interaction"<<endl
            //<< "      -c1 C    : The name of one chain that MUST be in interaction "<<endl
            //<< "      -c2 C    : To block to only One interaction" << endl
         << "      -n (1)   : Maximum number of treated interfaces (ranked by decreasing size)" << endl
         << endl
         << " May have some bugs on structures with more than two chains" << endl
         << "###########################################################################"<<endl
         << endl;

}


void IChemSwitch::detectPPI() const throw(ICMole::MoleExcept)
{

    // Recuperation des arguments primaires de detectppi
    const unsigned int InputSize = (unsigned int)Input_Values.size();
    //    cout << InputSize << endl;
    //cout << "Argument donné : " <<Input_Values.at(0) << " " << Input_Values.at(1) << " " << endl;
    string name = Input_Values.at(0);
    if (InputSize < 2) {throw MoleExcept(9010502,"IChem::detectPPi","Not enough parameters");}
    bool w_chain = false; bool w2_chain = false;
    bool w_ub = false;
    bool outputall = false;
    int nchain=1;
    const std::string RProt = Input_Values.at(1);


    MoleReader mread;
    Complex cplx;
    mread.loadNewFile(RProt);
    mread.loadInComplex(cplx,MoleType::PROTEIN);

    if (cplx.getNumMolecule()==0) throw MoleExcept(9010502,"IChem::detectPPi","No protein found in "+ RProt);
    //    cout << "Protein file : "<< RProt<<endl;

    // Recuperation des arguments secondaire
    string opt_name="",opt_val="";
    string chain="",chain2="";
    vector<string> vecChain;
    string outsvm="";
    for (std::map<std::string,std::vector<std::string> >::const_iterator  it = Opt_Values.begin(); it != Opt_Values.end(); it++)
    {
        opt_name = (*it).first;
        opt_val  = (*it).second[0];
        //  -rn -cn --values -sim -outInt -match -max -size --all_cliques -score
        if (opt_name.compare("-c")                 == 0) { w_ub=true; int pos=0;vecChain=(*it).second;}
        else if (opt_name.compare("-c1")           == 0) { w_chain=true;chain=opt_val.c_str();}
        else if (opt_name.compare("-c2")           == 0) { w2_chain=true;chain2=opt_val.c_str();}
        else if (opt_name.compare("-svm")          == 0) { outsvm = opt_val.c_str() ;}
        else if (opt_name.compare("-n")            == 0) { nchain = atof(opt_val.c_str()) ;}
        else if (opt_name.compare("--full")        == 0) { outputall = true;}
        else throw MoleExcept(9010505,"IChem::rundetectPPI","Unrecognized option : "+opt_name);
    }



    //% Retrouvé PDB_name et mole dans le reste du toolkit

    //    cout <<endl
    //    << "#####################################Begin" << endl;
    //    for (int boucle=0;boucle < vecChain.size();boucle ++){
    //        cout << "chain donné: " <<  vecChain[boucle] << endl;

    //    }
#ifdef ICHEM_DEBUG
    if (w_chain)    cout << "chaine sele: " << chain << endl;
    if (w2_chain)   cout << "chaine sele: " << chain2 << endl;
    for (size_t i=0; i< vecChain.size(); ++i) cout <<"Chain sele: "<< vecChain.at(i)<<endl;
#endif

    //    string PDB_name;


    PDBConvert PDConv;
    PDConv.loadPDB(RProt);
    PDConv.harmonizeSizeLine();
    PDConv.changeCSEtoCYS();
    PDConv.changeMSEtoMET();
    PDConv.moveHETATMtoend();
    PDConv.selAltAtm();
    PDConv.renumAtms();
    Complex complex;

    PDConv.toComplex(complex);
    Residu::loadRules();
    Molecule::loadRules();

    if (complex.getMole(MoleType::PROTEIN)==(Molecule*)NULL) throw MoleExcept(999999,"No protein found","");

    Molecule &protein=*complex.getMole(MoleType::PROTEIN);
    protein.setResiduTypes();
    protein.cleanUnwanted();

    AtomList failedAtom; ResiduList failedResidu;
    protein.createConnect(failedResidu, failedAtom,true);
    protein.matchTemplate(failedResidu, failedAtom);
//    protein.addHydrogen();


    Complex cp_sepa;
    /** Separation de la molecule complete, chaque chaine devient une molecule à part
     *
     */
    for (ItCChain itChainA = protein.firstChain();itChainA !=protein.lastChain() ;++itChainA)
    {
        //            if ((**itChainA).getName()=="XX")itChainA++;
        const Chain &chain = **itChainA;
        Molecule *Temp= new Molecule(MoleType::PROTEIN);
        ResiduList reslist;
        for (size_t iRes = 0; iRes < chain.getCountRes();++iRes)
        {
            Residu& res = chain.getResidu(iRes);

            reslist.push_back(&res);
        }
        // std::copy(chain.first(),chain.last(),reslist.begin());
        protein.moveResidu(reslist,*Temp);
        Temp->setName(chain.getName());
        //        Temp->checkMOL2();
        //        Temp->ringPerception();
        cp_sepa.addMolecule(Temp);

    }

    cp_sepa.genGrid(1.5);

    cp_sepa.genGrid(4.5);


    PPI *PPinter = new PPI(cp_sepa,name);

    //on repart de zero.
    //

    PPinter->detectPPi(vecChain,nchain,outputall);
    delete PPinter;
}
