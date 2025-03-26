#include "headers/ICTools/switch.h"
#include "headers/ICCalcs/interaction.h"
#include "headers/ICMole/similarity.h"
using namespace std;
using namespace ICMole;

void IChemSwitch::helpIFP() const
{
    cout << " IFP - Interaction FingerPrint"                                             <<endl
         << "    IFP protein ligand"<<endl
         << "    IFP protein ligand ligand_ref"<<endl
         << endl
         << "   [General options]"                                                       <<endl
         << "      -name   N (LIG)  Name of the fingerprint -Default: Name of the ligand"<<endl
         << "      --polar          Detect and output only polar interactions           "<<endl
         << "      --metal          Detect and output only metal interactions           "<<endl
         << "      --extended       Include within the fingeprint :"                     <<endl
         << "                          |--> Metal/Acceptor interaction"                  <<endl
         << "                          |--> Weak Hydrogen bonds "                        <<endl
         << "                          |--> PI-Cation interactions"                      <<endl
         << "     "<< endl
         << "    [testing options]" << endl
         << "      -D_Hb    N (3.5)  Hbond length                            (Angstroem)"<<endl
         << "      -D_Hyd   N (4.5)  Hydrophobic length                      (Angstroem)"<<endl
         << "      -D_Io    N (4.0)  Ionic length                            (Angstroem)"<<endl
         << "      -D_Me    N (2.8)  Metal/Acceptor length                   (Angstroem)"<<endl
         << "      -D_Ar    N (5.0)  Aromatic interaction length             (Angstroem)"<<endl
         << "      -D_Pic   N (5.0)  Pi cation interaction length            (Angstroem)"      <<endl
         << "      -a_H     N (Pi)   HBond angle                             (rad)"	   <<endl
         << "      -at_H    N (Pi/3) HBond tolerance angle                   (rad)"      <<endl
         << "      -a_ArFF  N (Pi)   Aromatic Face to Face interaction angle (rad)"      <<endl
         << "      -at_ArFF N (Pi/6) Aromatic Face to Face tolerance angle   (rad)"      <<endl
         << "      -a_ArEF  N (Pi/2) Aromatic Edge to Face interaction angle (rad)"      <<endl
         << "      -at_ArEF N (Pi/3) Aromatic Edge to Face tolerance angle   (rad)"      <<endl
         << "      -a_Pic   N (Pi)   Pi cation interaction angle             (rad)"      <<endl
         << "      -at_Pic  N (Pi/6) Pi cation tolerance angle               (rad)"      <<endl
//         << "      --newH   less permissive definitons of hydrophobic bonds "                    <<endl
         << "      --ligD            Print all possible ligand interactions "                    <<endl
         <<endl
        << "###########################################################################"<<endl
        <<endl;
}




void IChemSwitch::IFP() const throw(MoleExcept)
{
    const size_t InputSize = Input_Values.size();
    bool w_ref = false;
    if (InputSize == 3) w_ref=true;
    if (InputSize != 2 && w_ref==false)
        throw MoleExcept(9020101,
                         "IChem::BSACalc",
                         "Number of parameters must be 2 or 3");

    // USER INPUTS :
    const std::string& fProtein = Input_Values.at(0);
    const std::string& fLigand  = Input_Values.at(1);

    string fname="";
    bool polarOnly=false;
    bool metalOnly=false;
    bool extended=false;
    bool sol=true,cof=true,oldh=true, ligd=false, ifp_out=true;
    double dh=3.5,dhy=4.5,di=4.0,dm =2.8,da=4.0,dpi=4.0;
    double aH=M_PI,atH=M_PI/3,aA=M_PI,atA=M_PI/6,aAe=M_PI/2,atAe=M_PI/3,aPi=M_PI,atPi=M_PI/6;
    bool changedist=false;

    if (!Opt_Values.empty())
        for ( std::map<std::string,std::vector<std::string> >::const_iterator ItOption= Opt_Values.begin();
              ItOption!=Opt_Values.end();
              ++ItOption)
        {
            const std::string& key=(*ItOption).first;

            const std::vector<std::string>& val=(*ItOption).second;
            const std::string& value = val.at(0);

            if (key.compare("-name")==0) fname=val.at(0);
            else if (key.compare("--polar")   == 0) polarOnly=true;
            else if (key.compare("--metal")   == 0) metalOnly=true;
            else if (key.compare("--extended")== 0) extended=true;
            else if (key.compare("--solvent") == 0) sol=false;
            else if (key.compare("--cofactor")== 0) cof=false;
            else if (key.compare("--newH")    == 0) oldh=false;
            else if (key.compare("--ligD")    == 0) ligd=true;
            else if (key.compare("--z")       == 0) ifp_out=false;
            else if (key.compare("-D_Hb")     == 0) { changedist=true; dh  = atof(value.c_str());}
            else if (key.compare("-D_Hyd")    == 0) { changedist=true; dhy = atof(value.c_str());}
            else if (key.compare("-D_Io")     == 0) { changedist=true; di  = atof(value.c_str());}
            else if (key.compare("-D_Me")     == 0) { changedist=true; dm  = atof(value.c_str());}
            else if (key.compare("-D_Ar")     == 0) { changedist=true; da  = atof(value.c_str());}
            else if (key.compare("-D_Pic")    == 0) { changedist=true; dpi = atof(value.c_str());}
            else if (key.compare("-a_H")      == 0) { changedist=true; aH  = atof(value.c_str());}
            else if (key.compare("-at_H")     == 0) { changedist=true; atH = atof(value.c_str());}
            else if (key.compare("-a_ArFF")   == 0) { changedist=true; aA  = atof(value.c_str());}
            else if (key.compare("-at_ArFF")  == 0) { changedist=true; atA = atof(value.c_str());}
            else if (key.compare("-a_ArEF")   == 0) { changedist=true; aAe = atof(value.c_str());}
            else if (key.compare("-at_ArEF")  == 0) { changedist=true; atAe= atof(value.c_str());}
            else if (key.compare("-a_Pic")    == 0) { changedist=true; aPi = atof(value.c_str());}
            else if (key.compare("-at_Pic")   == 0) { changedist=true; atPi= atof(value.c_str());}
        }
//changedist=false;
    try{
        // PARAMETERS NEEDED FOR OUR WORK :
        Complex icomplex;
        MoleReader iread;


        if (sol){
            Residu::Rules[MoleType::PROTEIN][ResType::WATER]=MoleType::PROTEIN;

        }
        if (cof){
            Residu::Rules[MoleType::PROTEIN][ResType::COFACTOR] = MoleType::PROTEIN;
        }


        // STEP 1 - LOADING DATA :
        // 1.1 -KiB Loading protein :
        iread.loadNewFile(fProtein); // Just for the format 
        iread.get_format_file(); // Get the integer for Mol2
        iread.loadInComplex(icomplex,MoleType::PROTEIN); // ReadAtoms seems to read the file with no problems
        if (icomplex.getMole(MoleType::PROTEIN)==(Molecule*)NULL) // Icomplex cannot be nullpointer
            throw MoleExcept(9020102, "IChem::BSACalc", "No protein found in "+fProtein);
//        Molecule &protein = *icomplex.getMole(MoleType::PROTEIN);
//        for (ItCAtom itLA = protein.firstAtom(); itLA != protein.lastAtom();++itLA)
//        {
//            Atom &atomP=**itLA;
//        cout <<"PROTEIN ATOM \t"<<atomP.getIdentifier() //<< endl
//                <<"\t MOL2T: " << atomP.getMOL2Type() //<< endl
//               <<"\t props : "<< atomP.props.toString() << endl;

//        }
        // 1.15 Loading ref :
        if (w_ref){
            const std::string& fLigandref  = Input_Values.at(2);
            iread.loadNewFile(fLigandref);
            const size_t nRef=iread.getNumMolecules();
            iread.loadNewFile(fLigand);
            const size_t nLigand=iread.getNumMolecules();
            const size_t nLigandtoto = nRef+nLigand;

            iread.loadNewFile(fLigandref);

            Interactions ints(icomplex);
            if (changedist){
                ints.setDist_H(dh);
                ints.setDist_Hyd(dhy);
                ints.setDist_Ionic(di);
                ints.setDist_Metal(dm);
                ints.setDist_Arom(da);
                ints.setAngl_H(aH);
                ints.setAngl_Tol_H(atH);
                ints.setAngl_AromFF(aA);
                ints.setAngl_Tol_AromFF(atA);
                ints.setAngl_AromEF(aAe);
                ints.setAngl_Tol_AromEF(atAe);
                ints.setAngl_PICation(aPi);
                ints.setAngl_Tol_PICation(atPi);
                ints.setDist_PICation(dpi);
            }

            InterResults *interRes= new InterResults[nLigandtoto];
            std::vector<string> names;
            size_t currLig=0;
            size_t currRef=0;
            while (!iread.isEOF())
            {
                Molecule ligand;

                iread.loadNextMolecule(ligand,MoleType::LIGAND);
                names.push_back(ligand.getName());
                ligand.checkMOL2();
                ligand.ringPerception();
//                Interactions intsosef(icomplex);
//                InterResults osef;
//                intsosef.calcInteractions(*icomplex.getMole(MoleType::PROTEIN),osef,false,oldh,ligd);
                ints.calcInteractions(ligand,interRes[currLig],true,oldh,ligd);
                if      (!polarOnly && !extended && !metalOnly) ints.genIFP(interRes[currLig],0);
                else if ( polarOnly && !extended) ints.genIFP(interRes[currLig],1);
                else if (!polarOnly &&  extended) ints.genIFP(interRes[currLig],2);
                else if (metalOnly) ints.genIFP(interRes[currLig],4);
                else ints.genIFP(interRes[currLig],3);
                currLig++;
            }
            currRef=currLig;
            iread.loadNewFile(fLigand);
            while (!iread.isEOF())
            {
                Molecule ligand;

                iread.loadNextMolecule(ligand,MoleType::LIGAND);
                names.push_back(ligand.getName());
                ligand.checkMOL2();
                ligand.ringPerception();
                //                cout << "trait : " << ligand.getName() << endl;
                Interactions intsosef(icomplex);
                if (changedist){
                    intsosef.setDist_H(dh);
                    intsosef.setDist_Hyd(dhy);
                    intsosef.setDist_Ionic(di);
                    intsosef.setDist_Metal(dm);
                    intsosef.setDist_Arom(da);
                    intsosef.setAngl_H(aH);
                    intsosef.setAngl_Tol_H(atH);
                    intsosef.setAngl_AromFF(aA);
                    intsosef.setAngl_Tol_AromFF(atA);
                    intsosef.setAngl_AromEF(aAe);
                    intsosef.setAngl_Tol_AromEF(atAe);
                    ints.setAngl_PICation(aPi);
                    ints.setAngl_Tol_PICation(atPi);
                    ints.setDist_PICation(dpi);
                }

//                InterResults osef;
//                intsosef.calcInteractions(*icomplex.getMole(MoleType::PROTEIN),osef,false,oldh,ligd);
                ints.calcInteractions(ligand,interRes[currLig],true,oldh,ligd);
                if      (!polarOnly && !extended && !metalOnly) ints.genIFP(interRes[currLig],0);
                else if ( polarOnly && !extended) ints.genIFP(interRes[currLig],1);
                else if (!polarOnly &&  extended) ints.genIFP(interRes[currLig],2);
                else if (metalOnly) ints.genIFP(interRes[currLig],4);
                else ints.genIFP(interRes[currLig],3);

                currLig++;
            }
            if (ifp_out){
                for (size_t i=0;i< currLig;++i)
                {
                    cout << names.at(i)<<"\t"<<interRes[i].IFPString<<"\n"<<names.at(i)<<"\t"<<interRes[i].IFP.toString()<<"\n";
                }
            }
            Similarity sims(false);
            for (size_t i=0;i< currRef;++i)
            {
                sims.setRef(interRes[i].IFP);
                for (size_t j=i;j< currLig;++j)
                {
                    sims.setComp(interRes[j].IFP);
                    cout<< names.at(i)<<"\t"<<names.at(j)<<"\t"<<sims.Tanimoto()<<endl;
                }
            }

            delete[] interRes;
        }else {
            // 1.2 - Loading Ligand :
            iread.loadNewFile(fLigand); // Same thing for format mol2
            const size_t nLigand=iread.getNumMolecules(); // Get number of molecules (Essentially 1 ligand)
            if (nLigand==1)
            {
                iread.loadInComplex(icomplex,MoleType::LIGAND); // Load molecules => LoadNextMolecules which calls readMOL2Atom method

                if (icomplex.getMole(MoleType::LIGAND)==(Molecule*)NULL) // Cannot be empty
                    throw MoleExcept(9020102,
                                     "IChem::BSACalc",
                                     "No ligand found in "+fLigand);

                Molecule &ligand = *icomplex.getMole(MoleType::LIGAND);
                icomplex.genGrid(1.5);
                icomplex.genGrid(4.5);
                Interactions ints(icomplex);
                if (changedist){
                    ints.setDist_H(dh);
                    ints.setDist_Hyd(dhy);
                    ints.setDist_Ionic(di);
                    ints.setDist_Metal(dm);
                    ints.setDist_Arom(da);
                    ints.setAngl_H(aH);
                    ints.setAngl_Tol_H(atH);
                    ints.setAngl_AromFF(aA);
                    ints.setAngl_Tol_AromFF(atA);
                    ints.setAngl_AromEF(aAe);
                    ints.setAngl_Tol_AromEF(atAe);
                    ints.setAngl_PICation(aPi); 
                    ints.setAngl_Tol_PICation(atPi);
                    ints.setDist_PICation(dpi);
                }
                InterResults interRes;
                InterResults osef;
                ligand.checkMOL2();
                // ligand.ringPerception();
                Interactions intsosef(icomplex);
                if (changedist){
                    intsosef.setDist_H(dh);
                    intsosef.setDist_Hyd(dhy);
                    intsosef.setDist_Ionic(di);
                    intsosef.setDist_Metal(dm);
                    intsosef.setDist_Arom(da);
                    intsosef.setAngl_H(aH);
                    intsosef.setAngl_Tol_H(atH);
                    intsosef.setAngl_AromFF(aA);
                    intsosef.setAngl_Tol_AromFF(atA);
                    intsosef.setAngl_AromEF(aAe);
                    intsosef.setAngl_Tol_AromEF(atAe);
                    ints.setAngl_PICation(aPi);
                    ints.setAngl_Tol_PICation(atPi);
                    ints.setDist_PICation(dpi);
                }
//                intsosef.calcInteractions(*icomplex.getMole(MoleType::PROTEIN),osef,false,oldh,ligd);

                ints.calcInteractions(ligand,interRes,true,oldh,ligd);
                if      (!polarOnly && !extended && !metalOnly) ints.genIFP(interRes,0);
                else if ( polarOnly && !extended) ints.genIFP(interRes,1);
                else if (!polarOnly &&  extended) ints.genIFP(interRes,2);
                else if (metalOnly)               ints.genIFP(interRes,4);
                else ints.genIFP(interRes,3);
                cout << ints.toString(interRes)<<endl;
                if (ifp_out) {
                    cout << interRes.IFPString <<
                            "\n" <<
                            interRes.IFP.toString()<< "\n";
                }
            }
            else
            {
                icomplex.genGrid(1.5);
                icomplex.genGrid(4.5);
                Interactions ints(icomplex);
                if (changedist){
                    ints.setDist_H(dh);
                    ints.setDist_Hyd(dhy);
                    ints.setDist_Ionic(di);
                    ints.setDist_Metal(dm);
                    ints.setDist_Arom(da);
                    ints.setAngl_H(aH);
                    ints.setAngl_Tol_H(atH);
                    ints.setAngl_AromFF(aA);
                    ints.setAngl_Tol_AromFF(atA);
                    ints.setAngl_AromEF(aAe);
                    ints.setAngl_Tol_AromEF(atAe);
                    ints.setAngl_PICation(aPi);
                    ints.setAngl_Tol_PICation(atPi);
                    ints.setDist_PICation(dpi);
                }
                InterResults *interRes= new InterResults[nLigand+1];
                std::vector<string> names;
                size_t currLig=0;
                while (!iread.isEOF())
                {
                    Molecule ligand;

                    iread.loadNextMolecule(ligand,MoleType::LIGAND);
                    names.push_back(ligand.getName());
                   ligand.checkMOL2();
                   ligand.ringPerception();
                    Interactions intsosef(icomplex);
                    if (changedist){
                        intsosef.setDist_H(dh);
                        intsosef.setDist_Hyd(dhy);
                        intsosef.setDist_Ionic(di);
                        intsosef.setDist_Metal(dm);
                        intsosef.setDist_Arom(da);
                        intsosef.setAngl_H(aH);
                        intsosef.setAngl_Tol_H(atH);
                        intsosef.setAngl_AromFF(aA);
                        intsosef.setAngl_Tol_AromFF(atA);
                        intsosef.setAngl_AromEF(aAe);
                        intsosef.setAngl_Tol_AromEF(atAe);
                        ints.setAngl_PICation(aPi);
                        ints.setAngl_Tol_PICation(atPi);
                        ints.setDist_PICation(dpi);
                    }
//                    InterResults osef;
//                    intsosef.calcInteractions(*icomplex.getMole(MoleType::PROTEIN),osef,false,oldh,ligd);
                    ints.calcInteractions(ligand,interRes[currLig],true,oldh,ligd);
                    if      (!polarOnly && !extended && !metalOnly) ints.genIFP(interRes[currLig],0);
                    else if ( polarOnly && !extended) ints.genIFP(interRes[currLig],1);
                    else if (!polarOnly &&  extended) ints.genIFP(interRes[currLig],2);
                    else if (metalOnly)               ints.genIFP(interRes[currLig],4);
                    else ints.genIFP(interRes[currLig],3);

                    currLig++;
                }
                if (ifp_out) {
                for (size_t i=0;i< currLig;++i)
                    {
                        cout << names.at(i)<<"\t"<<interRes[i].IFPString<<"\n"<<names.at(i)<<"\t"<<interRes[i].IFP.toString()<<"\n";
                    }
                }
                delete[] interRes;
            }
        }
        /*
    static const int intToPos[4][NB_INTTYPE]=
        {{-1,3,4,5,6, 0,-1, 1, 2,-1,-1,-1,-1},
         {-1,0,1,2,3,-1,-1,-1,-1,-1,-1,-1,-1},
         {-1,3,4,5,6, 0,10, 1, 2, 9, 7, 8,-1},
         {-1,0,1,2,3,-1,7,-1,-1, 6, 4, 5,-1}};
    static const short length[4]={7,4,11,8};
    //     U ,H,H,I,I, H,M , A, A
    //     N ,B,B,O,O, Y,E , R, R
    //     D ,P,L,P,L, D,T , F, E
    //     E , , , , ,  ,  , F, F
    //     F
*/
    }
    catch (MoleExcept &e)
    {
        e.addTrace("IChem::BSACalc");
        throw;
    }



}