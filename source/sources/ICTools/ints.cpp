#include "headers/ICTools/switch.h"
#include "headers/ICCalcs/interaction.h"
#include "headers/ICPars/molereader.h"
#include "headers/ICMole/box.h"

using namespace std;
using namespace ICMole;

void IChemSwitch::helpints()      const
{
    cout
            << " INTERACTION GENERATOR : "                                                  <<endl
            << "    ints prot lig outfile"                                                  <<endl
            << "      -type      (CENT) Alter positionning output"                          <<endl
            << "                        Multiple values allowed, separated by space"        <<endl
            << "         PROT             InterPROT positionning "                          <<endl
            << "         LIG              InterLIG positionning "                           <<endl
            << "         CENT             Centered positionning "                           <<endl
            << "         MERG             Merged all 3 above"                               <<endl
            <<endl
           << "      -fgps    (STD)  Fingerprint format"                                   <<endl
           << "         STD              Standard   (1 0 21 0 0 3)"                        <<endl
           << "         SVM              SVM format (1:1 3:21 6:3)"                        <<endl
           << "         CMP              Compressed (1 [1 21 [2 3)"                        <<endl
           << "      --small           Compressed fingerprint "                            <<endl
           << endl
           << "   [General options]"                                                       <<endl
           << "      -name      (prot) Name of molecule in out file"                       <<endl
           << "      -logf             Name of log file "                                  <<endl
           << "      -D_Hb    N (3.5)  Hbond length                            (Angstroem)"<<endl
           << "      -D_Hyd   N (4.5)  Hydrophobic length                      (Angstroem)"<<endl
           << "      -D_Io    N (4.0)  Ionic length                            (Angstroem)"<<endl
           << "      -D_Me    N (2.8)  Metal/Acceptor length                   (Angstroem)"<<endl
           << "      -D_Ar    N (5.0)  Aromatic interaction length             (Angstroem)"<<endl
           << "      -D_Pic   N (5.0)  Pi cation interaction length            (Angstroem)"      <<endl
           << "      -a_H     N (Pi)   HBond angle                             (rad)"	   <<endl
           << "      -at_H    N (Pi/3) HBond tolerance angle                   (rad)"    <<endl
           << "      -a_ArFF  N (Pi)   Aromatic Face to Face interaction angle (rad)"       <<endl
           << "      -at_ArFF N (Pi/6) Aromatic Face to Face tolerance angle   (rad)"    <<endl
           << "      -a_ArEF  N (Pi/2) Aromatic Edge to Face interaction angle (rad)"    <<endl
           << "      -at_ArEF N (Pi/3) Aromatic Edge to Face tolerance angle   (rad)"    <<endl
           << "      -a_Pic   N (Pi)   Pi cation interaction angle             (rad)"      <<endl
           << "      -at_Pic  N (Pi/6) Pi cation tolerance angle               (rad)"      <<endl
           << "      --noMerge         Don't merge hydrophobic interactions"               <<endl
           << "      --newH            Less permissive definition of hydrophobic interaction"      <<endl
           << "      --stdout          Write results to standard output"      <<endl
           << endl
           << "###########################################################################"<<endl
           << endl;
}

void IChemSwitch::ints()      const throw(ICMole::MoleExcept)
{
    const unsigned int InputSize = (unsigned int)Input_Values.size();
    if (InputSize > 4|| InputSize<=1)throw MoleExcept(9010301,"IChem::runInts","Not enough parameters");
    const std::string protein = Input_Values.at(0);
    const std::string ligand  = Input_Values.at(1);
    std::string outfile;

    std::string logfile ="";
    std::string out_name = "";

    if (InputSize == 3) {
    out_name = Input_Values.at(2);
    }
    Fingerprint * fgp = (Fingerprint*)NULL;
    Complex cplx;
    MoleReader mread;
    //Interaction wth cofactor and water

    //% Mettre des arguments
    unsigned int fgp_format=0;
    bool out_prot=false,out_lig=false,out_cent=false;

    bool merge=true, out_merge=false,full=true, InterGen=true, oldh=true;
    bool sol = true, cof = true;
    bool enf =false, descripteur=false, plp_out=false, mono_prop=false, ifp_out=true;
    bool stdout = false;
    bool outints=false, outconv=false;
//    cout << "/!\ WARNING /!\ hydrophobic assignation modified" << endl;
    try{
        cplx.clear();
        if (verbose) cout <<"# READING MOLECULES"<<endl;
        mread.loadNewFile(protein); mread.loadInComplex(cplx,MoleType::PROTEIN);
        mread.loadNewFile(ligand);  mread.loadInComplex(cplx,MoleType::LIGAND);
//        mread.loadNewFile(protein); mread.loadInComplex(cplx,MoleType::NUCLEIC);

        if (cplx.getNumMolecule(MoleType::LIGAND)  ==0 || cplx.getMole(MoleType::LIGAND)->numAtom()==0) throw MoleExcept(9010302,"IChem::runInts","No ligand found in " + ligand);
        if (cplx.getNumMolecule(MoleType::PROTEIN) ==0 || cplx.getMole(MoleType::PROTEIN)->numAtom()==0) throw MoleExcept(9010303,"IChem::runInts","No protein found in "+ protein);

//        Residu::Rules[MoleType::PROTEIN][ResType::WATER]=MoleType::PROTEIN;
//        Residu::Rules[MoleType::PROTEIN][ResType::COFACTOR]=MoleType::PROTEIN;

        cplx.genGrid(1.5);
        cplx.genGrid(4.5);

        Interactions ins(cplx);

        if (verbose) cout <<"# SETTING PARAMETERS"<<endl;
        for (std::map<std::string,vector<std::string> >::const_iterator it = Opt_Values.begin(); it != Opt_Values.end(); it++)
        {
            const std::string &opt_name = (*it).first;
            const vector<string> & opt_val = (*it).second;
            const std::string& value = opt_val.at(0);

            ///% Tout mettre dans des variable et du coup lire lres fichier après, cela permet de mettre l'eau et le sovent en arguments



            //-type -name fgps -format -name -a_H -at_H -a_ArFF -at_ArFF -a_ArEF -at_ArEF --noMerge
            if (opt_name.compare("-d_Hb")          == 0) ins.setdist_H         (atof(value.c_str()));
            else if (opt_name.compare("-d_Hyd")    == 0) ins.setdist_Hyd       (atof(value.c_str()));
            else if (opt_name.compare("-d_Io")     == 0) ins.setdist_Ionic     (atof(value.c_str()));
            else if (opt_name.compare("-d_Me")     == 0) ins.setdist_Metal     (atof(value.c_str()));
            else if (opt_name.compare("-d_Ar")     == 0) ins.setdist_Arom      (atof(value.c_str()));
            else if (opt_name.compare("-D_Hb")     == 0) ins.setDist_H         (atof(value.c_str()));
            else if (opt_name.compare("-D_Hyd")    == 0) ins.setDist_Hyd       (atof(value.c_str()));
            else if (opt_name.compare("-D_Io")     == 0) ins.setDist_Ionic     (atof(value.c_str()));
            else if (opt_name.compare("-D_Me")     == 0) ins.setDist_Metal     (atof(value.c_str()));
            else if (opt_name.compare("-D_Ar")     == 0) ins.setDist_Arom      (atof(value.c_str()));
            else if (opt_name.compare("-a_H")      == 0) ins.setAngl_H         (atof(value.c_str()));
            else if (opt_name.compare("-at_H")     == 0) ins.setAngl_Tol_H     (atof(value.c_str()));
            else if (opt_name.compare("-d_Pic")    == 0) ins.setDist_PICation  (atof(value.c_str()));
            else if (opt_name.compare("-a_Pic")    == 0) ins.setAngl_PICation  (atof(value.c_str()));
            else if (opt_name.compare("-at_Pic")   == 0) ins.setAngl_Tol_PICation(atof(value.c_str()));
            else if (opt_name.compare("-a_ArFF")   == 0) ins.setAngl_AromFF    (atof(value.c_str()));
            else if (opt_name.compare("-at_ArFF")  == 0) ins.setAngl_Tol_AromFF(atof(value.c_str()));
            else if (opt_name.compare("-a_ArEF")   == 0) ins.setAngl_AromEF    (atof(value.c_str()));
            else if (opt_name.compare("-at_ArEF")  == 0) ins.setAngl_Tol_AromEF(atof(value.c_str()));
            else if (opt_name.compare("-name")     == 0) out_name=value;
            else if (opt_name.compare("-logf")     == 0) logfile =value;
            else if (opt_name.compare("--stdout")  == 0) stdout =true;
            else if (opt_name.compare("--newH")    == 0) oldh=false;
            else if (opt_name.compare("--noMerge") == 0) merge=false;
            else if (opt_name.compare("--solvent") == 0) sol=false;
            else if (opt_name.compare("--cofactor")== 0) cof=false;
            else if (opt_name.compare("--enf")     == 0) enf=true;
            else if (opt_name.compare("--small")   == 0) full=false;
            else if (opt_name.compare("--d")       == 0) descripteur=true;
            else if (opt_name.compare("--mono")    == 0) mono_prop=true;
            else if (opt_name.compare("--plp")     == 0) plp_out=true;
            else if (opt_name.compare("-fgps")     == 0)
            {
                InterGen = false;
                if (InputSize == 3){
                outfile = Input_Values.at(2);
                }
                else {
                    outfile = "fgps";

                }
                for(vector<string>::const_iterator itS = opt_val.begin(); itS != opt_val.end(); itS++)
                {
                    if (*itS == "STD") fgp_format+=1;
                    if (*itS == "SVM") fgp_format+=2;
                    if (*itS == "CMP") fgp_format+=4;
                }
            }
            else if (opt_name.compare("-type")     == 0)
            {

                for(vector<string>::const_iterator itS = opt_val.begin(); itS != opt_val.end(); itS++)
                {
                    cout <<  *itS << endl;
                    if (*itS == "PROT") out_prot=true;
                    if (*itS == "LIG")  out_lig=true;
                    if (*itS == "CENT") out_cent=true;
                    if (*itS == "MERG") out_merge=true;
                }
            }
            else throw MoleExcept(9010304,"IChem::runInts","Unrecognized option : "+opt_name);
        }
        if (verbose)
        {
            cout << "#  ||-->          H-Bond Distance threshold : "<< ins.getDist_H()         <<endl;
            cout << "#  ||-->     Hydrophobic Distance threshold : "<< ins.getDist_Hyd()       <<endl;
            cout << "#  ||-->           Ionic Distance threshold : "<< ins.getDist_Ionic()     <<endl;
            cout << "#  ||-->           Metal Distance threshold : "<< ins.getDist_Metal()     <<endl;
            cout << "#  ||--> Aromatic center Distance threshold : "<< ins.getDist_Arom()      <<endl;
            cout << "#  ||-->             H-Bond Angle threshold : "<< ins.getAngl_H()         <<endl;
            cout << "#  ||-->             H-Bond Angle tolerance : "<< ins.getAngl_Tol_H()     <<endl;
            cout << "#  ||-->       Face-to-Face Angle threshold : "<< ins.getAngl_AromFF()    <<endl;
            cout << "#  ||-->       Face-to-Face Angle tolerance : "<< ins.getAngl_Tol_AromFF()<<endl;
            cout << "#  ||-->       Edge-to-Face Angle threshold : "<< ins.getAngl_AromEF()    <<endl;
            cout << "#  ||-->       Edge-to-Face Angle tolerance : "<< ins.getAngl_Tol_AromEF()<<endl;
            cout << "# DETECTING INTERACTIONS"<<endl;
        }
        if (cof){
            Residu::Rules[MoleType::PROTEIN][ResType::COFACTOR]=MoleType::PROTEIN;
        }
        if (sol){
            Residu::Rules[MoleType::PROTEIN][ResType::WATER]=MoleType::PROTEIN;
        }

        InterResults results;
        results.listInters.reserve(10000);
        ins.calcInteractions(*cplx.getMole(MoleType::LIGAND),results,merge,oldh);
        //      cplx.getGrid(0).printInFile("gridints",cplx.getGrid(0).getAllBox(),false,true);

        if (plp_out) {

//            ofstream ofplp;
//            ofplp.open("../plp.txt",ios::out|ios::app);
//            Molecule& ligand = *cplx.getMole(MoleType::LIGAND);
//            double plp,plp_atm,plpA,plpB,plpC,plpD,plpE;
//            double dist;
////            Grid grid2 = cplx.genGrid(4.5);
////            grid2.rotateMole(ligand);
//            BoxList boxlist2;

//            for (ItCAtom itLA = ligand.firstAtom(); itLA != ligand.lastAtom();++itLA)
//            {
//                Atom &atomL=**itLA;
//                plp_atm = 0;
//                if (atomL.isHydrogen() || (atomL.getName()=="DuCy") || !atomL.isUsed()) continue;
////                if (atomL.getBox(&grid2) == (Box*)NULL) continue;
//                Box& boxL = *atomL.getBox(&cplx.getGrid(1));
//                boxlist2.clear();
//                cplx.getGrid(1).getAdjacency(boxlist2,&boxL,5.5,true);
//                for (ItBox itBAdj=boxlist2.begin();itBAdj != boxlist2.end();  ++itBAdj)
//                {

//                    for (ItCAtom itAA  = (*itBAdj)->firstAtom();
//                         itAA != (*itBAdj)->lastAtom();
//                         itAA++)
//                    {
//                        Atom &atomP = **itAA;
//                        if (&atomP.getParent()== &ligand || !atomP.isUsed() || atomP.isHydrogen())continue;
//                        dist = atomL.calcFixpos(atomP,5.5+0.1);
//                        plpA=0;plpB=0;plpC=0;plpD=0;plpE=0;
//                        if(      atomL.props.isMetal() && (atomP.props.isAcceptor() || atomP.props.isAnion())
//                           ||   (atomL.props.isAnion() || atomL.props.isAcceptor() ) && atomP.props.isMetal()){
//                            plpA = 1.4; plpB = 2.2; plpC = 3.1; plpD = 3.4;plpE = -5;
//                        }
//                        else if (   ((atomL.props.isAcceptor() || atomL.props.isAnion()) && (atomP.props.isDonor() || atomP.props.isCation()))
//                                    || ((atomL.props.isDonor() || atomL.props.isCation())  && (atomP.props.isAcceptor() || atomP.props.isAnion()))){
//                            plpA = 2.3; plpB = 2.6; plpC = 3.1; plpD = 3.4; plpE = -2;

//                        }
//                        else if (atomL.props.isApolar() && atomP.props.isApolar()){
//                            plpA = 3.4; plpB = 3.6; plpC = 4.5; plpD = 5.5; plpE = -0.4;
//                        }
//                        if (dist < plpA && plpA != 0){
//                            plp = 20*(plpA-dist)/plpA;
//                        }
//                        else if (dist >= plpA && dist < plpB){
//                            plp = plpE*(dist-plpA)/(plpB-plpA);
//                        }
//                        else if (dist >= plpB && dist < plpC){
//                            plp = plpE;
//                        }
//                        else if (dist >= plpC && dist < plpD){
//                            plp = plpE*(plpD-dist)/(plpD-plpC);
//                        }
//                        else if (dist >= plpD){
//                            plp = 0;
//                        }
//                        plp_atm += plp;
//                    }
//                }
//                cout << "plp " << plp_atm << endl;
//                atomL.setPartialCharge(plp_atm);
//                ofplp << atomL.getIdentifier() << "\t" << atomL.props.toString() << "\t" << plp_atm << endl;
//            }
//        ofplp.close();
        MoleWriter mwplp;
        if (stdout) {
            cout << ">>># ligandplp.mol2 #<<<" << endl;
            mwplp.printMOL2(cplx.getMole(MoleType::LIGAND));
        } else {
            mwplp.newFile("ligandplp.mol2");
            mwplp.writeMOL2(cplx.getMole(MoleType::LIGAND));
        }
        cplx.getMole(MoleType::LIGAND)->outplp();
        }

        if (InterGen)
        {
            if (!out_prot && !out_lig && !out_cent && !out_merge)out_merge=true;



            if (!logfile.empty())
            {
                outints = true;

                if (stdout) {
                    cout << ">>># "+logfile+" #<<<" << endl;
                    cout << ins.toString(results) << endl;
                } else {
                    ofstream of;
                    if (fopen(logfile.c_str(),"r")==NULL)
                    {
                        of.open(logfile.c_str(),ios::out|ios::app);
                    }
                    else
                    {
                        of.open(logfile.c_str(),ios::out|ios::app);
                    }
                    if (of.is_open())
                    {
                        of<<ins.toString(results);
                    }
                    of.close();
                }
            }

            if (out_name.empty())out_name=cplx.getMole(MoleType::PROTEIN)->getName();

            results.Ints.setName(out_name);
            MoleWriter mw;
            if (out_lig)
            {
                ins.interToMOL2(results,out_lig,false,false,false);
                if (stdout) {
                    cout << ">>># "+out_name+"_INTS_L.mol2 #<<<" << endl;
                    mw.printMOL2(&results.Ints);
                } else {
                    mw.newFile(out_name+"_INTS_L.mol2");
                    mw.writeMOL2(&results.Ints, outints);
                }

                if (descripteur){
                    ofstream ofs;
                    ofs.open("../../descripteurINTS.dsc",ios::out|ios::app);
                    for (ItCAtom itAtm= results.Ints.firstAtom(); itAtm != results.Ints.lastAtom(); ++itAtm)
                    {
                        Atom &atm = **itAtm;
                        for (ItCAtom itAtm2= itAtm+1; itAtm2 != results.Ints.lastAtom(); ++itAtm2)
                        {
                            Atom &atm2 = **itAtm2;
//                                                cout << atm.getName() << "-" <<atm2.getName() << "   ";
//                            cout << atm.getName() << " avec " << atm2.getName() << " " << atm.getFormalCharge() * atm2.getFormalCharge() << " d: "<< atm.fixpos.calcDist(atm2.fixpos) << endl;
                            ofs  << atm.getMOL2Type() << "-"<<  atm2.getMOL2Type() << " ";
                            ofs  << atm.fixpos.calcDist(atm2.fixpos) << endl;
                        }
                    }
                }

                if (enf){
                    ins.calcenfouiss(out_name+"_INTS_L");
                }
            }
            if (out_prot)
            {
                ins.interToMOL2(results,false,true,false,false);
                if (stdout) {
                    cout << ">>># "+out_name+"_INTS_P.mol2 #<<<" << endl;
                    mw.printMOL2(&results.Ints);
                } else {
                    mw.newFile(out_name+"_INTS_P.mol2");
                    mw.writeMOL2(&results.Ints, outints);
                }
                if (enf){
                    ins.calcenfouiss(out_name+"_INTS_P");
                }

            }
            if (out_cent)
            {
                ins.interToMOL2(results,false,false,true,false);
                if (stdout) {
                    cout << ">>># "+out_name+"_INTS_C.mol2 #<<<" << endl;
                    mw.printMOL2(&results.Ints);
                } else {
                    mw.newFile(out_name+"_INTS_C.mol2");
                    mw.writeMOL2(&results.Ints, outints);
                }
                if (enf){
                    ins.calcenfouiss(out_name+"_INTS_C");
                }
            }
            if (out_merge)
            {
                ins.interToMOL2(results,false,false,false,true);
                if (stdout) {
                    cout << ">>># "+out_name+"_INTS_M.mol2 #<<<" << endl;
                    mw.printMOL2(&results.Ints);
                } else {
                    mw.newFile(out_name+"_INTS_M.mol2");
                    mw.writeMOL2(&results.Ints, outints);
                }
                if (enf){
                    ins.calcenfouiss(out_name+"_INTS_M");
                }
            }
            /*
      if (verbose)cout << "# SAVING INTERACTIONS"<<endl;
      if (int_value == 0 && !out_merge) int_value=4;



      if (out_merge) ins.interMergedToMol2(outfile,out_name);
      if (int_value != 0) ins.interToMOL2File(outfile,int_value,out_name);*/
        }else{
            cout << "# SAVING FINGERPRINT"<<endl;
            if (results.listInters.size() < 3) throw MoleExcept(9010305,"IChemo::runInts","not enough interaction to generate fingerprint");
            if (out_name.empty())out_name=cplx.getMole(MoleType::PROTEIN)->getName();

            fgp = &ins.generateTriplets(results, full);

            if (fgp == (Fingerprint*) NULL) throw MoleExcept(9010306,"IChem::runInts","unable to generate fingerprint");
            if (fgp_format==0)fgp_format = 1;
            ostringstream oss;
            /*if (*itS == "STD") fgp_format+=1;
            if (*itS == "SVM") fgp_format+=2;
            if (*itS == "CMP") fgp_format+=4;*/


            ofstream of;
            ostringstream osFile;
//            cout <<  "here "<<outfile << endl;
//            cout << fgp_format << endl;
//            cout << fgp->toString() << endl;
            switch (fgp_format)
              {
              case 1:// STD

                osFile.str("");osFile<< outfile;         of.open(osFile.str().c_str(),ios::out|ios::app); if (out_name.length()) oss << out_name<<"\t"<< fgp->toStringbin()<<endl;         of << oss.str(); of.close();
                break;
              case 2:// SVM
                osFile.str("");osFile<< outfile;         of.open(osFile.str().c_str(),ios::out|ios::app); if (out_name.length()) oss << out_name<<"\t"<< fgp->toSVMString()<<endl;      of << oss.str(); of.close();
                break;
              case 3:// STD+SVM
                osFile.str("");osFile<< outfile<<"_STD"; of.open(osFile.str().c_str(),ios::out|ios::app); if (out_name.length()) oss << out_name<<"\t"<< fgp->toStringbin()<<endl;         of << oss.str(); of.close();
                osFile.str("");osFile<< outfile<<"_SVM"; of.open(osFile.str().c_str(),ios::out|ios::app); if (out_name.length()) oss << out_name<<"\t"<< fgp->toSVMString()<<endl;      of << oss.str(); of.close();
                break;
              case 4:// CMP
                osFile.str("");osFile<< outfile;         of.open(osFile.str().c_str(),ios::out|ios::app); if (out_name.length()) oss << out_name<<"\t"<< fgp->toCompressString()<<endl; of << oss.str(); of.close();
                break;
              case 5:// CMP+STD
                osFile.str("");osFile<< outfile<<"_STD"; of.open(osFile.str().c_str(),ios::out|ios::app); if (out_name.length()) oss << out_name<<"\t"<< fgp->toStringbin()<<endl;         of << oss.str(); of.close();
                osFile.str("");osFile<< outfile<<"_CMP"; of.open(osFile.str().c_str(),ios::out|ios::app); if (out_name.length()) oss << out_name<<"\t"<< fgp->toCompressString()<<endl; of << oss.str(); of.close();
                break;
              case 6: // SVM+CMP
                osFile.str("");osFile<< outfile<<"_SVM"; of.open(osFile.str().c_str(),ios::out|ios::app); if (out_name.length()) oss << out_name<<"\t"<< fgp->toSVMString()<<endl;      of << oss.str(); of.close();
                osFile.str("");osFile<< outfile<<"_CMP"; of.open(osFile.str().c_str(),ios::out|ios::app); if (out_name.length()) oss << out_name<<"\t"<< fgp->toCompressString()<<endl; of << oss.str(); of.close();
                break;
              case 7:// STD+SVM+CMP
                osFile.str("");osFile<< outfile<<"_STD"; of.open(osFile.str().c_str(),ios::out|ios::app); oss.str("");if (out_name.length()) oss << out_name<<"\t"<< fgp->toStringbin()<<endl;         of << oss.str(); of.close();
                osFile.str("");osFile<< outfile<<"_SVM"; of.open(osFile.str().c_str(),ios::out|ios::app); oss.str("");if (out_name.length()) oss << out_name<<"\t"<< fgp->toSVMString()<<endl;      of << oss.str(); of.close();
                osFile.str("");osFile<< outfile<<"_CMP"; of.open(osFile.str().c_str(),ios::out|ios::app); oss.str("");if (out_name.length()) oss << out_name<<"\t"<< fgp->toCompressString()<<endl; of << oss.str(); of.close();
                break;



              }
            delete fgp;

        }

        if (stdout) {
            cout << ">>># END #<<<" << endl;
        }

    }catch (MoleExcept &e)
    {

        throw;
    }


}
