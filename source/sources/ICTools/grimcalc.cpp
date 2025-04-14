#include "headers/ICCalcs/grim.h"
#include "headers/ICTools/switch.h"
#include "headers/ICMole/similarity.h"

using namespace std;
using namespace ICMole;


void IChemSwitch::helpgrim()      const
{
    cout << " GRIM - GRaph Interaction Matching :"                                       << endl
         << "    grim refProt refLig CompProt CompLig                                (1)"<< endl
         << "    grim refInts compInts                                               (2)"<< endl
         << "    grim refProt refFile dockFile                                       (3)"<< endl

         << endl
         << "   [Note] "<< endl
         << "       (1)  use --multim2 to use multimol2 ligand files"                    << endl
         << "       (3)  refFile & dockFile can be multimol2 files"                      << endl
         << "   [General options]"                                                       << endl
         << "   NOTE : INTERACTION GENERATION General Options also accessible"           << endl
         << "      -rn     N (Ref)   Reference name "                                    << endl
         << "      -cn     N (Comp)  Comparison name "                                   << endl
         << "      --values          Only output score and not alignment "               << endl
         << "      -sim    N (0)     Boolean telling whether the pair is similair or not"<< endl
         << "      -outInt   (MERG)  Output only one kind of interaction positionning "  << endl
         << "         MERG           All aligned interactions are outputed"              << endl
         << "         LIG            InterLIG positionning "                             << endl
         << "         CENT           Centered positionning "                             << endl
         << "         PROT           InterPROT positionning "                            << endl
         << "       NOTE : outInt useless when used with --values "                      << endl
         << "      -match  N (MERG)  Align only with a specific position "               << endl
         << "         MERG           Align with ALL interaction points"                  << endl
         << "         LIG            Align only with ligand interaction points"          << endl
         << "         PROT           Align only with protein interaction points"         << endl
         << "         CENT           Align only with centered interaction points"        << endl
         << "      -score  N (FCT)   Scoring method function"                            << endl
         << "         STD            Scored by decreasing SumCl and increasing RMSD"     << endl
         << "         FCT            Scored with scoring function"                       << endl
         << "      --newH            Less permissive definition of hydrophobe interaction "      <<endl
         << "   [Alignment options]"                                                       << endl
         << "      -max    N (1)     Maximal number of outputed cliques."                << endl
         << "      -size   N (3)     Minimal size of a clique."                          << endl
         << "      --all_cliques     Detect all cliques and not only maximal one"        << endl
         << "      -dsame  N (0.5)   Maximal difference distance between to similar point(LIG-LIG, ...)"        << endl
         << "      -dclose N (0.75)  Maximal difference distance between to close point (LIG-CENT or PROT-CENT)"<< endl
         << "      -dfar   N (1)     Maximal difference distance between to far point (LIG-CENT)"               << endl

         << endl
         << "###########################################################################"<< endl
         << endl;
}

void IChemSwitch::grim()      const throw(ICMole::MoleExcept)
{
    //    verbose = true;
    const size_t InputSize = (unsigned int)Input_Values.size();
    if (InputSize != 2 && InputSize != 4 && InputSize != 3)
    {throw MoleExcept(9050101,
                      "IChem::runGrim",
                      "Not enough parameters");}


    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////// ints OPTIONS ////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    Residu::Rules[MoleType::PROTEIN][ResType::WATER]=MoleType::PROTEIN;
    Residu::Rules[MoleType::PROTEIN][ResType::COFACTOR]=MoleType::PROTEIN;



    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////// TOOL OPTIONS ////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    // Maximal number of outputed cliques
    unsigned int  max_clique  =1;


    // Minimal size of a clique
    size_t  min_size_cl = 3;


    //Grim parameter
    float dif_dsame = 0.5;
    float dif_dclose = 0.75;
    float dif_dfar = 1;
    // With what type of interaction do we match :
    bool match_lig=false, match_prot=false, match_cent=false,match_merg=true;

    // What interaction type shall be output. 0=ALL 1=LIGAND 2=CENTERED 3=PROTEIN
    unsigned int  output_match=0;

    // Do we output values in an external file?
    bool          values      =false;

    // Do we search for all cliques (true) or only the maximal ones ?
    bool          all_cliques =false;

    // Scoring method
    bool          score_opt   =true;

    // multi mol2
    bool          multimol2  =false;
    bool          merge =true;

    // Reference name
    string        RName       ="";
    string        Rnamebu = "";

    // Comparison name
    string        CName       ="";
    string        bestName    ="";

    // Output name
    string        OName       ="GRIM";

    // Known similarity classification
    string        simils      ="";

    bool selected=false;

    // newH
    bool oldh=true;

    // verbose cut
    bool verbose2 = false;

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    /////////////////////////// PROGRAM OPTIONS ////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    InterResults reference;
    InterResults comparison;
    Complex complexR;
    Complex complexC;
    MoleWriter mwrite;
    MoleReader mread;

    ostringstream oss;
    vector<Fingerprint*> VFgp;
    // Similarity sim(false);












    //////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////
    /////////////////////// READING OPTION PARAMETER /////////////////////////////
    //////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////



    for (std::map<std::string,std::vector<std::string> >::const_iterator
         it = Opt_Values.begin();
         it != Opt_Values.end();
         it++)
    {
        const std::string &opt_name = (*it).first;
        const vector<string> & opt_val = (*it).second;
        const std::string& value = opt_val.at(0);
        switch (opt_val.size())
        {
        case 1:

            if      (opt_name.compare("-rn")            == 0) { RName=value;}
            else if (opt_name.compare("-cn")            == 0) { CName=value;}
            else if (opt_name.compare("-on")            == 0) { OName=value;}
            else if (opt_name.compare("--values")       == 0) { values=true;  }
            else if (opt_name.compare("-sim")           == 0) { simils=value;}
            else if (opt_name.compare("-max")           == 0) { max_clique = atoi(value.c_str());}
            else if (opt_name.compare("-size")          == 0) { min_size_cl = atoi(value.c_str());}
            else if (opt_name.compare("--all_cliques")  == 0) { all_cliques=true;}
            else if (opt_name.compare("--multim2")      == 0) { multimol2=true;}
            else if (opt_name.compare("--noMerge")      == 0) { merge=false;}
            else if (opt_name.compare("--newH")         == 0) { oldh=false;}
            else if (opt_name.compare("--v")            == 0) { verbose2=true;}
            else if (opt_name.compare("--vv")           == 0) { verbose=true;}
            else if (opt_name.compare("-dsame")         == 0) { dif_dsame = atof(value.c_str());}
            else if (opt_name.compare("-dclose")        == 0) { dif_dclose = atof(value.c_str());}
            else if (opt_name.compare("-dfar")          == 0) { dif_dfar = atof(value.c_str());}
            else if (opt_name.compare("-score")         == 0) {
                if      (value == "STD") score_opt=false;
                else if (value == "FCT") score_opt=true;
                else throw MoleExcept(9050102,
                                      "IChem::runGrim",
                                      "Unrecognized option for -score : "+value);
            } else if (opt_name.compare("-match")         == 0) {
                for (size_t i=0; i< opt_val.size();++i)
                {
                    if      (opt_val.at(i) == "MERG") match_merg = true;
                    else if (opt_val.at(i) == "LIG")  {match_lig=true;match_merg=false;}
                    else if (opt_val.at(i) == "PROT") {match_prot=true;match_merg=false;}
                    else if (opt_val.at(i) == "CENT") {match_cent=true;match_merg=false;}
                    else throw MoleExcept(9050105,
                                          "IChem::runGrim",
                                          "Unrecognized option for -match : "
                                          +opt_val.at(i));
                }
            }
            else throw MoleExcept(9050103,
                                  "IChem::runGrim",
                                  "Wrong number of values for param  "
                                  +opt_name+". Only one value expected");
            break;
        default:
            if (opt_name.compare("-outInt")        == 0) {
                for (size_t i=0; i< opt_val.size();++i)
                {
                    if      (opt_val.at(i) == "MERG") output_match += 8;
                    else if (opt_val.at(i) == "LIG")  output_match += 1;
                    else if (opt_val.at(i) == "PROT") output_match += 2;
                    else if (opt_val.at(i) == "CENT") output_match += 4;
                    else throw MoleExcept(9050104,
                                          "IChem::runGrim",
                                          "Unrecognized option for -outInt : "
                                          +opt_val.at(i));
                }
            }
            else if (opt_name.compare("-match")         == 0) {
                for (size_t i=0; i< opt_val.size();++i)
                {
                    if      (opt_val.at(i) == "MERG") match_merg = true;
                    else if (opt_val.at(i) == "LIG")  {match_lig=true;match_merg=false;}
                    else if (opt_val.at(i) == "PROT") {match_prot=true;match_merg=false;}
                    else if (opt_val.at(i) == "CENT") {match_cent=true;match_merg=false;}
                    else throw MoleExcept(9050105,
                                          "IChem::runGrim",
                                          "Unrecognized option for -match : "
                                          +opt_val.at(i));
                }
            }
            else throw MoleExcept(9010106,
                                  "IChem::runGrim",
                                  "Unrecognized option : "+opt_name);
        }
    }

    if (!multimol2){
        if (InputSize != 3){
            if (InputSize == 4)   // Case where we have protein/ligand file
            {
                const std::string& RProt = Input_Values.at(0);
                const std::string& RLig  = Input_Values.at(1);
                const std::string& CProt = Input_Values.at(2);
                const std::string& CLig  = Input_Values.at(3);

                mread.loadNewFile(RProt);mread.loadInComplex(complexR,MoleType::PROTEIN);
                mread.loadNewFile(RLig); mread.loadInComplex(complexR,MoleType::LIGAND);
                mread.loadNewFile(CProt);mread.loadInComplex(complexC,MoleType::PROTEIN);
                mread.loadNewFile(CLig); mread.loadInComplex(complexC,MoleType::LIGAND);

                complexR.genGrid(1.5);
                complexR.genGrid(4.5);
                complexC.genGrid(1.5);
                complexC.genGrid(4.5);

                if (complexR.getNumMolecule(MoleType::LIGAND)  ==0)
                    throw MoleExcept(9010107,
                                     "IChem::runGrim",
                                     "No ligand found as reference in : "+ RLig);
                if (complexR.getNumMolecule(MoleType::PROTEIN)  ==0)
                    throw MoleExcept(9010108,
                                     "IChem::runGrim",
                                     "No protein found as reference in "+ RProt);
                if (complexC.getNumMolecule(MoleType::LIGAND)  ==0)
                    throw MoleExcept(9010109,
                                     "IChem::runGrim",
                                     "No ligand found as reference in : "+ CLig);
                if (complexC.getNumMolecule(MoleType::PROTEIN)  ==0)
                    throw MoleExcept(9010110,
                                     "IChem::runGrim",
                                     "No protein found as reference in "+ CProt);


                if (RName.length()==0) RName=complexR.getMole(MoleType::PROTEIN)->getName();
                if (CName.length()==0) CName=complexC.getMole(MoleType::PROTEIN)->getName();

                Interactions interR(complexR);
                interR.calcInteractions(*complexR.getMole(MoleType::LIGAND),reference,merge,oldh);
                Interactions interC(complexC);
                interC.calcInteractions(*complexC.getMole(MoleType::LIGAND),comparison,merge,oldh);

                interR.interToMOL2(reference,match_lig,match_prot,match_cent,match_merg);
                //            interR.interToMOL2(reference,false,false,false,true);
                interC.interToMOL2(comparison,match_lig,match_prot,match_cent,match_merg);
                //            interC.interToMOL2(comparison,false,false,false,true);
                if (verbose){
                    cout << interR.toString(reference)<<endl;
                    cout << interC.toString(comparison)<<endl;
                }


            }
            else if (InputSize==2) // Case where we compare two ints file
            {
                const std::string& RInts = Input_Values.at(0);
                const std::string& CInts  = Input_Values.at(1);
                mread.loadNewFile(RInts); mread.loadNextMolecule(reference.Ints,MoleType::INT_MERG);
                mread.loadNewFile(CInts); mread.loadNextMolecule(comparison.Ints,MoleType::INT_MERG);


            }

            Grim grim(reference,comparison);
            grim.setMinCliqueSize(min_size_cl);
            grim.setDistClose(dif_dclose);
            grim.setDistSame(dif_dsame);
            grim.setdistfar(dif_dfar);
            if (!values)
            {

                mwrite.newFile("RefInts.mol2");
                mwrite.writeMOL2(const_cast<Molecule*>(&grim.getRefInterResult().Ints));
                mwrite.newFile("CompInts.mol2");
                mwrite.writeMOL2(const_cast<Molecule*>(&grim.getCompInterResult().Ints));

            }

            unsigned int interMatch=0;
            if (!match_merg)
            {
                if (match_lig)interMatch+=2;
                if (match_prot)interMatch+=4;
                if (match_cent)interMatch+=1;
                if (interMatch==7)interMatch=0;
            }

            int  nRLig=0,nRCent=0,nRProt=0,nCLig=0,nCCent=0,nCProt=0;
            for (size_t iRes=0; iRes < reference.Ints.numResidus();++iRes)
            {
                const Residu& residu = reference.Ints.getResidu(iRes);
                const std::string& name = residu.getName();
                const std::string pos = name.substr(2,1);
                if (pos == "L")nRLig+=residu.numAtom();
                else if (pos == "P")nRProt+=residu.numAtom();
                else if (pos == "C")nRCent+=residu.numAtom();
            }
            for (size_t iRes=0; iRes < comparison.Ints.numResidus();++iRes)
            {
                const Residu& residu =comparison.Ints.getResidu(iRes);
                const std::string& name =residu.getName();
                const std::string pos=name.substr(2,1);
                if (pos == "L")nCLig+=residu.numAtom();
                else if (pos == "P")nCProt+=residu.numAtom();
                else if (pos == "C")nCCent+=residu.numAtom();
            }







            grim.calcsCliques(interMatch);
            if (verbose){
                cout << "Number of pairs :         " <<grim.getResults().listPairs.size()  <<endl
                     << "Number of cliques found : " <<grim.getResults().cliques.size()<<endl
                     << "Num Int Ref :             " <<reference.Ints.numAtom()<<endl
                     << "Num Int comp:             " <<comparison.Ints.numAtom() <<endl;
            }
            char buff[500];
            ofstream ofsRES,ofsLOG;
            FILE *fp = NULL;
            fp=fopen("Grifp_res.csv","r");
            const bool isOpen= (fp==NULL)?false:true;
            if (fp!=NULL) fclose(fp);
            // If the result file does not exists, it will be created with an header
            ofsRES.open("Grifp_res.csv",ios::app|ios::out);
            if (!ofsRES.is_open()) throw MoleExcept(9010201,"IChem::calcGrim","Cannot open file ");
            if (!isOpen) ofsRES<<"NCli\tRef\tComp\tSimil\tLIG\tCENTER\tPROT\tSumCl\tRMSD\tRLig\tRCent\tRProt\tCLig\tCCent\tCProt\tGrSc\tRMSDAl\tNPol\n";
            // Otherwise, just opened

            if (grim.getResults().cliques.size()==0)
            {
                sprintf(buff,"%d\t%s\t%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%d\n",
                        0, RName.c_str(),CName.c_str(),simils.c_str(),0,0,0, 0.0,0.0,0,0,0,0,0,0,0.0,0.0,0);
                ofsRES<<buff;
            }


            if (!values)ofsLOG.open((OName+".log").c_str(),ios::out|ios::app);

            mwrite.setWithRotPos(true);
            Similarity sim(false);
            ostringstream ofss;
            unsigned int ncli=0;

            for (std::multimap<double, Results>::const_reverse_iterator
                 itClique=grim.firstScore();itClique != grim.lastScore();++itClique)
            {
                ++ncli;
                if(ncli> max_clique)break;
                const Results& result=(*itClique).second;
                Fingerprint *fgp = new Fingerprint(grim.getResults().listPairs.size());
                for (size_t i=0; i < result.clique.size(); i++) fgp->bitOn(result.clique.at(i));
                ostringstream oss;
                selected=true;
                if (VFgp.size())
                {
                    for (vector<Fingerprint*>::iterator itF = VFgp.begin(); itF != VFgp.end(); itF++)
                    {
                        sim.setRef(*fgp);
                        sim.setComp(*(*itF));

                        if (sim.Tanimoto() > 0.9){selected=false;break;}
                    }
                }
                if (!selected) {delete fgp;continue;}

                VFgp.push_back(fgp);
                if (verbose)
                {
                    cout << "##################"<<endl;
                    if (score_opt){
                        cout << " SCORE: " << result.score<<endl;
                    }
                    //cout << "CLIQUE SIZE : "<<result.clique.size()<<endl;
                    cout << " SumCl : " << result.Tc   <<endl;
                    cout << "RMSD  : " << result.rmsd <<endl;
                    cout << "LIG : " << result.NLig <<"\t PROT : " << result.NProt <<"\t CENTER : " << result.NCenter<<endl;
                }


                if (!values)
                {
                    if (InputSize == 4 )
                    { complexC.getMole(MoleType::PROTEIN)->rotateMolecule(result.matrix,result.TransMobil,result.TransRigid);
                        complexC.getMole(MoleType::LIGAND)->rotateMolecule(result.matrix,result.TransMobil,result.TransRigid);
                    }
                    comparison.Ints.rotateMolecule(result.matrix,result.TransMobil,result.TransRigid);
                    ofss.str("");
                    if (max_clique > 1)
                    {
                        ofss<<ncli<<"\t";
                    }

                    if (!values)
                    {
                        ofsLOG<<"#############"<<endl<<"#############"<<endl
                             << ofss.str()<<"SCORE\t"<<result.score<<endl
                             << ofss.str()<<"SumCl\t"<<result.Tc<<endl
                             << ofss.str()<<"RMSD\t"<<result.rmsd<<endl
                             << ofss.str()<<"LIG\t"<<result.NLig<<endl
                             << ofss.str()<<"CENT\t"<<result.NCenter<<endl
                             << ofss.str()<<"PROT\t"<<result.NProt<<endl;
                    }
                    if (max_clique > 1)ofsLOG<<"ID Clique\t";
                    if (!values)
                    {
                        if (InputSize == 4) ofsLOG <<"MATCH\tRef Ints Point\t\tComp Ints Point\tPair\tWeight\tdistance|| Reference Protein \t\t|| Reference Ligand \t\t|| Comparaison Protein \t\t|| Comparaison Ligand "<<endl;
                        else ofsLOG <<"MATCH\tRef Ints Point\t\tComp Ints Point\tPair\tWeight\tdistance "<<endl;
                        for (int i=0; i < (int)result.clique.size(); i++)
                        {
                            const Pair &Pr = grim.getResults().listPairs.at(result.clique.at(i));

                            ofsLOG <<ofss.str()
                                  <<"MATCH\t"<< Pr.ref->getIdentifier()
                                 <<"\t"<<Pr.comp->getIdentifier()
                                <<"\t"<<Pr.weight
                               <<"\t"<< Pr.ref->fixpos.calcDist(Pr.comp->fixpos);
//                            cout << pairRInt << " > " << (nRLig+nRProt) << endl;
//                            if (pairRInt> (nRLig+nRProt)) { pairRInt = pairRInt-(nRLig+nRProt);}
//                            else if (pairRInt> nRLig ) { pairRInt =   pairRInt-nRLig;}
//                            if (pairCInt> (nCLig+nCProt)) { pairCInt = pairCInt-(nCLig+nCProt);}
//                            else if (pairCInt> nCLig ) { pairCInt =   pairCInt-nCLig;}

//                                 //                              << "\t"<< Pr.ref->getNum()/3 << " to " << reference.listInters.size()
//                              << "\t"<< Pr.ref->getResiduName() ;
//                           ofsLOG << "\t"<< pairInt << " to " << reference.listInters.size()
                           if (InputSize == 4)
                           {
                               int pairCInt = Pr.comp->getTNum();
                               int pairRInt = Pr.ref->getTNum();
                               ofsLOG<< "\t|| "<< reference.listInters.at(pairRInt).Prot_Ref->getIdentifier()<< "\t||\t"<< reference.listInters.at(pairRInt).Lig_Ref->getIdentifier()
                            << "\t|| "<< comparison.listInters.at(pairCInt).Prot_Ref->getIdentifier()<< "\t||\t"<< comparison.listInters.at(pairCInt).Lig_Ref->getIdentifier()<<endl;

                           } else
                           {
                               ofsLOG << endl;
                           }
                        }
                    }
                    if (!values){
                        oss.str("");
                        if (InputSize == 4)
                        {

                            oss.str("");oss<<OName;if (max_clique > 1)oss<<ncli;oss <<"_prot.mol2";
                            mwrite.newFile(oss.str());
                            mwrite.writeMOL2(complexC.getMole(MoleType::PROTEIN));
                            oss.str("");oss<<OName;if (max_clique > 1)oss<<ncli;oss <<"_lig.mol2";
                            mwrite.newFile(oss.str());
                            mwrite.writeMOL2(complexC.getMole(MoleType::LIGAND));
                        }
                        oss.str("");oss<<OName;if (max_clique > 1)oss<<ncli;oss <<"_ints.mol2"; mwrite.newFile(oss.str());mwrite.writeMOL2(&comparison.Ints);
                    }


                }

                sprintf(buff,"%d\t%s\t%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%d\n",
                        ncli,
                        RName.c_str(),
                        CName.c_str(),
                        simils.c_str(),
                        result.NLig,
                        result.NCenter,
                        result.NProt,
                        result.Tc,
                        result.rmsd,
                        nRLig,
                        nRCent,
                        nRProt,
                        nCLig,
                        nCCent,
                        nCProt,
                        result.score,
                        result.devAl,
                        result.NPol);
                ofsRES<<buff;

                ncli++;
            }
            for (vector<Fingerprint*>::iterator itF = VFgp.begin(); itF != VFgp.end(); itF++)
            {
                delete *itF;
            }

        }
        else if (InputSize == 3){
            // Read of the reference

            const std::string& RProt = Input_Values.at(0);
            const std::string& RLig  = Input_Values.at(1);

            mread.loadNewFile(RProt);mread.loadInComplex(complexR,MoleType::PROTEIN);
            mread.loadNewFile(RLig);

            if (complexR.getNumMolecule(MoleType::PROTEIN)  == 0)
                throw MoleExcept(9010108,
                                 "IChem::runGrim",
                                 "No protein found as reference in "+ RProt);

            const size_t nLigandR=mread.getNumMolecules();


            if (nLigandR == 1 )
            {
                mread.loadInComplex(complexR,MoleType::LIGAND);

                complexR.genGrid(1.5);
                complexR.genGrid(4.5);


                if (complexR.getNumMolecule(MoleType::LIGAND)  == 0)
                    throw MoleExcept(9010107,
                                     "IChem::runGrim",
                                     "No ligand found as reference in : "+ RLig);



                if (RName.length() == 0) RName = complexR.getMole(MoleType::PROTEIN)->getName();

                Interactions interR(complexR);
                interR.calcInteractions(*complexR.getMole(MoleType::LIGAND),reference,merge,oldh);

                interR.interToMOL2(reference,match_lig,match_prot,match_cent,match_merg);
                if (verbose){
                    cout << interR.toString(reference)<<endl;
                }
                ///////////////////////////////////////////////

                // charge protein in comparison
                const std::string& CLig  = Input_Values.at(2);
                mread.loadNewFile(RProt); mread.loadInComplex(complexC,MoleType::PROTEIN);


                if (RName.length()==0) RName=complexC.getMole(MoleType::PROTEIN)->getName();

                mread.loadNewFile(CLig);
                const size_t nLigand=mread.getNumMolecules();

                if (nLigand == 1 )
                {
                    mread.loadInComplex(complexC,MoleType::LIGAND);

                    if (complexR.getMole(MoleType::LIGAND)==(Molecule*)NULL)
                        throw MoleExcept(9020102,
                                         "IChem::BSACalc",
                                         "No ligand found in "+CLig);

                    complexC.genGrid(1.5);
                    complexC.genGrid(4.5);

                    Interactions interC(complexC);
                    interC.calcInteractions(*complexC.getMole(MoleType::LIGAND),comparison,merge,oldh);
                    interC.interToMOL2(comparison,match_lig,match_prot,match_cent,match_merg);

                    Grim grim(reference,comparison);
                    grim.setMinCliqueSize(min_size_cl);
                    grim.setDistClose(dif_dclose);
                    grim.setDistSame(dif_dsame);
                    grim.setdistfar(dif_dfar);
                    if (!values)
                    {

                        mwrite.newFile("RefInts.mol2");
                        mwrite.writeMOL2(const_cast<Molecule*>(&grim.getRefInterResult().Ints));
                        mwrite.newFile("CompInts.mol2");
                        mwrite.writeMOL2(const_cast<Molecule*>(&grim.getCompInterResult().Ints));

                    }

                    unsigned int interMatch=0;
                    if (!match_merg)
                    {
                        if (match_lig)interMatch+=2;
                        if (match_prot)interMatch+=4;
                        if (match_cent)interMatch+=1;
                        if (interMatch==7)interMatch=0;
                    }


                    int  nRLig=0,nRCent=0,nRProt=0,nCLig=0,nCCent=0,nCProt=0;
                    for (size_t iRes=0; iRes < reference.Ints.numResidus();++iRes)
                    {
                        const Residu& residu = reference.Ints.getResidu(iRes);
                        const std::string& name = residu.getName();
                        const std::string pos = name.substr(2,1);
                        if (pos == "L")nRLig+=residu.numAtom();
                        else if (pos == "P")nRProt+=residu.numAtom();
                        else if (pos == "C")nRCent+=residu.numAtom();
                    }
                    for (size_t iRes=0; iRes < comparison.Ints.numResidus();++iRes)
                    {
                        const Residu& residu =comparison.Ints.getResidu(iRes);
                        const std::string& name =residu.getName();
                        const std::string pos=name.substr(2,1);
                        if (pos == "L")nCLig+=residu.numAtom();
                        else if (pos == "P")nCProt+=residu.numAtom();
                        else if (pos == "C")nCCent+=residu.numAtom();
                    }

                    grim.calcsCliques(interMatch);

                    char buff[500];
                    ofstream ofsRES,ofsLOG;
                    FILE *fp = NULL;
                    fp=fopen("Griimscreen.csv","r");
                    const bool isOpen= (fp==NULL)?false:true;
                    if (fp!=NULL) fclose(fp);
                    // If the result file does not exists, it will be created with an header
                    ofsRES.open("Grifp_res.csv",ios::app|ios::out);
                    if (!ofsRES.is_open()) throw MoleExcept(9010201,"IChem::calcGrim","Cannot open file ");
                    if (!isOpen) ofsRES<<"NCli\tRef\tComp\tSimil\tLIG\tCENTER\tPROT\tSumCl\tRMSD\tRLig\tRCent\tRProt\tCLig\tCCent\tCProt\tGrSc\tRMSDAl\tNPol\n";
                    // Otherwise, just opened

                    if (grim.getResults().cliques.size()==0)
                    {
                        sprintf(buff,"%d\t%s\t%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%d\n",
                                0, RName.c_str(),CName.c_str(),simils.c_str(),0,0,0, 0.0,0.0,0,0,0,0,0,0,0.0,0.0,0);
                        ofsRES<<buff;
                    }


                    if (!values)ofsLOG.open((OName+".log").c_str(),ios::out|ios::app);

                    mwrite.setWithRotPos(true);
                    Similarity sim(false);
                    ostringstream ofss;
                    unsigned int ncli=0;

                    for (std::multimap<double, Results>::const_reverse_iterator
                         itClique=grim.firstScore();itClique != grim.lastScore();++itClique)
                    {
                        ++ncli;
                        if(ncli> max_clique)break;
                        const Results& result=(*itClique).second;
                        Fingerprint *fgp = new Fingerprint(grim.getResults().listPairs.size());
                        for (size_t i=0; i < result.clique.size(); i++) fgp->bitOn(result.clique.at(i));
                        ostringstream oss;
                        selected=true;
                        if (VFgp.size())
                        {
                            for (vector<Fingerprint*>::iterator itF = VFgp.begin(); itF != VFgp.end(); itF++)
                            {
                                sim.setRef(*fgp);
                                sim.setComp(*(*itF));

                                if (sim.Tanimoto() > 0.9){selected=false;break;}
                            }
                        }
                        if (!selected) {delete fgp;continue;}

                        VFgp.push_back(fgp);
                        if (verbose)
                        {
                            cout << "##################"<<endl;
                            if (score_opt){
                                cout << " SCORE: " << result.score<<endl;
                            }
                            //cout << "CLIQUE SIZE : "<<result.clique.size()<<endl;
                            cout << " SumCl : " << result.Tc   <<endl;
                            cout << "RMSD  : " << result.rmsd <<endl;
                            cout << "LIG : " << result.NLig <<"\t PROT : " << result.NProt <<"\t CENTER : " << result.NCenter<<endl;
                        }


                        if (!values)
                        {
                            if (InputSize == 4 )
                            { complexC.getMole(MoleType::PROTEIN)->rotateMolecule(result.matrix,result.TransMobil,result.TransRigid);
                                complexC.getMole(MoleType::LIGAND)->rotateMolecule(result.matrix,result.TransMobil,result.TransRigid);
                            }
                            comparison.Ints.rotateMolecule(result.matrix,result.TransMobil,result.TransRigid);
                            ofss.str("");
                            if (max_clique > 1)
                            {
                                ofss<<ncli<<"\t";
                            }

                            if (!values)
                            {
                                ofsLOG<<"#############"<<endl<<"#############"<<endl
                                     << ofss.str()<<"SCORE\t"<<result.score<<endl
                                     << ofss.str()<<"SumCl\t"<<result.Tc<<endl
                                     << ofss.str()<<"RMSD\t"<<result.rmsd<<endl
                                     << ofss.str()<<"LIG\t"<<result.NLig<<endl
                                     << ofss.str()<<"CENT\t"<<result.NCenter<<endl
                                     << ofss.str()<<"PROT\t"<<result.NProt<<endl;
                            }
                            if (max_clique > 1)ofsLOG<<"ID Clique\t";
                            if (!values)
                            {
                                ofsLOG <<"MATCH\tRef Ints Point\tComp Ints Point\tPair Weight\tPair distance"<<endl;
                                for (int i=0; i < (int)result.clique.size(); i++)
                                {
                                    const Pair &Pr = grim.getResults().listPairs.at(result.clique.at(i));

                                    ofsLOG <<ofss.str()
                                          <<"MATCH\t"<< Pr.ref->getIdentifier()
                                         <<"\t"<<Pr.comp->getIdentifier()
                                        <<"\t"<<Pr.weight
                                       <<"\t"<< Pr.ref->fixpos.calcDist(Pr.comp->fixpos)<<endl;
                                }
                            }
                            if (!values){
                                oss.str("");
                                if (InputSize == 4)
                                {

                                    oss.str("");oss<<OName;if (max_clique > 1)oss<<ncli;oss <<"_prot.mol2";
                                    mwrite.newFile(oss.str());
                                    mwrite.writeMOL2(complexC.getMole(MoleType::PROTEIN));
                                    oss.str("");oss<<OName;if (max_clique > 1)oss<<ncli;oss <<"_lig.mol2";
                                    mwrite.newFile(oss.str());
                                    mwrite.writeMOL2(complexC.getMole(MoleType::LIGAND));
                                }
                                oss.str("");oss<<OName;if (max_clique > 1)oss<<ncli;oss <<"_ints.mol2"; mwrite.newFile(oss.str());mwrite.writeMOL2(&comparison.Ints);
                            }


                        }

                        sprintf(buff,"%d\t%s\t%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%d\n",
                                ncli,
                                RName.c_str(),
                                CName.c_str(),
                                simils.c_str(),
                                result.NLig,
                                result.NCenter,
                                result.NProt,
                                result.Tc,
                                result.rmsd,
                                nRLig,
                                nRCent,
                                nRProt,
                                nCLig,
                                nCCent,
                                nCProt,
                                result.score,
                                result.devAl,
                                result.NPol);
                        ofsRES<<buff;

                        ncli++;
                    }
                    for (vector<Fingerprint*>::iterator itF = VFgp.begin(); itF != VFgp.end(); itF++)
                    {
                        delete *itF;
                    }





                }else {

                    const Results *best_res;
                    double best_score = 0;
                    int best_lig = 0;

                    complexC.genGrid(1.5);
                    complexC.genGrid(4.5);
                    Interactions interC(complexC);
                    InterResults *interRes= new InterResults[nLigand+1];
                    std::vector<string> names;
                    size_t currLig=0;

                    while (!mread.isEOF())
                    {
                        try {
                            Molecule ligand;
                            mread.loadNextMolecule(ligand,MoleType::LIGAND);
                            names.push_back(ligand.getName());
                            CName= ligand.getName();
                            ligand.checkMOL2();
                            ligand.ringPerception();
                            interC.calcInteractions(ligand,interRes[currLig],merge,oldh);
                            interC.interToMOL2(interRes[currLig],match_lig,match_prot,match_cent,match_merg);
                            if (verbose){
                                cout << "comp  :" <<interC.toString(interRes[currLig]) << endl;
                            }
                            Grim grimscr(reference,interRes[currLig]);
                            grimscr.setMinCliqueSize(min_size_cl);

                            unsigned int interMatch=0;
                            if (!match_merg)
                            {
                                if (match_lig)interMatch+=2;
                                if (match_prot)interMatch+=4;
                                if (match_cent)interMatch+=1;
                                if (interMatch==7)interMatch=0;
                            }

                            int nRLig = 0, nRCent = 0, nRProt = 0, nCLig = 0, nCCent = 0, nCProt = 0;
                            for (size_t iRes = 0; iRes < reference.Ints.numResidus(); ++iRes)
                            {
                                const Residu &residu = reference.Ints.getResidu(iRes);
                                const std::string &name = residu.getName();
                                const std::string pos = name.substr(2, 1);
                                if (pos == "L")nRLig += residu.numAtom();
                                else if (pos == "P")nRProt += residu.numAtom();
                                else if (pos == "C")nRCent += residu.numAtom();
                            }
                            for (size_t iRes = 0; iRes < interRes[currLig].Ints.numResidus(); ++iRes)
                            {
                                const Residu &residu = interRes[currLig].Ints.getResidu(iRes);
                                const std::string &name = residu.getName();
                                const std::string pos = name.substr(2, 1);
                                if (pos == "L")nCLig += residu.numAtom();
                                else if (pos == "P")nCProt += residu.numAtom();
                                else if (pos == "C")nCCent += residu.numAtom();
                            }


                            cout << currLig << " : " << ligand.getName() << endl;
                            grimscr.calcsCliques(interMatch);


                            char buff[500];
                            ofstream ofsRES, ofsLOG;
                            FILE *fp = NULL;
                            fp = fopen("Grimscreen_res.tsv", "r");
                            const bool isOpen = (fp == NULL) ? false : true;
                            if (fp != NULL) fclose(fp);
                            // If the result file does not exists, it will be created with an header
                            ofsRES.open("Grimscreen_res.tsv", ios::app | ios::out);
                            if (!ofsRES.is_open()) throw MoleExcept(9010201, "IChem::calcGrimScr", "Cannot open file ");
                            if (!isOpen) ofsRES<< "NCli\tRef\tComp\tSimil\tLIG\tCENTER\tPROT\tSumCl\tRMSD\tRLig\tRCent\tRProt\tCLig\tCCent\tCProt\tGrSc\tRMSDAl\tNPol\n";
                            // Otherwise, just opened

                            if (grimscr.getResults().cliques.size() == 0)
                            {
                                sprintf(buff,"%d\t%s\t%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%d\n",
                                        0, RName.c_str(), CName.c_str(), simils.c_str(), 0, 0, 0, 0.0, 0.0, 0, 0, 0, 0,0, 0, 0.0, 0.0, 0);
                                ofsRES << buff;
                            }

                            mwrite.setWithRotPos(true);
                            Similarity sim(false);
                            ostringstream ofss;
                            unsigned int ncli = 0;

                            for (std::multimap<double, Results>::const_reverse_iterator
                                         itClique = grimscr.firstScore(); itClique != grimscr.lastScore(); ++itClique)
                            {
                                ++ncli;
                                if (ncli > max_clique)break;
                                const Results &result = (*itClique).second;
                                //                    Fingerprint *fgp = new Fingerprint(grimscr.getResults().listPairs.size());
                                //                    for (size_t i=0; i < result.clique.size(); i++) fgp->bitOn(result.clique.at(i));
                                //                    ostringstream oss;
                                //                    selected=true;
                                //                    if (VFgp.size())
                                //                    {
                                //                        for (vector<Fingerprint*>::iterator itF = VFgp.begin(); itF != VFgp.end(); itF++)
                                //                        {
                                //                            sim.setRef(*fgp);
                                //                            sim.setComp(*(*itF));

                                //                            if (sim.Tanimoto() > 0.9){selected=false;break;}
                                //                        }
                                //                    }
                                //                        if (!selected) {delete fgp;continue;}

                                //                    VFgp.push_back(fgp);

                                sprintf(buff,"%d\t%s\t%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%d\n",
                                        ncli,
                                        RName.c_str(),
                                        CName.c_str(),
                                        simils.c_str(),
                                        result.NLig,
                                        result.NCenter,
                                        result.NProt,
                                        result.Tc,
                                        result.rmsd,
                                        nRLig,
                                        nRCent,
                                        nRProt,
                                        nCLig,
                                        nCCent,
                                        nCProt,
                                        result.score,
                                        result.devAl,
                                        result.NPol);
                                ofsRES << buff;
                                if (result.score > best_score)
                                {
                                    if (verbose)
                                        cout << best_score << " < " << result.score << endl;
                                    best_score = result.score;
                                    best_res = &result;
                                    best_lig = currLig;
                                    bestName = CName;
                                }

                                ncli++;


                            }
                            //                for (vector<Fingerprint*>::iterator itF = VFgp.begin(); itF != VFgp.end(); itF++)
                            //                {
                            //                    delete *itF;
                            //                }


                            currLig++;
                        }catch(MoleExcept &e)
                        {
                            cerr << e.getCode()<<"\t"<<e.getSource()<<"\t"<< e.getData()<<endl;
                        }
                    }

//                    TODO: This comment block is to prevent the generation of the GrScreen.csv file which is NOT created with the released version
                    ofstream ofs("GrScreen.csv",ios::app|ios::out);
//
                    if (best_score > 0){
                        Grim grimscr(reference,interRes[best_lig]);
                        grimscr.setMinCliqueSize(min_size_cl);
                        unsigned int interMatch=0;
                        if (!match_merg)
                        {
                            if (match_lig)interMatch+=2;
                            if (match_prot)interMatch+=4;
                            if (match_cent)interMatch+=1;
                            if (interMatch==7)interMatch=0;
                        }
                        grimscr.calcsCliques(interMatch);

                        ofs<< bestName <<"\t"
                           << (*grimscr.firstScore()).second.score <<"\t"
                           << (*grimscr.firstScore()).second.clique.size()<<"\t"
                           << (*grimscr.firstScore()).second.NCenter<<"\t"
                           << (*grimscr.firstScore()).second.NLig <<"\t"
                           << (*grimscr.firstScore()).second.NProt <<"\t"
                           << (*grimscr.firstScore()).second.rmsd <<"\t"
                           << (*grimscr.firstScore()).second.Tc <<"\t"
                           << (*grimscr.firstScore()).second.devAl        <<"\t"
                           << (*grimscr.firstScore()).second.NPol         <<endl;
                    }
                    else
                    {
                        ofs << bestName <<"\t0\t0\t0\t0\t0\t0\t0\t0\t0"<<endl;
                    }
                }
            }
            else {

                complexR.genGrid(1.5);
                complexR.genGrid(4.5);
                if (RName.length() == 0) RName = complexR.getMole(MoleType::PROTEIN)->getName();

                while (!mread.isEOF())
                {
                    try {
                        Molecule ligand;
                        MoleReader mreadC;
                        mread.loadNextMolecule(ligand, MoleType::LIGAND);
                        CName = ligand.getName();
                        ligand.checkMOL2();
                        ligand.ringPerception();
                        Interactions interR(complexR);
                        interR.calcInteractions(ligand, reference, merge, oldh);
                        interR.interToMOL2(reference, match_lig, match_prot, match_cent, match_merg);
                        if (verbose) {
                            cout << interR.toString(reference) << endl;
                        }



                        ///////////////////////////////////////////////

                        // charge protein in comparison
                        const std::string &CLig = Input_Values.at(2);
                        mreadC.loadNewFile(RProt);mreadC.loadInComplex(complexC, MoleType::PROTEIN);


                        if (RName.length() == 0) RName = complexC.getMole(MoleType::PROTEIN)->getName();

                        mreadC.loadNewFile(CLig);
                        const size_t nLigand = mreadC.getNumMolecules();

                        if (nLigand == 1)
                        {
                            mreadC.loadInComplex(complexC, MoleType::LIGAND);

                            if (complexR.getMole(MoleType::LIGAND) == (Molecule *) NULL)
                                throw MoleExcept(9020102,
                                                 "IChem::BSACalc",
                                                 "No ligand found in " + CLig);

                            complexC.genGrid(1.5);
                            complexC.genGrid(4.5);

                            Interactions interC(complexC);
                            interC.calcInteractions(*complexR.getMole(MoleType::LIGAND), comparison, merge, oldh);
                            interC.interToMOL2(comparison, match_lig, match_prot, match_cent, match_merg);

                            Grim grim(reference, comparison);
                            grim.setMinCliqueSize(min_size_cl);
                            grim.setDistClose(dif_dclose);
                            grim.setDistSame(dif_dsame);
                            grim.setdistfar(dif_dfar);
                            if (!values)
                            {

                                mwrite.newFile("RefInts.mol2");
                                mwrite.writeMOL2(const_cast<Molecule *>(&grim.getRefInterResult().Ints));
                                mwrite.newFile("CompInts.mol2");
                                mwrite.writeMOL2(const_cast<Molecule *>(&grim.getCompInterResult().Ints));

                            }

                            unsigned int interMatch = 0;
                            if (!match_merg)
                            {
                                if (match_lig)interMatch += 2;
                                if (match_prot)interMatch += 4;
                                if (match_cent)interMatch += 1;
                                if (interMatch == 7)interMatch = 0;
                            }


                            int nRLig = 0, nRCent = 0, nRProt = 0, nCLig = 0, nCCent = 0, nCProt = 0;
                            for (size_t iRes = 0; iRes < reference.Ints.numResidus(); ++iRes)
                            {
                                const Residu &residu = reference.Ints.getResidu(iRes);
                                const std::string &name = residu.getName();
                                const std::string pos = name.substr(2, 1);
                                if (pos == "L")nRLig += residu.numAtom();
                                else if (pos == "P")nRProt += residu.numAtom();
                                else if (pos == "C")nRCent += residu.numAtom();
                            }
                            for (size_t iRes = 0; iRes < comparison.Ints.numResidus(); ++iRes)
                            {
                                const Residu &residu = comparison.Ints.getResidu(iRes);
                                const std::string &name = residu.getName();
                                const std::string pos = name.substr(2, 1);
                                if (pos == "L")nCLig += residu.numAtom();
                                else if (pos == "P")nCProt += residu.numAtom();
                                else if (pos == "C")nCCent += residu.numAtom();
                            }


                            grim.calcsCliques(interMatch);


                            char buff[500];
                            ofstream ofsRES, ofsLOG;
                            FILE *fp = NULL;
                            fp = fopen("Grimscreen.csv", "r");
                            const bool isOpen = (fp == NULL) ? false : true;
                            if (fp != NULL) fclose(fp);
                            // If the result file does not exists, it will be created with an header
                            ofsRES.open("Grifp_res.csv", ios::app | ios::out);
                            if (!ofsRES.is_open()) throw MoleExcept(9010201, "IChem::calcGrim", "Cannot open file ");
                            if (!isOpen) ofsRES<< "NCli\tRef\tComp\tSimil\tLIG\tCENTER\tPROT\tSumCl\tRMSD\tRLig\tRCent\tRProt\tCLig\tCCent\tCProt\tGrSc\tRMSDAl\tNPol\n";
                            // Otherwise, just opened

                            if (grim.getResults().cliques.size() == 0)
                            {
                                sprintf(buff,"%d\t%s\t%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%d\n",
                                        0, RName.c_str(), CName.c_str(), simils.c_str(), 0, 0, 0, 0.0, 0.0, 0, 0, 0, 0,0, 0, 0.0, 0.0, 0);
                                ofsRES << buff;
                            }


                            if (!values)ofsLOG.open((OName + ".log").c_str(), ios::out | ios::app);

                            mwrite.setWithRotPos(true);
                            Similarity sim(false);
                            ostringstream ofss;
                            unsigned int ncli = 0;

                            for (std::multimap<double, Results>::const_reverse_iterator
                                         itClique = grim.firstScore(); itClique != grim.lastScore(); ++itClique)
                            {
                                ++ncli;
                                if (ncli > max_clique)break;
                                const Results &result = (*itClique).second;
                                Fingerprint *fgp = new Fingerprint(grim.getResults().listPairs.size());
                                for (size_t i = 0; i < result.clique.size(); i++) fgp->bitOn(result.clique.at(i));
                                ostringstream oss;
                                selected = true;
                                if (VFgp.size())
                                {
                                    for (vector<Fingerprint *>::iterator itF = VFgp.begin(); itF != VFgp.end(); itF++)
                                    {
                                        sim.setRef(*fgp);
                                        sim.setComp(*(*itF));

                                        if (sim.Tanimoto() > 0.9) {selected = false;break;}
                                    }
                                }
                                if (!selected) {delete fgp;continue;}

                                VFgp.push_back(fgp);
                                if (verbose)
                                {
                                    cout << "##################" << endl;
                                    if (score_opt) {
                                        cout << " SCORE: " << result.score << endl;
                                    }
                                    //cout << "CLIQUE SIZE : "<<result.clique.size()<<endl;
                                    cout << " SumCl : " << result.Tc << endl;
                                    cout << "RMSD  : " << result.rmsd << endl;
                                    cout << "LIG : " << result.NLig << "\t PROT : " << result.NProt << "\t CENTER : "<< result.NCenter << endl;
                                }


                                if (!values) {
                                    if (InputSize == 4) {
                                        complexC.getMole(MoleType::PROTEIN)->rotateMolecule(result.matrix,
                                                                                            result.TransMobil,
                                                                                            result.TransRigid);
                                        complexC.getMole(MoleType::LIGAND)->rotateMolecule(result.matrix,
                                                                                           result.TransMobil,
                                                                                           result.TransRigid);
                                    }
                                    comparison.Ints.rotateMolecule(result.matrix, result.TransMobil, result.TransRigid);
                                    ofss.str("");
                                    if (max_clique > 1) {
                                        ofss << ncli << "\t";
                                    }

                                    if (!values) {
                                        ofsLOG << "#############" << endl << "#############" << endl
                                               << ofss.str() << "SCORE\t" << result.score << endl
                                               << ofss.str() << "SumCl\t" << result.Tc << endl
                                               << ofss.str() << "RMSD\t" << result.rmsd << endl
                                               << ofss.str() << "LIG\t" << result.NLig << endl
                                               << ofss.str() << "CENT\t" << result.NCenter << endl
                                               << ofss.str() << "PROT\t" << result.NProt << endl;
                                    }
                                    if (max_clique > 1)ofsLOG << "ID Clique\t";
                                    if (!values) {
                                        ofsLOG << "MATCH\tRef Ints Point\tComp Ints Point\tPair Weight\tPair distance"
                                               << endl;
                                        for (int i = 0; i < (int) result.clique.size(); i++) {
                                            const Pair &Pr = grim.getResults().listPairs.at(result.clique.at(i));

                                            ofsLOG << ofss.str()
                                                   << "MATCH\t" << Pr.ref->getIdentifier()
                                                   << "\t" << Pr.comp->getIdentifier()
                                                   << "\t" << Pr.weight
                                                   << "\t" << Pr.ref->fixpos.calcDist(Pr.comp->fixpos) << endl;
                                        }
                                    }
                                    if (!values) {
                                        oss.str("");
                                        if (InputSize == 4) {

                                            oss.str("");
                                            oss << OName;
                                            if (max_clique > 1)oss << ncli;
                                            oss << "_prot.mol2";
                                            mwrite.newFile(oss.str());
                                            mwrite.writeMOL2(complexC.getMole(MoleType::PROTEIN));
                                            oss.str("");
                                            oss << OName;
                                            if (max_clique > 1)oss << ncli;
                                            oss << "_lig.mol2";
                                            mwrite.newFile(oss.str());
                                            mwrite.writeMOL2(complexC.getMole(MoleType::LIGAND));
                                        }
                                        oss.str("");
                                        oss << OName;
                                        if (max_clique > 1)oss << ncli;
                                        oss << "_ints.mol2";
                                        mwrite.newFile(oss.str());
                                        mwrite.writeMOL2(&comparison.Ints);
                                    }


                                }

                                sprintf(buff,
                                        "%d\t%s\t%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%d\n",
                                        ncli,
                                        RName.c_str(),
                                        CName.c_str(),
                                        simils.c_str(),
                                        result.NLig,
                                        result.NCenter,
                                        result.NProt,
                                        result.Tc,
                                        result.rmsd,
                                        nRLig,
                                        nRCent,
                                        nRProt,
                                        nCLig,
                                        nCCent,
                                        nCProt,
                                        result.score,
                                        result.devAl,
                                        result.NPol);
                                ofsRES << buff;

                                ncli++;
                            }
                            for (vector<Fingerprint *>::iterator itF = VFgp.begin(); itF != VFgp.end(); itF++) {
                                delete *itF;
                            }


                        } else {

                            const Results *best_res;
                            double best_score = 0;
                            int best_lig = 0;

                            complexC.genGrid(1.5);
                            complexC.genGrid(4.5);
                            Interactions interC(complexC);
                            InterResults *interRes = new InterResults[nLigand + 1];
                            std::vector <string> names;
                            size_t currLig = 0;

                            while (!mreadC.isEOF()) {
                                Molecule ligand;
                                mreadC.loadNextMolecule(ligand, MoleType::LIGAND);
                                names.push_back(ligand.getName());
                                CName = ligand.getName();
                                ligand.checkMOL2();
                                ligand.ringPerception();
                                interC.calcInteractions(ligand, interRes[currLig], merge, oldh);
                                interC.interToMOL2(interRes[currLig], match_lig, match_prot, match_cent, match_merg);
                                if (verbose) {
                                    cout << "comp  :" << interC.toString(interRes[currLig]) << endl;
                                }
                                Grim grimscr(reference, interRes[currLig]);
                                grimscr.setMinCliqueSize(min_size_cl);

                                unsigned int interMatch = 0;
                                if (!match_merg) {
                                    if (match_lig)interMatch += 2;
                                    if (match_prot)interMatch += 4;
                                    if (match_cent)interMatch += 1;
                                    if (interMatch == 7)interMatch = 0;
                                }

                                int nRLig = 0, nRCent = 0, nRProt = 0, nCLig = 0, nCCent = 0, nCProt = 0;
                                for (size_t iRes = 0; iRes < reference.Ints.numResidus(); ++iRes) {
                                    const Residu &residu = reference.Ints.getResidu(iRes);
                                    const std::string &name = residu.getName();
                                    const std::string pos = name.substr(2, 1);
                                    if (pos == "L")nRLig += residu.numAtom();
                                    else if (pos == "P")nRProt += residu.numAtom();
                                    else if (pos == "C")nRCent += residu.numAtom();
                                }
                                for (size_t iRes = 0; iRes < interRes[currLig].Ints.numResidus(); ++iRes) {
                                    const Residu &residu = interRes[currLig].Ints.getResidu(iRes);
                                    const std::string &name = residu.getName();
                                    const std::string pos = name.substr(2, 1);
                                    if (pos == "L")nCLig += residu.numAtom();
                                    else if (pos == "P")nCProt += residu.numAtom();
                                    else if (pos == "C")nCCent += residu.numAtom();
                                }


                                cout << currLig << " : " << ligand.getName() << endl;
                                try {

                                    grimscr.calcsCliques(interMatch);
                                }
                                catch (MoleExcept &e) {
                                    cout << "Warning Graph not found in file" << endl;
                                    continue;
                                }


                                char buff[500];
                                ofstream ofsRES, ofsLOG;
                                FILE *fp = NULL;
                                fp = fopen("Grimscreen_res.tsv", "r");
                                const bool isOpen = (fp == NULL) ? false : true;
                                if (fp != NULL) fclose(fp);
                                // If the result file does not exists, it will be created with an header
                                ofsRES.open("Grimscreen_res.tsv", ios::app | ios::out);
                                if (!ofsRES.is_open())
                                    throw MoleExcept(9010201, "IChem::calcGrimScr", "Cannot open file ");
                                if (!isOpen)
                                    ofsRES
                                            << "NCli\tRef\tComp\tSimil\tLIG\tCENTER\tPROT\tSumCl\tRMSD\tRLig\tRCent\tRProt\tCLig\tCCent\tCProt\tGrSc\tRMSDAl\tNPol\n";
                                // Otherwise, just opened

                                if (grimscr.getResults().cliques.size() == 0) {
                                    sprintf(buff,
                                            "%d\t%s\t%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%d\n",
                                            0, RName.c_str(), CName.c_str(), simils.c_str(), 0, 0, 0, 0.0, 0.0, 0, 0, 0,
                                            0, 0, 0, 0.0, 0.0, 0);
                                    ofsRES << buff;
                                }

                                mwrite.setWithRotPos(true);
                                Similarity sim(false);
                                ostringstream ofss;
                                unsigned int ncli = 0;

                                for (std::multimap<double, Results>::const_reverse_iterator
                                             itClique = grimscr.firstScore();
                                     itClique != grimscr.lastScore(); ++itClique) {
                                    ++ncli;
                                    if (ncli > max_clique)break;
                                    const Results &result = (*itClique).second;
                                    //                    Fingerprint *fgp = new Fingerprint(grimscr.getResults().listPairs.size());
                                    //                    for (size_t i=0; i < result.clique.size(); i++) fgp->bitOn(result.clique.at(i));
                                    //                    ostringstream oss;
                                    //                    selected=true;
                                    //                    if (VFgp.size())
                                    //                    {
                                    //                        for (vector<Fingerprint*>::iterator itF = VFgp.begin(); itF != VFgp.end(); itF++)
                                    //                        {
                                    //                            sim.setRef(*fgp);
                                    //                            sim.setComp(*(*itF));

                                    //                            if (sim.Tanimoto() > 0.9){selected=false;break;}
                                    //                        }
                                    //                    }
                                    //                        if (!selected) {delete fgp;continue;}

                                    //                    VFgp.push_back(fgp);

                                    sprintf(buff,
                                            "%d\t%s\t%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%d\n",
                                            ncli,
                                            RName.c_str(),
                                            CName.c_str(),
                                            simils.c_str(),
                                            result.NLig,
                                            result.NCenter,
                                            result.NProt,
                                            result.Tc,
                                            result.rmsd,
                                            nRLig,
                                            nRCent,
                                            nRProt,
                                            nCLig,
                                            nCCent,
                                            nCProt,
                                            result.score,
                                            result.devAl,
                                            result.NPol);
                                    ofsRES << buff;
                                    if (result.score > best_score) {
                                        if (verbose)
                                            cout << best_score << " < " << result.score << endl;
                                        best_score = result.score;
                                        best_res = &result;
                                        best_lig = currLig;
                                        bestName = CName;
                                    }

                                    ncli++;


                                }
                                //                for (vector<Fingerprint*>::iterator itF = VFgp.begin(); itF != VFgp.end(); itF++)
                                //                {
                                //                    delete *itF;
                                //                }


                                currLig++;
                            }


                            ofstream ofs("GrScreen.csv", ios::app | ios::out);

                            if (best_score > 0) {
                                Grim grimscr(reference, interRes[best_lig]);
                                grimscr.setMinCliqueSize(min_size_cl);
                                unsigned int interMatch = 0;
                                if (!match_merg) {
                                    if (match_lig)interMatch += 2;
                                    if (match_prot)interMatch += 4;
                                    if (match_cent)interMatch += 1;
                                    if (interMatch == 7)interMatch = 0;
                                }
                                grimscr.calcsCliques(interMatch);

                                ofs << bestName << "\t"
                                    << (*grimscr.firstScore()).second.score << "\t"
                                    << (*grimscr.firstScore()).second.clique.size() << "\t"
                                    << (*grimscr.firstScore()).second.NCenter << "\t"
                                    << (*grimscr.firstScore()).second.NLig << "\t"
                                    << (*grimscr.firstScore()).second.NProt << "\t"
                                    << (*grimscr.firstScore()).second.rmsd << "\t"
                                    << (*grimscr.firstScore()).second.Tc << "\t"
                                    << (*grimscr.firstScore()).second.devAl << "\t"
                                    << (*grimscr.firstScore()).second.NPol << endl;
                            } else {
                                ofs << bestName << "\t0\t0\t0\t0\t0\t0\t0\t0\t0" << endl;
                            }
                        }
                    }catch(MoleExcept &e)
                    {
                        cerr << e.getCode()<<"\t"<<e.getSource()<<"\t"<< e.getData()<<endl;
                    }
                }
            }
        }
    }
    else
    {
        if (InputSize == 4)   // Case where we have protein/ligand multi file
        {

            const std::string& RProt = Input_Values.at(0);
            const std::string& RLig  = Input_Values.at(1);
            const std::string& CProt = Input_Values.at(2);
            const std::string& CLig  = Input_Values.at(3);

            mread.loadNewFile(RProt);mread.loadInComplex(complexR,MoleType::PROTEIN);
            mread.loadNewFile(RLig);

            if (complexR.getNumMolecule(MoleType::PROTEIN)  == 0)
                throw MoleExcept(9010108,
                                 "IChem::runGrim",
                                 "No protein found as reference in "+ RProt);

            const size_t nLigandR = mread.getNumMolecules();

            if (nLigandR == 0) {
                throw MoleExcept(9010108,
                                 "IChem::runGrim",
                                 "No ligand found as reference in "+ RLig);
            }


            complexR.genGrid(1.5);
            complexR.genGrid(4.5);
            if (RName.length() == 0) RName = complexR.getMole(MoleType::PROTEIN)->getName();
            Rnamebu = complexR.getMole(MoleType::PROTEIN)->getName();
            cout << "Docking poses : " << nLigandR << endl;
            int ligand_flag = 0;
            MoleReader mreadLC;
            MoleReader mreadPC;
            Complex complexMC;
            Molecule ProteinC;
            Molecule LigandC;
            while (!mread.isEOF())
            {

                Molecule ligand;
                InterResults referenceM;
                mread.loadNextMolecule(ligand,MoleType::LIGAND);
                RName = Rnamebu +"||"+ ligand.getName();
                ligand.checkMOL2();
                ligand.ringPerception();
                Interactions interR(complexR);
                interR.calcInteractions(ligand,referenceM,merge,oldh);
                interR.interToMOL2(referenceM,match_lig,match_prot,match_cent,match_merg);
                if (verbose){
                    cout << interR.toString(referenceM) << endl;
                }

                mreadLC.loadNewFile(CLig);
                mreadPC.loadNewFile(CProt);

                //                cout << "nombre de comp : " << mreadPC.getNumMolecules() << endl;

                ++ligand_flag;
                int comp_flag = 0 ;

                while (!mreadPC.isEOF())
                {


                    InterResults comparisonM;
                    if(verbose || verbose2){
                        cout << ligand_flag << " : " << ++comp_flag <<"\t"<< RName ;
                    }
                    mreadPC.loadNextMolecule(ProteinC,MoleType::PROTEIN);
                    mreadLC.loadNextMolecule(LigandC,MoleType::LIGAND);
                    CName= ProteinC.getName();

                    LigandC.checkMOL2();
                    LigandC.ringPerception();

                    complexMC.addMolecule(&ProteinC);
                    complexMC.genGrid(1.5);
                    complexMC.genGrid(4.5);
                    Interactions interC(complexMC);
                    interC.calcInteractions(LigandC,comparisonM,merge,oldh);
                    interC.interToMOL2(comparisonM,match_lig,match_prot,match_cent,match_merg);

                    Grim grimscr(referenceM,comparisonM);
                    grimscr.setMinCliqueSize(min_size_cl);
                    grimscr.setDistClose(dif_dclose);
                    grimscr.setDistSame(dif_dsame);
                    grimscr.setdistfar(dif_dfar);
                    unsigned int interMatch=0;
                    if (!match_merg)
                    {
                        if (match_lig)interMatch+=2;
                        if (match_prot)interMatch+=4;
                        if (match_cent)interMatch+=1;
                        if (interMatch==7)interMatch=0;
                    }

                    int  nRLig=0,nRCent=0,nRProt=0,nCLig=0,nCCent=0,nCProt=0;
                    for (size_t iRes=0; iRes < referenceM.Ints.numResidus();++iRes)
                    {
                        const Residu& residu = referenceM.Ints.getResidu(iRes);
                        const std::string& name = residu.getName();
                        const std::string pos = name.substr(2,1);
                        if (pos == "L")nRLig+=residu.numAtom();
                        else if (pos == "P")nRProt+=residu.numAtom();
                        else if (pos == "C")nRCent+=residu.numAtom();
                    }
                    for (size_t iRes=0; iRes < comparisonM.Ints.numResidus();++iRes)
                    {
                        const Residu& residu =comparisonM.Ints.getResidu(iRes);
                        const std::string& name =residu.getName();
                        const std::string pos=name.substr(2,1);
                        if (pos == "L")nCLig+=residu.numAtom();
                        else if (pos == "P")nCProt+=residu.numAtom();
                        else if (pos == "C")nCCent+=residu.numAtom();
                    }

                    if (!((comparisonM.CA != 0 && referenceM.CA != 0)||(comparisonM.CZ != 0 && referenceM.CZ != 0)||(comparisonM.N != 0 && referenceM.N != 0)||(comparisonM.NZ != 0 && referenceM.NZ != 0)||(comparisonM.O != 0 && referenceM.O != 0)||(comparisonM.OD1 != 0 && referenceM.OD1 != 0)||(comparisonM.Zn != 0 && referenceM.Zn != 0)))
                    {
                        if(verbose){
                            cout << "\tfail" << endl;
                        }
                        complexMC.deleteMole(&ProteinC);
                        char buff[500];
                        ofstream ofsRES;
                        FILE *fp = NULL;
                        fp=fopen("Grimscreen_res.tsv","r");
                        const bool isOpen= (fp==NULL)?false:true;
                        if (fp!=NULL) fclose(fp);
                        // If the result file does not exists, it will be created with an header
                        ofsRES.open("Grimscreen_res.tsv",ios::app|ios::out);
                        if (!ofsRES.is_open()) throw MoleExcept(9010201,"IChem::calcGrimScr","Cannot open file ");
                        if (!isOpen) ofsRES<<"NCli\tRef\tComp\tSimil\tLIG\tCENTER\tPROT\tSumCl\tRMSD\tRLig\tRCent\tRProt\tCLig\tCCent\tCProt\tGrSc\tRMSDAl\tNPol\n";
                        // Otherwise, just opened

                        if (grimscr.getResults().cliques.size()==0)
                        {
                            sprintf(buff,"%d\t%s\t%s\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n",
                                    0, RName.c_str(),CName.c_str());
                            ofsRES<<buff;
                        }
                        continue;
                    }

                    grimscr.calcsCliques(interMatch);

                    char buff[500];
                    ofstream ofsRES;
                    FILE *fp = NULL;
                    fp=fopen("Grimscreen_res.tsv","r");
                    const bool isOpen= (fp==NULL)?false:true;
                    if (fp!=NULL) fclose(fp);
                    // If the result file does not exists, it will be created with an header
                    ofsRES.open("Grimscreen_res.tsv",ios::app|ios::out);
                    if (!ofsRES.is_open()) throw MoleExcept(9010201,"IChem::calcGrimScr","Cannot open file ");
                    if (!isOpen) ofsRES<<"NCli\tRef\tComp\tSimil\tLIG\tCENTER\tPROT\tSumCl\tRMSD\tRLig\tRCent\tRProt\tCLig\tCCent\tCProt\tGrSc\tRMSDAl\tNPol\n";
                    // Otherwise, just opened

                    if (grimscr.getResults().cliques.size()==0)
                    {
                        sprintf(buff,"%d\t%s\t%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%d\n",
                                0, RName.c_str(),CName.c_str(),simils.c_str(),0,0,0, 0.0,0.0,0,0,0,0,0,0,0.0,0.0,0);
                        ofsRES<<buff;
                    }

                    mwrite.setWithRotPos(true);
                    Similarity sim(false);
                    ostringstream ofss;
                    unsigned int ncli=0;
                    for (std::multimap<double, Results>::const_reverse_iterator itClique=grimscr.firstScore();itClique != grimscr.lastScore();++itClique)
                    {
                        ++ncli;
                        if(ncli > max_clique)break;
                        const Results& result=(*itClique).second;
                        Fingerprint *fgp = new Fingerprint(grimscr.getResults().listPairs.size());
                        for (size_t i=0; i < result.clique.size(); i++) fgp->bitOn(result.clique.at(i));
                        ostringstream oss;
                        selected=true;
                        if (VFgp.size())
                        {
                            for (vector<Fingerprint*>::iterator itF = VFgp.begin(); itF != VFgp.end(); itF++)
                            {
                                sim.setRef(*fgp);
                                sim.setComp(*(*itF));

                                if (sim.Tanimoto() > 0.9){selected=false;break;}
                            }
                        }
                        if (!selected) {delete fgp;continue;}

                        VFgp.push_back(fgp);

                        sprintf(buff,"%d\t%s\t%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%d\n",
                                ncli,
                                RName.c_str(),
                                CName.c_str(),
                                simils.c_str(),
                                result.NLig,
                                result.NCenter,
                                result.NProt,
                                result.Tc,
                                result.rmsd,
                                nRLig,
                                nRCent,
                                nRProt,
                                nCLig,
                                nCCent,
                                nCProt,
                                result.score,
                                result.devAl,
                                result.NPol);
                        ofsRES<<buff;

                    }
                    if(verbose || verbose2){
                        cout << "\tok" << endl;
                    }
                    //                    complexMC.deleteMole(&ProteinC);
                }
            }
        }
        else{
            throw MoleExcept(9050101,
                             "IChem::runGrim",
                             "Wrong parameters");}
    }
}


