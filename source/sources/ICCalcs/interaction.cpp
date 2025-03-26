#include <iostream>
#include "headers/ICCalcs/interaction.h"
#include "headers/ICMole/complex.h"
#include "headers/ICMole/box.h"
#include "headers/ICMole/cycle.h"
#include "headers/ICCalcs/volsite.h"
#include "headers/ICPars/molereader.h"



//#define ICHEM_DEBUG
using namespace std;
using namespace ICMole;


double  Interactions::Dist_H;
double  Interactions::Dist_Hyd;
double  Interactions::Dist_Ionic;
double  Interactions::Dist_Metal;
double  Interactions::Dist_Arom;

double  Interactions::dist_H;
double  Interactions::dist_Hyd;
double  Interactions::dist_Ionic;
double  Interactions::dist_Metal;
double  Interactions::dist_Arom;

double  Interactions::Angl_H;
double  Interactions::AngT_H;
double  Interactions::Angl_AromFF;
double  Interactions::AngT_AromFF;
double  Interactions::Angl_AromEF;
double  Interactions::AngT_AromEF;

double Interactions::Dist_PiCation;
double Interactions::Angl_PiCation;
double Interactions::AngT_PiCation;


double Interactions::Dist_WHBond;
double Interactions::Angl_WHBond;
double Interactions::AngT_WHBond;


bool    Interactions::wInterType[NB_INTTYPE];

double *Interactions::path_dist;
unsigned int  Interactions::path_dist_size;
unsigned int  Interactions::size_triplets;
bool Interactions::Load_triplet=false;
int Interactions::vect_list[1000][1000];


Interactions::Interactions(ICMole::Complex &cp)
throw(ICMole::MoleExcept):complex(cp),grid(complex.genGrid(4.5))
{
    if (cp.getMole(MoleType::PROTEIN)== (Molecule*)NULL)
        throw MoleExcept(3020101,
                         "Interactions::Interactions",
                         "No protein found in complex");

    setStdRules();
} 


void Interactions::setStdRules()
{
    Dist_H       = 3.5;
    Dist_Hyd     = 4.5;
    Dist_Ionic   = 4.0;
    Dist_Metal   = 2.8;
    // Dist_Arom    = 4.0; // Changed to 5.0
    Dist_Arom    = 5.0;
    dist_H       = 2.3; // 2.5
    dist_Hyd     = 3.2;
    dist_Ionic   = 2.3; // 2.5
    dist_Metal   = 1.8;
    dist_Arom    = 3.2;
    // Dist_PiCation= 4.0; // Changed to 5.0
    Dist_PiCation= 5.0;
    Dist_WHBond  = 2.8;
    Angl_H       = M_PI;
    AngT_H       = M_PI/3;
    Angl_AromFF  = M_PI;
    AngT_AromFF  = M_PI/6;
    Angl_AromEF  = M_PI/2;
    AngT_AromEF  = M_PI/3;
    Angl_PiCation= M_PI;
    AngT_PiCation= M_PI/6;
    Angl_WHBond  = M_PI;
    AngT_WHBond  = M_PI/6;

    for (short i=0; i<NB_INTTYPE;++i)wInterType[i]=true;
    //wInterType[InterType::HBOND_PROT] =false;

}
struct resbest {
    Atom *atmP;
    double dist;
    unsigned short id;} ;


// void Interactions::displayAtomProperties(Molecule& ligand, Atom& atom, string molecule_name) const {

//         cout <<"LIGAND ATOM \t"<<atom.getIdentifier() 
//         << "\t MOL2T: " << atom.getMOL2Type() 
//         << "\t props : "<< atom.props.toString()
//         << "\t charge: " << atom.getFormalCharge() << endl;
// }  



// void Interactions::calcInteractions(Molecule& ligand, InterResults& interResult, bool wMerge, bool oldh, bool mono_prop, bool displayProperties, bool out_lig, bool stdoutt) const // , bool displayProperties
void Interactions::calcInteractions( Molecule& ligand, InterResults& interResult, bool wMerge, bool oldh, bool mono_prop, bool out_lig, bool stdout) const
{
    BoxList boxlist;
    double dist,angle;
    int clash = 0;
    int clash2 = 0;
    int NInter=0;
    double max_allowed_dist=4.5;
    double plp_max = 5.5;
    if (Dist_H > max_allowed_dist)max_allowed_dist=Dist_H;
    if (Dist_Hyd > max_allowed_dist)max_allowed_dist=Dist_Hyd;
    if (Dist_Ionic > max_allowed_dist)max_allowed_dist=Dist_Ionic;
    if (Dist_Metal > max_allowed_dist)max_allowed_dist=Dist_Metal;
    if (Dist_Arom > max_allowed_dist)max_allowed_dist=Dist_Arom;
    if (Dist_PiCation > max_allowed_dist)max_allowed_dist=Dist_PiCation;
    double min_allowed_dist = 1.5;
    double min_clash2 = 2.3;
    double plp,plp_atm,plpA,plpB,plpC,plpD,plpE;

    //    cout << "max allowed distance :" << max_allowed_dist << endl;
    //    cout << "distance pi cation : "<<Dist_PiCation << endl;
    //    cout << "angle pi cation : " << Angl_PiCation << " + " << AngT_PiCation << endl;

    map<Residu*,resbest> hydlist;
    map<Residu*,resbest> ::iterator itHydList;

    if (!complex.isMoleIn(&ligand)) grid.rotateMole(ligand);

    ofstream ofs;
    if (out_lig) {
        if (stdout) {
            cout << ">>># "+ligand.getName()+".ints #<<<" << endl;
            cout  << "#" <<  ligand.getName()<< endl;
        } else {
            string outname = ligand.getName() +".ints";
            ofs.open(outname.c_str(),ios::out|ios::app);
            ofs  << "#" <<  ligand.getName()<< endl;
        }
    }

    string molecule_name = ligand.getName();
    // if(displayProperties) {
    //     cout << "\n\n\t\t Molecule's name: " << molecule_name << "\n\n";
    // }

    for (ItCAtom itLA = ligand.firstAtom(); itLA != ligand.lastAtom();++itLA)
    {
        Atom &atomL=**itLA;
        bool hyd =false;
        plp_atm=0;

    // if(globalConfigOptions.displayProperties) {
    //     displayAtomProperties(ligand, atomL, molecule_name);
        // if(displayProperties) {
        //     displayAtomProperties(ligand, atomL, molecule_name);
        // }

    // }

// #ifdef ICHEM_DEBUG
    // cout << "displayProperties value in calcInteractions: " << displayProperties << "\n";
    // if(displayProperties) {
    //     cout <<"LIGAND ATOM \t"<<atomL.getIdentifier() //<< endl
    //     <<"\t MOL2T: " << atomL.getMOL2Type() //<< endl
    //     <<"\t props : "<< atomL.props.toString() //<< endl
    //     <<"\t charge: " << atomL.getFormalCharge() << endl;
    //     //  cout << " END" << endl;
    // }
// #endif

        if (atomL.isHydrogen() || (atomL.getName()=="DuCy")
                || !atomL.isUsed()) continue;
        if (atomL.getBox(&grid) == (Box*)NULL) continue;

        // out ligand possible link
        if (out_lig) {
            if (stdout) {
                cout  << atomL.getIdentifier() ;
            } else {
                ofs  << atomL.getIdentifier() ;
            }
        }
        int nbused = 0;
        Box& boxL = *atomL.getBox(&grid);
        boxlist.clear();
        grid.getAdjacency(boxlist,&boxL,plp_max*2,true);
        //        grid.printInFile(atomL.getIdentifier(),boxlist,true);

        hydlist.clear();;
        for (ItBox itBAdj=boxlist.begin();
             itBAdj != boxlist.end();
             ++itBAdj)
        {
            for (ItCAtom itAA  = (*itBAdj)->firstAtom();
                 itAA != (*itBAdj)->lastAtom();
                 itAA++)
            {
                Atom &atomP = **itAA;
                //                if (&atomP.getResidu()) //here

                if (&atomP.getParent()==&ligand
                        || !atomP.isUsed()
                        || atomP.isHydrogen())continue;
                dist = atomL.calcFixpos(atomP,plp_max+0.1);
                //                cout << atomP.props.toString() << "  "<< dist;
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////


                plp=0;plpA=0;plpB=0;plpC=0;plpD=0;plpE=0;

                if(  dist < 3.4  && ( atomL.props.isMetal() && (atomP.props.isAcceptor() || atomP.props.isAnion())
                         ||   (atomL.props.isAnion() || atomL.props.isAcceptor() ) && atomP.props.isMetal()
                         ||   (atomP.isMetallic() && atomL.getAtomicName() == "N"))){
                    plpA = 1.4; plpB = 2.2; plpC = 3.1; plpD = 3.4;plpE = -5;
                         atomL.props.setMetalA(true);
                }
                else if (   ((atomL.props.isAcceptor() || atomL.props.isAnion()) && (atomP.props.isDonor() || atomP.props.isCation()))
                            || ((atomL.props.isDonor() || atomL.props.isCation())  && (atomP.props.isAcceptor() || atomP.props.isAnion()))){
                    plpA = 2.3; plpB = 2.6; plpC = 3.1; plpD = 3.4; plpE = -2;

                }
                else if (atomL.props.isApolar() && atomP.props.isApolar()){
                    plpA = 3.4; plpB = 3.6; plpC = 4.5; plpD = 5.5; plpE = -0.4;
                }
                if (dist < plpA && plpA != 0){
                    plp = 20*(plpA-dist)/plpA;
                }
                else if (dist >= plpA && dist < plpB){
                    plp = plpE*(dist-plpA)/(plpB-plpA);
                }
                else if (dist >= plpB && dist < plpC){
                    plp = plpE;
                }
                else if (dist >= plpC && dist < plpD){
                    plp = plpE*(plpD-dist)/(plpD-plpC);
                }
                else if (dist >= plpD){
                    plp = 0;
                }
                plp_atm += plp;
                #ifdef ICHEM_DEBUG
                if (plp != 0 ){
                    cout << "PLP :\t" <<  atomL.getIdentifier() << " " <<atomL.props.toString() << "\t||\t" << atomP.getIdentifier() << " " << atomP.props.toString() << "\tdist\t"<< dist <<"\tplp " << plp << "\t" << plp_atm <<  endl;
                }

#endif

                if (dist > max_allowed_dist )
                {
                    continue;
                }
                if (dist < min_clash2){
                    clash2 ++;
                }if (dist < min_allowed_dist){
                    clash ++;
                    continue;
                }
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////






                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                if (atomL.props.isAcceptor())
                {
                    //                    cout <<atomL.toString() << endl;
                    // LIGAND ACCEPTOR <-> PROTEIN DONOR => HBOND_LIG
                    if (wInterType[InterType::HBOND_PROT]
                            && atomP.props.isDonor()
                            && dist <= Dist_H && dist >= dist_H)
                    {
                        for (size_t i=0; i< atomP.getNumBond();++i)
                        {
                            const Atom &atomP2 = atomP.getAtomLinked(i);

                            if (!atomP2.isHydrogen())     continue;

                            angle = atomP2.fixpos.calcAngle(atomP.fixpos,  atomL.fixpos);
                            //                            cout << "HBOND_PROT , angle : "<< angle << endl;

                            if (angle <= Angl_H-AngT_H
                                    || angle >= Angl_H+AngT_H)  { continue;}

#ifdef ICHEM_DEBUG
                            cout <<"HBLA\t"<<atomL.getIdentifier()
                                <<"\t"<< atomP.getIdentifier()
                               << "\t" << dist
                               << "\t"<< (angle*180/M_PI)
                               <<"\t"<<atomP2.getIdentifier()<<endl ;
#endif
                            if (out_lig) {
                                nbused++;
                            }
                            InterPoint IntP(NInter,
                                            &atomP,
                                            &atomL,
                                            (atomL.fixpos+atomP.fixpos)/2,
                                            InterType::HBOND_PROT,
                                            dist,
                                            angle);
                            NInter++;
                            interResult.listInters.push_back(IntP);
                            interResult.N++;

                        }
                    }

                    // LIGAND ACCEPTOR <-> PROTEIN METAL => METAL
                    //                    cout << wInterType[InterType::METAL] << atomP.props.isMetal() << Dist_Metal << endl;
                    if (wInterType[InterType::METAL]
                            && atomP.props.isMetal()
                            && dist <= Dist_Metal && dist >= dist_Metal)
                    {
                        if (out_lig) {nbused++;
                        }
                        InterPoint IntP(NInter,&atomP,&atomL,(atomL.fixpos+atomP.fixpos)/2, InterType::METAL, dist);
                        NInter++;
                        interResult.listInters.push_back(IntP);
                        interResult.Zn++;
#ifdef ICHEM_DEBUG
                        cout << "MEAC\t"<<atomL.getIdentifier()<<"\t"<<atomP.getIdentifier()<<" \tDist:"<<dist<<endl;
#endif
                    }

                    // LIGAND ACCEPTOR <-> PROTEIN WEAK DONOR => WEAK HBOND
                    if (wInterType[InterType::WHBOND_PROT]
                            && atomP.props.isweakDonor()
                            && dist <= Dist_WHBond)
                    {

                        for (size_t i=0; i< atomP.getNumBond();++i)
                        {
                            const Atom &atomP2 = atomP.getAtomLinked(i);
                            if (!atomP2.isHydrogen())     continue;

                            angle = atomP2.fixpos.calcAngle(atomP.fixpos,  atomL.fixpos);
                            if (angle <= Angl_WHBond-AngT_WHBond
                                    || angle >= Angl_WHBond+AngT_WHBond)  { continue;}
#ifdef ICHEM_DEBUG
                            cout <<"WHBLA\t"<<atomL.getIdentifier()
                                <<"\t"<< atomP.getIdentifier()
                               << "\t" << dist
                               << "\t"<< (angle*180/M_PI)
                               <<"\t"<<atomP2.getIdentifier()
                              <<endl ;
#endif
                            if (out_lig) {
                                nbused++;
                            }
                            InterPoint IntP(NInter,
                                            &atomP,
                                            &atomL,
                                            (atomL.fixpos+atomP.fixpos)/2,
                                            InterType::WHBOND_PROT,
                                            dist,
                                            angle);
                            NInter++;
                            interResult.listInters.push_back(IntP);
                        }
                    }
                }// END isAcceptor

                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////

                // LIGAND WEAK ACCEPTOR <-> PROTEIN WEAK DONOR OR DONOR => WEAK HBOND LIG
                if (atomL.props.isweakAcceptor()
                        && wInterType[InterType::WHBOND_PROT]
                        && dist < Dist_WHBond
                        && (atomP.props.isweakDonor() || atomP.props.isDonor()))
                {
                    for (size_t i=0; i< atomP.getNumBond();++i)
                    {
                        const Atom &atomP2 = atomP.getAtomLinked(i);
                        if (!atomP2.isHydrogen())     continue;

                        angle = atomP2.fixpos.calcAngle(atomP.fixpos,  atomL.fixpos);
                        if (angle <= Angl_WHBond-AngT_WHBond
                                || angle >= Angl_WHBond+AngT_WHBond)  { continue;}
#ifdef ICHEM_DEBUG
                        cout <<"WHBLA\t"<<atomL.getIdentifier()
                            <<"\t"<< atomP.getIdentifier()
                           << "\t" << dist
                           << "\t"<< (angle*180/M_PI)
                           <<"\t"<<atomP2.getIdentifier()
                          <<endl ;
#endif
                        if (out_lig) {
                            nbused++;
                        }
                        InterPoint IntP(NInter,
                                        &atomP,
                                        &atomL,
                                        (atomL.fixpos+atomP.fixpos)/2,
                                        InterType::WHBOND_PROT,
                                        dist,
                                        angle);
                        NInter++;
                        interResult.listInters.push_back(IntP);
                    }

                }




                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////

                if (atomL.props.isDonor() && mono_prop==false)
                {

                    // LIGAND DONOR <-> PROTEIN ACCEPTOR => HBOND_PROT
                    if (wInterType[InterType::HBOND_LIG]
                            && atomP.props.isAcceptor()
                            && dist <= Dist_H && dist >= dist_H)
                    {
                        for (size_t i=0; i< atomL.getNumBond();++i)
                        {
                            const Atom &atomL2 = atomL.getAtomLinked(i);
                            if (!atomL2.isHydrogen())     continue;

                            angle = atomL2.fixpos.calcAngle(atomP.fixpos,  atomL.fixpos);
                            if (angle <= Angl_H-AngT_H || angle >= Angl_H+AngT_H)  {continue;}

#ifdef ICHEM_DEBUG
                            cout <<"HBLA\t"<<atomL.getIdentifier()
                                <<"\t"<< atomP.getIdentifier()
                               << "\t" << dist
                               << "\t"<< (angle*180/M_PI)
                               <<"\t"<<atomL2.getIdentifier()<<endl ;
#endif
                            if (out_lig) {
                                nbused++;
                            }
                            InterPoint IntP(NInter,
                                            &atomP,
                                            &atomL,
                                            (atomL.fixpos+atomP.fixpos)/2,
                                            InterType::HBOND_LIG,
                                            dist,
                                            angle);
                            NInter++;
                            interResult.listInters.push_back(IntP);
                            interResult.O++;


                        }
                    }
                    // LIGAND DONOR <-> PROTEIN WEAK ACCEPTOR => WEAK HBOND PROT
                    if (wInterType[InterType::WHBOND_LIG]
                            && atomP.props.isweakAcceptor()
                            && dist <= Dist_WHBond)
                    {
                        for (size_t i=0; i< atomL.getNumBond();++i)
                        {
                            const Atom &atomL2 = atomL.getAtomLinked(i);
                            if (!atomL2.isHydrogen())     continue;

                            angle = atomL2.fixpos.calcAngle(atomP.fixpos,  atomL.fixpos);
                            if (angle <= Angl_WHBond-AngT_WHBond
                                    || angle >= Angl_WHBond+AngT_WHBond)  {continue;}

#ifdef ICHEM_DEBUG
                            cout <<"WHBLA\t"<<atomL.getIdentifier()
                                <<"\t"<< atomP.getIdentifier()
                               << "\t" << dist
                               << "\t"<< (angle*180/M_PI)
                               <<"\t"<<atomL2.getIdentifier()<<endl ;
#endif
                            if (out_lig) {
                                nbused++;
                            }
                            InterPoint IntP(NInter,
                                            &atomP,
                                            &atomL,
                                            (atomL.fixpos+atomP.fixpos)/2,
                                            InterType::WHBOND_LIG,
                                            dist,
                                            angle);
                            NInter++;
                            interResult.listInters.push_back(IntP);


                        }
                    }

                }// END isDonor

                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////

                if (atomL.props.isweakDonor())
                {

                    // LIGAND WEAK DONOR <-> PROTEIN  ACCEPTOR => WEAK HBOND
                    if (wInterType[InterType::WHBOND_LIG]
                            && (atomP.props.isAcceptor() || atomP.props.isweakAcceptor())
                            && dist <= Dist_WHBond)
                    {
                        for (size_t i=0; i< atomL.getNumBond();++i)
                        {
                            const Atom &atomL2 = atomL.getAtomLinked(i);
                            if (!atomL2.isHydrogen())     continue;

                            angle = atomL2.fixpos.calcAngle(atomP.fixpos,  atomL.fixpos);
                            if (angle <= Angl_WHBond-AngT_WHBond
                                    || angle >= Angl_WHBond+AngT_WHBond)  {continue;}

#ifdef ICHEM_DEBUG
                            cout <<"WHBLA\t"<<atomL.getIdentifier()
                                <<"\t"<< atomP.getIdentifier()
                               << "\t" << dist
                               << "\t"<< (angle*180/M_PI)
                               <<"\t"<<atomL2.getIdentifier()<<endl ;
#endif
                            if (out_lig) {
                                nbused++;
                            }
                            InterPoint IntP(NInter,
                                            &atomP,
                                            &atomL,
                                            (atomL.fixpos+atomP.fixpos)/2,
                                            InterType::WHBOND_LIG,
                                            dist,
                                            angle);
                            NInter++;
                            interResult.listInters.push_back(IntP);


                        }
                    }


                }// END isweakDonor



                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////

                if (wInterType[InterType::IONIC_PROT]
                        && atomL.props.isAnion()
                        && atomP.props.isCation()
                        && !atomP.props.isMetal()
                        && dist <= Dist_Ionic && dist >= dist_Ionic)
                {
                    if (out_lig) {
                        nbused++;
                    }
                    InterPoint IntP(NInter,
                                    &atomP,
                                    &atomL,
                                    (atomL.fixpos+atomP.fixpos)/2,
                                    InterType::IONIC_PROT,
                                    dist);
                    NInter++;
                    interResult.listInters.push_back(IntP);
                    interResult.NZ++;

#ifdef ICHEM_DEBUG
                    cout <<"ACCA\t"<<atomL.getIdentifier()
                        <<"\t"<< atomP.getIdentifier()
                       << "  \tDist=" << dist<<endl  ;
#endif
                }


                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                if (atomL.props.isCation() && mono_prop==false)
                {
                    if (wInterType[InterType::IONIC_LIG]
                            && atomP.props.isAnion()
                            && dist <= Dist_Ionic && dist >= dist_Ionic)
                    {
                        if (out_lig) {
                            nbused++;
                        }
                        InterPoint IntP(NInter,
                                        &atomP,
                                        &atomL,
                                        (atomL.fixpos+atomP.fixpos)/2,
                                        InterType::IONIC_LIG,
                                        dist);

                        NInter++;
                        interResult.listInters.push_back(IntP);
                        interResult.OD1++;
#ifdef ICHEM_DEBUG
                        cout <<"CAAC\t"<<atomL.getIdentifier()
                            <<"\t"<< atomP.getIdentifier()
                           << "  \tDist=" << dist<<endl  ;
#endif
                    }

                    if (wInterType[InterType::PICATION]
                            && atomP.getName()=="DuCy"
                            && dist < Dist_PiCation)
                    {
                        Cycle *cycleP = atomP.getParent().getCycleFromCenter(&atomP);
                        if (cycleP != (Cycle*)NULL)
                        {
                            cycleP->calcVector();
                            const double angle=cycleP->getCenter().fixpos.calcAngle(cycleP->getNormVector(),atomL.fixpos);
                            //                            cout << "Test pi cation dist :" << dist << " angle :" << angle << endl;
                            if (angle <= Angl_PiCation-AngT_PiCation
                                    || angle >= Angl_PiCation+AngT_PiCation)continue;
                            if (out_lig) {
                                nbused++;
                            }
                            InterPoint IntP(NInter,
                                            &cycleP->getCenter(),
                                            &atomL,
                                            (atomL.fixpos+cycleP->getCenter().fixpos)/2,
                                            InterType::PICATION,dist,angle);NInter++;
                            interResult.listInters.push_back(IntP);
                        }
                    }

                }
                // SPECIAL CASES :

                // N linked to sulphonamide in interaction with Metal

                if (wInterType[InterType::METAL]&&
                        atomL.isNitrogen()
                        &&atomP.props.isMetal()
                        &&dist < Dist_Metal)
                {
                    bool linkedtosulf=false;
                    for (size_t iAtmLink=0;iAtmLink < atomL.getNumBond();++iAtmLink)
                    {
                        const Atom &atml = atomL.getAtomLinked(iAtmLink);
                        if (atml.isSulfur()) linkedtosulf=true;
                    }
                    if (!linkedtosulf) continue;
                    if (out_lig) {
                        nbused++;
                    }
//                    cout << "mon ints metal " << atomL.props.toString() << endl;

                    InterPoint IntP(NInter,
                                    &atomP,
                                    &atomL,
                                    (atomL.fixpos+atomP.fixpos)/2,
                                    InterType::METAL,
                                    dist);
                    NInter++;
                    interResult.listInters.push_back(IntP);
                    interResult.Zn++;


#ifdef ICHEM_DEBUG
                    cout << "MEAC\t"<<atomL.getIdentifier()<<"\t"<<atomP.getIdentifier()<<" \tDist:"<<dist<<endl;
#endif

                }
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                if (!oldh){
                    if (wInterType[InterType::HYDROPHOBIC]
                            && atomL.props.isHydrophobic() && atomP.props.isHydrophobic()
                            && !(atomL.props.isAromatic() && atomP.props.isAromatic())
                            && dist <= Dist_Hyd && dist >= dist_Hyd)
                    {
                        AtomList listAtoms;
                        Box box = **itBAdj;
                        int nbhyd=1,nbatm=1;
                        grid.getAdjacentAtoms(listAtoms,box,4.5);
                        //                    cout << box.getId() << "  " << atomP.getIdentifier()      << endl;
                        //                    cout << "toto" << endl;
                        for (ItCAtom itAtm= listAtoms.begin(); itAtm != listAtoms.end(); ++itAtm)
                        {

                            Atom &atm = **itAtm;
                            if (atm.isHydrogen()) continue;
                            if (atm.props.isHydrophobic()){
                                nbhyd++;
                            }
                            nbatm++;
                            //                        cout << "\t\t" << atm.getIdentifier() << "\tdist : "<<  atomP.fixpos.calcDist(atm.fixpos) <<endl;
                        }
                        //                    cout << "tata " << endl;
                        //                        cout << atomP.getIdentifier() <<"\t\t" <<100*nbhyd/nbatm  << "\t|| "<<  nbhyd << "\t|| " << nbatm << endl;
                        if (nbhyd==0 || nbatm==0 )continue;
                        if ( 100*nbhyd/nbatm > 50){
                            //                        cout << "titi" << endl;
                            itHydList=hydlist.find(atomP.getResidu());
                            if (itHydList == hydlist.end())
                            {
                                resbest rbest;rbest.atmP=&atomP;rbest.dist=dist;
                                hydlist.insert(pair<Residu*,resbest>(atomP.getResidu(),rbest));
                            }
                            else
                            {
                                resbest &rbest = (*itHydList).second;
                                if ( dist < rbest.dist)
                                {
                                    rbest.dist=dist;
                                    rbest.atmP=&atomP;
                                }
                            }
                        }
                        //                    cout << "tutu" << endl;
                    }

                }else {
                    if (wInterType[InterType::HYDROPHOBIC]
                            && atomL.props.isHydrophobic() && atomP.props.isHydrophobic()
                            && !(atomL.props.isAromatic() && atomP.props.isAromatic())
                            && dist <= Dist_Hyd && dist >= dist_Hyd)
                    {

                        itHydList=hydlist.find(atomP.getResidu());
                        if (itHydList == hydlist.end())
                        {
                            resbest rbest;rbest.atmP=&atomP;rbest.dist=dist;
                            hydlist.insert(pair<Residu*,resbest>(atomP.getResidu(),rbest));
                        }
                        else
                        {
                            resbest &rbest = (*itHydList).second;
                            if ( dist < rbest.dist)
                            {
                                rbest.dist=dist;
                                rbest.atmP=&atomP;
                            }
                        }

                    }

                }
            }// END itAA

        }// END itBADj

        atomL.setPartialCharge(plp_atm);

        for (itHydList = hydlist.begin(); itHydList != hydlist.end(); ++itHydList)
        {
            if (out_lig) {
                nbused++;
            }

            //            cout << atomL.getResidu()->getIdentifier() << " avec " << (*itHydList).second.atmP->getResidu()->getIdentifier() << endl;

            if (atomL.getResidu()->getIdentifier() == (*itHydList).second.atmP->getResidu()->getIdentifier()){continue;}

            InterPoint IntP(NInter,
                            (*itHydList).second.atmP,
                            &atomL,
                            (atomL.fixpos+(*itHydList).second.atmP->fixpos)/2,
                            InterType::HYDROPHOBIC,
                            (*itHydList).second.dist);
            NInter++;
            interResult.listInters.push_back(IntP);
            interResult.CA++;
#ifdef ICHEM_DEBUG
            cout << "HYDR\t"<<atomL.getIdentifier()
                 <<"\t"<< (*itHydList).second.atmP->getIdentifier()
                <<"\t"<< (*itHydList).second.dist<<endl;
#endif
        }
        if (out_lig) {
            if (stdout) {
                cout << "\t| " << nbused << " | " << atomL.props.toString() << endl; ;
            } else {
                ofs << "\t| " << nbused << " | " << atomL.props.toString() << endl; ;
            }
        }
    }// END itLA*/




    if (!wInterType[InterType::AREDGEFACE]
            &&!wInterType[InterType::ARFACEFACE]
            &&!wInterType[InterType::PICATION])
    {
        if (wMerge) mergeInteractions(interResult);
        return;
    }

    AtomList toRotate;
#ifdef ICHEM_DEBUG
    cout << "AROMATIC INTERACTIONS"<<endl;
#endif
    CycleList cyclelist;
    bool possible =true;
    unsigned short EF=0;bool Ar=false;double distCENTER;



    for (ItCCycle itLC= ligand.firstCycle();
         itLC!= ligand.lastCycle();
         itLC++)
    {
        // Getting ligand cycle :
        Cycle& ligcycle=**itLC;

        // Not aromatic : continue;
        if (!ligcycle.isAromatic()) continue;

        // Getting atom center :
        Atom &atomL=ligcycle.getCenter();
#ifdef ICHEM_DEBUG
        cout << "LOOKING AT : "<<atomL.getIdentifier()
             <<"::"<<atomL.fixpos.toString()<<endl;
#endif
        // Checking that atom center is located in the cube :
        if (atomL.getBox(&grid) == (Box*)NULL)
        {
            // When not, we find its cube :
            toRotate.clear();
            toRotate.push_back(&atomL);
            grid.rotateAtoms(toRotate);
            // When we cannot find it, we ignore the cycle :
            if (atomL.getBox(&grid) == (Box*)NULL)
                continue;
        }

        // Fetching cube :
        Box& boxL = *atomL.getBox(&grid);

        ligcycle.calcVector();
        boxlist.clear();
        // Getting adjacent cube of this cube :
        grid.getAdjacency(boxlist,&boxL,max_allowed_dist*1.5,true);


        // Scanning every atom of every adjacent cube looking for Aromatic center:
        for (ItBox   itBAdj  = boxlist.begin();
             itBAdj != boxlist.end();
             ++itBAdj)
            for (ItCAtom itAA    = (*itBAdj)->firstAtom();
                 itAA   != (*itBAdj)->lastAtom();
                 itAA++)
            {
                Atom &atomP = **itAA;
                ItCAtom itAAA = ligcycle.first();
                Atom &atomcycle = **itAAA;
                //                cout << atomP.getResidu()->getIdentifier() <<" :: " <<  atomcycle.getResidu()->getIdentifier() << endl;
                if (atomP.getResidu()->getIdentifier() == atomcycle.getResidu()->getIdentifier() ) continue;
                //                    cout <<atomcycle.getResidu()->getIdentifier()  <<" :: " <<  atomP.getResidu()->getIdentifier() << endl;
                if (&atomP.getParent()==&ligand) continue;
                if (atomP.props.isCation())
                {
                    dist = atomL.calcFixpos(atomP,Dist_PiCation+0.1);
                    if (dist > Dist_PiCation)continue;// Cannot be aromatic either

                    dist = atomL.calcFixpos(atomP,Dist_PiCation+0.1);
                    if (dist > Dist_PiCation || dist < min_allowed_dist )continue;// Cannot be aromatic either

                    const double angle=ligcycle.getCenter().fixpos.calcAngle(ligcycle.getNormVector(),atomP.fixpos);
                    if ((angle >= Angl_PiCation-AngT_PiCation
                            && angle <= Angl_PiCation+AngT_PiCation)
                            ||(angle >= Angl_PiCation-AngT_PiCation-PI
                               && angle <= Angl_PiCation+AngT_PiCation-PI)
                            ){

                    InterPoint IntP(NInter,
                                    &atomP,
                                    &ligcycle.getCenter(),
                                    (atomP.fixpos+ligcycle.getCenter().fixpos)/2,
                                    InterType::PICATION,dist,angle);NInter++;
                    interResult.listInters.push_back(IntP);
                    }
                }


                if (atomP.getName()!="DuAr")continue;
                dist = atomL.calcFixpos(atomP,10.1);
                // Since we have both cycle centers here, we can check if they are
                // not too far away :
                if (dist > 10 || dist < min_allowed_dist) continue;

                // Getting the corresponding cycle from this atom center :
                Cycle* cycleP=atomP.getParent().getCycleFromCenter(&atomP);
                if (cycleP==(Cycle*)NULL)continue;
                // Saving it :
                cyclelist.push_back(cycleP);
            }

        // Sort and unique to avoid cycle redundance :

        SortUniqueVector(cyclelist);

        for (ItCCycle itC = cyclelist.begin()
             ; itC != cyclelist.end();itC++)
        {
            Cycle &cycleP = **itC;

            distCENTER = cycleP.getFixpos().calcDist(ligcycle.getFixpos());
#ifdef ICHEM_DEBUG
            cout << "  ||-->"<<cycleP.getCenter().fixpos.toString()
                 <<" " << ligcycle.getCenter().fixpos.toString()<<" " << distCENTER<<endl;
#endif
            // CASE AROMATIC CENTER TOO FAR => POSSIBLE HYDROPHOBIC INTERACTION :
            if (wInterType[InterType::HYDROPHOBIC] && distCENTER > Dist_Arom)
            {
                double bestdist=100;Atom* Hlig=(Atom*)NULL, *Hprot=(Atom*)NULL;
                // SCANNING LIGAND SIDE ATOMS IN CYCLE :
                for (ItCAtom itCLA = ligcycle.first();
                     itCLA!= ligcycle.end();
                     itCLA++)
                {
                    Atom& atomCLA = **itCLA;
                    if (!atomCLA.props.isHydrophobic()) continue;

                    // SCANNING PROTEIN SIDE ATOMS IN CYCLE :
                    for (ItCAtom itCLP = cycleP.first();
                         itCLP!= cycleP.end();
                         itCLP++)
                    {
                        Atom& atomCLP = **itCLP;
                        if (!atomCLP.props.isHydrophobic()) continue;

                        // GETTING CLOSEST ATOMS BETWEEN PROT/LIGAND
                        dist = atomCLA.calcFixpos(atomCLP);
                        if (dist < dist_H || bestdist > Dist_H) continue;
                        if (dist < bestdist)
                        {
                            Hlig=&atomCLA;
                            Hprot=&atomCLP;
                            bestdist=dist;
                        }
                    }

                }

#ifdef ICHEM_DEBUG
                cout << "HYD"<<endl;
#endif
                //                cout << " 2 : " << Hprot->getIdentifier() << " | " << Hlig->getIdentifier() << endl;
                // BELOW Dist_H : HYDROPHOBIC INTERACTION
                if (bestdist > Dist_H || bestdist < dist_H)continue;

                if (Hprot->getResidu()->getIdentifier() == Hlig->getResidu()->getIdentifier()) continue;

                InterPoint IntP(NInter,
                                Hprot,
                                Hlig,
                                (Hprot->fixpos+Hlig->fixpos)/2,
                                InterType::HYDROPHOBIC,
                                bestdist);
                NInter++;
                interResult.listInters.push_back(IntP);
                interResult.CA++;
                continue;
            }

            // CASE AROMATIC CENTER CLOSE :
            cycleP.calcVector();
            possible =true;
            EF=0; Ar=false;
            // SCANNING LIGAND ATOM VS PROTEIN ATOM :
            for (size_t itCLA =0; itCLA < ligcycle.getNumAtom(); itCLA++)
            {

                Atom& atomCLA =*ligcycle.getAtom(itCLA);
                if (!atomCLA.isCarbon())continue;
                const Coords& normalCLA = ligcycle.getNormVector(itCLA);

                for (size_t itCLP =0; itCLP < cycleP.getNumAtom(); itCLP++)

                {
                    Atom& atomCLP = *cycleP.getAtom(itCLP);
                    if (!atomCLP.isCarbon())continue;
                    dist = atomCLA.calcFixpos(atomCLP);

                    if (dist > Dist_Arom*3){possible=false;break;}
                    const Coords& normalCLP=cycleP.getNormVector(itCLP);
                    const double angle =
                            acos((atomCLA.fixpos.x-normalCLA.x)*(atomCLP.fixpos.x-normalCLP.x)
                                 +(atomCLA.fixpos.y-normalCLA.y)*(atomCLP.fixpos.y-normalCLP.y)
                                 +(atomCLA.fixpos.z-normalCLA.z)*(atomCLP.fixpos.z-normalCLP.z)
                                 );
                    if ((dist < Dist_Arom && dist > dist_Arom)
                            && ((angle >      Angl_AromFF-AngT_AromFF && angle <      Angl_AromFF+AngT_AromFF)
                                ||  (angle > M_PI+Angl_AromFF-AngT_AromFF && angle < M_PI+Angl_AromFF+AngT_AromFF)
                                ||  (angle >      Angl_AromFF-AngT_AromFF-M_PI && angle <      Angl_AromFF+AngT_AromFF-M_PI)))  Ar=true;
                    if ((dist < Dist_Arom && dist > dist_Arom)
                            && ((angle >      Angl_AromEF-AngT_AromEF && angle <      Angl_AromEF+AngT_AromEF)
                                || (angle > M_PI+Angl_AromEF-AngT_AromEF && angle < M_PI+Angl_AromEF+AngT_AromEF)
                                || (angle >      Angl_AromEF-AngT_AromEF-M_PI && angle <      Angl_AromEF+AngT_AromEF-M_PI)))  EF++;
                }
                if (!possible)break;
            }
            if (Ar || EF > 5){
                //                cout << "ligcenter " <<ligcycle.getCenter().props.toString() << endl;
                Box& boxL = *ligcycle.getCenter().getBox(&grid);
                boxlist.clear();
                grid.getAdjacency(boxlist,&boxL,max_allowed_dist*2,true);
                //        grid.printInFile(atomL.getIdentifier(),boxlist,true);

                hydlist.clear();;
                plp_atm =0;
                for (ItBox itBAdj=boxlist.begin();
                     itBAdj != boxlist.end();
                     ++itBAdj)
                {
                    for (ItCAtom itAA  = (*itBAdj)->firstAtom();
                         itAA != (*itBAdj)->lastAtom();
                         itAA++)
                    {
                        Atom &atomP = **itAA;
                        //                if (&atomP.getResidu()) //here
                        if (&atomP.getParent()==&ligand
                                || !atomP.isUsed()
                                || atomP.isHydrogen())continue;
                        dist = atomL.calcFixpos(atomP,max_allowed_dist+0.1);
                        if (dist > max_allowed_dist )
                        {
                            continue;
                        }
                        if (dist < min_clash2){
                            clash2 ++;
                        }if (dist < min_allowed_dist){
                            clash ++;
                            continue;
                        }
                        ////////////////////////////////////////////////////////////////////////////////
                        ////////////////////////////////////////////////////////////////////////////////


                        plpA=0;plpB=0;plpC=0;plpD=0;plpE=0;
                        if (atomP.props.isApolar()){
                            plpA = 3.4; plpB = 3.6; plpC = 4.5; plpD = 5.5; plpE = -0.4;
                        }
                        if (dist < plpA && plpA != 0){
                            plp = 20*(plpA-dist)/plpA;
                        }
                        else if (dist >= plpA && dist < plpB){
                            plp = plpE*(dist-plpA)/(plpB-plpA);
                        }
                        else if (dist >= plpB && dist < plpC){
                            plp = plpE;
                        }
                        else if (dist >= plpC && dist < plpD){
                            plp = plpE*(plpD-dist)/(plpD-plpC);
                        }
                        else if (dist >= plpD){
                            plp = 0;
                        }
                        plp_atm += plp;
                    }
                }

                ligcycle.getCenter().setPartialCharge(plp_atm);
            }
            if (wInterType[InterType::ARFACEFACE] && Ar)
            {

                InterPoint IntP(NInter,
                                &cycleP.getCenter(),
                                &ligcycle.getCenter(),
                                (cycleP.getCenter().fixpos+ligcycle.getCenter().fixpos)/2,
                                InterType::ARFACEFACE,distCENTER);NInter++;
                interResult.listInters.push_back(IntP);

                interResult.CZ++;
#ifdef ICHEM_DEBUG
                cout << "|||->AROMATIC FACE FACE"<<endl;
#endif
            }else if (wInterType[InterType::AREDGEFACE] && EF>5){
                InterPoint IntP(NInter,
                                &cycleP.getCenter(),
                                &ligcycle.getCenter(),
                                (cycleP.getCenter().fixpos+ligcycle.getCenter().fixpos)/2,
                                InterType::AREDGEFACE,distCENTER);NInter++;
                interResult.listInters.push_back(IntP);
                interResult.CZ++;
#ifdef ICHEM_DEBUG
                cout << "|||->AROMATIC EDGE FACE"<<endl;
#endif
            }


        }

    }

    if (wMerge){
        mergeInteractions(interResult);

    }
    ofs  << "Clash (<1.8): " <<  clash << "  Clash (<2.3): " <<  clash2 << endl;
}
void Interactions::calcInteractionsppi(InterResults& interResult, bool wMerge, bool wHydrogen, bool oldh) const
{

    BoxList boxlist;
    double dist,angle;

    int NInter=0;
    double max_allowed_dist=0;
    if (Dist_H > max_allowed_dist)max_allowed_dist=Dist_H;
    if (Dist_Hyd > max_allowed_dist)max_allowed_dist=Dist_Hyd;
    if (Dist_Ionic > max_allowed_dist)max_allowed_dist=Dist_Ionic;
    if (Dist_Metal > max_allowed_dist)max_allowed_dist=Dist_Metal;
    if (Dist_Arom > max_allowed_dist)max_allowed_dist=Dist_Arom;

    //    double min_allowed_dist = 1.8;
    //    double min_clash2 = 2.3;
    //    double dist_H = 2.5;
    //    double dist_Hyd = 3.2;
    double dist_Ionic = 2.5;
    //    double dist_Metal = 1.8;
    //    double dist_Arom = 3.2;
    map<Residu*,resbest> hydlist;
    map<Residu*,resbest> ::iterator itHydList;



    for (ItCAtom itLA = complex.firstAtom(); itLA != complex.lastAtom();++itLA)
    {
        Atom &atomL=**itLA;
#ifdef ICHEM_DEBUG
            cout <<"LIGAND ATOM \t"<<atomL.getIdentifier() << endl
            <<"\t MOL2T: " << atomL.getMOL2Type() << endl
            <<"\t props : "<< atomL.props.toString() << endl
            <<"\t charge: " << atomL.getFormalCharge() << endl;
            << " END" << endl;
        
#endif
        if (atomL.isHydrogen() || (atomL.getName()=="DuCy")
                || !atomL.isUsed()) continue;
        if (atomL.getBox(&grid) == (Box*)NULL) continue;

        Box& boxL = *atomL.getBox(&grid);
        boxlist.clear();
        grid.getAdjacency(boxlist,&boxL,max_allowed_dist*2,true);

        hydlist.clear();;

        for (ItBox itBAdj=boxlist.begin();
             itBAdj != boxlist.end();
             ++itBAdj)
        {
            for (ItCAtom itAA  = (*itBAdj)->firstAtom();
                 itAA != (*itBAdj)->lastAtom();
                 itAA++)
            {
                Atom &atomP = **itAA;
                if (&atomP.getParent()==&atomL.getParent()
                        || !atomP.isUsed()
                        || atomP.isHydrogen())continue;
                dist = atomL.calcFixpos(atomP,max_allowed_dist+0.1);
                if (dist > max_allowed_dist )
                {
                    continue;
                }

                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                if (atomL.props.isAcceptor())
                {
                    // LIGAND ACCEPTOR <-> PROTEIN DONOR => HBOND_LIG
                    if (wInterType[InterType::HBOND_LIG]
                            && atomP.props.isDonor()
                            && dist <= Dist_H)
                    {
                        //                        cout << "j'ai accepteur et donneur" << endl;
                        for (size_t i=0; i< atomP.getNumBond();++i)
                        {
                            const Atom &atomP2 = atomP.getAtomLinked(i);
#ifdef ICHEM_DEBUG
                            cout <<"PROT ATOM \t"<<atomP.getIdentifier() << endl
                                <<"\t MOL2T: " << atomP.getMOL2Type() << endl
                               <<"\t props : "<< atomP.props.toString() << endl
                              <<"\t charge: " << atomP.getFormalCharge() << endl
                             << " END" << endl;
#endif
                            if (!atomP2.isHydrogen())     continue;
                            angle = atomP2.fixpos.calcAngle(atomP.fixpos,  atomL.fixpos);

                            if (angle <= Angl_H-AngT_H
                                    || angle >= Angl_H+AngT_H)  { continue;}

#ifdef ICHEM_DEBUG
                            cout <<"HBLA\t"<<atomL.getIdentifier()
                                <<"\t"<< atomP.getIdentifier()
                               << "\t" << dist
                               << "\t"<< (angle*180/M_PI)
                               <<"\t"<<atomP2.getIdentifier()<<endl ;
#endif
                            InterPoint IntP(NInter,
                                            &atomP,
                                            &atomL,
                                            (atomL.fixpos+atomP.fixpos)/2,
                                            InterType::HBOND_PROT,
                                            dist,
                                            angle);
                            NInter++;
                            interResult.listInters.push_back(IntP);
                            interResult.N++;


                        }
                    }

                    // LIGAND ACCEPTOR <-> PROTEIN METAL => METAL
                    //                    cout << wInterType[InterType::METAL] << atomP.props.isMetal() << Dist_Metal << endl;
                    if (wInterType[InterType::METAL]
                            && atomP.props.isMetal()
                            && dist <= Dist_Metal)
                    {

                        InterPoint IntP(NInter,&atomP,&atomL,(atomL.fixpos+atomP.fixpos)/2, InterType::METAL, dist);
                        NInter++;
                        interResult.listInters.push_back(IntP);
                        interResult.Zn++;
#ifdef ICHEM_DEBUG
                        cout << "MEAC\t"<<atomL.getIdentifier()<<"\t"<<atomP.getIdentifier()<<" \tDist:"<<dist<<endl;
#endif
                    }

                }// END isAcceptor

                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////

                if (wInterType[InterType::IONIC_PROT]
                        && atomL.props.isAnion()
                        && atomP.props.isCation()
                        && !atomP.props.isMetal()
                        && dist <= Dist_Ionic && dist >= dist_Ionic )
                {
                    //                    cout << "un anion est proche d'un cation " << endl;

                    InterPoint IntP(NInter,
                                    &atomP,
                                    &atomL,
                                    (atomL.fixpos+atomP.fixpos)/2,
                                    InterType::IONIC_PROT,
                                    dist);
                    NInter++;
                    interResult.listInters.push_back(IntP);
                    interResult.NZ++;

#ifdef ICHEM_DEBUG
                    cout <<"ACCA\t"<<atomL.getIdentifier()
                        <<"\t"<< atomP.getIdentifier()
                       << "  \tDist=" << dist<<endl  ;
#endif
                }


                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////

                //                if (atomL.props.isCation() && !atomL.props.isMetal())
                //                {
                //                    if (wInterType[InterType::IONIC_LIG]
                //                            && atomP.props.isAnion()
                //                            && dist <= Dist_Ionic && dist >= dist_Ionic)
                //                    {
                //                        InterPoint IntP(NInter,
                //                                        &atomP,
                //                                        &atomL,
                //                                        (atomL.fixpos+atomP.fixpos)/2,
                //                                        InterType::IONIC_LIG,
                //                                        dist);

                //                        NInter++;
                //                        interResult.listInters.push_back(IntP);
                //                        interResult.OD1++;
                //#ifdef ICHEM_DEBUG
                //                        cout <<"CAAC\t"<<atomL.getIdentifier()
                //                            <<"\t"<< atomP.getIdentifier()
                //                           << "  \tDist=" << dist<<endl  ;
                //#endif
                //                    }


                //                }
                // SPECIAL CASES :

                // N linked to sulphonamide in interaction with Metal

                if (wInterType[InterType::METAL]&&
                        atomL.isNitrogen()
                        &&atomP.props.isMetal()
                        &&dist < Dist_Metal)
                {
                    bool linkedtosulf=false;
                    for (size_t iAtmLink=0;iAtmLink < atomL.getNumBond();++iAtmLink)
                    {
                        const Atom &atml = atomL.getAtomLinked(iAtmLink);
                        if (atml.isSulfur()) linkedtosulf=true;
                    }
                    if (!linkedtosulf) continue;
                    InterPoint IntP(NInter,
                                    &atomP,
                                    &atomL,
                                    (atomL.fixpos+atomP.fixpos)/2,
                                    InterType::METAL,
                                    dist);
                    NInter++;
                    interResult.listInters.push_back(IntP);
                    interResult.Zn++;


#ifdef ICHEM_DEBUG
                    cout << "MEAC\t"<<atomL.getIdentifier()<<"\t"<<atomP.getIdentifier()<<" \tDist:"<<dist<<endl;
#endif

                }
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////

                if (wInterType[InterType::HYDROPHOBIC]
                        && atomL.props.isHydrophobic() && atomP.props.isHydrophobic()
                        && !(atomL.props.isAromatic() && atomP.props.isAromatic())
                        && dist <= Dist_Hyd)
                {

                    itHydList=hydlist.find(atomP.getResidu());
                    if (itHydList == hydlist.end())
                    {
                        resbest rbest;rbest.atmP=&atomP;rbest.dist=dist;
                        hydlist.insert(pair<Residu*,resbest>(atomP.getResidu(),rbest));
                    }
                    else
                    {
                        resbest &rbest = (*itHydList).second;
                        if ( dist < rbest.dist)
                        {
                            rbest.dist=dist;
                            rbest.atmP=&atomP;
                        }
                    }

                }


            }// END itAA
        }// END itBADj

        for (itHydList = hydlist.begin(); itHydList != hydlist.end(); ++itHydList)
        {

            InterPoint IntP(NInter,
                            (*itHydList).second.atmP,
                            &atomL,
                            (atomL.fixpos+(*itHydList).second.atmP->fixpos)/2,
                            InterType::HYDROPHOBIC,
                            (*itHydList).second.dist);
            NInter++;
            interResult.listInters.push_back(IntP);
            interResult.CA++;
#ifdef ICHEM_DEBUG
            cout << "HYDR\t"<<atomL.getIdentifier()
                 <<"\t"<< (*itHydList).second.atmP->getIdentifier()
                <<"\t"<< (*itHydList).second.dist<<endl;
#endif
        }

    }// END itLA*/

    if (!wInterType[InterType::AREDGEFACE]
            &&!wInterType[InterType::ARFACEFACE]
            &&!wInterType[InterType::PICATION])
    {
        if (wMerge) mergeInteractions(interResult);
        return;
    }

    AtomList toRotate;
#ifdef ICHEM_DEBUG
    cout << "AROMATIC INTERACTIONS"<<endl;
#endif
    CycleList cyclelist;
    bool possible =true;
    unsigned short EF=0;bool Ar=false;double distCENTER;



    for (ItCMole itMole = complex.firstMole() ; itMole!=complex.lastMole();++itMole)
    {
        Molecule& chain = **itMole;

        if (chain.getName()=="XX" ){continue;}
        if (!chain.getAtom(0).isUsed()){continue;}
        if (!chain.getResidu(0).isUsed()){continue;}

        for (ItCCycle itLC= chain.firstCycle();itLC!= chain.lastCycle(); itLC++)
        {
            // Getting ligand cycle :
            Cycle& ligcycle=**itLC;

            // Not aromatic : continue;
            if (!ligcycle.isAromatic()) continue;
            for (ItCAtom itCLA = ligcycle.first();
                 itCLA!= ligcycle.end();
                 itCLA++)
            {
                Atom &atmC = **itCLA;
                if (!atmC.getResidu()->isUsed()){continue;}
                if (!atmC.getResidu()->getAtom(0).isUsed()){continue;}
            }

            // Getting atom center :
            Atom &atomL=ligcycle.getCenter();
#ifdef ICHEM_DEBUG
            cout << "LOOKING AT : "<<atomL.getIdentifier()
                 <<"::"<<atomL.fixpos.toString()<<endl;
#endif
            // Checking that atom center is located in the cube :
            if (atomL.getBox(&grid) == (Box*)NULL)
            {
                // When not, we find its cube :
                toRotate.clear();
                toRotate.push_back(&atomL);
                grid.rotateAtoms(toRotate);
                // When we cannot find it, we ignore the cycle :
                if (atomL.getBox(&grid) == (Box*)NULL)
                    continue;
            }

            // Fetching cube :
            Box& boxL = *atomL.getBox(&grid);

            ligcycle.calcVector();
            boxlist.clear();
            // Getting adjacent cube of this cube :
            grid.getAdjacency(boxlist,&boxL,max_allowed_dist*1.5,true);


            // Scanning every atom of every adjacent cube looking for Aromatic center:
            for (ItBox   itBAdj  = boxlist.begin();
                 itBAdj != boxlist.end();
                 ++itBAdj)
                for (ItCAtom itAA    = (*itBAdj)->firstAtom();
                     itAA   != (*itBAdj)->lastAtom();
                     itAA++)
                {
                    Atom &atomP = **itAA;
                    //                    cout << atomP.getIdentifier() <<" :: " << atomP.fixpos.toString() << endl;
                    if (&atomP.getParent()==&atomL.getParent()) continue;
                    if (atomP.props.isCation())
                    {
                        dist = atomL.calcFixpos(atomP,4.1);
                        if (dist > 4)continue;// Cannot be aromatic either
                        const double angle=ligcycle.getCenter().fixpos.calcAngle(ligcycle.getNormVector(),atomP.fixpos);
                        if ((angle >= Angl_PiCation-AngT_PiCation
                                && angle <= Angl_PiCation+AngT_PiCation)
                                ||(angle >= Angl_PiCation-AngT_PiCation-PI
                                   && angle <= Angl_PiCation+AngT_PiCation-PI)
                                ){
                        InterPoint IntP(NInter,
                                        &atomP,
                                        &ligcycle.getCenter(),
                                        (atomP.fixpos+ligcycle.getCenter().fixpos)/2,
                                        InterType::PICATION,dist,angle);NInter++;
                        interResult.listInters.push_back(IntP);
}

                    }


                    if (atomP.getName()!="DuAr")continue;
                    dist = atomL.calcFixpos(atomP,10.1);
                    // Since we have both cycle centers here, we can check if they are
                    // not too far away :
                    if (dist > 10) continue;

                    // Getting the corresponding cycle from this atom center :
                    Cycle* cycleP=atomP.getParent().getCycleFromCenter(&atomP);
                    if (cycleP==(Cycle*)NULL)continue;
                    // Saving it :
                    cyclelist.push_back(cycleP);
                }

            // Sort and unique to avoid cycle redundance :
            sortAndUnique(cyclelist);


            for (ItCCycle itC = cyclelist.begin()
                 ; itC != cyclelist.end();itC++)
            {
                Cycle &cycleP = **itC;
                distCENTER = cycleP.getFixpos().calcDist(ligcycle.getFixpos());
#ifdef ICHEM_DEBUG
                cout << "  ||-->"<<cycleP.getCenter().fixpos.toString()
                     <<" " << ligcycle.getCenter().fixpos.toString()<<" " << distCENTER<<endl;
#endif
                if (!cycleP.getAtom(0)->isUsed())  continue;
                // CASE AROMATIC CENTER TOO FAR => POSSIBLE HYDROPHOBIC INTERACTION :
                if (wInterType[InterType::HYDROPHOBIC] && distCENTER > Dist_Arom)
                {
                    double bestdist=100;Atom* Hlig=(Atom*)NULL, *Hprot=(Atom*)NULL;
                    // SCANNING LIGAND SIDE ATOMS IN CYCLE :
                    for (ItCAtom itCLA = ligcycle.first();
                         itCLA!= ligcycle.end();
                         itCLA++)
                    {
                        Atom& atomCLA = **itCLA;
                        if (!atomCLA.props.isHydrophobic())continue;

                        // SCANNING PROTEIN SIDE ATOMS IN CYCLE :
                        for (ItCAtom itCLP = cycleP.first();
                             itCLP!= cycleP.end();
                             itCLP++)
                        {
                            Atom& atomCLP = **itCLP;
                            if (!atomCLP.props.isHydrophobic())continue;

                            // GETTING CLOSEST ATOMS BETWEEN PROT/LIGAND
                            dist = atomCLA.calcFixpos(atomCLP);

                            if (dist < bestdist)
                            {
                                Hlig=&atomCLA;
                                Hprot=&atomCLP;
                                bestdist=dist;
                            }
                        }

                    }
                    // BELOW Dist_H : HYDROPHOBIC INTERACTION
                    if (bestdist >= Dist_Hyd)continue;
#ifdef ICHEM_DEBUG
                    cout << "HYD"<<endl;
#endif
                    InterPoint IntP(NInter,
                                    Hprot,
                                    Hlig,
                                    (Hprot->fixpos+Hlig->fixpos)/2,
                                    InterType::HYDROPHOBIC,
                                    bestdist);
                    NInter++;
                    interResult.listInters.push_back(IntP);
                    interResult.CA++;
                    continue;
                }

                // CASE AROMATIC CENTER CLOSE :
                cycleP.calcVector();
                possible = true;
                EF=0; Ar = false;
                // SCANNING LIGAND ATOM VS PROTEIN ATOM :
                for (size_t itCLA =0; itCLA < ligcycle.getNumAtom(); itCLA++)
                {

                    Atom& atomCLA =*ligcycle.getAtom(itCLA);
                    if (!atomCLA.isCarbon())continue;
                    const Coords& normalCLA = ligcycle.getNormVector(itCLA);

                    for (size_t itCLP =0; itCLP < cycleP.getNumAtom(); itCLP++)

                    {
                        Atom& atomCLP = *cycleP.getAtom(itCLP);
                        if (!atomCLP.isCarbon())continue;
                        dist = atomCLA.calcFixpos(atomCLP);

                        if (dist > Dist_Arom*3){possible=false;break;}
                        const Coords& normalCLP=cycleP.getNormVector(itCLP);
                        const double angle =
                                acos((atomCLA.fixpos.x-normalCLA.x)*(atomCLP.fixpos.x-normalCLP.x)
                                     +(atomCLA.fixpos.y-normalCLA.y)*(atomCLP.fixpos.y-normalCLP.y)
                                     +(atomCLA.fixpos.z-normalCLA.z)*(atomCLP.fixpos.z-normalCLP.z)
                                     );
                        if (dist < Dist_Arom
                                && ((angle >      Angl_AromFF-AngT_AromFF
                                     && angle <      Angl_AromFF+AngT_AromFF)
                                    ||  (angle > M_PI+Angl_AromFF-AngT_AromFF
                                         && angle < M_PI+Angl_AromFF+AngT_AromFF)
                                    ||  (angle >      Angl_AromFF-AngT_AromFF-M_PI
                                         && angle <      Angl_AromFF+AngT_AromFF-M_PI)))  Ar=true;
                        if (dist < Dist_Arom
                                && ((angle >      Angl_AromEF-AngT_AromEF
                                     && angle <      Angl_AromEF+AngT_AromEF)
                                    || (angle > M_PI+Angl_AromEF-AngT_AromEF
                                        && angle < M_PI+Angl_AromEF+AngT_AromEF)
                                    || (angle >      Angl_AromEF-AngT_AromEF-M_PI
                                        && angle <      Angl_AromEF+AngT_AromEF-M_PI)))  EF++;
                    }
                    if (!possible)break;
                }
                if (wInterType[InterType::ARFACEFACE] && Ar)
                {
                    InterPoint IntP(NInter,
                                    &cycleP.getCenter(),
                                    &ligcycle.getCenter(),
                                    (cycleP.getCenter().fixpos+ligcycle.getCenter().fixpos)/2,
                                    InterType::ARFACEFACE,distCENTER);NInter++;
                    interResult.listInters.push_back(IntP);

                    interResult.CZ++;
#ifdef ICHEM_DEBUG
                    cout << "|||->AROMATIC FACE FACE"<<endl;
#endif
                }
                else if (wInterType[InterType::AREDGEFACE] && EF>5){
                    InterPoint IntP(NInter,
                                    &cycleP.getCenter(),
                                    &ligcycle.getCenter(),
                                    (cycleP.getCenter().fixpos+ligcycle.getCenter().fixpos)/2,
                                    InterType::AREDGEFACE,distCENTER);NInter++;
                    interResult.listInters.push_back(IntP);
                    interResult.CZ++;
#ifdef ICHEM_DEBUG
                    cout << "|||->AROMATIC EDGE FACE"<<endl;
#endif
                }


            }

        }
    }
    if (wMerge){
        //        mergeInteractions(interResult);
        mergeSpeInts(interResult);
    }
}


/*
    Merge close proximity interactions within a molecule for optimization purposes.
    It looks for hydrophobic interactions between atoms in the molecule and combine them if they are close
*/
void Interactions::mergeInteractions(InterResults& interResult) const
{
    int NInter;
    Coords new_center;
    Atom *new_atmP = nullptr, *new_atmL = nullptr;
    bool modif = true ; double dist;
    
    while (modif)
    {
        //        cout << interResult.listInters.max_size() << endl;
        modif = false;
        // tmp vector for new interpoints
        std::vector<InterPoint> new_interpoints;
        for (size_t inter_i=0; inter_i < interResult.listInters.size();++inter_i)
        {
            InterPoint &intPi = interResult.listInters.at(inter_i);
            //            cout << interResult.listInters.size()  << " work on : " << intPi.point << endl;

            if (intPi.interaction != InterType::HYDROPHOBIC
                    ||  intPi.merged_to != -1) continue;

            const size_t sizeList=interResult.listInters.size();

            for (size_t inter_j=inter_i+1; inter_j < sizeList;++inter_j)
            {
                InterPoint &intPj = interResult.listInters.at(inter_j);

                if (intPj.interaction != InterType::HYDROPHOBIC
                        ||  intPj.merged_to != -1 || intPi.point==intPj.point) continue;
                dist = intPi.center.calcDist(intPj.center,1.1);
                if ( dist > 1 ) continue;
                //                if (interResult.listInters.size() == 2048) {
                //                    cout << "maximum size"  << endl;
                //                }

                modif=true;
                NInter = (int)interResult.listInters.size();
                new_center=(intPj.center+intPi.center)/2;

                if (new_center.calcDist(intPi.Prot_Ref->fixpos) <
                        new_center.calcDist(intPj.Prot_Ref->fixpos) )
                {
                    new_atmP=intPi.Prot_Ref;
                }else new_atmP=intPj.Prot_Ref;

                if (new_center.calcDist(intPi.Lig_Ref->fixpos) < new_center.calcDist(intPj.Lig_Ref->fixpos ))
                {
                    new_atmL=intPi.Lig_Ref;
                }else new_atmL=intPj.Lig_Ref;
                InterPoint IntP(NInter,
                                new_atmP,
                                new_atmL,
                                new_center,
                                InterType::HYDROPHOBIC,(intPi.dist<intPj.dist)? intPi.dist:intPj.dist);

                new_interpoints.push_back(IntP);
                if (intPi.merged_to != -1)
                {
                    for (size_t inter_k=0; inter_k < interResult.listInters.size();++inter_k)
                    {
                        InterPoint &intPk = interResult.listInters.at(inter_k);
                        if (intPk.merged_to == intPi.merged_to) intPk.merged_to=NInter;

                    }
                }
                intPi.merged_to=NInter;
                if (intPj.merged_to != -1)
                {
                    for (size_t inter_k=0; inter_k < interResult.listInters.size();++inter_k)
                    {
                        InterPoint &intPk = interResult.listInters.at(inter_k);
                        if (intPk.merged_to == intPj.merged_to) intPk.merged_to=NInter;

                    }
                }
                intPj.merged_to=NInter;
                NInter++;

            } //END inter_j
        } //END inter_i
        interResult.listInters.insert(interResult.listInters.end(), new_interpoints.begin(), new_interpoints.end());
    } //END WHILE
}



    //    double bestdist;InterPoint *toMove;
    //    for (size_t inter_i=0; inter_i < interResult.listInters.size();++inter_i)
    //    {
    //        InterPoint &intPi = interResult.listInters.at(inter_i);
    //        bestdist=100;
    //        if (intPi.interaction != InterType::IONIC_PROT
    //                &&intPi.interaction != InterType::IONIC_LIG) continue;

    //        for (size_t inter_j=inter_i+1; inter_j < interResult.listInters.size();++inter_j)
    //        {
    //            InterPoint &intPj = interResult.listInters.at(inter_j);
    //            if ((intPi.interaction == InterType::IONIC_PROT
    //                 && intPj.interaction != InterType::HBOND_PROT)
    //                    ||(intPi.interaction == InterType::IONIC_LIG
    //                       && intPj.interaction != InterType::HBOND_LIG) )continue;

    //            dist = intPi.center.calcDist(intPj.center,1.1);
    //            if (dist < bestdist) {bestdist = dist;toMove=&intPj;}


    //        }
    //        if (bestdist <1)
    //        {

    //#ifdef ICHEM_DEBUG
    //                cout << "#########################################"<<endl;
    //                if (intPi.interaction == InterType::IONIC_PROT)
    //                    cout << "MOVING IONIC PROT CATION TO H-BOND PROT DONOR POSITION"<<endl;
    //                else
    //                    cout << "MOVING IONIC LIG CATION TO H-BOND LIG DONOR POSITION"<<endl;
    //                cout << "Ionic pos: "
    //                     <<intPi.Lig_Ref->getIdentifier()
    //                    << "  "
    //                    << intPi.Prot_Ref->getIdentifier()
    //                    <<  intPi.center.toString()<<endl
    //                     << "HBond pos: "
    //                     <<toMove->Lig_Ref->getIdentifier()
    //                    << "  "
    //                    << toMove->Prot_Ref->getIdentifier()
    //                    << toMove->center.toString()<<endl
    //                       ;
    //#endif
    //            intPi.center = toMove->center;

    //        }
    //    }

// }

void Interactions::mergeSpeInts(InterResults& interResult)const
{
    int NInter;
    Coords new_center;
    Atom *new_atmP, *new_atmL;
    bool modif=true;double dist;

    while (modif)
    {
        modif =false;
        for (size_t inter_i=0; inter_i < interResult.listInters.size();++inter_i)
        {
            InterPoint &intPi = interResult.listInters.at(inter_i);
            if (intPi.interaction != InterType::HYDROPHOBIC
                    ||  intPi.merged_to != -1)continue;

            const size_t sizeList=interResult.listInters.size();
            for (size_t inter_j=inter_i+1; inter_j < sizeList;++inter_j)
            {
                InterPoint &intPj = interResult.listInters.at(inter_j);
                if (intPj.interaction != InterType::HYDROPHOBIC
                        ||  intPj.merged_to != -1 || intPi.point==intPj.point)continue;
                dist =intPi.center.calcDist(intPj.center,1.1);
                if (dist > 0.1 )continue;


                modif=true;
                NInter = (int)interResult.listInters.size();
                new_center=(intPj.center+intPi.center)/2;

                if (new_center.calcDist(intPi.Prot_Ref->fixpos) <
                        new_center.calcDist(intPj.Prot_Ref->fixpos) )
                {
                    new_atmP=intPi.Prot_Ref;
                }else new_atmP=intPj.Prot_Ref;

                if (new_center.calcDist(intPi.Lig_Ref->fixpos) < new_center.calcDist(intPj.Lig_Ref->fixpos ))
                {
                    new_atmL=intPi.Lig_Ref;
                }else new_atmL=intPj.Lig_Ref;
                InterPoint IntP(NInter,
                                new_atmP,
                                new_atmL,
                                new_center,
                                InterType::HYDROPHOBIC,(intPi.dist<intPj.dist)? intPi.dist:intPj.dist);

                interResult.listInters.push_back(IntP);
                if (intPi.merged_to != -1)
                {
                    for (size_t inter_k=0; inter_k < interResult.listInters.size();++inter_k)
                    {
                        InterPoint &intPk = interResult.listInters.at(inter_k);
                        if (intPk.merged_to == intPi.merged_to) intPk.merged_to=NInter;

                    }
                }
                intPi.merged_to=NInter;
                if (intPj.merged_to != -1)
                {
                    for (size_t inter_k=0; inter_k < interResult.listInters.size();++inter_k)
                    {
                        InterPoint &intPk = interResult.listInters.at(inter_k);
                        if (intPk.merged_to == intPj.merged_to) intPk.merged_to=NInter;

                    }
                }
                intPj.merged_to=NInter;
                NInter++;

            }//END inter_j
        }//END inter_i
    }//END WHILE

    double bestdist;InterPoint *toMove;
    for (size_t inter_i=0; inter_i < interResult.listInters.size();++inter_i)
    {
        InterPoint &intPi = interResult.listInters.at(inter_i);
        bestdist=100;
        if (intPi.interaction != InterType::IONIC_PROT
                &&intPi.interaction != InterType::IONIC_LIG) continue;

        for (size_t inter_j=inter_i+1; inter_j < interResult.listInters.size();++inter_j)
        {
            InterPoint &intPj = interResult.listInters.at(inter_j);
            if ((intPi.interaction == InterType::IONIC_PROT
                 && intPj.interaction != InterType::HBOND_PROT)
                    ||(intPi.interaction == InterType::IONIC_LIG
                       && intPj.interaction != InterType::HBOND_LIG) )continue;

            dist = intPi.center.calcDist(intPj.center,1.1);
            if (dist < bestdist) {bestdist = dist;toMove=&intPj;}


        }
        if (bestdist <1)
        {

#ifdef ICHEM_DEBUG
            cout << "#########################################"<<endl;
            if (intPi.interaction == InterType::IONIC_PROT)
                cout << "MOVING IONIC PROT CATION TO H-BOND PROT DONOR POSITION"<<endl;
            else
                cout << "MOVING IONIC LIG CATION TO H-BOND LIG DONOR POSITION"<<endl;
            cout << "Ionic pos: "
                 <<intPi.Lig_Ref->getIdentifier()
                << "  "
                << intPi.Prot_Ref->getIdentifier()
                <<  intPi.center.toString()<<endl
                 << "HBond pos: "
                 <<toMove->Lig_Ref->getIdentifier()
                << "  "
                << toMove->Prot_Ref->getIdentifier()
                << toMove->center.toString()<<endl
                   ;
#endif
            intPi.center = toMove->center;

        }
    }

}

// std::string Interactions::toString(const  InterResults& interResult)const
// {
//     ostringstream oss;

//     for (std::vector<InterPoint>::const_iterator it=interResult.listInters.begin();
//          it != interResult.listInters.end();it++)
//     {
//         const InterPoint& interpt=*it;
//         if (interpt.merged_to != -1)continue;
//         switch (interpt.interaction)
//         {
//         case InterType::HBOND_LIG:   oss << "HBond_LIG         ";break;
//         case InterType::HBOND_PROT:  oss << "HBond_PROT        ";break;
//         case InterType::IONIC_LIG:   oss << "Ionic_LIG         ";break;
//         case InterType::IONIC_PROT:  oss << "Ionic_PROT        ";break;
//         case InterType::METAL:       oss << "Metal_Acceptor    ";break;
//         case InterType::HYDROPHOBIC: oss << "Hydrophobic       ";break;
//         case InterType::ARFACEFACE:  oss << "Aromatic_Face/Face";break;
//         case InterType::AREDGEFACE:  oss << "Aromatic_Edge/Face";break;
//         case InterType::PICATION:    oss << "Pi/Cation         ";break;
//         case InterType::WHBOND_LIG:  oss << "Weak_HBond_LIG    ";break;
//         case InterType::WHBOND_PROT: oss << "Weak_HBond_PROT   ";break;
//         }
// //        if interpt.interaction != InterType::PICATION
//         if (interpt.interaction != InterType::AREDGEFACE
//                 &&interpt.interaction != InterType::ARFACEFACE
//                 &&interpt.interaction != InterType::PICATION)
//         {
//             oss << "\t"; oss.setf(ios::left);
//             oss << "|" << interpt.Prot_Ref->getIdentifier();
//             oss << "\t|" << interpt.Lig_Ref->getIdentifier();
//             oss << "\t|" << interpt.point;
//             oss << "\t|"; oss.width(7); oss.setf(ios::right); oss<< interpt.dist;
//             oss << "\t|";
//             if (interpt.angle != -100000)
//                 oss <<(interpt.angle*180/M_PI);
//             else oss<<"/";
//         }
//         else if (interpt.interaction== InterType::PICATION)
//         {
//             if (interpt.Prot_Ref->getName()=="DuAr")
//             {
//                 Cycle* cyc=interpt.Prot_Ref->getParent().getCycleFromCenter(interpt.Prot_Ref);

//                 oss << "\t"; oss.setf(ios::left);
//                 oss << "|"; oss.width(5);oss.setf(ios::left); oss<< cyc->getCenter().getName();
//                 oss << "|"; oss.width(5);oss.setf(ios::left); oss << "NULL";
//                 oss << "|"; oss << cyc->getAtom(0)->getResidu()->getIdentifier();

//                 oss << "\t|" << interpt.Lig_Ref->getIdentifier();
//                 oss << "\t|" << interpt.point;
//                 oss << "\t|"; oss.width(7); oss.setf(ios::right); oss << interpt.dist;
//                 oss << "\t|";
//                 if (interpt.angle != -100000)
//                     oss <<(interpt.angle*180/M_PI);
//                 else oss<<"/";


// //                oss<<"|";oss.width(7);oss<<"|";
// //                oss.width(5);oss.setf(ios::left); oss<< cyc->getCenter().getName();
// //                oss<<"|";
// //                oss.width(5); oss<< "NULL";
// //                oss <<"|";oss.width(9);oss<<cyc->getAtom(0)->getResidu()->getIdentifier();
// //                oss<< "|";oss.width(7);oss
// //                   <<interpt.Lig_Ref->getIdentifier();
// //                oss<< "|";oss.width(7);oss<<interpt.dist
// //                 <<oss.width(7);
// //                if (interpt.angle != -100000)
// //                    oss   <<(interpt.angle*180/M_PI);
// //                else oss<<"/";
//             }
//             else
//             {
//                 oss << "\t"; oss.setf(ios::left);
//                 oss << "|" << interpt.Prot_Ref->getIdentifier();
//                 Cycle* cyc = interpt.Lig_Ref->getParent().getCycleFromCenter(interpt.Lig_Ref);

//                 oss << "\t";
//                 oss << "|"; oss.setf(ios::left, ios::adjustfield); oss.width(5); oss << cyc->getCenter().getName(); // Without adjustfield, it does not work
//                 oss << "|"; oss.setf(ios::left); oss.width(5); oss << "0";
//                 oss << "|"; oss << "UNK1001-XX";
//                 oss << "\t|" << interpt.point;
//                 oss << "\t|"; oss.width(7); oss.setf(ios::right); oss<< interpt.dist;
//                 oss << "\t|";
//                 oss <<(interpt.angle*180/M_PI);

//             }
//         }
//         else
//         {   
//             // here
//             Cycle* cyc=interpt.Prot_Ref->getParent().getCycleFromCenter(interpt.Prot_Ref);
//             oss << "\t";
//             oss << "|"; oss.setf(ios::left, ios::adjustfield); oss.width(5); oss << cyc->getCenter().getName(); // Without adjustfield, it does not work
//             oss << "|"; oss.setf(ios::left); oss.width(5); oss << "NULL";
//             oss << "|"; oss << cyc->getAtom(0)->getResidu()->getIdentifier();

//             cyc=interpt.Lig_Ref->getParent().getCycleFromCenter(interpt.Lig_Ref);

//             oss << "\t";
//             oss << "|"; oss.setf(ios::left, ios::adjustfield); oss.width(5); oss << cyc->getCenter().getName(); // Without adjustfield, it does not work
//             oss << "|"; oss.setf(ios::left); oss.width(5); oss << "NULL";
//             oss << "|"; oss << cyc->getAtom(0)->getResidu()->getIdentifier();

//             oss << "\t|" << interpt.point;
//             oss << "\t|"; oss.width(7); oss.setf(ios::right); oss<< interpt.dist;
//             oss << "\t|/";

//         }
//         oss <<endl;
// //        for (std::vector<InterPoint>::const_iterator it2=interResult.listInters.begin();
// //             it2 != interResult.listInters.end();it2++)
// //        {
// //            const InterPoint& interpt2=*it2;
// //            if (interpt2.merged_to!= interpt.point) continue;
// //            oss<<"\t|-->"<<interpt2.Prot_Ref->getIdentifier()
// //              <<"\t"<<interpt2.Lig_Ref->getIdentifier()
// //             <<"\t"<<interpt2.dist
// //            <<"\t";
// //            if (interpt2.angle != -100000)
// //                oss   <<(interpt2.angle*180/M_PI);
// //            else oss<<"/";
// //            oss<<"\n";
// //        }
//     }
//     return oss.str();

// }


std::string Interactions::toString(const InterResults& interResult) const
{
    ostringstream oss;

    // Stores ligand cycle atoms per interaction index
    std::map<int, std::vector<std::string>> ligandCycleAtoms; 

    for (std::vector<InterPoint>::const_iterator it = interResult.listInters.begin();
         it != interResult.listInters.end(); it++)
    {
        const InterPoint& interpt = *it;
        if (interpt.merged_to != -1) continue;

        switch (interpt.interaction)
        {
        case InterType::HBOND_LIG:   oss << "HBond_LIG         "; break;
        case InterType::HBOND_PROT:  oss << "HBond_PROT        "; break;
        case InterType::IONIC_LIG:   oss << "Ionic_LIG         "; break;
        case InterType::IONIC_PROT:  oss << "Ionic_PROT        "; break;
        case InterType::METAL:       oss << "Metal_Acceptor    "; break;
        case InterType::HYDROPHOBIC: oss << "Hydrophobic       "; break;
        case InterType::ARFACEFACE:  oss << "Aromatic_Face/Face"; break;
        case InterType::AREDGEFACE:  oss << "Aromatic_Edge/Face"; break;
        case InterType::PICATION:    oss << "Pi/Cation         "; break;
        case InterType::WHBOND_LIG:  oss << "Weak_HBond_LIG    "; break;
        case InterType::WHBOND_PROT: oss << "Weak_HBond_PROT   "; break;
        }

        if (interpt.interaction != InterType::AREDGEFACE
            && interpt.interaction != InterType::ARFACEFACE
            && interpt.interaction != InterType::PICATION)
        {
            oss << "\t"; oss.setf(ios::left);
            oss << "|" << interpt.Prot_Ref->getIdentifier();
            oss << "\t|" << interpt.Lig_Ref->getIdentifier();
            oss << "\t|" << interpt.point;
            oss << "\t|"; oss.width(7); oss.setf(ios::right); oss << interpt.dist;
            oss << "\t|";
            if (interpt.angle != -100000)
                oss << (interpt.angle * 180 / M_PI);
            else
                oss << "/";
        }
        else if (interpt.interaction == InterType::PICATION)
        {
            if (interpt.Prot_Ref->getName() == "DuAr")
            {
                Cycle* cyc = interpt.Prot_Ref->getParent().getCycleFromCenter(interpt.Prot_Ref);
                oss << "\t"; oss.setf(ios::left);
                oss << "|"; oss.width(5); oss.setf(ios::left); oss << cyc->getCenter().getName();
                oss << "|"; oss.width(5); oss.setf(ios::left); oss << "NULL";
                oss << "|"; oss << cyc->getAtom(0)->getResidu()->getIdentifier();
                oss << "\t|" << interpt.Lig_Ref->getIdentifier();
                oss << "\t|" << interpt.point;
                oss << "\t|"; oss.width(7); oss.setf(ios::right); oss << interpt.dist;
                oss << "\t|";
                if (interpt.angle != -100000)
                    oss << (interpt.angle * 180 / M_PI);
                else
                    oss << "/";
            }
            else
            {
                oss << "\t"; oss.setf(ios::left);
                oss << "|" << interpt.Prot_Ref->getIdentifier();
                Cycle* cyc = interpt.Lig_Ref->getParent().getCycleFromCenter(interpt.Lig_Ref);
                oss << "\t";
                oss << "|"; oss.setf(ios::left, ios::adjustfield); oss.width(5); oss << cyc->getCenter().getName();
                oss << "|"; oss.setf(ios::left); oss.width(5); oss << "0";
                oss << "|"; oss << "UNK1001-XX";
                oss << "\t|" << interpt.point;
                oss << "\t|"; oss.width(7); oss.setf(ios::right); oss << interpt.dist;
                oss << "\t|";
                oss << (interpt.angle * 180 / M_PI);
            }
        }
        else
        {
            Cycle* cyc = interpt.Prot_Ref->getParent().getCycleFromCenter(interpt.Prot_Ref);
            oss << "\t";
            oss << "|"; oss.setf(ios::left, ios::adjustfield); oss.width(5); oss << cyc->getCenter().getName();
            oss << "|"; oss.setf(ios::left); oss.width(5); oss << "NULL";
            oss << "|"; oss << cyc->getAtom(0)->getResidu()->getIdentifier();

            cyc = interpt.Lig_Ref->getParent().getCycleFromCenter(interpt.Lig_Ref);
            oss << "\t";
            oss << "|"; oss.setf(ios::left, ios::adjustfield); oss.width(5); oss << cyc->getCenter().getName();
            oss << "|"; oss.setf(ios::left); oss.width(5); oss << "NULL";
            oss << "|"; oss << cyc->getAtom(0)->getResidu()->getIdentifier();

            oss << "\t|" << interpt.point;
            oss << "\t|"; oss.width(7); oss.setf(ios::right); oss << interpt.dist;
            oss << "\t|/";
        }

        oss << endl;

        // Store ligand cycle atoms for ARFACEFACE and AREDGEFACE interactions mapped by interaction index
        if (interpt.interaction == InterType::ARFACEFACE || interpt.interaction == InterType::AREDGEFACE)
        {
            Cycle* cyc = interpt.Lig_Ref->getParent().getCycleFromCenter(interpt.Lig_Ref);
            if (cyc != nullptr)
            {
                for (ItCAtom itA = cyc->first(); itA != cyc->end(); ++itA)
                {
                    ligandCycleAtoms[interpt.point].push_back((*itA)->getName());
                }
            }
        }
    }

    // Ligand cycle atoms printed, grouped by interaction index
    if (!ligandCycleAtoms.empty())
    {
        oss << "\n==== Ligand atoms involved in the aromatic cycle: ====\n";
        for (const auto& entry : ligandCycleAtoms)
        {
            oss << "Interaction " << entry.first << ": ";
            for (const auto& atom : entry.second)
            {
                oss << atom << " ";
            }
            oss << endl;
        }
    }
    return oss.str();
}





void Interactions::interToMOL2(  InterResults& interResult,
                                 const bool&Ligand,
                                 const bool&Protein,
                                 const bool&Center,
                                 const bool&Merge)const
{


    if (interResult.listInters.empty())return;



    bool exists=false;
    Atom *atml;

    const std::string AtmNames[NB_INTTYPE]= {"","N"   ,"O"  ,"NZ" ,"OD1","CA","Zn","CZ","CZ","","", "Zn"};
    const std::string mol2Names[NB_INTTYPE]={"","N.am","O.2","N.4","O.co2"  ,"C.3","Zn","C.ar","C.ar","","","O.co2"};
    const std::string AtomicName[NB_INTTYPE]={"","N"   ,"O"  ,"N" ,"O","C","Zn","C","C","","","O"};
    unsigned int NAtm=1;
    Molecule& mole=interResult.Ints;
    mole.clear();
    if (Merge) mole.setMoleType(MoleType::INT_MERG);
    else if (Ligand && !Protein && !Center) mole.setMoleType(MoleType::INT_LIG);
    else if (!Ligand && Protein && !Center) mole.setMoleType(MoleType::INT_PROT);
    else if (!Ligand && !Protein && Center) mole.setMoleType(MoleType::INT_CENT);
    else mole.setMoleType(MoleType::INT_MERG);


    if (Ligand || Merge)
    {



        Residu *listRes[NB_INTTYPE]={(Residu*)NULL,

                                     &mole.addResidu("Z",1,"SEL"),
                                     &mole.addResidu("Z",2,"ALL"),
                                     &mole.addResidu("Z",3,"LYL"),
                                     &mole.addResidu("Z",4,"ASL"),
                                     &mole.addResidu("Z",5,"GLL"),
                                     &mole.addResidu("Z",7,"ZIL"),
                                     &mole.addResidu("Z",6,"PHL"),
//                                     &mole.addResidu("Z",6,"PHL"),
                                     (Residu*)NULL,(Residu*)NULL,(Residu*)NULL};
        listRes[InterType::AREDGEFACE]=listRes[InterType::ARFACEFACE];

        for (size_t i=0; i< interResult.listInters.size();++i)
        {
            InterPoint &interP=interResult.listInters.at(i);
            if (interP.merged_to != -1)continue;
            exists=false;
            for (ItCAtom itC = mole.firstAtom();itC != mole.lastAtom();itC++)
            {
                Atom &atm = **itC;
                if (atm.fixpos.calcDist(interP.Lig_Ref->fixpos) > 0.005)continue;
                if (atm.getName()== AtmNames[interP.interaction]){ exists=true;break;}
            }
            if (exists)continue;
            if (listRes[interP.interaction] == (Residu*)NULL)continue;
            atml=&mole.addAtom(AtomicName[interP.interaction],
                    interP.Lig_Ref->fixpos,
                    AtmNames[interP.interaction],
                    mol2Names[interP.interaction],
                    listRes[interP.interaction]);
            atml->setFNum(NAtm);++NAtm;
            atml->setTNum(i);
            atml->setPartialCharge(interP.Lig_Ref->getPartialCharge());

        }

    }
    if (Protein || Merge)
    {

        Residu *listRes[NB_INTTYPE]={(Residu*)NULL,

                                     &mole.addResidu("Z",1,"SEP"),
                                     &mole.addResidu("Z",2,"ALP"),
                                     &mole.addResidu("Z",3,"LYP"),
                                     &mole.addResidu("Z",4,"ASP"),
                                     &mole.addResidu("Z",5,"GLP"),
                                     &mole.addResidu("Z",7,"ZIP"),
                                     &mole.addResidu("Z",6,"PHP"),
                                     //                                     &mole.addResidu("Z",6,"PHP"),
                                     (Residu*)NULL,(Residu*)NULL,(Residu*)NULL};
        listRes[InterType::AREDGEFACE]=listRes[InterType::ARFACEFACE];

        for (size_t i=0; i< interResult.listInters.size();++i)
        {
            InterPoint &interP=interResult.listInters.at(i);
            if (interP.merged_to != -1)continue;
            exists=false;
            for (ItCAtom itC = mole.firstAtom();itC != mole.lastAtom();itC++)
            {
                Atom &atm = **itC;
                if (atm.fixpos.calcDist(interP.Prot_Ref->fixpos) > 0.005)continue;
                if (atm.getName()== AtmNames[interP.interaction]){ exists=true;break;}
            }
            if (exists)continue;
            if (listRes[interP.interaction] == (Residu*)NULL)continue;
            atml=&mole.addAtom(AtomicName[interP.interaction],
                    interP.Prot_Ref->fixpos,
                    AtmNames[interP.interaction],
                    mol2Names[interP.interaction],
                    listRes[interP.interaction]);
            atml->setFNum(NAtm);++NAtm;
            atml->setTNum(i);

        }

    }
    if (Center || Merge)
    {


        Residu *listRes[NB_INTTYPE]={(Residu*)NULL,

                                     &mole.addResidu("Z",1,"SEC"),
                                     &mole.addResidu("Z",2,"ALC"),
                                     &mole.addResidu("Z",3,"LYC"),
                                     &mole.addResidu("Z",4,"ASC"),
                                     &mole.addResidu("Z",5,"GLC"),
                                     &mole.addResidu("Z",7,"ZIC"),
                                     &mole.addResidu("Z",6,"PHC"),
                                     //                                     &mole.addResidu("Z",6,"PHC"),
                                     (Residu*)NULL,(Residu*)NULL,(Residu*)NULL};
        listRes[InterType::AREDGEFACE]=listRes[InterType::ARFACEFACE];

        for (size_t i=0; i< interResult.listInters.size();++i)
        {
            exists = false;
            InterPoint &interP=interResult.listInters.at(i);
            if (interP.merged_to != -1)continue;
            if (listRes[interP.interaction] == (Residu*)NULL)continue;
            for (ItCAtom itC = mole.firstAtom();itC != mole.lastAtom();itC++)
            {
                Atom &atm = **itC;
                if (atm.fixpos.calcDist(interP.center) > 0.005)continue;
                if (atm.getName()== AtmNames[interP.interaction]){
                    exists=true;
//                    cout << "superposé Center ";
                    atm.setTNum(interP.point);
//                    for (int flo:atm.getTNum())
//                        cout << flo << " ";
//                    cout << endl;
                    break;

                }
            }
            if (exists)continue;
            atml=&mole.addAtom(AtomicName[interP.interaction],
                    interP.center,
                    AtmNames[interP.interaction],
                    mol2Names[interP.interaction],
                    listRes[interP.interaction]);
            atml->setFNum(NAtm);++NAtm;
//            atml->setTNum(i);
            atml->setTNum(interP.point);
        }
    }




}

void Interactions::calcenfouiss(string mol_name) // calcul l'enfouissement des interfaces
{

    cout << "calcul enfouiss"<< endl;
    VolSite volsite(complex);
    volsite.proj(55-20,false);
    MoleReader mread2;

    string PDB_name =     complex.getMole(0)->getName();
    string enffile = "../../ints.enfs"; //*$

    ofstream offenf;
    if (fopen(enffile.c_str(),"r")==NULL)
    {
        offenf.open(enffile.c_str(),ios::out|ios::app);
    }
    else
    {
        offenf.open(enffile.c_str(),ios::out|ios::app);
    }

    //        cout << A << "-" <<  B << " avec "  << intersize[i]->size << endl;
    //        cout << PDB_name+"_ints_"+A+"_"+B+".mol2" << endl;
    Molecule test;
    mread2.loadNewFile(mol_name+".mol2");
    mread2.loadNextMolecule(test,MoleType::LIGAND);


    //        cout << "ma proteine " << test.toString() << endl;
    //        return;

    vector<double> FGPCount;
    vector<vector<double> > FGPBuried;
    for (size_t j=0;j<8;++j)
    {
        FGPCount.push_back(0);
        if (j==0)continue;
        vector<double> tmp;
        for (size_t j=0;j<10;++j)tmp.push_back(0); FGPBuried.push_back(tmp);
    }
    for (size_t iAtm=0; iAtm < test.numAtom();++iAtm)
    {
        const Atom& atm= test.getAtom(iAtm);
        Box* box = volsite.getGrid().getBox(atm.fixpos);
        //cout << atm.getName()<<"\t"<<box->getId()<< " " << atm.fixpos.calcDist(box->fixpos)<< " " << volsite.getProjValue(*box)<<endl;
        double pos= floor((volsite.getProjValue(*box)-30)/10);
        if (offenf.is_open())
        {
            offenf  << atm.getName()<< " " << volsite.getProjValue(*box)  << endl;
        }

        //            cout << volsite.getProjValue(*box) << "\t col : " << pos+1 << endl;
        //pas de distinction protein-lignadd
        //passer tous les <5 a <7
        if (atm.getName()=="CA")      {FGPBuried[0][pos]+=1;FGPCount[1]++;FGPCount[0]++;}
        else if (atm.getName()=="O")  {FGPBuried[1][pos]+=1;FGPCount[2]++;FGPCount[0]++;}
        else if (atm.getName()=="OD1"){FGPBuried[2][pos]+=1;FGPCount[3]++;FGPCount[0]++;}
        else if (atm.getName()=="N")  {FGPBuried[1][pos]+=1;FGPCount[2]++;FGPCount[0]++;}
        else if (atm.getName()=="NZ") {FGPBuried[2][pos]+=1;FGPCount[3]++;FGPCount[0]++;}
        else if (atm.getName()=="CZ") {FGPBuried[3][pos]+=1;FGPCount[4]++;FGPCount[0]++;}
        else if (atm.getName()=="Zn") {FGPBuried[4][pos]+=1;FGPCount[5]++;FGPCount[0]++;}
    }
    for (size_t i=0; i<4;++i)
    {

        for (size_t j=0;j<10;++j)
        {
            if (FGPBuried[i][j]==0)continue;
            //                    cout << i << ":" << j << "valeur " << FGPBuried[i][j]<< endl;
            FGPBuried[i][j]=FGPBuried[i][j];

        }
        FGPCount[i+1]=FGPCount[i+1]/FGPCount[0]*100.0;

    }

    //    cout <<complex.ge<<"\t";
    ofstream ofscore;
    ofscore.open("../../descriptor.sre", ios::out|ios::app);
    for (size_t i=0;i<=4;++i) {
        cout << FGPCount[i]<< " ";
        ofscore << FGPCount[i]<< " ";
    }
    for (size_t i=0; i<4;++i) // 4 sans metal 5 avec metal
    {

        for (size_t j=0;j<10;++j)
        {

            cout << FGPBuried[i][j]<< " ";
            ofscore <<FGPBuried[i][j]<< " ";
        }

    }


    cout <<"#"<< PDB_name << endl;

    ofscore <<"#"<< PDB_name<< endl;
    ofscore.close();
    offenf.close();

    return;
}


void Interactions::loadTripletsComb(const bool full)
{
    ostringstream oss;
    ifstream Finput;
    oss << ICMole::get_IChem_LIB_Path()<<"/datas/PH4_list_6_7"<< ((full)? "F":"M")<<".csv";
    Finput.open(oss.str().c_str(),ios::in);
    if (!Finput.is_open()) {cerr << "ERROR WHILE OPENING " << oss.str()<<endl;return;}


    const unsigned int _pos=7;
    const unsigned int dist_size=6;
    map<char,int> Prop_To_Int;
    string ligne;
    Prop_To_Int.insert(pair<char,int>('A',InterType::HBOND_PROT));
    Prop_To_Int.insert(pair<char,int>('B',InterType::HBOND_LIG));
    Prop_To_Int.insert(pair<char,int>('C',InterType::IONIC_PROT));
    Prop_To_Int.insert(pair<char,int>('D',InterType::IONIC_LIG));
    Prop_To_Int.insert(pair<char,int>('E',InterType::AREDGEFACE));
    Prop_To_Int.insert(pair<char,int>('F',InterType::HYDROPHOBIC));
    Prop_To_Int.insert(pair<char,int>('G',InterType::METAL));

//    HBOND_PROT=1,  /*!< HBond Protein donor        */
//    HBOND_LIG=2,  /*!< HBond Ligand donor         */
//    IONIC_PROT=3,  /*!< Ionic Protein cation       */
//    IONIC_LIG=4,  /*!< Ionic Ligand cation        */
//    HYDROPHOBIC=5,  /*!< Hydrophobic                */
//    METAL=6,  /*!< Metal acceptor interaction */
//    ARFACEFACE=7,  /*!< Aromatic face to face      */
//    AREDGEFACE=8,  /*!< Aromatic edge to face      */
//    PICATION=9,  /*!< PI cation                  */
//    WHBOND_LIG=10,
//    WHBOND_PROT=11,

    path_dist = new double[7];
    path_dist[0]=0;
    path_dist[1]=4;
    path_dist[2]=6;
    path_dist[3]=9;
    path_dist[4]=13;
    path_dist[5]=17;
    path_dist[6]=32;

    path_dist_size=6;




    //  if (debug){
    //      cout << "DEBUG ******************************************************"<<endl
    //           << "DEBUG ********** LOADING TRIPLETS COMBINAISONS *************"<<endl
    //           << "DEBUG ******************************************************"<<endl;
    //      cout << oss.str() <<endl;
    //      cout << "DEBUG Number of dist_size : " << dist_size <<endl
    //           << "DEBUG Number of interaction groups : " << _pos <<endl;
    //    }

    for (unsigned int i=0; i<1000; i++)
        for (unsigned int j=0; j<1000; j++) vect_list[i][j]=-1;


    std::getline( Finput, ligne ); size_triplets=atoi(ligne.c_str());

    //  if (debug)    cout << "Number of possible triplets : "<< (size_triplets+1) <<endl;
    int ni=0,ni2=0,pos=0;

    while (!Finput.eof())
    {
        std::getline( Finput, ligne );
        if (ligne.length() < 8) continue; // Idk why this line is here. But without it, "core dumped" issues occurred
        ni = atoi(ligne.substr(0,3).c_str());
        oss.str("");
        oss << Prop_To_Int[ligne[4]]<<Prop_To_Int[ligne[5]]<<Prop_To_Int[ligne[6]];
        ni2 = atoi(oss.str().c_str());
        pos = atoi(ligne.substr(8).c_str());
        vect_list[ni][ni2]=pos;
        //      if (debug)cout << ni << " " << ni2 << " " << pos <<endl;
        if (ligne[4] == 'E')
        {
            oss.str("");
            oss << InterType::ARFACEFACE<<Prop_To_Int[ligne[5]]<<Prop_To_Int[ligne[6]];
            ni2 = atoi(oss.str().c_str());
            vect_list[ni][ni2]=pos;
            //          if (debug)cout << ni << " " << ni2 << " " << pos <<endl;
        }
        if (ligne[5] == 'E')
        {
            oss.str("");
            oss << Prop_To_Int[ligne[4]]<<InterType::ARFACEFACE<<Prop_To_Int[ligne[6]];
            ni2 = atoi(oss.str().c_str());
            vect_list[ni][ni2]=pos;
            //          if (debug)       cout << ni << " " << ni2 << " " << pos <<endl;
        }
        if (ligne[6] == 'E')
        {
            oss.str("");
            oss << Prop_To_Int[ligne[4]]<<Prop_To_Int[ligne[5]]<<InterType::ARFACEFACE;
            ni2 = atoi(oss.str().c_str());
            vect_list[ni][ni2]=pos;
            //          if (debug)     cout << ni << " " << ni2 << " " << pos <<endl;
        }
    }
    Finput.close();
    Interactions::Load_triplet=true;


}

Fingerprint& Interactions::generateTriplets(InterResults& interResult,const bool full)
{




    if (!Load_triplet) loadTripletsComb(full);

    Fingerprint *fgp = new Fingerprint(size_triplets+1,true);
    double dist1, dist2, dist3;
    double dist11, dist12, dist13;
    const unsigned int size = interResult.listInters.size();
    if (size < 3) { return *fgp;}
    double TableDist[size][size];
    unsigned int IV1=0, IV2=0,IV3=0;

    // Compute the pairwise IPA distances
    for (std::vector<InterPoint>::const_iterator pt1=interResult.listInters.begin();pt1 != interResult.listInters.end();pt1++){
        //  for (ItPInter pt1 = interResult.listInters.begin(); pt1 != interResult.listInters.end(); pt1++){
        IV2=IV1;
        for (std::vector<InterPoint>::const_iterator pt2=pt1 ; pt2 != interResult.listInters.end(); pt2++)	{

            dist1 = getCenter(interResult,pt1->point).calcDist(getCenter(interResult,pt2->point));
            TableDist[IV1][IV2]=dist1;
            TableDist[IV2][IV1]=dist1;
            IV2++;

        }
        IV1++;
    }

    //cout << "#############  BEGIN "<<endl;
    int num=0,num1=0;
    IV1=0; IV2=0;IV3=0;
    int intdist1,intdist2,intdist3;
    unsigned int intP1, intP2,intP3;
    unsigned int intP11, intP12,intP13;
    int nd=0, np=0;
    for (std::vector<InterPoint>::const_iterator pt1=interResult.listInters.begin();pt1 != interResult.listInters.end();pt1++){
        IV2=IV1+1;intP1 = pt1->interaction;
        // TODO: Correct the bug with pt2 = pt1 + 1
        for (std::vector<InterPoint>::const_iterator pt2=pt1 ; pt2 != interResult.listInters.end(); pt2++)	{
            IV3=IV2+1;intP2 = pt2->interaction;
            for (std::vector<InterPoint>::const_iterator pt3 = pt2+1; pt3 != interResult.listInters.end(); pt3++) {
                //  if (!(pt1->interaction == InterType::HYDROPHOBIC|| pt2->interaction == InterType::HYDROPHOBIC || pt3->interaction == InterType::HYDROPHOBIC))continue;
                //  if (pt1->interaction != InterType::HYDROPHOBIC && pt2->interaction != InterType::HYDROPHOBIC && pt3->interaction != InterType::HYDROPHOBIC)continue;
                // cout << "ST"<<endl;
                num1++;
                intP3 = pt3->interaction;
                dist1 = TableDist[IV1][IV2];
                dist2 = TableDist[IV2][IV3];
                dist3 = TableDist[IV3][IV1];


                //                cout << dist1<< " " << dist2 << " " << dist3 << " " << pt1->interaction << " " << pt2->interaction << " " <<pt3->interaction<<endl;
                if (dist1 > dist2+dist3){cerr << "TRIANGULAR INEGALITY NOT VERIFIED !! " <<endl;continue;}
                if (dist2 > dist1+dist3){cerr << "TRIANGULAR INEGALITY NOT VERIFIED !! " <<endl;continue;}
                if (dist3 > dist1+dist2){cerr << "TRIANGULAR INEGALITY NOT VERIFIED !! " <<endl;continue;}


                dist11=dist1;
                intP11=intP1;

                if (dist2 > dist1) {
                    dist11 = dist2;
                    intP11=intP2;
                    dist12=dist1;
                    intP12=intP1;
                } else {
                    dist12=dist2;
                    intP12=intP2;
                }

                if (dist3 > dist11){
                    dist13=dist12;
                    intP13=intP12;
                    dist12=dist11;
                    intP12=intP11;
                    dist11=dist3;
                    intP11=intP3;
                } else if (dist3 > dist12) {
                    dist13=dist12;
                    intP13=intP12;
                    dist12=dist3;
                    intP12=intP3;
                } else {
                    dist13=dist3;
                    intP13=intP3;
                }

                if (dist13 == 0) dist13 = 0.001;
                intdist1 = -1; intdist2=-1;intdist3=-1;
                for (unsigned int pos = 0 ; pos < path_dist_size; pos++) if (dist11 >= path_dist[pos] && dist11 < path_dist[pos+1])   intdist1=(pos+1);
                for (unsigned int pos = 0 ; pos < path_dist_size; pos++) if (dist12 >= path_dist[pos] && dist12 < path_dist[pos+1])   intdist2=(pos+1);
                for (unsigned int pos = 0 ; pos < path_dist_size; pos++) if (dist13 >= path_dist[pos] && dist13 < path_dist[pos+1])   intdist3=(pos+1);

                if (intdist1 == -1 || intdist2 == -1 || intdist3 == -1) continue;
                nd = intdist1*100+intdist2*10+intdist3;
                np = intP11*100+intP12*10+intP13;
                //                if (debug)   cout << nd << " " << np << " " << vect_list[nd][np]<<endl;

               if (vect_list[nd][np] != -1) {
                    fgp->pushOne(vect_list[nd][np]);
                }


                num++;
                IV3++;
            }
            IV2++;
        }
        IV1++;
    }

    //    if (debug) cout << num << " " << num1<<endl;// <<" " << fgp->__ToString()<<endl;

    return *fgp;
}



void Interactions::genIFP(InterResults& interResult, const unsigned int& fgpType) const {
    static const int intToPos[5][NB_INTTYPE]= {
        {-1,3,4,5,6, 0,-1, 1, 2,-1,-1,-1,-1,-1},
        {-1,0,1,2,3,-1,4,-1,-1,-1,-1,-1,-1,-1},
        {-1,3,4,5,6, 0,8, 1, 2, 7, -1, -1,-1,-1},
        {-1,0,1,2,3,-1,7,-1,-1, 6, 4, 5,-1,-1},
        {-1,-1,-1,-1,-1,-1,0,-1,-1,-1,-1,-1,-1}
    };

    static const short length[5]={7,5,9,8,1};
    //     U ,H,H,I,I, H,M , A, A
    //     N ,B,B,O,O, Y,E , R, R
    //     D ,P,L,P,L, D,T , F, E
    //     E , , , , ,  ,  , F, F
    //     F
    
    map<int,Residu*> NumtoRes;
    for (ItCRes itR = complex.firstResidu();itR != complex.lastResidu();++itR) 
    {
        Residu *res =*itR;
        if (!res->isUsed() || res->getParent()->getMoleType()==MoleType::LIGAND) continue;
        if (Residu::Rules[res->getParent()->getMoleType()][res->getResType()] == MoleType::UNDEFINED || Molecule::Rules[res->getParent()->getMoleType()] == MoleType::UNDEFINED)
            continue;
        NumtoRes.insert(pair<int,Residu*>(res->getNum(),res));
    }



    multimap<Residu*,InterPoint*> listRes;

    for (size_t i=0; i< interResult.listInters.size();++i)
    {
        InterPoint& interP = interResult.listInters.at(i);
        if (interP.merged_to != -1)continue;
        Residu *res =interP.Prot_Ref->getResidu();
        if (interP.interaction==InterType::AREDGEFACE || interP.interaction==InterType::ARFACEFACE)
        {
            Cycle* cyc=interP.Prot_Ref->getParent().getCycleFromCenter(interP.Prot_Ref);
            res=cyc->getAtom(0)->getResidu();

        }
        listRes.insert(pair<Residu*,InterPoint*>(res,&interP));

    }

    interResult.IFP=Fingerprint(NumtoRes.size()*length[fgpType]);
    int NRes=0;
    ostringstream oss;
    for (map<int,Residu*>::iterator it =NumtoRes.begin(); it != NumtoRes.end();++it ) {

        std::pair <std::multimap<Residu*,InterPoint*>::iterator, std::multimap<Residu*,InterPoint*>::iterator> ret;
        ret = listRes.equal_range((*it).second);

        for (std::multimap<Residu*,InterPoint*>::iterator it2=ret.first; it2!=ret.second; ++it2) {
            if (intToPos[fgpType][(*it2).second->interaction] == -1)
                continue;

            interResult.IFP.bitOn(NRes*length[fgpType]+intToPos[fgpType][(*it2).second->interaction]);
        }
        interResult.IFPString+="|";
        oss.str("");
        if ((*it).second->getResType() == ResType::STD_AA)
        {
            bool found=false;
            for (size_t pData=0;pData <NBAA;++pData)
            {

                if (AAcid[pData].name==(*it).second->getName())
                {
                    oss << AAcid[pData].code << (*it).second->getFNum();
                    found=true;
                    break;
                }
            }
            if (!found){
                oss << (*it).second->getName() << (*it).second->getFNum();
            }
        }
        else {      oss <<(*it).second->getName()<<(*it).second->getFNum();
        }
        //        cout << oss.str() << endl;
        interResult.IFPString+=(*it).second->getChainName();
        for (int ni=0; ni < length[fgpType]-(int)oss.str().length()-2;ni++)
            interResult.IFPString += " ";
        interResult.IFPString += oss.str();


        NRes++;
    }
}







