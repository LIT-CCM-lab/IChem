#include <iostream>
#include "headers/ICMole/molecule.h"
#include "headers/ICMole/bond.h"
#include "headers/ICMole/vertex.h"
#include "headers/ICMole/edge.h"
#include "headers/ICMole/cycle.h"
#include "headers/ICMole/complex.h"
#include "headers/ICMole/box.h"
//#define ICHEM_DEBUG 1
#ifdef ICHEM_DEBUG
#include "headers/ICPars/molewriter.h"
#endif
#include "headers/ICPars/molewriter.h"
using namespace std;
using namespace ICMole;






struct NewHAtom
{
  Atom* atmA;
  Coords HPos;
};

vector<NewHAtom> addPlanar2H(Atom& atm)
{
    vector<NewHAtom> tmpListNewH;
    NewHAtom NHA;NHA.atmA=&atm;
  const Atom &atmL = atm.getAtomLinked(0);
  Coords axe= atm.fixpos-atmL.fixpos;axe.normalize();
  Coords center(axe);center*=0.49;center+=atm.fixpos;
  Coords atmCheckAngl;
  if (&atmL.getAtomLinked(0) != &atm)
  {
   atmCheckAngl= atmL.getAtomLinked(0).fixpos-atmL.fixpos;//-atmL.fixpos+atm.fixpos;
  }
  else
  {
      atmCheckAngl= atmL.getAtomLinked(1).fixpos-atmL.fixpos;//-atmL.fixpos+atm.fixpos;
  }
  //atmCheckAngl.normalize();
 // cout << atmCheckAngl.toString()<<endl;
  atmCheckAngl+=atm.fixpos;
 Coords A2(atmL.fixpos);A2.x+=0.1;

  Coords vect_x=center.getNormal(atmL.fixpos,A2);
  Coords vect_N=center-atm.fixpos;vect_N.normalize();
  Coords vect_U=vect_x-atm.fixpos;vect_U.normalize();
  const double scal=vect_N.x*vect_U.x+vect_N.y*vect_U.y+vect_N.z*vect_U.z;
  Coords vectoriel(vect_N.y*vect_U.z-vect_N.z*vect_U.y,
                   vect_N.z*vect_U.x-vect_N.x*vect_U.z,
                   vect_N.x*vect_U.y-vect_N.y*vect_U.x);
  //double diff
   double best_angl=170*M_PI/180, best_i=-10;
  for (double i=0; i<=2*M_PI;i+=0.25)
  {
       Coords v =vect_U*cos(i)+vect_N*(1-cos(i))*scal+vectoriel*sin(i)+atm.fixpos;
       Coords v2=vect_U*cos(i+M_PI)+vect_N*(1-cos(i+M_PI))*scal+vectoriel*sin(i+M_PI)+atm.fixpos;
       double angl1=atm.fixpos.calcAngle(v,atmCheckAngl);
       double angl2=atm.fixpos.calcAngle(v2,atmCheckAngl);
       if (angl1 > best_angl)
       {
       best_angl=angl1; best_i=i;
       continue;
       }
       if (angl2 > best_angl)
       {
       best_angl=angl2; best_i=i;
       continue;
       }

   }
  if (best_i>-10)
  {
      Coords v =vect_U*cos(best_i)+vect_N*(1-cos(best_i))*scal+vectoriel*sin(best_i)+atm.fixpos;
      Coords v2=vect_U*cos(best_i+M_PI)+vect_N*(1-cos(best_i+M_PI))*scal+vectoriel*sin(best_i+M_PI)+atm.fixpos;

       NHA.HPos=v; tmpListNewH.push_back(NHA);
        NHA.HPos=v2; tmpListNewH.push_back(NHA);
  }
  return tmpListNewH;

}

void Molecule::addHydrogen(const std::string &verbose)
{

    const bool wVerbose= !verbose.empty();
    if (wVerbose) cout << verbose <<"START\t\n";
    AtomList FAILED_ATOM;BondList FAILED_BOND;
    AtomList newHydrogens;

{
    const size_t nAtoms = Atoms.size();
    for (size_t iAtm=0;iAtm < nAtoms; ++iAtm)
    {
        Atom& atomRef = *Atoms.at(iAtm);
        if (atomRef.isHydrogen())continue;
        const std::string& MOL2 = atomRef.getMOL2Type();
        const size_t nLink      = atomRef.getNumBond();


        // Take a tetrahedric carbon. If you put 2 over 4 bonds on look at them
        // on a plan, the two other bonds will be above and below this first plan.
        // Moreover, the two missing atoms will also be in the plan formed by
        // the reference atom and the geometrical center of the two bonds
        // So if you have like here only two bonds, you can describe them onto
        // an plan to find the two other bonds.
        // So let's be bond1, bond2 the two bonds in the plan.
        // Let's be atomC1 and atomC2 the two atoms involved in the two bonds
        // The normal vector is the vector perpendicular to the plane formed by
        // the two bonds. This normal vector is x_axis and will be parallel to
        // the axe between the two missing bonds.
        // midVect is here to represent the vector between the reference atom
        // and the geometrical center.
        // Then we just need a multiplication to find the position of the two atoms.

        if ((MOL2 == "C.3"||MOL2=="N.4") && nLink >=2&& nLink != 4)
        {

            Bond& bond1 =*atomRef.getBond(0);
            Bond& bond2 =*atomRef.getBond(1);
            const Atom& atomC1= bond1.getOtherAtom(atomRef);
            const Atom& atomC2= bond2.getOtherAtom(atomRef);
            const Coords x_axis= atomRef.fixpos.getNormal(atomC1.fixpos,atomC2.fixpos);
            const Coords diffC1(atomC1.fixpos-atomRef.fixpos,true);
            const Coords diffC2(atomC2.fixpos-atomRef.fixpos,true);
            const Coords midVect((diffC1+diffC2)/-2);
            const Coords newH1(atomRef.fixpos-(x_axis-atomRef.fixpos)*0.81+midVect*0.99);
            const Coords newH2(atomRef.fixpos+(x_axis-atomRef.fixpos)*0.81+midVect*0.99);
            if (nLink==2)
            {

                Atom &atomH1=addAtom("H",newH1,"H","H",atomRef.getResidu());
                addBond(atomH1,atomRef,BondType::SINGLE).setBondType(BondType::SINGLE);
                newHydrogens.push_back(&atomH1);

                Atom &atomH2=addAtom("H",newH2,"H","H",atomRef.getResidu());
                addBond(atomH2,atomRef,BondType::SINGLE).setBondType(BondType::SINGLE);
                newHydrogens.push_back(&atomH2);

            }
            else if (nLink==3)
            {
                const Atom& atomC3= atomRef.getAtomLinked(2);
                const double ang1=atomRef.fixpos.calcAngle(atomC3.fixpos,newH1);
                const double ang2=atomRef.fixpos.calcAngle(atomC3.fixpos,newH2);

                if (ang1> ang2)
                {
                    Atom &atomH1=addAtom("H",newH1,"H","H",atomRef.getResidu());
                    addBond(atomH1,atomRef,BondType::SINGLE).setBondType(BondType::SINGLE);
                    newHydrogens.push_back(&atomH1);
                }
                else
                {
                    Atom &atomH2=addAtom("H",newH2,"H","H",atomRef.getResidu());
                    addBond(atomH2,atomRef,BondType::SINGLE).setBondType(BondType::SINGLE);
                    newHydrogens.push_back(&atomH2);
                }

            }// END ELSE IF nLink C.3

        }// END C.3

        // Here it's simpler, since all three bonds are on the same plan,
        // and each bond is the opposite vector of the geometrical center of the two others.
        else if ((MOL2 == "C.ar"|| MOL2=="N.am" || MOL2=="N.pl3" || MOL2=="C.2") && nLink == 2)
        {
            Bond& bond1 =*atomRef.getBond(0);
            Bond& bond2 =*atomRef.getBond(1);
            const Atom& atomC1= bond1.getOtherAtom(atomRef);
            const Atom& atomC2= bond2.getOtherAtom(atomRef);
            const Coords diffC1(atomC1.fixpos-atomRef.fixpos,true);
            const Coords diffC2(atomC2.fixpos-atomRef.fixpos,true);
            const Coords midVect((diffC1+diffC2)/-2,true);
            Atom& atomH1=addAtom("H",atomRef.fixpos+midVect,"H","H",atomRef.getResidu());
            addBond(atomH1,atomRef,BondType::SINGLE).setBondType(BondType::SINGLE);
            newHydrogens.push_back(&atomH1);
        }
        else if ((MOL2 == "N.pl3" || MOL2=="N.am") && nLink==1)
        {

            vector<NewHAtom>tmpListNewH;
            tmpListNewH=    addPlanar2H(atomRef);
            for (size_t i=0; i<tmpListNewH.size();i++)
            {
                Atom& atomH1=addAtom("H",tmpListNewH.at(i).HPos,"H","H",atomRef.getResidu());
                addBond(atomH1,atomRef,BondType::SINGLE).setBondType(BondType::SINGLE);
                newHydrogens.push_back(&atomH1);
            }

        }


    }// END iAtm
}//END STEP1

checkMOL2();
    bool useGrid=false;
    Grid *grid=(Grid*)NULL;
    if (complex != (Complex*)NULL)
    {
        grid=&complex->getGrid(0);
        useGrid=true;
        if (!newHydrogens.empty())
        {
            grid->rotateAtoms(newHydrogens);
            newHydrogens.clear();
        }

    }


//STEP2
    {
    const size_t nAtoms = Atoms.size();
    AtomList atmAccH, atmClose,adjacentList;
    for (size_t iAtom=0;iAtom < nAtoms;++iAtom)
    {
        Atom &atomRef = *Atoms.at(iAtom);

        if (!(atomRef.getMOL2Type() == "C.3"||
            atomRef.getMOL2Type() == "N.4"))continue;
        if (atomRef.getNumBond()  != 1) continue;

        // Getting close atoms for clashes
        atmAccH.clear();atmClose.clear();
        if (!useGrid || atomRef.getBox(grid)==(Box*)NULL)
        {
            for (size_t iAtm2=0;iAtm2 < nAtoms; ++iAtm2)
            {
                Atom &atm2=*Atoms.at(iAtm2);
                if (&atm2!= &atomRef
                        && atm2.fixpos.calcDist(atomRef.fixpos,3.6)< 3.5)
                {

                    if (atm2.props.isAcceptor())atmAccH.push_back(&atm2);
                    atmClose.push_back(&atm2);
                }

            }// END FOR iAtm2
        }
        else if (atomRef.getBox(grid)!=(Box*)NULL)
        {
            grid->getAdjacentAtoms(adjacentList,atomRef,3.51);
              for (ItCAtom itA = adjacentList.begin(); itA != adjacentList.end();++itA)
              {
                  Atom &atm2 = **itA;
                  if (&atm2!= &atomRef)
                  {
                      if (atm2.props.isAcceptor() || atm2.props.isweakAcceptor())atmAccH.push_back(&atm2);
                      atmClose.push_back(&atm2);
                  }
              }
        }
        // The onyl atom linked
        const Atom &atmL = atomRef.getAtomLinked(0);

        // First atom in beta of C.3 so we can check that the angle is 70° between new hydrogen, C.3 and this atom
        Coords atmCheckAngl;
        if (&atmL.getAtomLinked(0) != &atomRef)
            atmCheckAngl= atmL.getAtomLinked(0).fixpos;
        else
            atmCheckAngl= atmL.getAtomLinked(1).fixpos;


        Coords axe= atomRef.fixpos-atmL.fixpos;axe.normalize();
        Coords center(axe);center*=0.4;center+=atomRef.fixpos;
        Coords A2(atmL.fixpos);A2.x+=0.1;

        Coords vect_x=center.getNormal(atmL.fixpos,A2);
        Coords vect_N=center-atomRef.fixpos;vect_N.normalize();
        Coords vect_U=vect_x-atomRef.fixpos;vect_U.normalize();
        const double scal=vect_N.x*vect_U.x+vect_N.y*vect_U.y+vect_N.z*vect_U.z;
        Coords vectoriel(vect_N.y*vect_U.z-vect_N.z*vect_U.y,
                         vect_N.z*vect_U.x-vect_N.x*vect_U.z,
                         vect_N.x*vect_U.y-vect_N.y*vect_U.x);
double best_diff=100, best_i=0;
        for (double i=0; i<=2*M_PI/3;i+=0.02)
        {
            Coords v =vect_U*cos(i)+vect_N*(1-cos(i))*scal+vectoriel*sin(i)+atomRef.fixpos;
            Coords v2=vect_U*cos(i+2*M_PI/3)+vect_N*(1-cos(i+2*M_PI/3))*scal+vectoriel*sin(i+2*M_PI/3)+atomRef.fixpos;
            Coords v3=vect_U*cos(i+4*M_PI/3)+vect_N*(1-cos(i+4*M_PI/3))*scal+vectoriel*sin(i+4*M_PI/3)+atomRef.fixpos;
           // cout << atomRef.fixpos.calcAngle(v,atmCheckAngl)*180/M_PI<< " " << atomRef.fixpos.calcAngle(v2,atmCheckAngl)*180/M_PI<< " " <<atomRef.fixpos.calcAngle(v3,atmCheckAngl)*180/M_PI<<endl;
           // cout << atomRef.fixpos.calcAngle(v,atmCheckAngl)<< " " << atomRef.fixpos.calcAngle(v2,atmCheckAngl)<< " " <<atomRef.fixpos.calcAngle(v3,atmCheckAngl)<<endl;
            const double angl1=atomRef.fixpos.calcAngle(v,atmCheckAngl);
            const double angl2=atomRef.fixpos.calcAngle(v2,atmCheckAngl);
            const double angl3=atomRef.fixpos.calcAngle(v3,atmCheckAngl);

            double diffs=0;
            if (fabs(angl1-angl2)< 0.2)diffs = (angl1+angl2)/2;
            else if (fabs(angl1-angl3)< 0.2)diffs = (angl1+angl3)/2;
            else if (fabs(angl2-angl3)< 0.2)diffs = (angl2+angl3)/2;
            //cout << diffs<<endl;


            if (fabs(diffs-M_PI/2) < best_diff) {

                best_diff=fabs(diffs-M_PI/2);best_i=i;
            }
        }
        //cout << atomRef.getIdentifier()<<"\t"<<best_i<<"\t"<< best_diff<<endl;

        if ((best_diff < 100))
        {
            Coords v =vect_U*cos(best_i)+vect_N*(1-cos(best_i))*scal+vectoriel*sin(best_i)+atomRef.fixpos;
            Coords v2=vect_U*cos(best_i+2*M_PI/3)+vect_N*(1-cos(best_i+2*M_PI/3))*scal+vectoriel*sin(best_i+2*M_PI/3)+atomRef.fixpos;
            Coords v3=vect_U*cos(best_i+4*M_PI/3)+vect_N*(1-cos(best_i+4*M_PI/3))*scal+vectoriel*sin(best_i+4*M_PI/3)+atomRef.fixpos;

            Atom &atomH1=addAtom("H",v,"H","H",atomRef.getResidu());
            addBond(atomH1,atomRef,BondType::SINGLE);
            newHydrogens.push_back(&atomH1);

            Atom &atomH2=addAtom("H",v2,"H","H",atomRef.getResidu());
            addBond(atomH2,atomRef,BondType::SINGLE);
            newHydrogens.push_back(&atomH2);

            Atom &atomH3=addAtom("H",v3,"H","H",atomRef.getResidu());
            addBond(atomH3,atomRef,BondType::SINGLE);
            newHydrogens.push_back(&atomH3);

        }
        else
        {

            FAILED_ATOM.push_back(&atomRef);
            if (wVerbose)cerr << verbose<<"ERROR\tUnable to find hydrogen for"<<atomRef.getIdentifier()<<endl;
        }

    } // END FOR iAtom


}   //END STEP2




    AtomList AccHAtm,DonHAtm;

    vector<string> ListAtmAH,ListAtmDH;
    ListAtmAH.push_back("O.2");
    ListAtmAH.push_back("O.3");
    ListAtmAH.push_back("O.co2");
    ListAtmAH.push_back("S.3");
    ListAtmAH.push_back("N.ar");
    ListAtmAH.push_back("N.2");
    ListAtmAH.push_back("O.spc");
    ListAtmAH.push_back("O.t3p");
    ListAtmDH.push_back("O.3");
    ListAtmDH.push_back("N.pl3");
    ListAtmDH.push_back("N.4");
    ListAtmDH.push_back("N.ar");
    ListAtmDH.push_back("O.spc");
    ListAtmDH.push_back("O.t3p");

try{
    checkMOL2();
    if (useGrid)
    {
        grid->rotateAtoms(newHydrogens);
        newHydrogens.clear();
    }

}
    catch (MoleExcept &e)
    {
        cerr << verbose<<"ERROR\t"<<e.getCode()<<"\t"
             << e.getData()<<"\t"
             << e.getTrace()<<endl;
    }

     vector<NewHAtom> tmpListNewH;
     string ResName="";
     for (size_t iRes = 0; iRes < Residues.size();iRes++)
     {
         Residu &res = *Residues.at(iRes);
         ResName=res.getName();
         if (ResName != "ASN" && ResName != "GLN"&& ResName != "HIS")continue;
         if (ResName == "ASN" || ResName == "GLN")
         {

             Atom* Nitrogen=(Atom*)NULL;
             Atom* Oxygen  =(Atom*)NULL;
             Atom* NH1     =(Atom*)NULL;
             Atom* NH2     =(Atom*)NULL;
             bool clashOx=false,clashNi=false;
             for (ItCAtom itCA = res.firstCAtom();itCA != res.lastCAtom();++itCA)
             {
                 Atom &atm = **itCA;

                   if (atm.getName()=="OE1" || atm.getName()=="OD1") Oxygen=&atm;
                   else if (atm.getName()=="NE2" || atm.getName()=="ND2")
                   {
                       Nitrogen=&atm;
                       for (size_t nNbd=0; nNbd < atm.getNumBond();nNbd++)
                       {
                            Atom &atmLink = atm.getBond(nNbd)->getOtherAtom(atm);
                           if (atmLink.getName()!="H")continue;

                           if (NH1==(Atom*)NULL)NH1=&atmLink;
                           else NH2=&atmLink;

                       }
                   }
             }
             if (Nitrogen==(Atom*)NULL || Oxygen==(Atom*)NULL || NH1==(Atom*)NULL||NH2==(Atom*)NULL)continue;
             for (size_t iAtm2 = 0; iAtm2 < Atoms.size();++iAtm2)
             {
                 Atom &atm2 = *Atoms.at(iAtm2);
                 if (&atm2==Oxygen || &atm2==Nitrogen || &atm2==NH1 ||&atm2==NH2)continue;
                 if (atm2.isNitrogen() && atm2.fixpos.calcDist(Nitrogen->fixpos,3.1)<3)clashNi=true;
                 else if (atm2.isOxygen() && atm2.fixpos.calcDist(Oxygen->fixpos,3.1)<3)clashOx=true;
                 if (clashNi && clashOx)break;
             }

             if (!(clashNi && clashOx)) continue;
             delAtom(NH1);delAtom(NH2);
              Nitrogen->fixpos.swap(Oxygen->fixpos);
              tmpListNewH.clear();
              tmpListNewH=  addPlanar2H(*Nitrogen);
              for (size_t i=0; i<tmpListNewH.size();i++)
              {
                  Atom &atomH1=addAtom("H",tmpListNewH.at(i).HPos,"H","H",Nitrogen->getResidu());
                          newHydrogens.push_back(&atomH1);
                  addBond(&atomH1,Nitrogen,BondType::SINGLE);
              }
             if (wVerbose) cout << verbose<<"MODIFS\t"<<  Nitrogen->getIdentifier()<<"<->"<<Oxygen->getIdentifier()<<endl;

         }
         else if (res.getName()=="HIS")
         {

              Atom* NE2=(Atom*)NULL;Atom* HNE2=(Atom*)NULL;
              Atom* ND1=(Atom*)NULL;Atom* HND1=(Atom*)NULL;Atom*Ccommon=(Atom*)NULL;
              bool toSwitch=false;
              bool clashND=false,clashNE=false;
              for (ItCAtom itCA = res.firstCAtom();itCA != res.lastCAtom();++itCA)
              {
                  Atom &atm = **itCA;
                    if (atm.getName()=="ND1") ND1=&atm;
                    else if (atm.getName()=="NE2") NE2=&atm;
                    else if (atm.getName()=="H" && atm.getNumBond()==1)
                    {
                        if (atm.getAtomLinked(0).getName()=="ND1")HND1=&atm;
                        else if (atm.getAtomLinked(0).getName()=="NE2")HNE2=&atm;

                    }

              }
              if (NE2==(Atom*)NULL || ND1==(Atom*)NULL || (HND1==(Atom*)NULL && HNE2==(Atom*)NULL))continue;
              for (size_t iNBd = 0; iNBd < NE2->getNumBond();++iNBd)
              {
                   Atom& atmLink = const_cast<Atom&>(NE2->getAtomLinked(iNBd));
                  if (atmLink.hasBondWith(*ND1))Ccommon=&atmLink;
              }
              if (Ccommon==(Atom*)NULL)continue;
              for (size_t iAtm2 = 0; iAtm2 < Atoms.size();++iAtm2)
              {
                  Atom &atm2 = *Atoms.at(iAtm2);
                  if (atm2.isMetallic())
                  {
                      if (HND1!=(Atom*)NULL &&
                              atm2.fixpos.calcDist(ND1->fixpos,3.1) < 3 &&
                              fabs(HND1->fixpos.calcAngle(ND1->fixpos,atm2.fixpos)-M_PI) <M_PI/6) {toSwitch=true;break;}
                      if (HNE2!=(Atom*)NULL &&
                              atm2.fixpos.calcDist(NE2->fixpos,3.1) < 3 &&
                              fabs(HNE2->fixpos.calcAngle(NE2->fixpos,atm2.fixpos)-M_PI) <M_PI/6) {toSwitch=true;break;}

                  }
                  if (atm2.props.isDonor() || atm2.props.isweakDonor())
                  {
                      if (HND1!=(Atom*)NULL && atm2.fixpos.calcDist(ND1->fixpos,3.1) < 3) clashND=true;
                      if (HNE2!=(Atom*)NULL && atm2.fixpos.calcDist(NE2->fixpos,3.1) < 3) clashNE=true;
                  }
                  if (atm2.props.isAcceptor() || atm2.props.isweakAcceptor())
                  {
                      if (HND1==(Atom*)NULL && atm2.fixpos.calcDist(ND1->fixpos,3.1) < 3) clashND=true;
                      if (HNE2==(Atom*)NULL && atm2.fixpos.calcDist(NE2->fixpos,3.1) < 3) clashNE=true;
                  }
              }
              if (!toSwitch && (!clashND && !clashNE))continue;
              if (wVerbose)cout <<verbose<<"MODIFS\t"<<NE2->getIdentifier()<<"<->"<<ND1->getIdentifier()<<endl;
              NE2->setName("ND1");
              ND1->setName("NE2");
            //  CLC.clear();

              if (HND1!=(Atom*)NULL)
              {
                  delAtom(HND1);
                  Bond& bond1 =*NE2->getBond(0);
                  Bond& bond2 =*NE2->getBond(1);
                  const Atom& atomC1= bond1.getOtherAtom(NE2);
                  const Atom& atomC2= bond2.getOtherAtom(NE2);
                  const Coords diffC1(atomC1.fixpos-NE2->fixpos,true);
                  const Coords diffC2(atomC2.fixpos-NE2->fixpos,true);
                  const Coords midVect((diffC1+diffC2)/-2);
                  Atom& atomH1=addAtom("H",NE2->fixpos+midVect*2,"H","H",NE2->getResidu());
                  addBond(&atomH1,NE2,BondType::SINGLE).setBondType(BondType::SINGLE);
                  newHydrogens.push_back(&atomH1);
                  const_cast<Bond*>(NE2->getBondWith(*Ccommon))->setBondType(BondType::SINGLE);
                  const_cast<Bond*>(ND1->getBondWith(*Ccommon))->setBondType(BondType::DOUBLE);
                  NE2->setMOL2Type("N.pl3");
                  ND1->setMOL2Type("N.2");
              }
              else if (HNE2!=(Atom*)NULL)
              {
                  delAtom(HNE2);
                  Bond& bond1 =*ND1->getBond(0);
                  Bond& bond2 =*ND1->getBond(1);
                  const Atom& atomC1= bond1.getOtherAtom(ND1);
                  const Atom& atomC2= bond2.getOtherAtom(ND1);
                  const Coords diffC1(atomC1.fixpos-ND1->fixpos,true);
                  const Coords diffC2(atomC2.fixpos-ND1->fixpos,true);
                  const Coords midVect((diffC1+diffC2)/-2);
                  Atom& atomH1=addAtom("H",ND1->fixpos+midVect*2,"H","H",ND1->getResidu());
                  newHydrogens.push_back(&atomH1);
                  addBond(&atomH1,ND1,BondType::SINGLE).setBondType(BondType::SINGLE);
                  const_cast<Bond*>(ND1->getBondWith(*Ccommon))->setBondType(BondType::SINGLE);
                  const_cast<Bond*>(NE2->getBondWith(*Ccommon))->setBondType(BondType::DOUBLE);
                  ND1->setMOL2Type("N.pl3");
                  NE2->setMOL2Type("N.2");


              }




         }
     }



     try{
         checkMOL2();
         if (useGrid)
         {
             grid->rotateAtoms(newHydrogens);
             newHydrogens.clear();
         }

     }
         catch (MoleExcept &e)
         {
             cerr << verbose<<"ERROR\t"<<e.getCode()<<"\t"
                  << e.getData()<<"\t"
                  << e.getTrace()<<endl;
         }


    const size_t nAtoms = Atoms.size();
    for (size_t iAtom=0;iAtom < nAtoms;++iAtom)
    {
        const std::string& mol2 = Atoms.at(iAtom)->getMOL2Type();
        if (find(ListAtmAH.begin(),ListAtmAH.end(),mol2) != ListAtmAH.end()) AccHAtm.push_back(Atoms.at(iAtom));
        if (find(ListAtmDH.begin(),ListAtmDH.end(),mol2) != ListAtmDH.end()) DonHAtm.push_back(Atoms.at(iAtom));

      }

    //  cout << mole->toString()<<endl;








const double minAngl=104.9*M_PI/180, maxAngl=105.1*M_PI/180;


for (ItAtom itD = DonHAtm.begin(); itD != DonHAtm.end(); itD++){
    Atom &AtmD = **itD;
    const bool one_Bond=(AtmD.getNumBond()== 1)? true:false;
    if (!one_Bond)continue;
    double best_Dist=100;
    Coords best_coord(0,0,0);
multimap<double,Coords> bDist;

    Coords tmpIfFail;
    AtomList adjacentList;
    if (useGrid && AtmD.getBox(grid)!=(Box*)NULL)
    {
       grid->getAdjacentAtoms(adjacentList,AtmD,4.1);

    }

    for (ItAtom itA = AccHAtm.begin(); itA != AccHAtm.end(); itA++)

    {

        Atom &AtmA = **itA;

        if (!adjacentList.empty() && find(adjacentList.begin(),adjacentList.end(),&AtmA)==adjacentList.end())continue;

        const double dist =AtmA.fixpos.calcDist(AtmD.fixpos,4.1);
        if (dist >= 4 || dist == 0) continue;



        const Atom &atmL = AtmD.getAtomLinked(0);
        // First atom in beta of C.3 so we can check that the angle is 70° between new hydrogen, C.3 and this atom
        Coords atmCheckAngl;
        if (&atmL != &AtmD)
        {
            atmCheckAngl= atmL.getAtomLinked(0).fixpos;//-atmL.fixpos+atm->fixpos;
        }
        else
        {
            atmCheckAngl= atmL.getAtomLinked(1).fixpos;//-atmL.fixpos+atm->fixpos;
        }

        const Coords axe(AtmD.fixpos-atmL.fixpos,true);
        Coords center(axe);center*=0.2678;center+=AtmD.fixpos;
        Coords A2(atmL.fixpos);A2.x+=0.1;

        const Coords vect_x=center.getNormal(atmL.fixpos,A2);
        const Coords vect_N(center-AtmD.fixpos,true);
        const Coords vect_U(vect_x-AtmD.fixpos,true);
        const double scal=vect_N.x*vect_U.x+vect_N.y*vect_U.y+vect_N.z*vect_U.z;
        Coords vectoriel(vect_N.y*vect_U.z-vect_N.z*vect_U.y,
                         vect_N.z*vect_U.x-vect_N.x*vect_U.z,
                         vect_N.x*vect_U.y-vect_N.y*vect_U.x);
        best_Dist=100;
        best_coord.setCoords(0,0,0);

        Coords H =vect_U*cos(1.3)+vect_N*(1-cos(1.3))*scal+vectoriel*sin(1.3)+AtmD.fixpos;
       tmpIfFail=H;

        for (double i=0; i<=2*M_PI;i+=0.1)
        {
            Coords H =vect_U*cos(i)+vect_N*(1-cos(i))*scal+vectoriel*sin(i)+AtmD.fixpos;

            const double angl_HDDl = (AtmD.fixpos.calcAngle(atmL.fixpos,H));

            if (angl_HDDl >= maxAngl || angl_HDDl <= minAngl)continue;

            const double anglAHD = (H.calcAngle(AtmA.fixpos,AtmD.fixpos))*180/PI_CONST;         if (anglAHD <=130) {continue;}
            const double alpha = 34 - H.calcDist(AtmA.fixpos)*sin( 1*sin(angl_HDDl)/dist);
            const double angl_ADDl=((AtmD.fixpos.calcAngle(atmL.fixpos,AtmA.fixpos))*180/PI_CONST);if (angl_ADDl > 109.5+alpha || angl_ADDl < 109.5-alpha) {continue;}
            const double curr_dist=H.calcDist(AtmA.fixpos,best_Dist+0.1);

                bool allowed=true;
            for (size_t iLink=0;iLink < AtmA.getNumBond();++iLink)
            {
                 Atom &atmLinked=AtmA.getBond(iLink)->getOtherAtom(AtmA);
                if (atmLinked.fixpos.calcDist(H) <  (1.50+atmLinked.getVdWRadius())/2){allowed=false;break;}
            }
            if (!allowed)continue;

            if (curr_dist < best_Dist) {

                best_Dist=curr_dist; best_coord=H;}

        }// END FOR i
        if (best_Dist != 100)
        {

            bDist.insert(pair<double,Coords> (best_Dist,best_coord));
        }
    }// END IT ACC

    if (!bDist.empty())
    {
        bool distOk=true;

        multimap<double,Coords>::iterator it=bDist.begin();
        for (size_t iAtmH=0;iAtmH<Atoms.size();iAtmH++)
        {
            const Atom& atmH=*Atoms.at(iAtmH);
            if (!atmH.isHydrogen())continue;
            if (atmH.fixpos.calcDist((*it).second,1.1) < 1){distOk=false;break;}
        }

        if (distOk)
        {
            Atom& atmXH=addAtom("H",(*(it)).second,"H","H",AtmD.getResidu());
            addBond(atmXH,AtmD,BondType::SINGLE);
            newHydrogens.push_back(&atmXH);
           //cout << atmXH.toString();
        }
        else if (bDist.size()>1)
        {

            it=bDist.begin(); it++;
            Atom &atomH1=addAtom("H",(*it).second,"H","H",AtmD.getResidu());
            addBond(atomH1,AtmD,BondType::SINGLE);
            newHydrogens.push_back(&atomH1);
        }
        else
        {

            Atom &atomH1=addAtom("H",(*(it)).second,"H","H",AtmD.getResidu());
            newHydrogens.push_back(&atomH1);
            FAILED_BOND.push_back(&addBond(atomH1,AtmD,BondType::SINGLE));

        }

    }else
    {
        if (tmpIfFail == Coords(0,0,0))
        {
            for (double x=-1;x<=1;x+=0.1){
            for (double y=-1;y<=1;y+=0.1)
            for (double z=-1;z<=1;z+=0.1)
            {
                tmpIfFail.setCoords(x,y,z);
                tmpIfFail.normalize();
                tmpIfFail+=AtmD.fixpos;
                if (tmpIfFail.calcDist(AtmD.fixpos,1.2) > 1.1 || tmpIfFail.calcDist(AtmD.fixpos,1.1)<1)continue;
                double angle=AtmD.fixpos.calcAngle(tmpIfFail,AtmD.getAtomLinked(0).fixpos);
                if (angle < 104*M_PI/180 || angle > 106*M_PI/180)continue;
                break;
            }}

        }
        if (tmpIfFail != Coords(0,0,0))
        {

            Atom &atomH1=addAtom("H",tmpIfFail,"H","H",AtmD.getResidu());
            addBond(atomH1,AtmD,BondType::SINGLE);
            newHydrogens.push_back(&atomH1);
        }
        else if (wVerbose) cerr <<verbose<<"ERRORS\tUnable to find Hydrogen position for "<< AtmD.getIdentifier()<<endl;
    }





}// END IT DON H


try{
    checkMOL2();
    if (useGrid)
    {
        grid->rotateAtoms(newHydrogens);
        newHydrogens.clear();
    }

}
    catch (MoleExcept &e)
    {
        cerr << verbose<<"ERROR\t"<<e.getCode()<<"\t"
             << e.getData()<<"\t"
             << e.getTrace()<<endl;
    }



//vector<NewHAtom> tmpListNewH;
//string ResName="";
for (size_t iRes = 0; iRes < Residues.size();iRes++)
{
    Residu &res = *Residues.at(iRes);
    ResName=res.getName();
    if (ResName != "ASN" && ResName != "GLN"&& ResName != "HIS")continue;
    if (ResName == "ASN" || ResName == "GLN")
    {

        Atom* Nitrogen=(Atom*)NULL;
        Atom* Oxygen  =(Atom*)NULL;
        Atom* NH1     =(Atom*)NULL;
        Atom* NH2     =(Atom*)NULL;
        bool clashOx=false,clashNi=false;
        for (ItCAtom itCA = res.firstCAtom();itCA != res.lastCAtom();++itCA)
        {
            Atom &atm = **itCA;

              if (atm.getName()=="OE1" || atm.getName()=="OD1") Oxygen=&atm;
              else if (atm.getName()=="NE2" || atm.getName()=="ND2")
              {
                  Nitrogen=&atm;
                  for (size_t nNbd=0; nNbd < atm.getNumBond();nNbd++)
                  {
                       Atom &atmLink = atm.getBond(nNbd)->getOtherAtom(atm);
                      if (atmLink.getName()!="H")continue;

                      if (NH1==(Atom*)NULL)NH1=&atmLink;
                      else NH2=&atmLink;

                  }
              }
        }
        if (Nitrogen==(Atom*)NULL || Oxygen==(Atom*)NULL || NH1==(Atom*)NULL||NH2==(Atom*)NULL)continue;
        for (size_t iAtm2 = 0; iAtm2 < Atoms.size();++iAtm2)
        {
            Atom &atm2 = *Atoms.at(iAtm2);
            if (&atm2==Oxygen || &atm2==Nitrogen || &atm2==NH1 ||&atm2==NH2)continue;
            if (atm2.isNitrogen() && atm2.fixpos.calcDist(Nitrogen->fixpos,3.1)<3)clashNi=true;
            else if (atm2.isOxygen() && atm2.fixpos.calcDist(Oxygen->fixpos,3.1)<3)clashOx=true;
            if (clashNi && clashOx)break;
        }

        if (!(clashNi && clashOx)) continue;
        delAtom(NH1);delAtom(NH2);
         Nitrogen->fixpos.swap(Oxygen->fixpos);
         tmpListNewH.clear();
         tmpListNewH=  addPlanar2H(*Nitrogen);
         for (size_t i=0; i<tmpListNewH.size();i++)
             addBond(&addAtom("H",tmpListNewH.at(i).HPos,"H","H",Nitrogen->getResidu()),Nitrogen,BondType::SINGLE);
         if (wVerbose)cout << verbose<<"MODIFS\t"<<  Nitrogen->getIdentifier()<<"<->"<<Oxygen->getIdentifier()<<endl;

    }
    else if (res.getName()=="HIS")
    {
         Atom* NE2=(Atom*)NULL;Atom* HNE2=(Atom*)NULL;
         Atom* ND1=(Atom*)NULL;Atom* HND1=(Atom*)NULL;Atom*Ccommon=(Atom*)NULL;
         bool toSwitch=false;
         bool clashND=false,clashNE=false;
         for (ItCAtom itCA = res.firstCAtom();itCA != res.lastCAtom();++itCA)
         {
             Atom &atm = **itCA;
               if (atm.getName()=="ND1") ND1=&atm;
               else if (atm.getName()=="NE2") NE2=&atm;
               else if (atm.getName()=="H" && atm.getNumBond()==1)
               {
                   if (atm.getAtomLinked(0).getName()=="ND1")HND1=&atm;
                   else if (atm.getAtomLinked(0).getName()=="NE2")HNE2=&atm;

               }

         }
         if (NE2==(Atom*)NULL || ND1==(Atom*)NULL || (HND1==(Atom*)NULL && HNE2==(Atom*)NULL))continue;
         for (size_t iNBd = 0; iNBd < NE2->getNumBond();++iNBd)
         {
              Atom& atmLink = const_cast<Atom&>(NE2->getAtomLinked(iNBd));
             if (atmLink.hasBondWith(*ND1))Ccommon=&atmLink;
         }
         if (Ccommon==(Atom*)NULL)continue;
         for (size_t iAtm2 = 0; iAtm2 < Atoms.size();++iAtm2)
         {
             Atom &atm2 = *Atoms.at(iAtm2);
             if (atm2.getResidu() == &res)continue;
             if (atm2.isMetallic())
             {
                 if (HND1!=(Atom*)NULL &&
                         atm2.fixpos.calcDist(ND1->fixpos,3.1) < 3 &&
                         fabs(HND1->fixpos.calcAngle(ND1->fixpos,atm2.fixpos)-M_PI) <M_PI/6) {toSwitch=true;break;}
                 if (HNE2!=(Atom*)NULL &&
                         atm2.fixpos.calcDist(NE2->fixpos,3.1) < 3 &&
                         fabs(HNE2->fixpos.calcAngle(NE2->fixpos,atm2.fixpos)-M_PI) <M_PI/6) {toSwitch=true;break;}

             }
             if (atm2.props.isDonor() || atm2.props.isweakDonor())
             {
                 if (HND1!=(Atom*)NULL && atm2.fixpos.calcDist(ND1->fixpos,3.1) < 3)
                 {
                     for (size_t iBd=0;iBd < atm2.getNumBond();++iBd)
                     {
                         if (!atm2.getAtomLinked(iBd).isHydrogen())continue;
                         if (HND1->fixpos.calcDist(atm2.getAtomLinked(iBd).fixpos) < 1.5)
                         {
                             clashND=true;
                         }
                     }

                 }
                 if (HNE2!=(Atom*)NULL && atm2.fixpos.calcDist(NE2->fixpos,3.1) < 3)
                 {
                     for (size_t iBd=0;iBd < atm2.getNumBond();++iBd)
                     {
                         if (!atm2.getAtomLinked(iBd).isHydrogen())continue;
                         if (HNE2->fixpos.calcDist(atm2.getAtomLinked(iBd).fixpos) < 1.5)clashNE=true;
                     }
                 }
             }
          /*   if (atm2.props.isAcceptor() || atm2.props.isweakAcceptor())
             {
                 if (HND1==(Atom*)NULL && atm2.fixpos.calcDist(ND1->fixpos,3.1) < 3) clashND=true;
                 if (HNE2==(Atom*)NULL && atm2.fixpos.calcDist(NE2->fixpos,3.1) < 3) clashNE=true;
             }*/
         }
         if (!toSwitch && (!clashND && !clashNE))continue;
         if (wVerbose) cout <<verbose<<"MODIFS\t"<<NE2->getIdentifier()<<"<->"<<ND1->getIdentifier()<<endl;
         NE2->setName("ND1");
         ND1->setName("NE2");
       //  CLC.clear();

         if (HND1!=(Atom*)NULL)
         {
             delAtom(HND1);
             Bond& bond1 =*NE2->getBond(0);
             Bond& bond2 =*NE2->getBond(1);
             const Atom& atomC1= bond1.getOtherAtom(NE2);
             const Atom& atomC2= bond2.getOtherAtom(NE2);
             const Coords diffC1(atomC1.fixpos-NE2->fixpos,true);
             const Coords diffC2(atomC2.fixpos-NE2->fixpos,true);
             const Coords midVect((diffC1+diffC2)/-2);
             Atom& atomH1=addAtom("H",NE2->fixpos+midVect*2,"H","H",NE2->getResidu());
             addBond(&atomH1,NE2,BondType::SINGLE).setBondType(BondType::SINGLE);
             const_cast<Bond*>(NE2->getBondWith(*Ccommon))->setBondType(BondType::SINGLE);
             const_cast<Bond*>(ND1->getBondWith(*Ccommon))->setBondType(BondType::DOUBLE);
             NE2->setMOL2Type("N.pl3");
             ND1->setMOL2Type("N.2");

         }
         else if (HNE2!=(Atom*)NULL)
         {
             delAtom(HNE2);
             Bond& bond1 =*ND1->getBond(0);
             Bond& bond2 =*ND1->getBond(1);
             const Atom& atomC1= bond1.getOtherAtom(ND1);
             const Atom& atomC2= bond2.getOtherAtom(ND1);
             const Coords diffC1(atomC1.fixpos-ND1->fixpos,true);
             const Coords diffC2(atomC2.fixpos-ND1->fixpos,true);
             const Coords midVect((diffC1+diffC2)/-2);
             Atom& atomH1=addAtom("H",ND1->fixpos+midVect*2,"H","H",ND1->getResidu());
             addBond(&atomH1,ND1,BondType::SINGLE).setBondType(BondType::SINGLE);
             const_cast<Bond*>(ND1->getBondWith(*Ccommon))->setBondType(BondType::SINGLE);
             const_cast<Bond*>(NE2->getBondWith(*Ccommon))->setBondType(BondType::DOUBLE);
             ND1->setMOL2Type("N.pl3");
             NE2->setMOL2Type("N.2");


         }




    }
}


checkMOL2();


// SETTING HYDROGEN NAMES && CHECKING CLASHES :
unsigned int NPos=0;
  AtomList adjacentList;
for (size_t iAtm = 0 ; iAtm < Atoms.size();++iAtm)
{

    Atom& atm = *Atoms.at(iAtm);

    if (atm.isMetallic())
    {
        if (atm.getMOL2Type()=="Mg"|| atm.getMOL2Type()=="Ca"||atm.getMOL2Type()=="Zn"||atm.getMOL2Type()=="Mn"||atm.getMOL2Type()=="Cu") atm.setFormalCharge(2);
        else if (atm.getMOL2Type()=="Na") atm.setFormalCharge(1);
        else if (atm.getMOL2Type()=="Al") atm.setFormalCharge(3);
        else if (atm.getMOL2Type()=="Fe")
        {
            for (size_t iAtm2 = 0 ; iAtm2 < Atoms.size();++iAtm2)
            {
                Atom& atm2 = *Atoms.at(iAtm2);
                if (!atm2.props.isAcceptor() && !atm2.props.isDonor())continue;
                if (atm2.fixpos.calcDist(atm.fixpos,3.1)> 3)continue;

            }
        }
        if (useGrid && atm.getBox(grid)!=(Box*)NULL)
        {
           grid->getAdjacentAtoms(adjacentList,atm,3.1);
           for (ItAtom itA2 = adjacentList.begin();
                       itA2!= adjacentList.end();
                ++itA2)
           {
               Atom& atm2 = **itA2;
               const double distD_M=atm2.fixpos.calcDist(atm.fixpos,3.6);
               if (distD_M > 3.5)continue;
               for (size_t i=0; i<atm2.getNumBond();i++)
               {
                   Atom &H=const_cast<Atom&>(atm2.getAtomLinked(i));
                   if (!H.isHydrogen())continue;
                   if (H.fixpos.calcDist(atm.fixpos,distD_M+0.1) > distD_M)continue;
                   delAtom(&H);
                   if (atm2.getMOL2Type()=="N.4") {atm2.setMOL2Type("N.3");}
                   else if (atm2.getMOL2Type()=="N.ar"){atm2.setMOL2Type("N.pl3");
                   atm2.setFormalCharge(-1);}
                   else atm2.setFormalCharge(-1);
                   if (wVerbose)cout << verbose<<"|MODIFS\t"<< atm2.getIdentifier()<<" DELETE HYDROGEN CLOSE TO METAL\n";
               }
           }

        }
        else
        {
            for (size_t iAtm2 = 0 ; iAtm2 < Atoms.size();++iAtm2)
            {
                Atom& atm2 = *Atoms.at(iAtm2);
                if (!atm2.props.isDonor())continue;
                double distD_M=atm2.fixpos.calcDist(atm.fixpos,3.6);
                if (distD_M > 3.5)continue;
                for (size_t i=0; i<atm2.getNumBond();i++)
                {
                    Atom &H=const_cast<Atom&>(atm2.getAtomLinked(i));
                    if (!H.isHydrogen())continue;
                    if (H.fixpos.calcDist(atm.fixpos,distD_M+0.1) > distD_M)continue;
                    delAtom(&H);
                    if (atm2.getMOL2Type()=="N.4") {atm2.setMOL2Type("N.3");}
                    else if (atm2.getMOL2Type()=="N.ar"){atm2.setMOL2Type("N.pl3");
                    atm2.setFormalCharge(-1);}
                    else atm2.setFormalCharge(-1);
                    if (wVerbose)cout << verbose<<"MODIFS\t"<< atm2.getIdentifier()<<" DELETE HYDROGEN CLOSE TO METAL\n";
                }

            }
        }
    }
    if (atm.getMOL2Type()=="C.cat")
    {
        unsigned int N_N=0;
        for (unsigned int i=0;i < atm.getNumBond();i++) if (atm.getAtomLinked(i).isNitrogen())N_N++;
        if (N_N==3) atm.setFormalCharge(+1);

    }
    if (atm.getMOL2Type()=="H")
    {
        if (atm.getNumBond()!=1)
        {
            cerr << atm.toString()<<" has not bond"<<endl;
            continue;
        }
        if (useGrid && atm.getAtomLinked(0).getBox(grid)!=(Box*)NULL)
        {
           grid->getAdjacentAtoms(adjacentList,atm.getAtomLinked(0),3);
           for (ItAtom itA2 = adjacentList.begin();
                       itA2!= adjacentList.end();
                ++itA2)
           {
               Atom& atm2 = **itA2;
               if (&atm2==&atm || atm2.getNum() <= atm.getNum())continue;
               if (atm2.fixpos.calcDist(atm.fixpos,1.1) < 1 && !atm2.hasBondWith(atm))
               {
                    ItBond it=std::find(FAILED_BOND.begin(),FAILED_BOND.end(),atm2.getBondWith(atm));
                    if (it == FAILED_BOND.end())
                  if (wVerbose) cerr << verbose<<"ERROR\tHydrogen clash "<<atm.getIdentifier()<<"\t"<<atm2.getIdentifier()<<"\t"<<atm2.fixpos.calcDist(atm.fixpos)<<"\n";
               }
           }

        }
        else
        for (size_t iAtm2 = 0 ; iAtm2 < Atoms.size();++iAtm2)
        {
            Atom& atm2 = *Atoms.at(iAtm2);
            if (&atm2==&atm || atm2.getNum() <= atm.getNum())continue;
            if (atm2.fixpos.calcDist(atm.fixpos,1.1) < 1 && !atm2.hasBondWith(atm))
            {
                 ItBond it=std::find(FAILED_BOND.begin(),FAILED_BOND.end(),atm2.getBondWith(atm));
                 if (it == FAILED_BOND.end())
                if (wVerbose)cerr << verbose<<"ERROR\tHydrogen clash "<<atm.getIdentifier()<<"\t"<<atm2.getIdentifier()<<"\t"<<atm2.fixpos.calcDist(atm.fixpos)<<"\n";
            }
        }
    }
    else
    {

        for (size_t i=0; i <NBHN;i++)
        {
            if (AtomName[i].resName!=atm.getResidu()->getName())continue;
            if (AtomName[i].atmName!=atm.getName())continue;
            NPos=0;
            for (size_t iBd = 0; iBd < atm.getNumBond(); iBd++)

            {
                Atom &atmL=atm.getBond(iBd)->getOtherAtom(atm);
                if (atmL.getName() != "H") continue;
                switch (NPos)
                {
                 case 0:if (AtomName[i].H1 != "") atmL.setName(AtomName[i].H1); else if (wVerbose) cerr <<verbose<<"ERRORS\tUnknown hydrogen name "<< atm.getName()<<" " << atm.getResidu()->getName()<<" " << atm.getResidu()->getChainName()<< " " <<atmL.getIdentifier()<<" " << NPos<<endl;break;
                 case 1:if (AtomName[i].H2 != "") atmL.setName(AtomName[i].H2); else if (wVerbose)cerr <<verbose<<"ERRORS\tUnknown hydrogen name "<< atm.getName()<<" " << atm.getResidu()->getName()<<" " << atm.getResidu()->getChainName()<< " " <<atmL.getIdentifier()<<" " << NPos<<endl;break;
                 case 2:if (AtomName[i].H3 != "") atmL.setName(AtomName[i].H3); else if (wVerbose)cerr <<verbose<<"ERRORS\tUnknown hydrogen name "<< atm.getName()<<" " << atm.getResidu()->getName()<<" " << atm.getResidu()->getChainName()<< " " <<atmL.getIdentifier()<<" " << NPos<<endl;break;
                }
                NPos++;
            }
            break;
        }
    }
    atm.assignBits();
}
for (size_t iBd = 0; iBd < Bonds.size();++iBd) Bonds.at(iBd)->assignBits();

if (wVerbose) cout << verbose <<"END\t\n";


}// END addHydrogen



void Molecule::setUse(const bool& use)
{
    const size_t nChain=Chains.size();
    for (size_t iCh=0; iCh < nChain;++iCh)
    {
        Chain &chain= *Chains.at(iCh);
        chain.setUse(use);
    }
}
