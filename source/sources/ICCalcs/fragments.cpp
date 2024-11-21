#include <iostream>
#include "headers/ICCalcs/fragments.h"
#include "headers/ICMole/cycle.h"
#include "headers/ICMole/bond.h"
#include "headers/ICPars/molewriter.h"
#include "headers/ICMole/complex.h"
using namespace std;
using namespace ICMole;

FragGenerator::FragGenerator(Complex& complex):
   complex(complex)
{
listFragments=(Fragment*)NULL;NCL=0;
Interactions inters(complex);
InterResults intersR;
inters.calcInteractions(*complex.getMole(MoleType::LIGAND),intersR);
                        inters.interToMOL2(intersR,false,false,false,true);
                        MoleWriter mw("xx.mol2");
                        mw.writeMOL2(&intersR.Ints);
}

FragGenerator::~FragGenerator()
{
    if (listFragments != (Fragment*)NULL)
    delete[] listFragments;
}


bool isSubsti(const Bond& bd,const Atom&look)
{
    vector<const Atom*> AL,checked,newAL;
    checked.push_back(&look);
    vector<const Atom*>::const_iterator ALit;AL.push_back(&bd.getOtherAtom(look));
    //cout << bd.getAtom1().getName()<<"<->"<<bd.getAtom2().getName()<<" => "<<look.getName()<<endl;
    while(!AL.empty())
    {
        newAL.clear();

        for (vector<const Atom*>::const_iterator itA = AL.begin(); itA != AL.end(); itA++)
        {
            checked.push_back(*itA);
            const Atom &atomAL=**itA;
            for (size_t i=0; i < atomAL.getNumBond(); ++i)
            {
                const Atom &atm2 = atomAL.getAtomLinked(i);
                ALit = find(checked.begin(), checked.end(), &atm2);
                if (ALit != checked.end())continue;
                if (atm2.props.isInRing() || atm2.getMOL2Type()=="Du") return false;
                newAL.push_back(&atm2);
            }
        }
        AL=newAL;

    }
    return true;

}

Bond* toBreak(BondList &BDL, Bond& bd,const Atom &look,unsigned int &NBondLength)
{

    if (look.props.isHydrophobic() && find(BDL.begin(),BDL.end(),&bd) != BDL.end())
    {
        return &bd;}

    vector<const Atom*> AL,checked,newAL;
    checked.push_back(&bd.getOtherAtom(look));
    vector<const Atom*>::const_iterator ALit;AL.push_back(&look);
    //cout << "T::"<<bd.getAtom1().getName()<<"<->"<<bd.getAtom2().getName()<<" => "<<look.getName()<<endl;
    while(!AL.empty())
    {
        newAL.clear();
        NBondLength++;
        for (vector<const Atom*>::const_iterator itA = AL.begin(); itA != AL.end(); itA++)
        {
            checked.push_back(*itA);
            const Atom &atomAL=**itA;
            for (size_t i=0; i < atomAL.getNumBond(); ++i)
            {
                Bond &Bd = *atomAL.getBond(i);
                if (!Bd.isUsed())continue;
                const Atom &atm2 = atomAL.getAtomLinked(i);
               // cout << "||T::"<<Bd.getAtom1().getName()<<"<->"<<Bd.getAtom2().getName()<<" => "<<atm2.getName()<<endl;
                if (find(checked.begin(), checked.end(), &atm2)!= checked.end())continue;
                if (atm2.isHydrogen())continue;
                if (atm2.props.isHydrophobic() && find(BDL.begin(),BDL.end(),&Bd) != BDL.end()){//cout << "TODEL:2X"<<endl;
                    return &Bd;}
                if (atm2.props.isInRing())
                {
                    if (atm2.props.isHydrophobic() == false) return (Bond*)NULL;
                  //  cout << "TODEL:3X"<<endl;
                    return &Bd;
                }

                ALit = find(checked.begin(), checked.end(), &atm2);
                if (ALit != checked.end())continue;
                newAL.push_back(&atm2);
            }
        }
        AL=newAL;

    }

    return (Bond*)NULL;
}



void FragGenerator::fragmentLigand()
{
    Molecule &ligand = *complex.getMole(MoleType::LIGAND);
    ligand.ringPerception();

    ///////////////////////////////////////////////////////////////////////
    //////////////////// STEP1 - Merging cycles ////////////////////
    ///////////////////////////////////////////////////////////////////////
    map<Cycle*,unsigned int> FusedRing;
    map<Atom*,Cycle*> AtomtoCycle;
    map<Atom*,Cycle*>::iterator acpos1, acpos2;
    BondList PossBreak;
    //cout << ligand.numCycles()<<endl;
    for (size_t i=0; i< ligand.numCycles();++i)
    {

        Cycle &cy1 = **(ligand.firstCycle()+i);
        bool hasBeenMerged=false;
        for (unsigned int j=i+1; j<ligand.numCycles(); j++)
        {
            Cycle &cy2 = **(ligand.firstCycle()+j);
            vector<const Atom*> AList;// Contains list of fused atoms between cy1 && cy2
            bool merged=false;
            for (ItCAtom itAcy1 = cy1.first(); itAcy1 != cy1.end(); ++itAcy1)
            {
                const Atom& atom1 = **itAcy1;
                for (ItCAtom itAcy2 = cy2.first(); itAcy2 != cy2.end(); ++itAcy2)
                {
                    const Atom& atom2 = **itAcy2;
                    if (&atom1 == &atom2){merged=true;AList.push_back(&atom1);break;}

                }
                if (merged)break;
            }
            if (merged)
            {
                hasBeenMerged=true;
                for (ItCAtom itAcy2 = cy2.first(); itAcy2 != cy2.end(); ++itAcy2)
                {
                    Atom& atom2 = **itAcy2;
                    vector<const Atom*>::const_iterator itA = find(AList.begin(),AList.end(),&atom2);// We check that the atom of cy2 is not a fused atom
                    if (itA == AList.end()) cy1.addAtom(&atom2);// so we can add it
                }

                j=i;
                ligand.delCycle(&cy2);
            }
        }
        if (hasBeenMerged) FusedRing.insert(pair<Cycle*,unsigned int> (&cy1,2));
        else              FusedRing.insert(pair<Cycle*,unsigned int> (&cy1,1));
        for (ItCAtom itAcy1 = cy1.first(); itAcy1 != cy1.end(); ++itAcy1)
        {
            Atom& atom1 = **itAcy1;
            AtomtoCycle.insert(pair<Atom*,Cycle*>(&atom1,&cy1));
        }
    }

    {
        bool isLinkO=false;
        for (ItCAtom itA=ligand.firstAtom();itA != ligand.lastAtom();++itA)
        {
            Atom &atm = **itA;
            if (atm.isSulfur())
            {


                for (size_t i=0; i < atm.getNumBond();++i)
                {
                    if (atm.getAtomLinked(i).isOxygen())isLinkO=true;
                }
                if (isLinkO) atm.props.setHydrophobic(false);
            }
        }
    }



    ///////////////////////////////////////////////////////////////////////
    //////////////////// STEP1.4 - Bond ruels ////////////////////
    ///////////////////////////////////////////////////////////////////////
    for (size_t ibd= 0; ibd < ligand.numBonds();++ibd)
    {
        Bond& bond = const_cast<Bond&>(ligand.getBond(ibd));
        if (bond.getBondType() != BondType::SINGLE)continue;
        if (bond.getAtom1().props.isInRing() && bond.getAtom2().props.isInRing())
        {
            PossBreak.push_back(&bond);continue;
        }
        if (bond.getAtom1().isHydrogen()
                || bond.getAtom2().isHydrogen())   continue;

        if ((!bond.getAtom1().props.isHydrophobic() && !bond.getAtom1().props.isInRing()) &&
                (!bond.getAtom2().props.isHydrophobic() && !bond.getAtom2().props.isInRing()))continue;

        if ((bond.getAtom1().props.isInRing() && !bond.getAtom2().props.isHydrophobic())
                ||(bond.getAtom2().props.isInRing() && !bond.getAtom1().props.isHydrophobic()))continue;
      /*  cout << bond.getAtom1().getName()<<"::"<<bond.getAtom1().props.toString()<<"\t"
             << bond.getAtom2().getName()<<"::"<<bond.getAtom2().props.toString()<<"\t"<<endl;*/
        PossBreak.push_back(&bond);
    }



    ostringstream oZname;unsigned int NDel=0;
    unsigned int NBondLength=0;

    for (size_t i=0; i< ligand.numCycles();++i)
    {

        Cycle &cyLOOK = **(ligand.firstCycle()+i);
       // cout << "####################################"<<endl;

        cyLOOK.uniqueAtoms();
        for (ItCAtom itAcy1 = cyLOOK.first(); itAcy1 != cyLOOK.end(); ++itAcy1)
        {


        //cout <<  "Number of bonds linked to this cycle : "<< bdLinkedToCycle.size()<<endl;
            const Atom& atomLOOK = **itAcy1;
          //  cout << "#########"<< atomLOOK.getIdentifier()<<"#########"<<endl;
            for (size_t iBd=0; iBd < atomLOOK.getNumBond();++iBd)
            {
                Bond &bond =*atomLOOK.getBond(iBd);
                if (!bond.isUsed())continue;
                Atom&atomBD=bond.getOtherAtom(atomLOOK);
                if (atomBD.getMOL2Type()=="Du" || atomBD.isHydrogen())continue;
               // cout << "LINKED TO : "<< atomBD.getIdentifier()<<"::"<<atomBD.props.toString()<<endl;
                if (atomBD.props.isInRing() && AtomtoCycle.at(&atomBD)==&cyLOOK)continue;


                bool in=true,sa=false;
                BondList PossBreak_G=PossBreak;
                Bond *bdL=(Bond*)NULL;
                do
                {
                    in=true;sa=false;
                    NBondLength=1;
                    bdL=(Bond*)NULL;
                    bdL = toBreak(PossBreak_G,bond,atomBD,NBondLength);

                    if (bdL == (Bond*)NULL || find(PossBreak_G.begin(),PossBreak_G.end(), bdL) == PossBreak_G.end()){//cout << "NULL"<<endl;
                        break;}
                    else
                    {
                        if (isSubsti(bond,atomLOOK) && NBondLength <= 3){
                        //    cout << "SUBSTITUENT"<<endl;
                            PossBreak_G.erase(find(PossBreak_G.begin(),PossBreak_G.end(),bdL));
                            continue;}

                        if (bdL->getAtom1().props.isInRing() && bdL->getAtom2().props.isInRing())
                        {in=false;

                        }
                        else if (bdL == &bond) {in =false; sa=true;}

                        acpos1 = AtomtoCycle.find(&bdL->getAtom1());
                        acpos2= AtomtoCycle.find(&bdL->getAtom2());
                        if (acpos1 != AtomtoCycle.end() && acpos2 != AtomtoCycle.end())
                        {
                            cout << "CYCLES : "<< (*acpos1).second << " " <<(*acpos2).second<<endl;
                            if ((*acpos1).second == (*acpos2).second)continue;
                        }
                        NDel++;
                        oZname.str("");
                        oZname<< "Z"<<NDel;
                     //  cout << bdL->toString()<<endl;
                     //  cout << "DELETION OF : "<< bdL->getAtom1().getIdentifier()<<"<->"<<bdL->getAtom2().getIdentifier()<<endl;
                    //   Bond& bd1=ligand.addBond(ligand.addAtom("",bdL->getAtom1().fixpos,oZname.str(),"Du",bdL->getAtom1().getResidu()),bdL->getAtom2(),BondType::SINGLE);
                    //   Bond &bd2=ligand.addBond(ligand.addAtom("",bdL->getAtom2().fixpos,oZname.str(),"Du",bdL->getAtom2().getResidu()),bdL->getAtom1(),BondType::SINGLE);
                     //  cout << "NEW BD:"<<bd1.getAtom1().getIdentifier()<<"<->"<<bd1.getAtom2().getIdentifier()<<endl;
                     // cout << "NEW BD:"<<bd2.getAtom1().getIdentifier()<<"<->"<<bd2.getAtom2().getIdentifier()<<endl;
                       bdL->setUse(false,false);
                      /* ostringstream ossT;ossT<<"frag"<<tmpi<<".mol2";
                       tmpi++;
                       MoleWriter mw(ossT.str());mw.setOnlyUsed(true);mw.writeMOL2(&ligand);*/
                       anchorBond.push_back(bdL);
                        if(sa) iBd=0;
                        //}
                    }
                }while(in);

            }


        }

   }





    while(ligand.numCycles() != 0)
    {
        ligand.delCycle(*ligand.firstCycle());
    }
/*MoleWriter mw("outX.mol2");
mw.setOnlyUsed(true);
mw.writeMOL2(&ligand);*/



    vector<const Atom*> al,new_al,checked;
    for (size_t i=0; i < ligand.numAtom();++i)
    {

        NCLUS.insert(pair<unsigned int,unsigned int>(ligand.getAtom(i).getNum(),0));
    }

    ligand.checkMOL2();
    for (size_t i=0; i < ligand.numAtom();++i)
    {
        const Atom &atm = ligand.getAtom(i);
       // cout << atm.getIdentifier()<<" " << atm.getMOL2Type()<<endl;
       // if (atm.getMOL2Type()!="Du")continue;

        if (NCLUS.at(atm.getNum()) != 0)continue;

        checked.clear();al.clear();
        checked.push_back(&atm);
        NCLUS.at(atm.getNum())=NCL;
        for (size_t iBd=0; iBd < atm.getNumBond();++iBd)
        {
            if (!atm.getBond(iBd)->isUsed())continue;
            al.push_back(&atm.getAtomLinked(iBd));

        }

        while(!al.empty())
        {
            new_al.clear();
            for (vector<const Atom*>::const_iterator itA = al.begin(); itA != al.end(); itA++)
            {
                const Atom &atomA=**itA;
                checked.push_back(&atomA);
               // cout << "ADDING : "<< atomA.getIdentifier()<<" TO "<<NCL<<endl;
                NCLUS.at((*itA)->getNum())=NCL;
                for (size_t iBd=0; iBd < atomA.getNumBond();++iBd)
                {
                    if (!atomA.getBond(iBd)->isUsed())continue;
                    const Atom& atm2=atomA.getAtomLinked(iBd);
                    if (find(checked.begin(), checked.end(), &atm2) != checked.end())continue;
                    new_al.push_back(&atm2);
                }
            }
            al=new_al;
        }
        NCL++;

    }
}

void FragGenerator::calcInteractions()
{
    Molecule &ligand = *complex.getMole(MoleType::LIGAND);
    for (size_t iAtm=0; iAtm< ligand.numAtom();++iAtm)
    {
        Atom &atm = const_cast<Atom&>(ligand.getAtom(iAtm));
        atm.props.clear();;
    }
    ligand.ringPerception();
    Interactions inters(complex);


 listFragments=new Fragment[NCL];
 bool inCycle=false;
 for (unsigned int i=0; i< NCL;++i)
 {
     listFragments[i].wCycle=false;
listFragments[i].id=i+1;
inCycle=false;
    for (size_t iAtm=0; iAtm< ligand.numAtom();++iAtm)
    {
        Atom &atm = const_cast<Atom&>(ligand.getAtom(iAtm));
        if (NCLUS.at(atm.getNum()) != i){atm.setUse(false,false,false);continue;}
        atm.setUse(true,false,false);
        if (atm.props.isInRing())inCycle=true;
        atm.checkMOL2type();
        listFragments[i].atoms.push_back(&atm);
        for (size_t iBd=0; iBd < atm.getNumBond();++iBd)
        {
            Bond &bond =*atm.getBond(iBd);
            if (bond.isUsed())continue;
            listFragments[i].anchor.push_back(&bond);
            ItBond pos = find(anchorBond.begin(),anchorBond.end(),&bond);
            listFragments[i].anchorID.push_back(std::distance(anchorBond.begin(),pos)+1);
        }
    }

    cout << "########## FRAGMENT " << listFragments[i].id<< " ##########"<<endl;
    for (ItAtom itA = listFragments[i].atoms.begin(); itA != listFragments[i].atoms.end();++itA)
    {
        Atom &atm = **itA;
        cout << atm.getIdentifier()<<endl;
    }
    for (size_t k=0; k< listFragments[i].anchor.size();++k)
    {
        cout << listFragments[i].anchor.at(k)->toString()<<" " << listFragments[i].anchorID.at(k)<<endl;
    }
    if (!inCycle)
    {
continue;
    }
    listFragments[i].wCycle=true;
    inters.calcInteractions(ligand,listFragments[i].inters);
    inters.interToMOL2(listFragments[i].inters,false,false,false,true);
    cout << inters.toString(listFragments[i].inters);

 }



}



void FragGenerator::saveFragments(const string &name)
{
    ostringstream oss;
    MoleWriter mw;
    Molecule &ligand = *complex.getMole(MoleType::LIGAND);
    AtomList anchors;
    const size_t nAtms=ligand.numAtom();
    unsigned short curr_FRAG=1;
    unsigned short NHeavy=0,NLink=0,NRing=0,NArRing=0;
    cout << name<<"|FRAG\tID\tNHeavy\tNLink\tNRing\tNArRing\tTOT\tCA\tCZ\tO\tOD1\tN\tNZ\tZN"<<endl;

    for (unsigned int i=1; i< NCL;++i)
    {
        anchors.clear();
        if (!listFragments[i].wCycle)continue;

        NHeavy=0;NLink=0;
        for (size_t iAtm=0; iAtm< nAtms;++iAtm)
        {
            Atom &atm = const_cast<Atom&>(ligand.getAtom(iAtm));
            if (NCLUS.at(atm.getNum()) != i){atm.setUse(false,false,false);continue;}
            atm.setUse(true,false,false);
            if (!atm.isHydrogen())NHeavy++;

            for (size_t k=0; k< listFragments[i].anchor.size();++k)
            {   Bond & bd=*listFragments[i].anchor.at(k);
                if (find(atm.first(),atm.last(),&bd) != atm.last())
                {
                    Atom& atmL=bd.getOtherAtom(atm);
                    oss.str("");oss<<"Z"<<listFragments[i].anchorID.at(k);
                    Atom& Z=ligand.addAtom("Du",atmL.fixpos,oss.str(),"Du",atm.getResidu());
                    //cout << Z.toString();
                    ligand.addBond(Z,atm,BondType::SINGLE);
                    anchors.push_back(&Z);
                    NHeavy--;
                    NLink++;
                }

            }
        }
        NRing=0;
        NArRing=0;
        for (ItCCycle itCC=ligand.firstCycle();itCC != ligand.lastCycle();++itCC)
        {
            bool allIn=true;
            const Cycle& cycle = **itCC;
            for (ItCAtom itA = cycle.first(); itA != cycle.end();++itA)
            {
                const Atom& atm =**itA;
                if (atm.isUsed()) continue;
                allIn=false;
                break;
            }
            if (!allIn)continue;
            NRing++;
            if (cycle.isAromatic())NArRing++;

        }


        oss.str("");
        if (!name.empty()) oss <<name<<"_FRAG_"<<  curr_FRAG;
        else oss <<"FRAG_"<<  curr_FRAG;
mw.newFile(oss.str()+"_MOLE.mol2");
mw.setOnlyUsed(true);
mw.writeMOL2(&ligand);
        for (ItAtom itA = anchors.begin(); itA != anchors.end();++itA)
        {
            ligand.delAtom(*itA,false);
        }
        oss.str("");
        if (!name.empty())oss<<name<<"_FRAG_"<< curr_FRAG<<"_INTS.mol2";
        else oss<<"FRAG_"<< curr_FRAG<<"_INTS.mol2";
        mw.setOnlyUsed(false);
       // cout << listFragments[i].inters.Ints.toString()<<endl;
        mw.newFile(oss.str());
        mw.writeMOL2(&listFragments[i].inters.Ints);



        unsigned short CA,CZ,O,OD1,N,NZ,Zn;
        CA=0;CZ=0;O=0;OD1=0;N=0;NZ=0;Zn=0;
        for (size_t iAtm=0;iAtm < listFragments[i].inters.Ints.numAtom();++iAtm)
        {
            const Atom& atm =  listFragments[i].inters.Ints.getAtom(iAtm);
            if (atm.getResiduName().substr(2,1) != "C")continue;
            if (atm.getName()=="CA")CA++;
            else if (atm.getName()=="CZ")CZ++;
            else if (atm.getName()=="O")O++;
            else if (atm.getName()=="OD1")OD1++;
            else if (atm.getName()=="N")N++;
            else if (atm.getName()=="NZ")NZ++;
            else if (atm.getName()=="Zn")Zn++;

        }
        cout << name<<"|FRAG\t"<<curr_FRAG<<"\t"<<NHeavy<<"\t"<<NLink<<"\t"
             << NRing<<"\t"<<NArRing<<"\t"
             << (CA+CZ+O+OD1+N+NZ+Zn)<<"\t"
             << CA<<"\t"
             << CZ<<"\t"
             << O<<"\t"
             << OD1<<"\t"
             << N<<"\t"
             << NZ<<"\t"
             << Zn<<"\n";

curr_FRAG++;




    }
}

