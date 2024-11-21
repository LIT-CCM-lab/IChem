#include "headers/ICMole/complex.h"
#include "headers/ICMole/bond.h"
#include "headers/ICMole/cycle.h"
#include "headers/ICPars/molewriter.h"
#include "headers/ICCalcs/volsite.h"
#include "headers/ICMole/box.h"
using namespace std;
using namespace ICMole;

//#define ICHEM_DEBUG
Complex::Complex()
{


    for (size_t i=0;i < NB_MOLETYPE;++i)
    {

        lastmolecule[i]=(Molecule*)NULL;
    }
    uptoGrid=false;
    for (unsigned int i=0; i < N_GRID_CPLX; i++) boxStep[i]=-1;
}

/**
 * @brief Complex::~Complex
 * Delete all molecules within it
 */
Complex::~Complex()
{

    for (size_t i=0;i < NB_MOLETYPE;++i)
    {
        for (ItMole it = molecules[i].begin(); it != molecules[i].end();++it)
            delete (*it);
    }


}


void Complex::clearGrid()
{

    for (unsigned int i=0; i < N_GRID_CPLX; i++)
    {
        if (boxStep[i] != -1) grid[i].clear();
        boxStep[i]=-1;
    }
}



/**
 * @brief Complex::addMolecule
 * @param molecule: molecule to add
 *
 * Add the given molecule to this complex. From now own,
 * the complex owe the molecule. Still, all atoms, bonds, residus, chains
 * will be stored (but not own) within the complex.
 *
 */
void Complex::addMolecule(Molecule *const molecule) throw(MoleExcept)
{
    allmolecules.push_back(molecule);
    molecules[molecule->getMoleType()].push_back(molecule);
    lastmolecule[molecule->getMoleType()]=molecule;
    molecule->setComplex(this);
    for (ItCAtom  itA=molecule->firstAtom()  ;itA!=molecule->lastAtom()  ;++itA){  atomlist.push_back(*itA);  }
    for (ItCRes   itR=molecule->firstResidu();itR!=molecule->lastResidu();++itR){residulist.push_back(*itR);  }
    for (ItCChain itC=molecule->firstChain() ;itC!=molecule->lastChain() ;++itC){ chainlist.push_back(*itC);  }
    for (ItCBond  itB=molecule->firstBond()  ;itB!=molecule->lastBond()  ;++itB){  bondlist.push_back(*itB);  }
    uptoGrid=false;
}









/**
  * @brief Complex::getNumMolecule
  * @param MoleType : Type of molecule you want the count of
  * @return Number of molecuels of this molecular type
  */
size_t Complex::getNumMolecule(const unsigned int &MoleType) const throw(MoleExcept)
{
    if (MoleType >NB_MOLETYPE )
        throw MoleExcept(1100101,
                         "Complex::getNumMolecule"
                         "Given MoleType is above the maximum");
    return molecules[MoleType].size();
}




/**
 * @brief Complex::getNumMolecule
 * @return The number of molecules within this complex
 */
size_t Complex::getNumMolecule() const
{
    return allmolecules.size();

}





/**
 * @brief Complex::clear
 * Delete all molecules within this complex
 */
void  Complex::clear() throw (MoleExcept)
{
    for (ItMole it = allmolecules.begin(); it != allmolecules.end();++it)
        delete (*it);

    for (size_t i=0;i < NB_MOLETYPE;++i)
    {
        molecules[i].clear();
        lastmolecule[i]=(Molecule*)NULL;
    }
    atomlist.clear();
    bondlist.clear();
    residulist.clear();
    chainlist.clear();
    uptoGrid=false;
}


void  Complex::deleteMole(Molecule *const molecule) throw(MoleExcept)
{
    if (molecule == (Molecule*)NULL)
        throw MoleExcept(1100701,
                         "Complex::updateMoleType"
                         "No molecule given");

    ItMole it = find(allmolecules.begin(),
                     allmolecules.end(),
                     molecule);
    if (it == allmolecules.end())
        throw MoleExcept(1100702,
                         "Complex::updateMoleType"
                         "Given molecule is not part of this complex");
    allmolecules.erase(it);
    it = find(molecules[molecule->getMoleType()].begin(),
            molecules[molecule->getMoleType()].end(),
            molecule);
    if (it == molecules[molecule->getMoleType()].end())
    {

        for (size_t i=0;i < NB_MOLETYPE;++i)
        {
            it = find(molecules[i].begin(),
                      molecules[i].end(),
                      molecule);
            if (it != molecules[i].end())molecules[i].erase(it);
        }

    }
    else molecules[molecule->getMoleType()].erase(it);
    //    delete molecule;
    uptoGrid=false;
    loadData();
}


void  Complex::removeMole(Molecule *const molecule) throw(MoleExcept)
{
    if (molecule == (Molecule*)NULL)
        throw MoleExcept(1100701,
                         "Complex::updateMoleType"
                         "No molecule given");

    ItMole it = find(allmolecules.begin(),
                     allmolecules.end(),
                     molecule);
    if (it == allmolecules.end())
        throw MoleExcept(1100702,
                         "Complex::updateMoleType"
                         "Given molecule is not part of this complex");
    allmolecules.erase(it);
    it = find(molecules[molecule->getMoleType()].begin(),
            molecules[molecule->getMoleType()].end(),
            molecule);
    if (it == molecules[molecule->getMoleType()].end())
    {

        for (size_t i=0;i < NB_MOLETYPE;++i)
        {
            it = find(molecules[i].begin(),
                      molecules[i].end(),
                      molecule);
            if (it != molecules[i].end())molecules[i].erase(it);
        }

    }
    else molecules[molecule->getMoleType()].erase(it);
    uptoGrid=false;
    loadData();


}

bool Complex::isMoleIn(Molecule* const molecule) const
{
    const MoleList &listM= molecules[molecule->getMoleType()];
    ItCMole it = find(listM.begin(),
                      listM.end(),
                      molecule);
    if (it == molecules[molecule->getMoleType()].end()) return false;
    return true;

}
bool Complex::isMoleIn(const Molecule&  molecule) const
{
    const MoleList &listM= molecules[molecule.getMoleType()];
    ItCMole it = find(listM.begin(),
                      listM.end(),
                      &molecule);
    if (it == molecules[molecule.getMoleType()].end()) return false;
    return true;

}


void Complex::loadData(const bool& only_used)
{
    size_t mtype=0;
    atomlist.clear();
    bondlist.clear();
    residulist.clear();
    chainlist.clear();
    for (;mtype < NB_MOLETYPE; mtype++)
        for (ItCMole it = molecules[mtype].begin();
             it !=molecules[mtype].end();
             it++)
        {
            const Molecule &mole=**it;
            for (ItCAtom itA=mole.firstAtom();
                 itA!=mole.lastAtom();
                 ++itA)
            {
                Atom &atm=**itA;
                if (only_used && !atm.isUsed() ) continue;
                atomlist.push_back(&atm);
            }
            for (ItCRes itR=mole.firstResidu();
                 itR!=mole.lastResidu();
                 ++itR)
            {
                Residu &res=**itR;
                if (only_used && !res.isUsed() ) continue;
                residulist.push_back(&res);
            }
            for (ItCBond itB=mole.firstBond();
                 itB!=mole.lastBond();
                 ++itB)
            {
                Bond &bond=**itB;
                if (only_used && !bond.isUsed() ) continue;
                bondlist.push_back(&bond);
            }
            for (ItCChain itC=mole.firstChain();
                 itC!=mole.lastChain();
                 ++itC)
            {
                Chain &chain=**itC;
                if (only_used && !chain.isUsed() ) continue;
                chainlist.push_back(&chain);
            }
        }
}



void Complex::addAtom( Atom* const atom) throw(MoleExcept)
{
    if (atom == (Atom*)NULL) throw MoleExcept(1100301,
                                              "Complex::addAtom"
                                              "No Atom given");
    atomlist.push_back(atom);
    uptoGrid=false;
}


void Complex::addBond( Bond* const bond) throw(MoleExcept)
{
    if (bond == (Bond*)NULL) throw MoleExcept(1100401,
                                              "Complex::addBond"
                                              "No Bond given");
    bondlist.push_back(bond);
    uptoGrid=false;
}

void Complex::addResidu( Residu* const residu) throw(MoleExcept)
{
    if (residu == (Residu*)NULL) throw MoleExcept(1100501,
                                                  "Complex::addResidu"
                                                  "No Residu given");
    residulist.push_back(residu);
    uptoGrid=false;
}


void Complex::addChain( Chain* const chain) throw(MoleExcept)
{
    if (chain == (Chain*)NULL) throw MoleExcept(1100601,
                                                "Complex::addChain"
                                                "No Chain given");
    chainlist.push_back(chain);
    uptoGrid=false;
}


void Complex::deleteAtom(const Atom* const atom) throw(MoleExcept)
{
    if (atom == (Atom*)NULL) throw MoleExcept(1100301,
                                              "Complex::addAtom"
                                              "No Atom given");
    const ItAtom pos=find(atomlist.begin(),atomlist.end(),atom);
    if (pos == atomlist.end())
        throw MoleExcept(1100302,
                         "Complex::addAtom"
                         "Atom not in complex");
    atomlist.erase(pos);
    uptoGrid=false;
}


void Complex::deleteBond(const Bond * const bond) throw(MoleExcept)
{
    if (bond == (Bond*)NULL) throw MoleExcept(1100301,
                                              "Complex::addAtom"
                                              "No Atom given");
    const ItBond pos=find(bondlist.begin(),bondlist.end(),bond);
    if (pos ==bondlist.end())
        throw MoleExcept(1100302,
                         "Complex::addAtom"
                         "Atom not in complex");
    bondlist.erase(pos);
    uptoGrid=false;
}


void Complex::deleteResidu(const Residu* const residu) throw(MoleExcept)
{
    if (residu == (Residu*)NULL) throw MoleExcept(1100301,
                                                  "Complex::addAtom"
                                                  "No Atom given");
    const ItRes pos=find(residulist.begin(),residulist.end(),residu);
    if (pos ==residulist.end())
        throw MoleExcept(1100302,
                         "Complex::addAtom"
                         "Atom not in complex");
    residulist.erase(pos);
    uptoGrid=false;
}
void Complex::deleteChain (const Chain * const chain ) throw(MoleExcept)
{
    if (chain == (Chain*)NULL) throw MoleExcept(1100301,
                                                "Complex::addAtom"
                                                "No Atom given");
    const ItChain pos=find(chainlist.begin(),chainlist.end(),chain);
    if (pos ==chainlist.end())
        throw MoleExcept(1100302,
                         "Complex::addAtom"
                         "Atom not in complex");
    chainlist.erase(pos);
    uptoGrid=false;
}

void Complex::updateMoleType(Molecule* const molecule,const unsigned int & type)
{
    if (molecule == (Molecule*)NULL)
        throw MoleExcept(1100701,
                         "Complex::updateMoleType"
                         "No molecule given");

    ItMole it = find(allmolecules.begin(),
                     allmolecules.end(),
                     molecule);
    if (it == allmolecules.end())
        throw MoleExcept(1100702,
                         "Complex::updateMoleType"
                         "Given molecule is not part of this complex");
    it = find(molecules[molecule->getMoleType()].begin(),
            molecules[molecule->getMoleType()].end(),
            molecule);
    if (it == molecules[molecule->getMoleType()].end())
    {
        if (Moleaccess!= Levels::NONE)
        {
            throw MoleExcept(1100703,
                             "Complex::updateMoleType"
                             "Given molecule is not associated with the good MoleType");
        }else
        {
            for (size_t i=0;i < NB_MOLETYPE;++i)
            {
                it = find(molecules[i].begin(),
                          molecules[i].end(),
                          molecule);
                if (it != molecules[i].end())molecules[i].erase(it);
            }
        }
    }
    else molecules[molecule->getMoleType()].erase(it);
    molecules[type].push_back(molecule);
}






Grid& Complex::genGrid(const double& boxValue,const bool& wHydrogen) throw(MoleExcept)
{
    // Position in the grid Array :
    unsigned int gr_pos=0;

    for (gr_pos=0;gr_pos<=N_GRID_CPLX; gr_pos++)
    {
        if (gr_pos== N_GRID_CPLX)
            throw MoleExcept(1100801,
                             "Complex::genGrid",
                             "Too much grid generated");
        if (boxStep[gr_pos] == -1)break;
        if (boxStep[gr_pos]== boxValue) return grid[gr_pos];
        if( fmod(boxValue,boxStep[gr_pos]) != 0)
            throw MoleExcept(1100802,
                             "Complex::genGrid",
                             "Given boxStep value is not a factor of existing grid");

    }

    if (gr_pos==0)
    {
        grid[gr_pos].createComplexBoxes(*this,boxValue);
        grid[gr_pos].setWHydrogen(wHydrogen);
        AtomList toRotate;
        for (ItCMole itM =allmolecules.begin(); itM != allmolecules.end(); ++itM)
        {
            grid[gr_pos].colorCube(**itM);

            for (ItCCycle itLC= (*itM)->firstCycle();itLC!= (*itM)->lastCycle();itLC++)
            {
                Cycle& ligcycle=**itLC;
                Atom &atomL=ligcycle.getCenter();
                toRotate.push_back(&atomL);
            }
        }

        grid[gr_pos].rotateAtoms(toRotate);

    }
    else grid[gr_pos].createParentGrid(grid[gr_pos-1],boxValue);
    boxStep[gr_pos]=boxValue;
    uptoGrid=true;
    return grid[gr_pos];

}





int Complex::getMaxGrid() const
{
    for (int gr_pos=0;gr_pos<=N_GRID_CPLX; gr_pos++)
    {
        if (boxStep[gr_pos] != -1)
            return gr_pos;
    }
    return -1;

}




void Complex::splitMoleculeInToComplex(Molecule& molecule,const std::string &verbose) throw(ICMole::MoleExcept)
{
    const bool wVerbose= !verbose.empty();

    Molecule *newMolecule=(Molecule*)NULL;
    if (wVerbose) cout << verbose<<"START\t"<<endl;
    // Moving water molecules out
    size_t nlig=0;
    for (size_t iRes = 0; iRes <molecule.numResidus();++iRes)
    {
        Residu &residu = molecule.getResidu(iRes);
        const size_t nResAtm = residu.numAtom();
        const unsigned int& ResType=residu.getResType();
        newMolecule=(Molecule*)NULL;

        try
        {
            if (ResType == ResType::WATER)
            {
                newMolecule= new Molecule(3,2,MoleType::WATER);
                molecule.moveResidu(residu,*newMolecule);
                if (wVerbose) cout << verbose<<"MODIFS\t"<<residu.getIdentifier()<<"\tWATER\t"<< residu.getNInterResidu()<<endl;
                // cout << watermole->toString()<<endl;
                addMolecule(newMolecule);
                --iRes;
            }
            else if ((ResType == ResType::COFACTOR ||
                      ResType == ResType::UNWANTED ||
                      ResType == ResType::LIGAND)
                     && nResAtm > 1 && residu.getNInterResidu()==0)
            {

                if (wVerbose) cout << verbose<<"MODIFS\t"<<residu.getIdentifier()<<"\tCOFACTOR\t"<< residu.getNInterResidu()<<endl;

                newMolecule = new Molecule(nResAtm,10,(ResType==ResType::COFACTOR)?MoleType::COFACTOR:MoleType::LIGAND);
                newMolecule->setName(residu.getName());
                molecule.moveResidu(residu,*newMolecule);
                addMolecule(newMolecule);
                ++nlig;

                --iRes;
            }


        }catch (MoleExcept &e)
        {
            if (newMolecule != (Molecule*)NULL) delete newMolecule;
            if (wVerbose)cerr << verbose<<"ERRORS|\t"<< e.getCode()<<"\t"<<e.getData()<<endl;
            else cerr << e.getCode()<<"\n"
                      << e.getData()<<"\n"
                      << e.getSource()<<"\n"
                      << e.getTrace()<<endl;
        }
    }//END iRes



    // Scanning chains to see if we have peptides :
    double N_RES[NB_RESTYPE];
    const size_t nChain= molecule.numChains();
    vector<Residu*> reslist;
    for (size_t iChain=0; iChain<nChain;++iChain)
    {
        for (size_t j=0; j < NB_RESTYPE;j++) N_RES[j]=0;
        const Chain &chain= molecule.getChain(iChain);
        const size_t nRes =chain.getCountRes();
        size_t nAtm=0;
        reslist.clear();;
        for (size_t iRes = 0; iRes < nRes;++iRes)
        {
            Residu& res = chain.getResidu(iRes);
            nAtm+=res.numAtom();
            N_RES[res.getResType()]++;
            reslist.push_back(&res);
        }

        // Check for peptide:
        if (nRes <= 8 && N_RES[ResType::STD_AA]+N_RES[ResType::MOD_AA]>= (double)nRes/2)
        {
            newMolecule = new Molecule(nAtm,10,MoleType::LIGAND);
            molecule.moveResidu(reslist,*newMolecule);
            addMolecule(newMolecule);
            ++nlig;
            ostringstream osname;
            for (size_t iRes = 0; iRes < nRes;++iRes)
            {
                Residu& res = chain.getResidu(iRes);
                osname<<"_"<<res.getName();
            }
            newMolecule->setName(osname.str());
        }
    }// END iChain




    multimap<Residu*,Residu*> mapres;
    vector<set<const Residu*> > clustres;

    unsigned int ncount=0;
    const size_t nBond = molecule.numBonds();
    for (size_t iBd = 0; iBd < nBond; ++iBd)
    {
        const Bond& bond = molecule.getBond(iBd);
        const Atom& atom1=bond.getAtom1();
        const Atom& atom2=bond.getAtom2();
        if (atom1.getFNum()< atom2.getFNum() && atom1.getResidu() != atom2.getResidu())
        {
            mapres.insert(pair<Residu*,Residu*>(atom1.getResidu(),atom2.getResidu()));
            ncount++;
        }
    }// END iBd



    set<const Residu*> newset;
    newset.insert(&molecule.getResidu(0)); clustres.push_back(newset);
    bool foundR1=false, foundR2=false;
    for (multimap<Residu*,Residu*>::iterator it = mapres.begin(); it != mapres.end(); it++)
    {

        foundR1=false; foundR2=false;
        Residu *r1 = (*it).first, *r2 = (*it).second;
        for (vector<set<const Residu*> >::iterator Rit = clustres.begin(); Rit != clustres.end(); Rit++)
        {

            for (set<const Residu*>::iterator Sit = (*Rit).begin(); Sit != (*Rit).end(); Sit++)
            {
                if (*Sit == r1) foundR1=true;
                if (*Sit == r2) foundR2=true;
            }
            if (foundR1 && foundR2)break;
            if (foundR1) { (*Rit).insert(r2);break;}
            if (foundR2) { (*Rit).insert(r1);break;}
        }
        if (!foundR1 && !foundR2)
        {
            set<const Residu*> newset;
            newset.insert(r1);newset.insert(r2); clustres.push_back(newset);
        }
    }
    // Merging clusters :
    foundR2=true;
    while(foundR2)
    {
        foundR2=false;
        for (unsigned int i=0; i < clustres.size(); i++){
            for (unsigned int j=i+1; j < clustres.size(); j++)
            {
                foundR1=false;
                for (set<const Residu*>::const_iterator Rit = clustres.at(i).begin(); Rit != clustres.at(i).end(); Rit++)
                {
                    for (set<const Residu*>::const_iterator Rit2 = clustres.at(j).begin(); Rit2 != clustres.at(j).end(); Rit2++)
                    {
                        if (*Rit2 == *Rit) {foundR1=true;break;}
                    }
                    if (foundR1) break;
                }
                if (foundR1){
                    clustres.at(i).insert(clustres.at(j).begin(),clustres.at(j).end());
                    // for (set<Residu*>::iterator Rit2 = clustres.at(j).begin(); Rit2 != clustres.at(j).end(); Rit2++)
                    //{
                    //    clustres.at(i).insert(*Rit2);
                    // }

                    clustres.erase(clustres.begin()+j);

                    foundR2=true;
                    break;
                }

            }
            if (foundR2)break;
        }
    }


    //

    for (vector<set<const Residu*> >::iterator Rit = clustres.begin(); Rit != clustres.end(); Rit++)
    {
        for (size_t i=0; i < NB_RESTYPE;i++) N_RES[i]=0;
        //if (verbose)cout << PDB_name<<"|SPLITM|STATIS\tNEW SET : "<<endl;
        for (set<const Residu*>::iterator Sit = (*Rit).begin(); Sit != (*Rit).end(); Sit++)
        {
            N_RES[(*Sit)->getResType()]++;
            // if (verbose)cout << PDB_name<<"|SPLITM|STATIS\t"<<(*Sit)->getName()<< " ("<<(*Sit)->getChain()->getName()<<"-"<<(*Sit)->getMapNum()<<")"<<endl;
        }
        // if (verbose)cout <<endl;
        const double RSize = (*Rit).size();

        if (RSize>8 )continue;
        if (N_RES[ResType::SUGAR]/RSize >= 0.5 || N_RES[ResType::UNWANTED]/RSize >= 0.5)
        {
            for (set<const Residu*>::iterator Sit = (*Rit).begin(); Sit != (*Rit).end(); Sit++)
            {
                molecule.delResidu(*const_cast<Residu*>(*Sit));
            }
        }
    }
    if (molecule.numResidus() <=2)
    {
        //     cout << PDB_name<<"|SPLITM|WARNIN\tNo molecule left"<<endl;
        throw MoleExcept(999999,"","");
    }

    if(wVerbose)  cout <<verbose<<"STATIS\t" << nlig<<"\tNumber of ligand found "<<"\n"
                      << verbose<<"END\t"<<endl;
    if (nlig==0) throw MoleExcept(2040201,"MOLEConvert::splitmoleincomplex","No Ligand found");

}



std::string Complex::toString() const
{
    ostringstream oss;
    for (size_t i=0; i < NB_MOLETYPE;++i)
    {
        switch (i)
        {
        case MoleType::LIGAND: oss<<"LIGAND";break;
        case MoleType::COFACTOR: oss<<"COFACTOR";break;
        case MoleType::PROTEIN: oss<<"PROTEIN";break;
        case MoleType::SITE: oss<<"SITE";break;
        case MoleType::OTHER: oss<<"OTHER";break;
        case MoleType::WATER: oss<<"WATER";break;
        case MoleType::CAV_ALL: oss<<"CAV_ALL";break;
        case MoleType::CAV_4: oss<<"CAV_4";break;
        case MoleType::CAV_6: oss<<"CAV_6";break;
        case MoleType::CAV_8: oss<<"CAV_8";break;
        case MoleType::CAV_12: oss<<"CAV_12";break;
        case MoleType::INT_LIG: oss<<"INT_LIG";break;
        case MoleType::INT_CENT: oss<<"INT_CENT";break;
        case MoleType::INT_PROT: oss<<"INT_PROT";break;
        case MoleType::INT_MERG: oss<<"INT_MERG";break;
        case MoleType::UNDEFINED: oss<<"UNDEFINED";break;
        case MoleType::NUCLEIC: oss<<"NUCLEIC";break;
        }
        oss << "\t"<<molecules[i].size()<<"\n";
    }
    return oss.str();

}


double calcHScore(const Coords& H,
                  const Coords& AtmCD,
                  const Coords &AtmCA,
                  const Coords &AtmD_L,
                  const double& dist,
                  const bool &water)
{
    const double distH_D = H.calcDist(AtmCD,1.2);
    const double distH_A = H.calcDist(AtmCA,dist+0.1);
    const double angl_HDDl = (AtmCD.calcAngle(AtmD_L,H));
    const double anglAHD = (H.calcAngle(AtmCA,AtmCD))*180/PI;
    if (anglAHD <=130) { return 0;}
    const double alpha = 34 - distH_A*sin( distH_D*sin(angl_HDDl)/dist);
    const double angl_ADDl=((AtmCD.calcAngle(AtmD_L,AtmCA))*180/PI);
    if (angl_ADDl > 109.5+alpha || angl_ADDl < 109.5-alpha) {return 0;}
    const double diffAHD=180-anglAHD;
    if (diffAHD <= 30) {if (water) return 1+(diffAHD)/50; else return 1;}
    else if (diffAHD <= 80) return 1-(diffAHD-30)/50;
    else return 0;

}


struct bestOrientH
{
    double x,y,z,i;
    Coords H1,H2;

};
void Complex::selectWater(const std::string& verbose) throw(MoleExcept)
{
    const bool wVerbose= !verbose.empty();
    if (wVerbose) cout << verbose<<"START"<<endl;
    unsigned short NMod=0;
    // STEP 1 - FINDING WATER RESIDUES:
    const size_t nMole=allmolecules.size();
    ResiduList listWater,waterDel;
    for (size_t iMole =0; iMole < nMole;++iMole)
    {
        Molecule& molecule=*allmolecules.at(iMole);

        const size_t nResidu = molecule.numResidus();
        for (size_t iRes =0; iRes <nResidu;++iRes)
        {
            Residu& residu = molecule.getResidu(iRes);
            if (residu.getResType() != ResType::WATER)continue;
            listWater.push_back(&residu);
        }
    }
    if (wVerbose)  cout << verbose<<"STATIS\t"<<listWater.size()<<"\tNumber of water inital molecules "<<endl;
    Grid& grid=genGrid(1.5);
    bool found=false;

    //STEP 2 ???

    AtomList atmList;
    vector<const Atom*> accHList,donHList,clashList;
    Box* bx;
    for (size_t i=0; i< listWater.size();++i)
    {
        Residu & water = *listWater.at(i);
        bx=(Box*)NULL;
        Atom *tmpAtom=(Atom*)NULL;
        for (size_t iA=0;iA < water.numAtom();++iA)
        {
            Atom &atomO = water.getAtom(iA);
            if (!atomO.isOxygen())continue;
            tmpAtom=&atomO;
        }
        if (tmpAtom ==(Atom*)NULL)continue;
        Atom &atomO= *tmpAtom;
        for (size_t iAx=0; iAx< atomO.getNumBond();++iAx)
        {
            Atom& atmH=atomO.getBond(iAx)->getOtherAtom(atomO);

            water.getParent()->delAtom(&atmH);
        }
        found=false;
        //        accHList.clear();
        //        donHList.clear();
        //        clashList.clear();


        if (atomO.getBox(&grid)==(Box*)NULL)
        {
            bx=grid.getBox(atomO.fixpos);
            if (bx == (Box*)NULL)
            {
                waterDel.push_back(&water);
                continue;
            }

        }
        else bx=atomO.getBox(&grid);
        grid.getAdjacentAtoms(atmList,*bx,6.6);

        for (ItCAtom itCA= atmList.begin();itCA != atmList.end();++itCA)
        {
            const Atom& atom = **itCA;

            if (atom.getParent().getMoleType()==MoleType::LIGAND
                    || atom.getParent().getMoleType() == MoleType::COFACTOR)
            {
                found=true;break;}
        }

        if (!found) {
            waterDel.push_back(&water);

            continue;
        }

        accHList.clear();
        donHList.clear();
        clashList.clear();
        grid.getAdjacentAtoms(atmList,atomO,4);
        for (ItCAtom itCA= atmList.begin();itCA != atmList.end();++itCA)
        {
            const Atom& atom = **itCA;
            if (atom.getResidu()->getResType()==ResType::WATER
                    ||atom.getResidu()->getResType()==ResType::LIGAND)continue;
            if (atom.fixpos.calcDist(atomO.fixpos) <= 3.5)
            {
                if (atom.props.isDonor())   {  donHList.push_back(&atom);}
                if (atom.props.isAcceptor()){  accHList.push_back(&atom);}
                if (!atom.props.isDonor() &&
                        !atom.props.isAcceptor())  clashList.push_back(&atom);
            }
            if (atom.getParent().getMoleType()==MoleType::LIGAND
                    ||atom.getParent().getMoleType() ==MoleType::COFACTOR)
            {

                found=true;break;
            }
        }

        const size_t accSize=accHList.size();
        const size_t donSize=donHList.size();
        double distADs[accSize], distDDs[donSize];
        multimap<double,bestOrientH> bestResults;
        for (size_t i=0; i < accSize; i++)
        {
            distADs[i]=accHList.at(i)->fixpos.calcDist(atomO.fixpos);
        }
        for (size_t i=0; i < donSize; i++)
        {
            distDDs[i]=donHList.at(i)->fixpos.calcDist(atomO.fixpos);
        }
        Coords H1;
        double dist, score1,score2,maxscore, sumscore,curr_i_score1,curr_i_score2;
        bestOrientH currBOH1,currBOH2;currBOH1.i=-10000;currBOH2.i=-10000;

        for (double x=-1;x<=1;x+=0.25)
            for (double y=-1;y<=1;y+=0.25)
                for (double z=-1;z<=1;z+=0.25)
                {
                    H1.setCoords(atomO.fixpos.x+x,atomO.fixpos.y+y,atomO.fixpos.z+z);
                    dist = H1.calcDist(atomO.fixpos);
                    if (dist>1.05 || dist <0.95 )continue;
                    Coords axe= atomO.fixpos-H1;axe.normalize();
                    Coords center(axe);center*=0.2678;center+=atomO.fixpos;
                    Coords A2(H1);A2.x+=0.1;
                    Coords vect_x=center.getNormal(H1,A2);
                    Coords vect_N=center-atomO.fixpos;vect_N.normalize();
                    Coords vect_U=vect_x-atomO.fixpos;vect_U.normalize();
                    const double scal=vect_N.x*vect_U.x+vect_N.y*vect_U.y+vect_N.z*vect_U.z;
                    Coords vectoriel(vect_N.y*vect_U.z-vect_N.z*vect_U.y,
                                     vect_N.z*vect_U.x-vect_N.x*vect_U.z,
                                     vect_N.x*vect_U.y-vect_N.y*vect_U.x);
                    curr_i_score1=0;
                    curr_i_score2=0;


                    for (double i=0; i<=2*M_PI;i+=0.1)
                    {
                        sumscore=0;
                        Coords H2 =vect_U*cos(i)+
                                vect_N*(1-cos(i))*scal+
                                vectoriel*sin(i)+
                                atomO.fixpos;
                        for (size_t i_acc=0; i_acc < accSize; i_acc++)
                        {
                            score1=calcHScore(H1,
                                              atomO.fixpos,
                                              accHList.at(i_acc)->fixpos,
                                              H2,
                                              distADs[i_acc],
                                              true);
                            score2=calcHScore(H2,
                                              atomO.fixpos,
                                              accHList.at(i_acc)->fixpos,
                                              H1,
                                              distADs[i_acc],
                                              true);
                            if (score1>score2)maxscore=score1; else maxscore=score2;

                            sumscore+=maxscore;
                        }
                        for (size_t i_don=0; i_don < donSize; i_don++)
                        {
                            const Atom &donH=*donHList.at(i_don);
                            for (size_t iBd=0; iBd < donH.getNumBond();++iBd)
                            {
                                const Coords& othAtm =donH.getAtomLinked(iBd).fixpos;
                                if (othAtm.calcDist(H1) < 2) sumscore=-1000;
                                if (othAtm.calcDist(H2) < 2) sumscore=-1000;

                            }
                            if (sumscore == -1000)break;
                        }
                        if (sumscore != -1000)
                            for (vector<const Atom*>::const_iterator
                                 itA  = clashList.begin();
                                 itA != clashList.end();
                                 itA++)
                            {
                                if (H1.calcDist((*itA)->fixpos) < 1
                                        || H2.calcDist((*itA)->fixpos) < 1) {sumscore=-1000;break;}
                            }
                        if (sumscore != -1000 && sumscore >=curr_i_score1)
                        {
                            curr_i_score2=curr_i_score1;
                            currBOH2=currBOH1;
                            curr_i_score1=sumscore;


                            currBOH1.i=i;currBOH1.x=x;currBOH1.y=y;currBOH1.z=z;
                            currBOH1.H1=H1;currBOH1.H2=H2;

                        }


                    }
                    if (currBOH1.i != -10000 && curr_i_score1>=2)
                    {
                        bestResults.insert(pair<double,bestOrientH>(curr_i_score1,currBOH1));
                    }
                    if (currBOH2.i != -10000 && curr_i_score2 >= 2)
                    {
                        bestResults.insert(pair<double,bestOrientH>(curr_i_score2,currBOH2));
                    }


                } // END X Y Z

        if (bestResults.empty()){waterDel.push_back(&water);continue;}
        for (multimap<double,bestOrientH>::reverse_iterator
             itr= bestResults.rbegin();
             itr != bestResults.rend();
             itr++)
        {
            score1=(*itr).first;
            bestOrientH &boh = (*itr).second;
            Molecule &parent=*water.getParent();
            parent.addBond(parent.addAtom("H",boh.H1,"H","H",&water),atomO,BondType::SINGLE);
            parent.addBond(parent.addAtom("H",boh.H2,"H","H",&water),atomO,BondType::SINGLE);
            // NMod++;
            if (wVerbose) cout << verbose<<"SELECT\t"<< water.getIdentifier()<<"\t"<<score1<<endl;
            NMod++;
            break;
        }


        /* water.getParent()->renumAtom();
        MoleWriter mw(water.getIdentifier()+".mol2");
        mw.writeMOL2(water.getParent());*/






        //  break;
    }







    for (ItCRes itD = waterDel.begin(); itD != waterDel.end(); itD++)
    {
        Residu &residu=**itD;
        //  cout << "WATER DELETION : "<<residu.getIdentifier()<<endl;
        if (residu.getParent()->numResidus()==1)
        {
            deleteMole(residu.getParent());
            continue;
        }
        else
        {
            residu.getParent()->delResidu(residu);
        }
    }


    if (wVerbose) cout << verbose<<"STATIS\t"<<NMod<<"\tNumber of selected waters"<<endl
                       << verbose<<"END\t"<<endl;
}





template <class IT,class BA> double Complex::calcMoleSmallDist(const IT& ref, const BA& comp, const double& best_dist, const bool& wHydrogen) const
{
    double bestDist=1000;
    const size_t nRefAtm=ref.numAtom();
    const size_t nCompAtm=comp.numAtom();

    for (size_t iRef = 0; iRef < nRefAtm; ++iRef)
    {
        const Atom& atomRef = ref.getAtom(iRef);
        if (atomRef.isHydrogen() && !wHydrogen)continue;
        for (size_t iComp = 0; iComp < nCompAtm; ++iComp)
        {
            const Atom& atomComp = comp.getAtom(iComp);
            if (atomComp.isHydrogen() && !wHydrogen)continue;
            const double dist=atomRef.fixpos.calcDist(atomComp.fixpos,bestDist+0.1);
            if (dist < bestDist)
            {
                bestDist = dist;
                if (bestDist <= best_dist) return bestDist;
            }

        }
    }
    return bestDist;
}


template<class BA> bool Complex::calcSmallDist(const Molecule& ref, const BA& comp, const double& best_dist, const bool& wHydrogen)
{
    double bestDist=1000;
    const size_t nRefAtm=ref.numAtom();
    const size_t nCompAtm=comp.numAtom();
    AtomList listadjacent,listfulladjacent;
    bool useGrid=false;
    Grid *grid=(Grid*)NULL;

    grid=&getGrid(0);
    useGrid=true;




    for (size_t iRef = 0; iRef < nRefAtm; ++iRef)
    {
        const Atom& atomRef = ref.getAtom(iRef);
        if (atomRef.isHydrogen() && !wHydrogen)continue;
        if (useGrid && atomRef.getBox(grid) !=(Box*)NULL)
        {
            grid->getAdjacentAtoms(listadjacent,atomRef,best_dist,wHydrogen);
            listfulladjacent.insert(listfulladjacent.end(),listadjacent.begin(),listadjacent.end());
        }
        else
        {
            for (size_t iComp = 0; iComp < nCompAtm; ++iComp)
            {
                const Atom& atomComp = comp.getAtom(iComp);
                if (atomComp.isHydrogen() && !wHydrogen)continue;
                const double dist=atomRef.fixpos.calcDist(atomComp.fixpos,bestDist+0.1);
                if (dist < bestDist)
                {
                    bestDist = dist;
                    if (bestDist <= best_dist) return true;
                }

            }
        }
    }
    sortAndUnique(listfulladjacent);
    for (size_t iComp = 0; iComp < nCompAtm; ++iComp)
    {
        const Atom& atomComp = comp.getAtom(iComp);
        if (atomComp.isHydrogen() && !wHydrogen)continue;
        if (find(listfulladjacent.begin(),listfulladjacent.end(),&atomComp)==listfulladjacent.end())continue;
        return true;
    }

    return false;
}







void Complex::selectLigand(const std::string& PDB_DIR, const ResiduList& failedRes,const AtomList& failedAtm, const std::string &verbose, const bool &wUnDrugg)
{
    const bool wVerbose= !verbose.empty();
    if (wVerbose) cout << verbose<<"START"<<endl;


    // Listing all cofactor and ligand:
    MoleList List_Poten_Lig;
    Molecule *tempMole=(Molecule*)NULL;
    vector<std::string> chainnames;
    const size_t nMole = allmolecules.size();

    // We first list ligand molecule:
    for (size_t iMole=0; iMole < nMole;++iMole)
    {
        Molecule &mole = *allmolecules.at(iMole);
        for (size_t iCh=0; iCh < mole.numChains();++iCh) chainnames.push_back(mole.getChain(iCh).getName());
        if (mole.getMoleType()==MoleType::PROTEIN) tempMole=&mole;
        if (mole.getMoleType() != MoleType::LIGAND) continue;
        List_Poten_Lig.push_back(&mole);


    }
    // AND THEN COFACTORS !!!
    for (size_t iMole=0; iMole < nMole;++iMole)
    {
        Molecule &mole = *allmolecules.at(iMole);
        if (mole.getMoleType() != MoleType::COFACTOR) continue;
        List_Poten_Lig.push_back(&mole);
    }
    if (tempMole==(Molecule*)NULL) return;

    Molecule& protein=*tempMole;


    // Since the structure is exploded into multiple molecules,
    // to be sure that we don't miss any chain, we have listed them above
    // and now we make the unique
    sortAndUnique(chainnames);
    const size_t nChain = chainnames.size();
    vector<bool> chainsallowed(nChain,false);




    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    //////////////////////////// INITIALIZE VARIABLES ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    bool            ResAllowed=false,MoleAllowed=false,
            isCof=false;
    unsigned int    NSITE=0;       // Iterator for ligand ID.

    ostringstream   ligOutput;
    map<string,double> PercPres;

    size_t NCavAll=0;
    size_t bestICav=100000;
    double bestCavRecov=0,
            bestLigRecov=0;

    unsigned int
            STDAA=0,                            // Number of Standard amino acid within binding site
            MODAA=0,                            // Number of Modified amino acid within binding site
            WAT=0,                            // Number of water molecules within binding site
            ION=0,                            // Number of ions within binding site
            COF=0,                            // Number of cofactors within binding site
            UNK=0;                            // Number of unknown within binding site
    ostringstream   ossCOF,
            ossION;
    VolSiteResult* volligand= (VolSiteResult*)NULL;
    VolSiteResult*  cavALL= (VolSiteResult*)NULL;
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////// SCANNING EACH LIGAND ////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    for (ItMole itL = List_Poten_Lig.begin(); itL != List_Poten_Lig.end(); itL++)
    {
        try
        {

            // STEP 1 - Defining the ligand
            //1.1 - Getting ligand reference:
            Molecule &ligand =**itL;
            // cout << "LOOK AT : "<<ligand.getResidu(0).getIdentifier()<<endl;
            const bool isPeptide=(ligand.numResidus() >1)?true:false;
            //1.2 - Is the ligand a cofactor
            if (ligand.getMoleType() == MoleType::COFACTOR) isCof=true; else isCof=false;

            //1.3 - Ligand ID

            ligOutput.str("");
            ligand.setUse(true);
            volligand= (VolSiteResult*)NULL;
            cavALL= (VolSiteResult*)NULL;
            NCavAll=0;
            bestICav=100000;
            bestCavRecov=0;
            bestLigRecov=0;


            //1.4 - Getting ligand geometrical center
            const Coords LigPos = ligand.getFixpos();

            //1.5 IF COFACTOR : Checking that no ligand is within the same site - Otherwise skipping it
            if (isCof)
            {
                //cout << "  ||--> COFACTOR "<<endl;
                bool check=true;
                for (ItMole itL2 = List_Poten_Lig.begin(); itL2 != List_Poten_Lig.end(); itL2++)
                {
                    const Molecule &molecule2=**itL2;
                    if (&ligand== &molecule2)continue;
                    // if (molecule2.getMoleType() == MoleType::LIGAND && calcMoleSmallDist(ligand,molecule2, 6.5)<= 6.5 ) {check=false;break;}
                    if (molecule2.getMoleType() == MoleType::LIGAND && calcSmallDist(ligand,molecule2, 6.5)) {check=false;break;}
                }
                if (!check)
                {
                    //  cout << "   ||--> Ligand within the same site => stop process "<<endl;
                    NSITE--;
                    continue;
                }
            }






            // STEP 2 - CREATING THE PROTEIN :

            PercPres.clear();


            for (size_t iCh=0;iCh<nChain;++iCh)
                chainsallowed.at(iCh)=false;

            // All molecules can potentially be in protein.mol2:
            for (size_t iMole=0; iMole < nMole;++iMole)
            {
                Molecule &mole = *allmolecules.at(iMole);
                if (mole.getMoleType()!=MoleType::PROTEIN
                        &&mole.getMoleType()!=MoleType::COFACTOR)continue;
                // So we consider them as usable:
                mole.setUse(true);
                MoleAllowed=false;
                for (size_t iChain=0; iChain < mole.numChains(); ++iChain)
                {
                    ResAllowed=false;
                    Chain& chain= mole.getChain(iChain);

                    const size_t nRes = chain.getCountRes();
                    for (size_t iRes = 0; iRes < nRes;++iRes)
                    {
                        Residu &residu = chain.getResidu(iRes);
                        if (residu.getAtom(0).fixpos.calcDist(ligand.getAtom(0).fixpos) > 50)
                        {
                            continue;
                        }
                        const double dist =calcMoleSmallDist(ligand,residu,6.5);
                        if (dist < 6.5){ ResAllowed=true;    break;}
                    }
                    if (!ResAllowed) { continue;}
                    chainsallowed.at(getPos(chainnames,chain.getName()))=true;
                }


            }
            if (isPeptide)
            {
                for (size_t iCh=0; iCh < ligand.numChains();++iCh)
                {
                    const std::string &n=ligand.getChain(iCh).getName();

                    chainsallowed.at(getPos(chainnames,n))=false;
                }
            }
            for (size_t iCh=0;iCh<nChain;++iCh)
            {
                if (chainsallowed.at(iCh)) PercPres.insert(pair<string,double>(chainnames.at(iCh),0));
            }


            /*                 for (size_t iCh=0;iCh<nChain;++iCh)
                 {
                     cout << chainnames.at(iCh)<<"\t"<< chainsallowed.at(iCh)<<endl;
                 }*/
            MoleList keptedwater;
            for (size_t iMole=0; iMole < nMole;++iMole)
            {
                Molecule &mole = *allmolecules.at(iMole);
                mole.clearSets();
                for (size_t iChain=0; iChain < mole.numChains(); ++iChain)
                {
                    Chain& chain= mole.getChain(iChain);
                    if (!chainsallowed.at(getPos(chainnames,chain.getName())))
                        chain.setUse(false);
                }
                if (mole.getMoleType()==MoleType::WATER)
                {
                    //cout << mole.getResidu(0).getIdentifier()<< " " << calcMoleSmallDist(ligand,mole,6.5)<< " " <<mole.getChain(0).getName()<< " " << chainsallowed.at(getPos(chainnames,mole.getChain(0).getName()))<< " ";
                    const double dist=calcMoleSmallDist(ligand,mole,6.5);
                    mole.setUse(false);
                    /* if (dist >=6.5)
                    {mole.setUse(false);}
                    else*/
                    if (!mole.getChain(0).isUsed() && dist <6.5)
                    {
                        AtomSet HYDESet;
                        HYDESet.name="HYDE_BINDING_SITE";
                        HYDESet.subtype="<user>";
                        HYDESet.type="STATIC";
                        HYDESet.list.clear();
                        for (ItCAtom itA = mole.firstAtom();itA != mole.lastAtom();++itA)
                        {
                            Atom& atm=**itA;
                            if (atm.isOxygen()) HYDESet.list.push_back(&atm);
                        }
                        mole.addAtomSet(HYDESet);
                        keptedwater.push_back(&mole);
                    }
                    //mole.setUse(true);

                }
                else if (mole.getMoleType()==MoleType::PROTEIN)
                {
                    const size_t nRes = mole.numResidus();
                    for (size_t iRes=0;iRes < nRes;++iRes)
                    {
                        Residu& residu = mole.getResidu(iRes);
                        const double dist=calcMoleSmallDist(ligand,residu,6.5);
                        if (residu.getAtom(0).fixpos.calcDist(ligand.getAtom(0).fixpos) > 50)
                        {
                            residu.setUse(false);
                        }
                        else if (residu.getResType()==ResType::METAL &&
                                 dist >=6.5)residu.setUse(false);
                    }
                }
                else if (mole.getMoleType() == MoleType::COFACTOR)
                {
                    const double dist=calcMoleSmallDist(ligand,mole,6.5);
                    if (dist >=6.5)
                    {mole.setUse(false);}
                    else if (!mole.getChain(0).isUsed() && dist <6.5) mole.setUse(true);

                }
                else if (mole.getMoleType() == MoleType::LIGAND)
                {
                    mole.setUse(false);
                }
            }

            //


            for (size_t iMole=0; iMole < nMole;++iMole)
            {
                Molecule &mole = *allmolecules.at(iMole);

                AtomSet HYDESet;
                HYDESet.name="HYDE_BINDING_SITE";
                HYDESet.subtype="<user>";
                HYDESet.type="STATIC";
                HYDESet.list.clear();
                const size_t nRes = mole.numResidus();
                for (size_t iRes=0;iRes < nRes;++iRes)
                {
                    Residu& residu = mole.getResidu(iRes);
                    if (!residu.isUsed())continue;
                    const bool isSite=(calcMoleSmallDist(ligand,residu,6.5)>=6.5)?false:true;
                    for (size_t iAtm=0; iAtm < residu.numAtom();++iAtm)
                    {
                        Atom& atm=residu.getAtom(iAtm);
                        if (atm.isCarbon())
                        {
                            if (&mole==&ligand)
                            {
                                atm.setColor(255,127,0);
                                mole.addAtomColor(&atm);
                            }
                            else if (isSite)
                            {

                                atm.setColor(0,255,0);
                            }
                            else atm.setColor(255,0,255);
                            mole.addAtomColor(&atm);
                        }
                        if (isSite && !atm.isHydrogen() && &mole != &ligand)HYDESet.list.push_back(&atm);
                    }

                }
                if (HYDESet.list.empty())continue;
                mole.addAtomSet(HYDESet);
            }

            //protein.addAtomSet(HYDESet);
            unsigned short nFailedResLig=0,
                    nFailedAtmLig=0,
                    nFailedResSit=0,
                    nFailedAtmSit=0;


            NSITE++;
            ostringstream osLigName;
            osLigName.str("");
            if (!verbose.empty()) osLigName<< verbose.substr(0,4)<<"_";
            osLigName<<NSITE;
            double MWeight=0, NHeavy=0;
            for (size_t i=0;i< ligand.numResidus();++i)
            {
                Residu& residu=ligand.getResidu(i);
                if (find(failedRes.begin(),failedRes.end(),&residu) != failedRes.end())
                {
                    nFailedResLig++;
                }
                const std::string& name=residu.getName();
                osLigName<< "_"<<name;

                MWeight+=residu.getWeight();
                for (size_t iAtm=0;iAtm < residu.numAtom();++iAtm)
                {
                    Atom &atom=residu.getAtom(iAtm);
                    if (atom.isHydrogen())continue;
                    if (find(failedAtm.begin(),failedAtm.end(),&atom) != failedAtm.end() || atom.getMOL2Type()=="Du")
                    {
                        nFailedAtmLig++;
                    }
                    NHeavy++;
                }
            }

            ligand.setName(osLigName.str());
            ///// STEP - GENERATING CAVITIES :
            ligand.setUse(true);
            //Molecule::Rules[MoleType::WATER]=0;

            try{
                if (isCof) ligand.setMoleType(MoleType::LIGAND);
                VolSite volsite(*this,ligand);
#ifdef ICHEM_DEBUG
                volsite.grid.printInFile(osLigName.str()+"_INI",volsite.grid.AllBoxes,false,true);
#endif
                const double proj=55;
                volsite.proj(proj-20);
#ifdef ICHEM_DEBUG
                volsite.grid.printInFile(osLigName.str()+"_PROJ",volsite.grid.AllBoxes,false,true);
#endif
                volsite.treatments();
#ifdef ICHEM_DEBUG
                volsite.grid.printInFile(osLigName.str()+"_TREAT",volsite.grid.AllBoxes,false,true);
#endif
                volsite.clusterCAV(proj);
#ifdef ICHEM_DEBUG
                double dist;

                volsite.output_cavs_lig(ligand,osLigName.str()+"CAVSELS",true,dist,false,1);
#endif

                if (volsite.getNumCav()==0) continue;
                volligand = new VolSiteResult[5];

                volsite.structCavLig(volligand,ligand,ligand.getName(),true);
                ligand.setUse(false);
                bool nocav=false;
                for (size_t i=0;i<5;++i)
                {
                    if (volligand[i].nbcube<10)nocav=true;
                }
                if (volligand[4].nbcube>30 && nocav)nocav=false;
                if (nocav)
                {
                    delete[]volligand;
                    volligand=(VolSiteResult*)NULL;
                    cout << "no cav found"<<endl;
                    if (isCof) ligand.setMoleType(MoleType::COFACTOR);
                    continue;
                }
                VolSite volFull(*this,1.5);
#ifdef ICHEM_DEBUG
                volFull.grid.printInFile(osLigName.str()+"_FINI",volsite.grid.AllBoxes,false,true);
#endif
                volFull.proj(proj-20);
#ifdef ICHEM_DEBUG
                volFull.grid.printInFile(osLigName.str()+"_FPROJ",volsite.grid.AllBoxes,false,true);
#endif
                volFull.treatments();
#ifdef ICHEM_DEBUG
                volFull.grid.printInFile(osLigName.str()+"_FTREAT",volsite.grid.AllBoxes,false,true);
#endif
                volFull.clusterCAV(proj);
                NCavAll=volFull.getNumCav();
                cavALL=new VolSiteResult[NCavAll];
                volFull.structCavAll(cavALL,ligand.getName()+"_CAVALL.mol2",true);


            }catch (MoleExcept &e)
            {
                cerr << e.getCode()<<endl<<e.getData()<<e.getSource()<<e.getTrace()<<endl;
                if (volligand != (VolSiteResult*)NULL)
                {
                    delete[]volligand;
                    volligand=(VolSiteResult*)NULL;
                }
                if (cavALL != (VolSiteResult*)NULL)
                {
                    delete[]cavALL;
                    cavALL = (VolSiteResult*)NULL;
                }
                if (isCof) ligand.setMoleType(MoleType::COFACTOR);
                continue;
            }

            ////// STEP - Seeking if cavity is druggable or not :
            ///

            // Seeking FULL cavity closest to ligand :
            double currLigRecov=0;


            bool found;

            for (size_t iCav=0; iCav < NCavAll;++iCav)
            {
                const Molecule& cavity= cavALL[iCav].moleRes;
                currLigRecov=0;
                for (size_t iLig=0;iLig < ligand.numAtom();++iLig)
                {
                    Atom& ligAtm=const_cast<Atom&>(ligand.getAtom(iLig));
                    if (ligAtm.isHydrogen())continue;
                    for (size_t iCav=0;iCav < cavity.numAtom();++iCav)
                    {
                        const Atom& cavAtm=cavity.getAtom(iCav);
                        found=false;
                        if (ligAtm.fixpos.calcDist(cavAtm.fixpos) > ligAtm.getVdWRadius())continue;
                        found=true;break;

                    }
                    if (!found)continue;
                    currLigRecov++;
                }

                if (currLigRecov < bestLigRecov) continue;
                bestLigRecov=currLigRecov;
                bestICav=iCav;


                for (size_t iCav=0;iCav < cavity.numAtom();++iCav)
                {
                    const Atom& cavAtm=cavity.getAtom(iCav);
                    found=false;
                    for (size_t iLig=0;iLig < ligand.numAtom();++iLig)
                    {
                        Atom& ligAtm=const_cast<Atom&>(ligand.getAtom(iLig));
                        if (ligAtm.isHydrogen())continue;
                        if (ligAtm.fixpos.calcDist(cavAtm.fixpos) > ligAtm.getVdWRadius())continue;
                        found=true;

                    }
                    if (!found)continue;
                    bestCavRecov++;
                }


            }



            if (bestICav==100000)
            {
                NSITE--;
                if (volligand != (VolSiteResult*)NULL)
                {
                    delete[]volligand;
                    volligand=(VolSiteResult*)NULL;
                }
                if (cavALL != (VolSiteResult*)NULL)
                {
                    delete[]cavALL;
                    cavALL = (VolSiteResult*)NULL;
                }
                if (isCof) ligand.setMoleType(MoleType::COFACTOR);
                continue;
            }
            VolSiteResult& volResALL=cavALL[bestICav];




            if (!wUnDrugg && volResALL.drugg < 0)
            {
                NSITE--;
                if (volligand != (VolSiteResult*)NULL)
                {
                    delete[]volligand;
                    volligand=(VolSiteResult*)NULL;
                }
                if (cavALL != (VolSiteResult*)NULL)
                {
                    delete[]cavALL;
                    cavALL = (VolSiteResult*)NULL;
                }
                if (isCof) ligand.setMoleType(MoleType::COFACTOR);
                continue;
            }


            for (ItMole it = keptedwater.begin();it != keptedwater.end();++it)
                (*it)->setUse(true);

            MoleWriter mw;
            mw.setOnlyUsed(true);
            ligand.setUse(true);

            mw.newFile(PDB_DIR+"/"+ligand.getName()+".mol2");
            ligand.renumAtom();
            ligand.renumBond();
            ligand.renumResidu();
            mw.writeMOL2(&ligand);
            ligand.clearSets();
            ligand.setUse(false);
            protein.renumAtom();
            protein.renumBond();
            protein.renumResidu();

            mw.newFile(PDB_DIR+"/"+ligand.getName()+"_PROT.mol2");
            mw.writeMOL2(*this,ligand.getName()+"_PROT");


            STDAA=0;MODAA=0;WAT=0;ION=0;COF=0;UNK=0;ossCOF.str("");ossION.str("");
            double BFactor=0;
            double NA=0;
            size_t NSiteSubStruct=0;
            for (size_t iMole=0; iMole < nMole;++iMole)
            {
                Molecule &mole = *allmolecules.at(iMole);
                const size_t nRes = mole.numResidus();
                for (size_t iRes=0;iRes < nRes;++iRes)
                {
                    Residu& residu = mole.getResidu(iRes);
                    if (!residu.isUsed())continue;
                    if (calcMoleSmallDist(ligand,residu,6.5)>=6.5)
                    {
                        residu.setUse(false);
                        continue;
                    }
                    if (find(failedRes.begin(),failedRes.end(),&residu) != failedRes.end())
                    {
                        nFailedResSit++;
                    }
                    if (residu.getResType()!=ResType::WATER)
                    {
                        PercPres.at(residu.getChainName())++; // And we calculate presence percentages
                    }
                    NSiteSubStruct++;
                    switch (residu.getResType())
                    {
                    case ResType::STD_AA:   STDAA++;break;
                    case ResType::MOD_AA:   MODAA++;break;
                    case ResType::COFACTOR: COF++;ossCOF<<residu.getName()<<" ";break;
                    case ResType::METAL:    ION++;ossION<<residu.getName()<<" ";break;
                    case ResType::WATER:    WAT++;break;
                    default: UNK++;
                    }
                    for (size_t iAtm=0;iAtm < residu.numAtom();++iAtm)
                    {
                        const Atom&atom = residu.getAtom(iAtm);
                        if (find(failedAtm.begin(),failedAtm.end(),&atom) != failedAtm.end() || atom.getMOL2Type()=="Du")
                        {
                            nFailedAtmSit++;
                        }
                        if (atom.isHydrogen())continue;

                        NA++;
                        BFactor +=atom.getBFactor();
                    }


                }

            }

            if (NA > 0)
            BFactor/=NA;
            else BFactor=0;


            protein.renumAtom();
            protein.renumBond();
            protein.renumResidu();
            mw.newFile(PDB_DIR+"/"+ligand.getName()+"_SITE.mol2");
            mw.writeMOL2(*this,ligand.getName()+"_SITE");


            Molecule& cavity= cavALL[bestICav].moleRes;
            cavity.setName(ligand.getName()+"_FULL");
            mw.newFile(PDB_DIR+"/"+ligand.getName()+"_FULL.mol2");
            mw.writeMOL2(&cavity);


            Grid gridBSA;
            gridBSA.createComplexBoxes(*this,ligand,0.5,Coords(20,20,20));
            double BSA,volLig;
            gridBSA.calcBSA(*this,ligand,BSA,volLig);

            cout << verbose<<"SITE\t\n"
                 << verbose<<"HEADER\tBSA\tBFactor\tMW\tNATM\tHETNAME\n"
                 << verbose<<"VALUES\t"<<BSA<<"\t"<<BFactor<<"\t"<<MWeight<<"\t"<<NHeavy<<"\t"<<ligand.getName()<<"\n"
                 << verbose<<"HEADER\tNRES\tSTDAA\tMODAA\tNWAT\tION\tCOF\tCOFList\tIONList\tCHAINS\n"
                 << verbose<<"VALUES\t"<<NSiteSubStruct<<"\t"<<STDAA<<"\t"<<MODAA<<"\t"<< WAT <<"\t"<< ION<<"\t"<<COF<<"\t"<<ossCOF.str() <<"\t"<< ossION.str() <<"\t";
            for (size_t iCh=0;iCh<nChain;++iCh) if (chainsallowed.at(iCh)) cout << chainnames.at(iCh)<<" ";
            cout <<"\n";
            cout << verbose<<"ERRORS\tRESLIG\tRESSIT\tATMLIG\tATMSIT\n"
                 << verbose<<"VALUES\t"<< nFailedResLig<<"\t"<<nFailedResSit<<"\t"<<nFailedAtmLig<<"\t"<<nFailedAtmSit<<"\n"
                 << verbose<<"HEADER\tX\tY\tZ\n"
                 << verbose<<"VALUES\t"<<LigPos.x <<"\t"<< LigPos.y     <<"\t" << LigPos.z     <<"\n"
                 << verbose<<"HEADER\tPERC_PRES\t"<<PercPres.size()<<"\n";

            for (map<string,double>::iterator it = PercPres.begin() ; it != PercPres.end(); it++)
            {
                cout <<verbose<<"VALUES\t"<< (*it).first<< "\t"<<
                       (double)((*it).second/(double)(NSiteSubStruct-WAT))<<"\n";
            }

            cout << verbose << "CAVITY\tSize\tNPts\tVolume\tCA\tCZ\tO\tOG\tOD1\tN\tNZ\tDU\tDrugg\tLigRecov\tCavRecov"<<endl;
            cout << verbose<<"VALUES\t"
                 <<"FULL\t"
                <<volResALL.nbcube<<"\t"
               <<volResALL.volume<<"\t"
              <<volResALL.CA<<"\t"
             <<volResALL.CZ<<"\t"
            <<volResALL.O<<"\t"
            <<volResALL.OG<<"\t"
            <<volResALL.OD1<<"\t"
            <<volResALL.N<<"\t"
            <<volResALL.NZ<<"\t"
            <<volResALL.DU<<"\t"
            <<volResALL.drugg<<"\t"
            <<(bestLigRecov/NHeavy*100)<<"\t"
            <<((bestCavRecov/(double)volResALL.nbcube)*100)<<endl;
            for (size_t i=0;i<5;++i)
            {
                VolSiteResult& volRes =volligand[i];
                cout << verbose<<"VALUES\t";
                mw.newFile(PDB_DIR+"/"+volRes.moleRes.getName()+".mol2");
                mw.writeMOL2(&volRes.moleRes);
                switch (volRes.moleRes.getMoleType())
                {

                case MoleType::CAV_4:cout << "4";break;
                case MoleType::CAV_6:cout << "6";break;
                case MoleType::CAV_8:cout << "8";break;
                case MoleType::CAV_12:cout << "12";break;
                case MoleType::CAV_ALL:cout << "ALL";break;

                }
                cout <<"\t"<<volRes.nbcube<<"\t"
                    << volRes.volume<<"\t"
                    <<volRes.CA<<"\t"
                   <<volRes.CZ<<"\t"
                  <<volRes.O<<"\t"
                 <<volRes.OG<<"\t"
                <<volRes.OD1<<"\t"
                <<volRes.N<<"\t"
                <<volRes.NZ<<"\t"
                <<volRes.DU<<"\t"
                <<volRes.drugg<<"\t"
                <<((volRes.ligRecov/NHeavy)*100)<<"\t"
                <<((volRes.cavRecov/(double)volRes.nbcube)*100)<<endl;
            }
            cout << verbose<<"ENDSITE\t\n";


            if (volligand != (VolSiteResult*)NULL)
            {
                delete[]volligand;
                volligand=(VolSiteResult*)NULL;
            }
            if (cavALL != (VolSiteResult*)NULL)
            {
                delete[]cavALL;
                cavALL = (VolSiteResult*)NULL;
            }
            if (isCof) ligand.setMoleType(MoleType::COFACTOR);


        }
        catch (MoleExcept &e)
        {
            cout << e.getData()<<endl;

        }

    }
    if (wVerbose) cout << verbose<<"END"<<endl;


}

void Complex::moveMoleToComplex(Molecule* const molecule,Complex &complex)

{

    if (molecule == (Molecule*)NULL)

        throw MoleExcept(1100701,

                         "Complex::updateMoleType"

                         "No molecule given");

    ItMole it = find(allmolecules.begin(),

                     allmolecules.end(),

                     molecule);

    if (it == allmolecules.end())

        throw MoleExcept(1100702,

                         "Complex::updateMoleType"

                         "Given molecule is not part of this complex");

    allmolecules.erase(it);

    it = find(molecules[molecule->getMoleType()].begin(),

            molecules[molecule->getMoleType()].end(),

            molecule);

    if (it == molecules[molecule->getMoleType()].end())

    {

        for (size_t i=0;i < NB_MOLETYPE;++i)

        {

            it = find(molecules[i].begin(),

                      molecules[i].end(),

                      molecule);

            if (it != molecules[i].end())molecules[i].erase(it);

        }

    }

    else molecules[molecule->getMoleType()].erase(it);

    complex.addMolecule(molecule);

}
