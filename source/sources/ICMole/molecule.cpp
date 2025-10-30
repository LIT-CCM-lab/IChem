#include <iostream>
#include <stdexcept>
#include "headers/ICMole/bond.h"
#include "headers/ICMole/molecule.h"
#include "headers/ICMole/cycle.h"
#include "headers/ICMole/complex.h"

//#define ICHEM_DEBUG
using namespace std;
using namespace ICMole;


 unsigned int Molecule::Rules[NB_MOLETYPE];
 bool Molecule::rules_loaded=false;


void Molecule::loadRules()
{
    if (rules_loaded)return;
    rules_loaded=true;
    for (unsigned int I=0;I < NB_MOLETYPE;I++)
    {Rules[I]=0;

    }
    Rules[MoleType::PROTEIN] =MoleType::PROTEIN;
    Rules[MoleType::SITE]    =MoleType::PROTEIN;
    Rules[MoleType::LIGAND]  =MoleType::LIGAND;
    Rules[MoleType::WATER]   =MoleType::PROTEIN;
    Rules[MoleType::COFACTOR]=MoleType::LIGAND;
    Rules[MoleType::NUCLEIC] =MoleType::PROTEIN;
}





/*! \fn Molecule::Molecule(const unsigned int NAtom, const unsigned int NBonds)
  * \brief Constructor
  * \param NAtom : Number of Atom of this Molecule. Can be 0. It is used to allocate memory and makes the code faster. Default:50
  * \param NBonds : Number of Bonds of this Molecule. Can be 0. It is used to allocate memory and makes the code faster. Default:100
  */

Molecule::Molecule(const unsigned int &NAtom,
                   const unsigned int &NBonds,
                   const unsigned int MoleType):
    Atoms(0),
    maxNumAtom(0),
    maxNumBond(0),
    maxNumResidu(0),
    typemol(MoleType),
    bary_check(false),
    name(""),
    dummyChain(this,"XX"),
    unknownRes(this, 1000, 1000,"UNK",dummyChain,ResType::UNKNOWN),
    aromaticRes(this, 1001, 1001,"UNK",dummyChain,ResType::UNKNOWN),
    cyclicRes(this, 1002, 1002,"UNK",dummyChain,ResType::UNKNOWN),
    complex((Complex*)NULL),
    fixpos(0,0,0),
    rotpos(0,0,0),
    extremH(0,0,0),
    extremL(0,0,0)

{
    Atoms.reserve(NAtom);
    Bonds.reserve(NBonds);
}


Molecule::~Molecule()
{
    for(size_t i = 0; i < Cycles.size(); ++i) delete Cycles[i];

    for(size_t i = 0; i < Residues.size(); ++i)
    {
        delete Residues[i];
    }

    for(size_t i = 0; i < Chains.size(); ++i)
    {
        delete Chains[i];
    }

    for(size_t i = 0; i < Bonds.size(); ++i) delete Bonds[i];

    for(size_t i = 0; i < Atoms.size(); ++i) delete Atoms[i];

}






////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// COORDS ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

const Coords& Molecule::getFixpos(const bool& withH)
{
    if (bary_check) return fixpos;
    fixpos.clear();
    rotpos.clear();
    bary_check=true;
    double n=0;
    for (ItCAtom itA = Atoms.begin();
                 itA!= Atoms.end();
               ++itA)
    {
        if (!(*itA)->isUsed()||  (!withH && (*itA)->isHydrogen()))continue;
        n++;
        fixpos+=(*itA)->fixpos;
        rotpos+=(*itA)->rotpos;
    }
if (n ==0) throw MoleExcept(1071401, "Molecule::fixpos", "Division by 0");
    fixpos/=n;
    rotpos/=n;
    return fixpos;
}

const Coords& Molecule::getRotpos(const bool& withH)
{
    if (bary_check) return rotpos;
    bary_check=true;
    fixpos.clear();
    rotpos.clear();
    double n=0;
    for (ItCAtom itA = Atoms.begin();
                 itA!= Atoms.end();
               ++itA)
    {
        if (!(*itA)->isUsed()|| (!withH && (*itA)->isHydrogen()))continue;
        n++;
        fixpos+=(*itA)->fixpos;
        rotpos+=(*itA)->rotpos;
    }
    fixpos/=n;
    rotpos/=n;
    return rotpos;
}


Coords Molecule::getalpha(const int residu, const bool& withH){
    int n = 0;
    for (ItCAtom itA = Atoms.begin(); itA!= Atoms.end(); ++itA) {
        const Atom& atm = **itA;
        if (!(*itA)->isUsed()|| (!withH && (*itA)->isHydrogen())) continue;
        if (atm.getName() == "CA") {
          if (n == residu) {
            return atm.fixpos;
          }
          n++;
        }
    }
    //    return Coords(0,0,0);
}

void Molecule::getextrem(const bool& withH){
    int n =0;
    for (ItCAtom itA = Atoms.begin();
                 itA!= Atoms.end();
               ++itA)
    {
        if (!(*itA)->isUsed()|| (!withH && (*itA)->isHydrogen()))continue;
        n++;
        if ((*itA)->fixpos.x > extremH.x) {
            extremH.x = (*itA)->fixpos.x;
        }else if ((*itA)->fixpos.x < extremH.x){
            extremL.x = (*itA)->fixpos.x;
        }if ((*itA)->fixpos.y > extremH.y) {
            extremH.y = (*itA)->fixpos.y;
        }else if ((*itA)->fixpos.y < extremH.y){
            extremL.y = (*itA)->fixpos.y;
        }        if ((*itA)->fixpos.z > extremH.z) {
            extremH.z = (*itA)->fixpos.z;
        }else if ((*itA)->fixpos.z < extremH.z){
            extremL.z = (*itA)->fixpos.z;
        }

    }
}

/*! \fn   Coords Molecule::calcInertialMoment(bool With_H,bool rot)
      * \param rot (default: false). If true, apply the calculation on the mobile coordinates
      * \return Coordinates of the inertial moment
      */
Coords Molecule::calcInertialMoment(const bool& withH,
                                    const bool& onRotpos)
{
  Coords moment;
  bary_check=true;
  fixpos.clear();
  rotpos.clear();
  double n=0;
  for (ItCAtom itA = Atoms.begin();
               itA!= Atoms.end();
             ++itA)
  {
      const Atom& atm = **itA;
      if (!(*itA)->isUsed() || (!withH && (*itA)->isHydrogen()))continue;
      n++;
      fixpos+=(*itA)->fixpos;
      rotpos+=(*itA)->rotpos;
  }
  if (n==0) return Coords(0,0,0);
  fixpos/=n;
  rotpos/=n;

  double na=1;
  for (ItCAtom itA = Atoms.begin(); itA != Atoms.end(); itA++)
    {
      Atom &atm = **itA;
      if ((atm.isHydrogen() && !withH) || !atm.isUsed()) continue;
      if (onRotpos) moment += (atm.rotpos+rotpos*(-1.0))*(atm.rotpos+rotpos*(-1.0))*atm.getWeight();
      else          moment += (atm.fixpos+fixpos*(-1.0))*(atm.fixpos+fixpos*(-1.0))*atm.getWeight();
      na++;

    }
//  moment.setCoords(moment.x,moment.y,moment.z);
  if (na <=1) {return Coords(0,0,0);}
  else {moment/=(na-1);}
//  moment.setCoords(sqrt(moment.x),sqrt(moment.y),sqrt(moment.z));
  if (onRotpos) moment+=rotpos;
//  else moment+=fixpos;

  return moment;
}


void Molecule::rotateMolecule(const std::vector<double>& matrix,
                              const Coords& transMobil,
                              const Coords& transRigid
                              )
{
    if (matrix.size() < 9) return;
    Coords Center(getFixpos());
    Coords Translate(Center-transMobil);
    Translate.setCoords(Translate.x*matrix[0]+Translate.y*matrix[1]+Translate.z*matrix[2],
                        Translate.x*matrix[3]+Translate.y*matrix[4]+Translate.z*matrix[5],
                        Translate.x*matrix[6]+Translate.y*matrix[7]+Translate.z*matrix[8]);



    for (ItCAtom itA = Atoms.begin(); itA != Atoms.end(); ++itA)
        {
        Atom *atm = *itA;
            atm->rotpos.setCoords( atm->fixpos.x*matrix[0]+atm->fixpos.y*matrix[1]+atm->fixpos.z*matrix[2],
                                   atm->fixpos.x*matrix[3]+atm->fixpos.y*matrix[4]+atm->fixpos.z*matrix[5],
                                   atm->fixpos.x*matrix[6]+atm->fixpos.y*matrix[7]+atm->fixpos.z*matrix[8]);
        }

bary_check=false;
    Center = getRotpos();

    for (ItCAtom itA = Atoms.begin(); itA != Atoms.end(); ++itA)
        {
        Atom *atm = *itA;
        atm->rotpos+=(Translate+transRigid-Center);
    }

}















///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// ATOMS //////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Atom& Molecule::addAtom(
               const std::string& atomicName,
               const std::string& name,
               const std::string& mol2t,
               Residu*            residu)

{
    if (residu != (Residu*)NULL && residu->molecule != this)
        throw MoleExcept(1070101,
                         "Molecule::addAtom",
                         "The given residu is not in this molecule: "
                         +residu->identifier);

    if (residu == (Residu*)NULL) residu = &unknownRes;
    Atom *atom=(Atom*)NULL;
    try{

      atom = new Atom(atomicName,this,residu,name,mol2t);
      atom->num=maxNumAtom;
      if (atom->fNum==0)atom->fNum=maxNumAtom;
      maxNumAtom++;
      residu->atoms.push_back(atom);
      Atoms.push_back(atom);
      bary_check=false;
      if (complex != (Complex*)NULL) complex->addAtom(atom);
    }
      catch (std::bad_alloc& ba)
      {
        Moleaccess=Levels::FATAL;
        cerr << ba.what()<<endl;
        throw MoleExcept(1070102,
                         "Molecule::addAtom",
                         "Bad allocation exception - exiting");
      }
      return *atom;
}


Atom& Molecule::addAtom(
        const std::string& atmName,
        const double& x,
        const double& y,
        const double& z,
        const std::string& name,
        const std::string& mol2t,
        Residu*  residu) 
{

    if (residu != (Residu*)NULL && residu->molecule != this)
        throw MoleExcept(1070201,
                         "Molecule::addAtom",
                         "The given residu is not in this molecule: "
                         +residu->identifier);

    if (residu == (Residu*)NULL) residu = &unknownRes;
  Atom *atom=(Atom*)NULL;
  try{

    atom = new Atom(atmName,this,residu,x,y,z,name,mol2t);
    atom->num=maxNumAtom;
    if (atom->fNum==0)atom->fNum=maxNumAtom;
    maxNumAtom++;
    residu->atoms.push_back(atom);
    Atoms.push_back(atom);
    if (complex != (Complex*)NULL) complex->addAtom(atom);
    bary_check=false;

  }
  catch (MoleExcept &e)
  {
    e.addTrace("Molecule::addAtom()");
    switch (e.getCode())
      {
      case 1050402: // No corresponding atom found during loadAtomicData()
        throw;break;

      case 1050401: // Atom is not biologically relevant
        switch (ICMole::Moleaccess)
          {
          case Levels::STRICT: throw;break;
          case Levels::NOTICE: cerr << "NOTICE : "<< e.getData()<<endl;break;
          case Levels::WARNING: cerr<< "WARNING : "<< e.getData()<<endl;break;
          }
        throw;
      }
  }
  catch (std::bad_alloc& ba)
  {
    Moleaccess=Levels::FATAL;
    cerr << ba.what()<<endl;
    throw MoleExcept(1070202,
                     "Molecule::addAtom",
                     "Bad allocation exception - exiting");
  }
  return *atom;

}






Atom& Molecule::addAtom(
        const std::string& atomicName,
        const Coords &coords,
        const std::string& name,
        const std::string& mol2t,
        Residu*  residu)  
{
    if (residu != (Residu*)NULL && residu->molecule != this)
        throw MoleExcept(1070201,
                         "Molecule::addAtom",
                         "The given residu is not in this molecule: "
                         +residu->identifier);

    if (residu == (Residu*)NULL) residu = &unknownRes;
  Atom *atom=(Atom*)NULL;
  try{

    atom = new Atom(atomicName,this,residu,coords,name,mol2t);
    atom->num=maxNumAtom;
    if (atom->fNum==0)atom->fNum=maxNumAtom;
    maxNumAtom++;
    residu->atoms.push_back(atom);
    Atoms.push_back(atom);
    if (complex != (Complex*)NULL) complex->addAtom(atom);
    bary_check=false;

  }
  catch (MoleExcept &e)
  {
    e.addTrace("Molecule::addAtom()");
    switch (e.getCode())
      {
    case 1070403:
      case 1050402: // No corresponding atom found during loadAtomicData()
        throw;break;

      case 1050401: // Atom is not biologically relevant
        switch (ICMole::Moleaccess)
          {
          case Levels::STRICT: throw;break;
          case Levels::NOTICE: cerr << "NOTICE : "<< e.getData()<<endl;break;
          case Levels::WARNING: cerr<< "WARNING : "<< e.getData()<<endl;break;
          }
        throw;
      }
  }
  catch (std::bad_alloc& ba)
  {
    Moleaccess=Levels::FATAL;
    cerr << ba.what()<<endl;
    throw MoleExcept(1070402,
                     "Molecule::addAtom",
                     "Bad allocation exception - exiting");
  }
  return *Atoms.at(Atoms.size()-1);
}

/*! \fn void Molecule::addAtoms(const unsigned int N)
 * \brief Creates N Atoms within this Molecule
 * \param N : Number of Atom to be created
 * \throw 1040102 : When bad allocation appends
 */

void Molecule::addAtoms(const unsigned int &N,
                        Residu* const residu) 
{
    if (N==0)return;
    if (N+Atoms.size() > Atoms.capacity()) Atoms.reserve(N+Atoms.size());
    if (residu != (Residu*)NULL && residu->molecule != this)
        throw MoleExcept(1070201,
                         "Molecule::addAtom",
                         "The given residu is not in this molecule: "
                         +residu->identifier);
    Atom *atom= (Atom*)NULL;
    bary_check=false;
    try
    {
    for (unsigned int i=0; i<N;i++)
    {
        atom= new Atom("Du",this,residu,"Du","Du");
        atom->num=maxNumAtom;
        Atoms.push_back(atom);
        if (complex != (Complex*)NULL) complex->addAtom(atom);
        maxNumAtom++;
        atom=(Atom*)NULL;
    }
    }
    catch (std::bad_alloc& ba)
    {
        if (atom != (Atom*)NULL)delete atom;
        string error("Bad allocation exception - exiting : ");
        error+=ba.what();
        throw MoleExcept(1040102,"Molecule::addAtoms",error);
    }
}



/*! \fn void Molecule::delAtom(Atom* const Atom, const bool with_maxnum)
  * \brief Remove the given Atom from the Molecule.
  * It will first search in the Molecule if this Atom is associated to this Molecule, then delete the Atom
  * It also call the cleanBond() function of the Atom and therefore delete all the associated Bonds of this Molecule
  * Also update the maxNumVe value
  * \throw 1040301 - No Atom given
  * \throw 1040302 - If the Molecule doesn't have the given Atom
  * \throw Exceptions from Molecule::delBond() && Molecule::cleanBond()
  * \param Atom : Atom to delete
  * \param with_maxnum : Update the maximum number for Atom (strong recommended)
  * \warning Will call cleanBond() function even if the given Atom is not in the Molecule
  */
void Molecule::delAtom( Atom *const atom, const bool& with_maxnum) 
{
    try{
        if (atom == (Atom*)NULL)
            throw MoleExcept(1040301,
                             "Molecule::delAtom",
                             "No Atom given");

// Searching the position of the Atom within the Atom list :
        const ItAtom itatom=find(Atoms.begin(), Atoms.end(),atom);
        if (itatom == Atoms.end())
            throw MoleExcept(1040302,
                             "Molecule::delAtom",
                             "This Molecule doesn't contains the following Atom : "+atom->toString());

// Deleting all Bonds related :
        // REMOVING BONDS :
        BondList& links= atom->links;
        const size_t linkSize=links.size();
        if (linkSize)
        for (size_t nlink=0;nlink <linkSize ;++nlink)
        {
            // Here the bond deletion imply a deleting in links vector
            // Therefore we must fetch the first element in vector and not the
            // nlink th element.
            const Bond* const bond=links.at(0);
            if (bond == (Bond*)NULL)
                throw MoleExcept(1040403,
                                 "Molecule::delAtoms",
                                 "No Bond given");
#ifdef ICHEM_DEBUG
                cout << "DEBUG|ATOM|DELATOMS|Deleting bond : " <<bond->toString()<<endl;
#endif
            // Check Bond existence
            const ItBond itPED=find(Bonds.begin(), Bonds.end(),bond); // Deletion if it finds in Bonds the corresponding Bond
            if (itPED == Bonds.end())
                throw MoleExcept(1040404,
                                 "Molecule::delAtoms",
                                 "This Molecule doesn't contains the following Bond : "+bond->toString());

            // Deleting the Bond from Molecule :
            if (complex != (Complex*)NULL) complex->deleteBond(links.at(0));
            delete bond;
            Bonds.erase(itPED);

        }
        links.clear();
        const unsigned int numatm= atom->getNum();

// Deleting cycle in which it can be involved with:
        if (atom->inCycle)
        {
            for (ItCycle it =Cycles.begin(); it != Cycles.end(); it++)
            {
                Cycle* cycle=*it;

                for (ItCAtom itA = cycle->first();itA != cycle->end();itA++)
                {
                    if (*itA != atom)continue;
                    Cycles.erase(it);
                    it=Cycles.begin();break;
                }
            }
        }
// Delete residu when no atoms in it:
        if (atom->getResidu() != (Residu*)NULL)
        {
            const Residu &residu = *atom->getResidu();
            if (residu.numAtom()==1)
            {
                const ItRes itR = find(Residues.begin(), Residues.end(), &residu);
                if (itR == Residues.end())
                    throw MoleExcept(1040303,
                                     "Molecule::delAtom",
                                     "Residu of the given atom is not part of this molecule");

                // Deleting chain when no residu left except this one:
                if (residu.getChain()->getCountRes()==1)
                {
                    const ItChain itC = find(Chains.begin(),Chains.end(),residu.getChain());
                    if (itC == Chains.end())
                        throw MoleExcept(1040304,
                                         "Molecule::delAtom",
                                         "Chain of the given atom is not part of this molecule");
                    if (complex != (Complex*)NULL) complex->deleteChain(residu.getChain());
                    delete *itC;
                    Chains.erase(itC);
                }
                if (complex != (Complex*)NULL) complex->deleteResidu(atom->getResidu());
                delete *itR;
                Residues.erase(itR);


            }
            atom->getResidu()->delAtom(*atom);

            // Delete chain when no residu in it :

        }

// Deleting the Bond from Molecule :
        delete *itatom;
        Atoms.erase(itatom);
        bary_check=false;
        if (complex != (Complex*)NULL) complex->deleteAtom(atom);
        if (!with_maxnum)return;

// Updating the maxNumVe - Unecessary if the Atom num is not the maxNumVe
        if (numatm+1 != maxNumAtom) return;
        bool zer=false;
        maxNumAtom=1;
        for(ItAtom it2 = Atoms.begin(); it2 != Atoms.end(); it2++)
        {
            if ((*it2)->getNum() < maxNumAtom)continue;
             maxNumAtom= (*it2)->getNum();
             if (maxNumAtom==1) zer=true;
        }
        if (maxNumAtom >1|| zer==true) maxNumAtom++;
    }
    catch (MoleExcept &e)
    {
        e.addTrace("Molecule::delAtom()");
        switch (e.getCode())
        {
        case 1040701:
        case 1040702:
        case 1040202:
        case 1040101:
            switch (Moleaccess)
            {
            case Levels::NOTICE:
            case Levels::WARNING:cerr <<e.getCode()<<"\t"
                                      <<e.getData()<<" in "
                                      <<e.getSource()<<endl;break;
            case Levels::STRICT:
            case Levels::FATAL:cerr <<e.getCode()<<"\t"
                                    <<e.getData()<<" in "
                                    <<e.getSource()<<endl;
                                e.setCode(999999);throw;
            }
        }

    }

}


void Molecule::delAtoms( AtomList& AtomsDel) 
{
    try{

        if (AtomsDel.empty())return;
        const size_t sizeAtm = AtomsDel.size();
        for (size_t nAtom = 0; nAtom < sizeAtm; nAtom++)
        {
        Atom* atom = AtomsDel.at(nAtom);
#ifdef ICHEM_DEBUG
                cout << "DEBUG|ATOM|DELATOMS|Deleting atom : "<< atom->getIdentifier()<<endl;

#endif
        if (atom == (Atom*)NULL)
            throw MoleExcept(1040401,
                             "Molecule::delAtoms",
                             "No Atom given");


         ItAtom itATM=find(Atoms.begin(), Atoms.end(),atom); // Deletion if it finds in Atoms the corresponding Atom
        if (itATM== Atoms.end())
            throw MoleExcept(1040402,
                             "Molecule::delAtoms",
                             "This Molecule doesn't contains the following Atom : "+atom->toString());
        // REMOVING BONDS :
        BondList& links= atom->links;
        const size_t linkSize=links.size();
        if (linkSize)
        for (size_t nlink=0;nlink <linkSize ;++nlink)
        {
            // Here the bond deletion imply a deleting in links vector
            // Therefore we must fetch the first element in vector and not the
            // nlink th element.
            const Bond* const bond=links.at(0);
            if (bond == (Bond*)NULL)
                throw MoleExcept(1040403,
                                 "Molecule::delAtoms",
                                 "No Bond given");
#ifdef ICHEM_DEBUG
                cout << "DEBUG|ATOM|DELATOMS|Deleting bond : " <<bond->toString()<<endl;
#endif
            // Check Bond existence
            const ItBond itPED=find(Bonds.begin(), Bonds.end(),bond); // Deletion if it finds in Bonds the corresponding Bond
            if (itPED == Bonds.end())
                throw MoleExcept(1040404,
                                 "Molecule::delAtoms",
                                 "This Molecule doesn't contains the following Bond : "+bond->toString());

            // Deleting the Bond from Molecule :
             if (complex != (Complex*)NULL) complex->deleteBond(links.at(0));
            delete bond;
            Bonds.erase(itPED);

        }
        links.clear();
        //REMOVING ATOM FROM RESIDU :
        if (atom->getResidu() != (Residu*)NULL)
        {
            Residu &residu = *atom->getResidu();
            residu.delAtom(*atom);
            if (residu.numAtom()==0)
            {
                ItRes itR = find(Residues.begin(), Residues.end(), &residu);
                if (itR == Residues.end())
                    throw MoleExcept(1040405,
                                     "Molecule::delAtoms"
                                     "Residu is not part of this protein");
#ifdef ICHEM_DEBUG
                cout << "DEBUG|ATOM|DELATOMS|Deleting residu : "
                     << residu.getIdentifier()<<endl;
#endif
                Chain &chain = *residu.getChain();
                if (chain.getCountRes()==0)
                {
                    ItChain itC = find(Chains.begin(), Chains.end(), &chain);
                    if (itC == Chains.end())
                        throw MoleExcept(1040406,
                                         "Molecule::delAtoms"
                                         "Chain is not part of this protein");
                    if (complex != (Complex*)NULL) complex->deleteChain(residu.getChain());
                    delete *itC;
                    Chains.erase(itC);
                }
                if (complex != (Complex*)NULL) complex->deleteResidu(atom->getResidu());
                delete &residu;
                Residues.erase(itR);
            }
        }
if (complex != (Complex*)NULL) complex->deleteAtom(atom);
        delete *itATM;
        Atoms.erase(itATM);

        }



        // Updating the maxNumVe - Unecessary if the Atom num is not the maxNumVe
        bool zer=false;
        maxNumAtom=1;
        bary_check=false;
        for(ItCAtom it2 = Atoms.begin(); it2 != Atoms.end(); it2++)
        {
            if ((*it2)->getNum() < maxNumAtom) continue;
             maxNumAtom= (*it2)->getNum();if (maxNumAtom==1) zer=true;
        }
        if (maxNumAtom >1|| zer==true) maxNumAtom++;

        maxNumBond=1;
         zer=false;
        for(ItCBond it2 = Bonds.begin(); it2 != Bonds.end(); it2++)
          {
            if ((*it2)->getNum() >= maxNumBond) {maxNumBond= (*it2)->getNum();if (maxNumBond==1) zer=true;}
          }
        if (maxNumBond >1|| zer==true) maxNumBond++;
    }
    catch (MoleExcept &e)
    {
        e.addTrace("Molecule::delAtom()");
        switch (e.getCode())
        {
        case 1040701:
        case 1040702:
        case 1040202:
        case 1040101:
            switch (Moleaccess)
            {
            case Levels::NOTICE:
            case Levels::WARNING:cerr <<e.getCode()<<"\t"<<e.getData()<<" in "<<e.getSource()<<endl;break;
            case Levels::STRICT:
            case Levels::FATAL:cerr <<e.getCode()<<"\t"<<e.getData()<<" in "<<e.getSource()<<endl; e.setCode(999999);throw;
            }
        }

    }

}




/*! \fn void Molecule::renumAtom()
  *\brief Renum all Atoms of the Molecule
  *
  * Scan all Atom of the Molecule and renumerotate all. Call setNum() function for all Atomes
  */
void Molecule::renumAtom()
{
  maxNumAtom=1;
  for (ItRes itR = Residues.begin(); itR != Residues.end();++itR)
  {
      Residu &residu =**itR;
      for (ItAtom it = residu.firstAtom(); it != residu.lastAtom();++it)
      {
          ((*it))->setNum(maxNumAtom);
          maxNumAtom++;
      }
  }
 /* for (ItAtom it = Atoms.begin(); it != Atoms.end(); it++)
    {
      ((*it))->setNum(maxNumAtom);
      maxNumAtom++;
    }*/

}


/*! \fn Atom& Molecule::getAtom(const unsigned int n,bool pos)
 *  \brief Return the Atom associated with the number n.
 *  \param pos: When pos=true (default) : n corresponds to the position within the \
 *  list of Atoms. Otherwise if the id thought getNum() function that will be checked to n
 *  \param n : Position
 *  \return The corresponding Atom
 *  \throw 1040801- When pos is true and n is above the number of Atom in the Molecule
 *
 */
const Atom& Molecule::getAtom(const unsigned int &n, const bool &pos) const 
{
  if (pos)
    {
      if (n >= Atoms.size()) throw MoleExcept(1040501,"Molecule::getAtom"," given value above the number of Atoms within the Molecule");
      return *Atoms.at(n);
    }
  for (size_t i=0; i<Atoms.size();i++) if (Atoms.at(i)->getNum() == n) return *Atoms.at(i);

  throw MoleExcept(1040502,"Molecule::getAtom"," No Atom found with the given id");

}





/**
 * @brief Molecule::reserveAtomSize
 * @param NAtm : Number of atom that we expect
 * Requests that the atom array capacity be at least enough to contain
 * NAtm elements.
 */
void Molecule::reserveAtomSize(const unsigned int &NAtm)
{

    Atoms.reserve(NAtm+Atoms.size());

}











///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// BondS ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*!
  * \fn Bond& Molecule::addBond(const Atom &ve1,const Atom &ve2) 
  * \param ve1 : First Atom of the Bond
  * \param ve2 : Second Atom of the Bond
  * \throw MoleExcept code 1040601 when ve1 and ve2 are the same Atom. Returns the last added Bond
  * \throw MoleExcept code 1040602 when ve1 are not part of this Molecule. Returns the last added Bond
  * \throw MoleExcept code 1040603 when ve2 are not part of this Molecule. Returns the last added Bond
  * \throw MoleExcept code 1040604 when memory allocation of the Bond failed : Fatal error
  * \return Reference of the newly created Bond
  * Add an Bond to Bonds list. <br/>
  * Each time a bond is added to a molecule, the Bond associated to the bond is added to the Molecule associated to the molecule
  */
Bond& Molecule::addBond(Atom &Atom1, Atom &Atom2, const unsigned int &BondType) 
{
    Bond *bd=(Bond*)NULL;
    try
    {
        if (&Atom1 == &Atom2)         throw MoleExcept(1040601, "Molecule::addBond", "Atoms given for the creation of the Bond are the same");
        if (&Atom1.getParent() != this) throw MoleExcept(1040602, "Molecule::addBond", "First Atom is not part of this Molecule");
        if (&Atom2.getParent() != this) throw MoleExcept(1040603, "Molecule::addBond", "Second Atom is not part of this Molecule");
        bd = new Bond(Atom1,Atom2,maxNumBond,BondType,*this);

        maxNumBond++;
        Bonds.push_back(bd);

        if (typemol==MoleType::PROTEIN &&
              Atom1.getResidu()!= Atom2.getResidu()
             && Atom1.getResidu()!= (Residu*)NULL
             && Atom2.getResidu()!= (Residu*)NULL)
        {
            //cout << "A"<<Atom1.getIdentifier()<<"\t"<<Atom2.getIdentifier()<<endl;
            Atom1.getResidu()->InterRes.push_back(Atom2.getResidu());
            Atom2.getResidu()->InterRes.push_back(Atom1.getResidu());
        }
        if (complex != (Complex*)NULL) complex->addBond(bd);
        return *bd;
    }
    catch(MoleExcept &e)
    {
        if (bd!=(Bond*)NULL) delete bd;
        e.addTrace("Bond::addBond");
        throw;
    }
    catch (std::bad_alloc& ba)
    {
        Moleaccess=Levels::FATAL;
        string error("Bad allocation exception - exiting :"); error+=ba.what();
        throw MoleExcept(1040604,"Molecule::addBond", error);
    }
}






Bond&  Molecule::addBond( Atom *const Atom1,
                       Atom *const Atom2,
                       const unsigned int &BondType) 
{
    Bond *bd=(Bond*)NULL;
    try
    {
        if (Atom1 == (Atom*)NULL)    throw MoleExcept(1040601, "Molecule::addBond", "Atom 1 is null");
        if (Atom2 == (Atom*)NULL)    throw MoleExcept(1040602, "Molecule::addBond", "Atom 2 is null");
        if (Atom1 == Atom2)          throw MoleExcept(1040601, "Molecule::addBond", "Atoms given for the creation of the Bond are the same");
        if (Atom1->molecule != this) throw MoleExcept(1040602, "Molecule::addBond", "First Atom is not part of this Molecule");
        if (Atom2->molecule != this) throw MoleExcept(1040603, "Molecule::addBond", "Second Atom is not part of this Molecule");
        bd = new Bond(*Atom1,*Atom2,maxNumBond,BondType,*this);
        if (typemol==MoleType::PROTEIN &&
              Atom1->getResidu()!= Atom2->getResidu()
             && Atom1->getResidu()!= (Residu*)NULL
             && Atom2->getResidu()!= (Residu*)NULL)
        {
            Atom1->getResidu()->InterRes.push_back(Atom2->getResidu());
            Atom2->getResidu()->InterRes.push_back(Atom1->getResidu());
        }
        if (complex != (Complex*)NULL) complex->addBond(bd);
        maxNumBond++;
        Bonds.push_back(bd);
        return *bd;
    }
    catch(MoleExcept &e)
    {
        if (bd!=(Bond*)NULL) delete bd;
        e.addTrace("Bond::addBond");
        throw;
    }
    catch (std::bad_alloc& ba)
    {
        Moleaccess=Levels::FATAL;
        string error("Bad allocation exception - exiting :"); error+=ba.what();
        throw MoleExcept(1040604,"Molecule::addBond", error);
    }
}






//Atom::cleanBond
void Molecule::delBond(const Bond* const bond) 
{
    if (bond == (Bond*)NULL) throw MoleExcept(1040701,"Molecule::delBond","No Bond given");
    const unsigned int NumBond= bond->getNum();
// Check Bond existence
    const ItBond itPED=find(Bonds.begin(), Bonds.end(),bond); // Deletion if it finds in Bonds the corresponding Bond
    if (itPED == Bonds.end())throw MoleExcept(1040702,"Molecule::delBond","This Molecule doesn't contains the following Bond : "+bond->toString());
     if (complex != (Complex*)NULL) complex->deleteBond(bond);
// Deleting the Bond from Molecule :
    delete *itPED;
    Bonds.erase(itPED);


// Updating the maxNumBond - Unecessary if the Bond num is not the maxNumBond
    if (NumBond+1 != maxNumBond) return;
    maxNumBond=1;
    bool zer=false;
    for(ItCBond it2 = Bonds.begin(); it2 != Bonds.end(); it2++)
      {
        if ((*it2)->getNum() >= maxNumBond) {maxNumBond= (*it2)->getNum();if (maxNumBond==1) zer=true;}
      }
    if (maxNumBond >1|| zer==true) maxNumBond++;
}


void Molecule::cleanBond()
{
    for(size_t i = 0; i < Bonds.size(); ++i)
    {
        if (complex != (Complex*)NULL) complex->deleteBond(Bonds[i]);
        // Deleting the Bond from Molecule :
        Residu* const res1=Bonds[i]->getAtom1().getResidu();
        Residu* const res2=Bonds[i]->getAtom2().getResidu();
        if (res1 != (Residu*)NULL
          &&res2 != (Residu*)NULL
          &&res1!=res2)
        {
            const ItRes pos1=find(res1->InterRes.begin(),res1->InterRes.end(),res2);
            if (pos1 != res1->InterRes.end())res1->InterRes.erase(pos1);
            const ItRes pos2=find(res2->InterRes.begin(),res2->InterRes.end(),res1);
            if (pos2 != res2->InterRes.end())res2->InterRes.erase(pos2);
        }


        delete Bonds[i];
    }
    Bonds.clear();
}

void Molecule::delBond(const BondList& EdList) 
{
    for (ItCBond itEd = EdList.begin() ; itEd != EdList.end(); itEd++)
    {
        const Bond* const bond=*itEd;
        if (bond == (Bond*)NULL) throw MoleExcept(1040801,"Molecule::delBond","No Bond given");
        // Check Bond existence
        const ItBond itPED=find(Bonds.begin(), Bonds.end(),bond); // Deletion if it finds in Bonds the corresponding Bond
        if (itPED == Bonds.end())throw MoleExcept(1040802,"Molecule::delBond","This Molecule doesn't contains the following Bond : "+bond->toString());
        if (complex != (Complex*)NULL) complex->deleteBond(bond);
        // Deleting the Bond from Molecule :
        Residu* const res1=bond->getAtom1().getResidu();
        Residu* const res2=bond->getAtom2().getResidu();
    if (res1 != (Residu*)NULL
      &&res2 != (Residu*)NULL
      &&res1!=res2)
    {
        const ItRes pos1=find(res1->InterRes.begin(),res1->InterRes.end(),res2);
        if (pos1 != res1->InterRes.end())res1->InterRes.erase(pos1);
        const ItRes pos2=find(res2->InterRes.begin(),res2->InterRes.end(),res1);
        if (pos2 != res2->InterRes.end())res2->InterRes.erase(pos2);
    }
        delete *itPED;
        Bonds.erase(itPED);

    }

// Updating the maxNumBond - Unecessary if the Bond num is not the maxNumBond
    maxNumBond=1;
    bool zer=false;
    for(ItCBond it2 = Bonds.begin(); it2 != Bonds.end(); it2++)
      {
        if ((*it2)->getNum() >= maxNumBond) {maxNumBond= (*it2)->getNum();if (maxNumBond==1) zer=true;}
      }
    if (maxNumBond >1|| zer==true) maxNumBond++;
}

/*! \fn void   Molecule::delBond  (const Bond&  ed) 
  * \brief delete the given Bond from this Molecule
  * \throw 1040301 - Given Bond not part of this Molecule
  * \throw 1020101 when the given Bond is not found in one of its two Atoms
  * \param ed : Bond to delete
  * Remove the given Bond from the Molecule.
  * It will first search in the Molecule if this Bond is associated to this Molecule, then delete the Bond
  * The deletion of the Bond automatically call delBond() function of Atoms involved.
  * Also update the maxNumBond value
  */
void  Molecule::delBond(const       Bond &  bond) 
{
  const unsigned int NumBond= bond.getNum();

  // Check Bond existence
  const ItBond itPED=find(Bonds.begin(), Bonds.end(),&bond); // Deletion if it finds in Bonds the corresponding Bond
  if (itPED == Bonds.end())throw MoleExcept(1040901,"Molecule::delBond","This Molecule doesn't contains the following Bond : "+bond.toString());
if (complex != (Complex*)NULL) complex->deleteBond(&bond);

    Residu* const res1=bond.getAtom1().getResidu();
    Residu* const res2=bond.getAtom2().getResidu();
if (res1 != (Residu*)NULL
  &&res2 != (Residu*)NULL
  &&res1!=res2)
{
    const ItRes pos1=find(res1->InterRes.begin(),res1->InterRes.end(),res2);
    if (pos1 != res1->InterRes.end())res1->InterRes.erase(pos1);
    const ItRes pos2=find(res2->InterRes.begin(),res2->InterRes.end(),res1);
    if (pos2 != res2->InterRes.end())res2->InterRes.erase(pos2);
}
  delete *itPED;
  Bonds.erase(itPED);

  // Updating the maxNumBond - Unecessary if the Bond num is not the maxNumBond
  if (NumBond+1 != maxNumBond) return;
  maxNumBond=1;
  bool zer=false;
  for(ItCBond it2 = Bonds.begin(); it2 != Bonds.end(); it2++)
    {
      if ((*it2)->getNum() >= maxNumBond) {maxNumBond= (*it2)->getNum();if (maxNumBond==1) zer=true;}
    }
  if (maxNumBond >1|| zer==true) maxNumBond++;
}




/** \fn void Molecule::renumBond()
  \brief Renum all Bonds of the Molecule

  Scan all Bonds of the Molecule and renumerotate all. Call setNum() function for all Bonds
  */
void Molecule::renumBond()
{
  maxNumBond=1;
  for (ItBond it= Bonds.begin(); it != Bonds.end(); it++){
      ((*it))->num=maxNumBond;
      maxNumBond++;
    }
}

void Molecule::renumResidu()
{
  maxNumResidu=1;
  for (ItRes it= Residues.begin(); it != Residues.end(); it++){
      ((*it))->num=maxNumResidu;
      maxNumResidu++;
    }
}

/**
 * @brief Molecule::reserveBondSize
 * @param NBond : Number of bond to add
 * Requests that the bond array capacity be at least enough to contain
 * NBond elements.
 */
void Molecule::reserveBondSize(const unsigned int &NBond)
{

    Bonds.reserve(NBond+Bonds.size());

}



/**
 * @brief Molecule::getBond
 * @param n:Position in the bond array
 * @param pos : When pos=true (default) : n corresponds to the position within the \
 *  list of Atoms. Otherwise if the id thought getNum() function that will be checked to n
 * @return The corresponding Bond
 */
const Bond& Molecule::getBond(const size_t &n, const bool &pos) const 
{
  if (pos)
    {
      if (n >= Bonds.size()) throw MoleExcept(1040501,"Molecule::getBond"," given value above the number of Atoms within the Molecule");
      return *Bonds.at(n);
    }
  for (size_t i=0; i<Bonds.size();i++) if (Bonds.at(i)->getNum() == n) return *Bonds.at(i);

  throw MoleExcept(1040502,"Molecule::getBond"," No Atom found with the given id");

}








////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// RESIDUS ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Molecule::reserveResiduSize
 * @param NResidu : Number of residu to contain
 * Requests that the residu capacity be at least enough to contain
 * NResidu elements.
 */
void Molecule::reserveResiduSize(const unsigned int &NResidu)
{
    Residues.reserve(Residues.size()+NResidu);
}



Residu& Molecule::getResidu(const size_t &pos) 
{
  if (pos >= Residues.size()) throw MoleExcept(1071301,"Molecule::getResidu","Position is above the number of residu within the molecule");
  else return *Residues.at(pos);
}

Chain& Molecule::getChain(const size_t pos) 
{
    if (pos >= Chains.size()) throw MoleExcept(1071401,"Molecule::getChain","Position is above the number of chains within the molecule");
    else return *Chains.at(pos);
}



/**
 * @brief Molecule::addResidu
 * @param chain : Name of the chain related to this residu
 * @param resid : residu ID as given by the file
 * @param resname : residu name
 * @return The newly created residu when it does not exist
 *  or the existing residu
 */
Residu& Molecule::addResidu(const std::string& chain,
                            const  int& resid,
                            const std::string& resname)
{
  // Aim to check if the given chain already exists or not
  Chain *check_chain=(Chain*)NULL;

  // Scanning chains
  if (!Chains.empty())
  for (ItCChain itC = Chains.begin(); itC != Chains.end(); itC++)
  {
      if ((*itC)->getName() != chain)continue;
      check_chain=*itC;
      // Since same chain name
      // Scan residu to look whether the residu already exists :
      for (ItCRes itCR = (*itC)->firstC(); itCR != (*itC)->lastC(); itCR++)
      {
          if ((*itCR)->getFNum() == resid &&
              (*itCR)->getName() == resname) return **itCR;
      }
  }

// The given residu doesn't exists, we create it :
try
  {
  // First by creating the chain when not existing
  if (check_chain == (Chain*)NULL)
    {
      check_chain = new Chain(this,chain);
      if (complex != (Complex*)NULL) complex->addChain(check_chain);
      Chains.push_back(check_chain);

    }
  // Then the residu
  Residu *res = new Residu(this,
                           maxNumResidu,
                           resid,
                           resname,
                           *check_chain,
                           ResType::UNKNOWN);
  if (complex != (Complex*)NULL) complex->addResidu(res);
  Residues.push_back(res);maxNumResidu++;
  return *res;
  }
  catch (std::bad_alloc& ba)
  {
    Moleaccess=Levels::FATAL;
    cerr << ba.what()<<endl;
    throw MoleExcept(1071001,"Molecule::addResidu", "Bad allocation exception - exiting");
  }

}


void Molecule::delResidu(Residu& residu)
{
    ItRes itt = std::find(Residues.begin(),Residues.end(),&residu);

    if (&residu == &unknownRes) throw MoleExcept(1071102,"Molecule::delResidu","Temporary residu cannot be deleted");
    if (&residu == &aromaticRes) throw MoleExcept(1071103,"Molecule::delResidu","Cyclic residu cannot be deleted");
    if (&residu == &cyclicRes) throw MoleExcept(1071103,"Molecule::delResidu","Cyclic residu cannot be deleted");
    if (itt == Residues.end())throw MoleExcept(1071104,"Molecule::delResidu","No residu in the molecule");
    if (residu.getChain()->getCountRes() == 1)
      { Chain * const chain=residu.getChain();
        ItChain itC=std::find(Chains.begin(),Chains.end(), chain);
        if (itC != Chains.end()) Chains.erase(itC);
      }
    else residu.getChain()->delResidu(&residu);

    AtomList atomToDel(residu.atoms);

    delAtoms(atomToDel);





}


void Molecule::setResToAtom(Atom *const atm, Residu* const res)
{
  if (atm == (Atom*)NULL)return;
  if (res == (Residu*)NULL)return;
  if (&atm->getParent() != this)
    {

      atm->molecule=this;

    }
  if (atm->getResidu() != (Residu*)NULL)
    {
      atm->getResidu()->delAtom(*atm);

    }
  res->addAtom(*atm);
  atm->residu=res;

}




Cycle* Molecule::getCycleFromCenter(Atom* const atom)
{
    if (atom == (Atom*)NULL)return (Cycle*)NULL;
    for (ItCycle itC = Cycles.begin();
         itC != Cycles.end();
         itC++)
    {
         Cycle& cycle = **itC;
        if (&cycle.getCenter()== atom) return &cycle;
    }
    return (Cycle*)NULL;
}



void Molecule::delCycle(Cycle*const cy) 
{
  if (cy == (Cycle*)NULL) throw MoleExcept(1071201,"Molecule::delCycle","No cycle given");
  ItCycle itt = std::find(Cycles.begin(),Cycles.end(),cy);
  if (itt == Cycles.end()) throw MoleExcept(1071202,"Molecule::delCycle","No cycle in this molecule");
    delete cy;
  Cycles.erase(itt);

}




/**
 * @brief Molecule::moveResidu
 * @param residu
 * @param molecule
 * @param delBonds
 * @warning : Can cause issues when switching residu from a molecule to another molecule
 * of another complex than the reference one
 */
void Molecule::moveResidu(Residu& residu,
                          Molecule& molecule,
                          const bool &delBonds) 
{

    // First check that the given residu is indeed in this molecule:
     ItRes respos = find(Residues.begin(), Residues.end(), &residu);
    if (respos == Residues.end())
        throw MoleExcept(1071301,
                         "Molecule::moveResidu",
                         "Given residu ("+residu.getIdentifier()+") is not part of this molecule "+name);
    // Getting the list of atoms:
    const AtomList listAtoms = residu.atoms;


    BondList listBonds, // List of bonds that are in this residu:
             failedBond;// List of bonds that are between this residu and another one
    CycleList listCycle;
    // Checking that all atoms are within this molecule:
    // And getting all bonds
    for (ItCAtom itA = listAtoms.begin();
                 itA!= listAtoms.end();
               ++itA)
    {
        const Atom& atom = **itA;


        if (find(Atoms.begin(),Atoms.end(),*itA) == Atoms.end())
        {
            throw MoleExcept(1071302,
                             "Molecule::moveResidu",
                             "Atom ("+(*itA)->getIdentifier()+") is not part of this molecule "+name);
        }
        const size_t nBonds = atom.getNumBond();
        // Now fetching all bonds that are within this residu
        for (size_t iBd=0; iBd < nBonds;++iBd)
        {
             Bond &bond = *atom.getBond(iBd);
            if (bond.getOtherAtom(atom).getResidu() != &residu)
            {
                failedBond.push_back(&bond);
                continue;
            }

                listBonds.push_back(&bond);

        }// END iBD
        // Fetching cycles:
        if (atom.props.isInRing())
        {
            for (ItCCycle itCC = Cycles.begin();
                          itCC != Cycles.end();
                          itCC++)
            {
                Cycle& cycle = **itCC;
                if (cycle.hasAtom(atom) &&
                    find(listCycle.begin(),
                         listCycle.end(),
                         &cycle)==listCycle.end())
                    listCycle.push_back(&cycle);
            }
        }

    }// END ItCAtom
    sortAndUnique(listBonds);
    // Checking that there is no inter-residu bonds:
    if (!delBonds && !failedBond.empty())
    {
        throw MoleExcept(1071303,
                         "Molecule::moveResidu",
                         "A bond has been detected between two residues without possibility of deletion in molecule "+name);
    }
    else if (!failedBond.empty())// Or we delete them all
    {
        for (ItCBond itB = failedBond.begin(); itB != failedBond.end();++itB)
        {
            Bond &bond = **itB;
            delBond(&bond);
        }
    }

    // Now moving all atoms into new molecule:
    for (size_t i =0;i< listAtoms.size();++i)
    {
        Atom &atom=*listAtoms.at(i);

      molecule.Atoms.push_back(&atom);
      const ItAtom itpos= find(Atoms.begin(),Atoms.end(),&atom);
      Atoms.erase(itpos);
      atom.molecule=&molecule;
    }

    // Now moving all bonds into new molecule:
    if (!listBonds.empty())
    for (size_t i=0; i< listBonds.size();++i)
    {
        Bond &bond=*listBonds.at(i);

        ItBond itBB= find(Bonds.begin(),Bonds.end(),&bond);
        if (itBB == Bonds.end())
        {

            throw MoleExcept(1071303,
                             "Molecule::moveResidu",
                             "Unable to find bond in molecule "+name+"\n"+bond.toString());
        }

            Bonds.erase(itBB);
            molecule.Bonds.push_back(&bond);
        bond.molecule=&molecule;

    }
    // Now moving all cycle into new molecule:
    if (!listCycle.empty())
    for (size_t i=0; i< listCycle.size();++i)
    {
        molecule.Cycles.push_back(listCycle.at(i));
        Cycles.erase(find(Cycles.begin(),Cycles.end(),listCycle.at(i)));
    }
    // Now moving residu:
    molecule.Residues.push_back(&residu);
    residu.molecule=&molecule;
    // Now updating chains :
    Chain *mole_chain = (Chain*)NULL;
    Chain& former_chain = *residu.chain;
    for (ItChain itcc=molecule.Chains.begin(); itcc != molecule.Chains.end();++itcc)
    {
        Chain &chain = **itcc;
        if (chain.getName()==residu.chain->getName())mole_chain=&chain;
    }

    if (mole_chain != (Chain*)NULL) mole_chain->addResidu(&residu);
    else
    {
        Chain *check_chain = new Chain(molecule,residu.chain->getName());
        molecule.Chains.push_back(check_chain);
        check_chain->addResidu(&residu);
    }
    // Final check, see if former chain has residues
    if (former_chain.getCountRes()==0)
    {
        const ItChain itc=find(Chains.begin(),Chains.end(),&former_chain);
        delete &former_chain;
        Chains.erase(itc);
    }

Residues.erase(respos);


}











void Molecule::moveResidu(ResiduList& residuList,
                          Molecule& molecule,
                          const bool &delBonds) 
{


    AtomList listAtoms;
       // First check that the given residu is indeed in this molecule:
    for (ItRes ipos = residuList.begin(); ipos != residuList.end();++ipos)
    {
        Residu &residu=**ipos;
        const ItRes respos=find(Residues.begin(), Residues.end(), &residu);
        if (respos == Residues.end())
            throw MoleExcept(1071301,
                             "Molecule::moveResidu",
                             "Given residu ("+residu.getIdentifier()+") is not part of this molecule "+name);
        listAtoms.insert(listAtoms.end(),residu.firstAtom(),residu.lastAtom());

    }


    BondList listBonds, // List of bonds that are in this residu:
             failedBond;// List of bonds that are between this residu and another one
    CycleList listCycle;

    // Checking that all atoms are within this molecule:
    // And getting all bonds
    for (ItCAtom itA = listAtoms.begin();
                 itA!= listAtoms.end();
               ++itA)
    {
        const Atom& atom = **itA;


        if (find(Atoms.begin(),Atoms.end(),&atom) == Atoms.end())
        {
            throw MoleExcept(1071302,
                             "Molecule::moveResidu",
                             "Atom ("+(*itA)->getIdentifier()+") is not part of this molecule "+name);
        }

        const size_t nBonds = atom.getNumBond();
        // Now fetching all bonds that are within this residu
        for (size_t iBd=0; iBd < nBonds;++iBd)
        {
             Bond &bond = *atom.getBond(iBd);
            if (find(residuList.begin(),
                     residuList.end(),
                     bond.getOtherAtom(atom).getResidu())==residuList.end())
            {
                failedBond.push_back(&bond);
                continue;
            }
            listBonds.push_back(&bond);

        }// END iBD
        // Fetching cycles:
        if (atom.props.isInRing())
        {
            for (ItCCycle itCC = Cycles.begin();
                          itCC != Cycles.end();
                          itCC++)
            {
                Cycle& cycle = **itCC;
                if (cycle.hasAtom(atom) &&
                    find(listCycle.begin(),
                         listCycle.end(),
                         &cycle)==listCycle.end())
                    listCycle.push_back(&cycle);
            }
        }

    }// END ItCAtom
     sortAndUnique(listBonds);
    // Checking that there is no inter-residu bonds:
    if (!delBonds && !failedBond.empty())
    {
        throw MoleExcept(1071303,
                         "Molecule::moveResidu",
                         "A bond has been detected between two residues without possibility of deletion in molecule "+name);
    }
    else if (!failedBond.empty())// Or we delete them all
    {
        for (ItCBond itB = failedBond.begin(); itB != failedBond.end();++itB)
        {
            Bond &bond = **itB;
            delBond(&bond);
        }
    }

    // Now moving all atoms into new molecule:
    for (size_t i =0;i< listAtoms.size();++i)
    {
        Atom &atom=*listAtoms.at(i);
      molecule.Atoms.push_back(&atom);
      const ItAtom itpos=find(Atoms.begin(),Atoms.end(),&atom);
      Atoms.erase(itpos);
      atom.molecule=&molecule;
    }

    // Now moving all bonds into new molecule:
    if (!listBonds.empty())
    for (size_t i=0; i< listBonds.size();++i)
    {
        Bond &bond=*listBonds.at(i);
        molecule.Bonds.push_back(&bond);
        ItBond itBB= find(Bonds.begin(),Bonds.end(),&bond);
        if (itBB == Bonds.end())
        {

            throw MoleExcept(1071303,
                             "Molecule::moveResidu",
                             "Unable to find bond in molecule "+name+"\n"+bond.toString());
        }

            Bonds.erase(itBB);
        bond.molecule=&molecule;

    }
    // Now moving all cycle into new molecule:
    if (!listCycle.empty())
    for (size_t i=0; i< listCycle.size();++i)
    {
        molecule.Cycles.push_back(listCycle.at(i));
        Cycles.erase(find(Cycles.begin(),Cycles.end(),listCycle.at(i)));
    }
    ChainList chainsList;
    // Now moving residu:
    molecule.Residues.insert(molecule.Residues.end(),residuList.begin(),residuList.end());
    for (ItRes ipos = residuList.begin(); ipos != residuList.end();++ipos)
    {
        Residu &residu=**ipos;
        if (find(chainsList.begin(),chainsList.end(),residu.getChain())==chainsList.end())chainsList.push_back(residu.getChain());
    residu.molecule=&molecule;
    }

    map<Chain*,Chain*> oldToNewChain;
    // Now updating chains :
    for (ItChain itCh=chainsList.begin();itCh != chainsList.end();itCh++)
    {
        Chain *mole_chain = (Chain*)NULL;
        Chain& former_chain = **itCh;
        // Checking if the new molecule has a chain with the same name:
        for (ItChain itcc  = molecule.Chains.begin();
                     itcc != molecule.Chains.end();
                   ++itcc)
        {
            Chain &chain = **itcc;
            if (chain.getName()==former_chain.getName())
            {
                mole_chain=&chain;
                break;
            }
        }

        // No chain with same name => Create one:
        if (mole_chain == (Chain*)NULL)
        {
            mole_chain = new Chain(molecule,former_chain.getName());
            molecule.Chains.push_back(mole_chain);
        }

        // Making the mapping:
        oldToNewChain.insert(pair<Chain*,Chain*>(&former_chain,mole_chain));

        // Final check, see if former chain has residues
        if (former_chain.getCountRes()==0)
        {
            const ItChain itc=find(Chains.begin(),Chains.end(),&former_chain);
            delete &former_chain;
            Chains.erase(itc);
        }
    }

    // Now updating residu chain
    for (size_t ipos = 0; ipos < residuList.size();++ipos)
    {
        Residu &residu=*residuList.at(ipos);
        Chain *chain = oldToNewChain.at(residu.chain);
        oldToNewChain.at(residu.chain)->residus.push_back(&residu);
        residu.chain=chain;

        Residues.erase(find(Residues.begin(),Residues.end(),&residu));
    }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// GENERICS //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////




/**
 * @brief Molecule::toString
 * @return Description of the whole molecule, including atom, chain, residu,
 * bond and cycles
 */
std::string Molecule::toString(const bool &onlyUse) const
{
  ostringstream oss;oss << "########## Molecule "<< this<< "##########"<<endl;

  for (ItCAtom it = Atoms.begin(); it != Atoms.end(); it++)
  {

      if ((onlyUse&& (*it)->isUsed())||!onlyUse)oss << (*it)->toString()<<endl;
  }
  for (ItCBond it = Bonds.begin(); it != Bonds.end(); it++)
      if ((onlyUse&& (*it)->isUsed())||!onlyUse)oss <<"Bond-"<< (*it)->toString()<<endl;
  for (ItCRes it = Residues.begin(); it != Residues.end(); it++)
     if ((onlyUse&& (*it)->isUsed())||!onlyUse) oss << (*it)->toString()<<endl;
  for (ItCCycle it = Cycles.begin(); it != Cycles.end(); it++)
      oss << (*it)->toString()<<endl;
  return oss.str();
}


void Molecule::setMoleType(const unsigned int & type)
{

    if (complex != (Complex*)NULL) complex->updateMoleType(this,type);
    typemol=type;
}







std::string Molecule::getLongMoletype() const
{
  switch (typemol)
    {
    case MoleType::COFACTOR:   return "COFACTOR";
    case MoleType::LIGAND:     return "LIGAND";
    case MoleType::OTHER :     return "OTHER";
    case MoleType::PROTEIN :   return "PROTEIN";
    case MoleType::WATER :     return "WATER";
    case MoleType::SITE:       return "SITE";
    case MoleType::CAV_ALL:    return "CAVITY ALL";
    case MoleType::CAV_4:      return "CAVITY 4";
    case MoleType::CAV_6:      return "CAVITY 6";
    case MoleType::CAV_8:      return "CAVITY 8";
    case MoleType::CAV_12:     return "CAVITY 12";
    case MoleType::INT_LIG:    return "Ligand-Centered interaction";
    case MoleType::INT_CENT:   return "Center-Centered interaction";
    case MoleType::INT_PROT:   return "Protein-Centered interaction";
    case MoleType::INT_MERG:   return "Full interaction";
    case MoleType::NUCLEIC:    return "NUCLEIC";
    default:                   return "UNKNOWN";
    }
}



void Molecule::outplp() const{
    double plp = 0;
    for (ItCAtom it = Atoms.begin(); it != Atoms.end(); it++)
    {
        Atom &atom= **(it);
        plp += atom.getPartialCharge();

    }
    cout << "plp de la molecule : "  << plp << endl;

}

