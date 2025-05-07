#include <iterator>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "headers/ICMole/atom.h"
#include "headers/ICMole/molecule.h"
#include "headers/ICMole/bond.h"
#include "headers/ICMole/box.h"

/*!
\todo Add (Atom*)NULL as initializer of Atom
\todo Check in which case we need to look at the first letter of atomicName
*/

using namespace ICMole;
using namespace std;

bool Atom::prep_loaded=false;
std::vector<std::string> Atom::Hyd_Db;

std::vector<std::string> Atom::Hyd_DbS;

std::vector<std::string> Atom::A_Atm;

std::vector<std::string> Atom::Met_list;    
void Atom::loadData()
{

    Hyd_Db.push_back("N.4");
    Hyd_Db.push_back("N.2");
    Hyd_Db.push_back("O.2");
    Hyd_Db.push_back("N.ar");
    Hyd_Db.push_back("N.1");
    Hyd_Db.push_back("N.am");
//    Hyd_Db.push_back("O.3");
//    Hyd_Db.push_back("O.co2");
//    Hyd_Db.push_back("O.1");
    Hyd_DbS.push_back("N.am");
    Hyd_DbS.push_back("O.3");
    Hyd_DbS.push_back("O.2");
    Hyd_DbS.push_back("N.2");
    Hyd_DbS.push_back("N.3");
    Hyd_DbS.push_back("N.pl3");

    A_Atm.push_back("N.pl3");
    A_Atm.push_back("N.am");
    A_Atm.push_back("N.4");
    Met_list.push_back("Ca");
    Met_list.push_back("Cd");
    Met_list.push_back("Co");
    Met_list.push_back("Cu");
    Met_list.push_back("Fe");
    Met_list.push_back("Mg");
    Met_list.push_back("Mn");
    Met_list.push_back("Ni");
    Met_list.push_back("Zn");
    Met_list.push_back("Na");
    Met_list.push_back("CA");
    Met_list.push_back("CD");
    Met_list.push_back("CO");
    Met_list.push_back("CU");
    Met_list.push_back("FE");
    Met_list.push_back("MG");
    Met_list.push_back("MN");
    Met_list.push_back("NI");
    Met_list.push_back("ZN");
    Met_list.push_back("NA");
    prep_loaded=true;
//    plp_Alist.push_back("D");
//    plp_Alist.push_back("Neg");
//    plp_Dlist.push_back("A");
//    plp_Alist.push_back("Pos");
//    plp_Hlist.push_back("H");
//    plp_Hlist.push_back("Ar");
//    plp_Mlist.push_back("A");
//    plp_Mlist.push_back("Neg");
//    plp_Mlist.push_back("AD");

}







////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// CONSTRUCTORS /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

Atom::Atom():num(0),           fNum(0),                atomicNum(0) ,
    partial_charge(0),formal_charge(0),       BFactor(0) ,     data_loaded(false),
    inUse(true),inCycle(false),
    name(""),       atomicName(""), mol2type(""), bits(""),
    molecule((Molecule*)NULL), residu((Residu*)NULL),            fixpos(0,0,0),   rotpos(0,0,0)
{
    if (!prep_loaded) loadData();
    links.reserve(6);
    atomlinked.reserve(6);
#ifdef GUI_DESIGN
    color[0] = 0;color[1] =0;color[2] =0;
    visible=true;
    showtype=ShowType::SMALLSPHERE;
#endif
}

Atom::Atom(const std::string& atomicName,
           Molecule *const parent,
           Residu *const res,
           const std::string& name,
           const std::string& mol2t) throw(MoleExcept)
    :num(0),           fNum(0),                atomicNum(0) ,
      partial_charge(0),formal_charge(0),       BFactor(0) ,     data_loaded(false),
      inUse(true),inCycle(false),
      name(name),       atomicName(atomicName), mol2type(mol2t), bits(""),
      molecule(parent), residu(res),            fixpos(0,0,0),   rotpos(0,0,0)
{
    if (!prep_loaded) loadData();
    if (parent == (Molecule*)NULL)
        throw MoleExcept(1050101,
                         "Atom::Atom",
                         "No parent molecule given");
    if (res == (Residu*)NULL)
        throw MoleExcept(1050102,
                         "Atom::Atom",
                         "No residu given");
    try
    {
        links.reserve(6);
        atomlinked.reserve(6);

#ifdef GUI_DESIGN
        color[0] = 0;color[1] =0;color[2] =0;
        visible=true;
        showtype=ShowType::SMALLSPHERE;
#endif

        if (!atomicName.empty()) loadAtomicData();
    }
    catch (MoleExcept &e)
    {
        e.addTrace("Atom::Atom");

    }
}



Atom::Atom(const std::string& atomicName,
           Molecule *const parent,
           Residu *const res,
           const double& x,
           const double& y,
           const double& z,
           const std::string& name,
           const std::string& mol2t) throw(MoleExcept)
    :num(0),           fNum(0),                atomicNum(0) ,
      partial_charge(0),formal_charge(0),       BFactor(0) ,     data_loaded(false),
      inUse(true),inCycle(false),
      name(name),       atomicName(atomicName), mol2type(mol2t), bits(""),
      molecule(parent), residu(res),            fixpos(x,y,z),   rotpos(0,0,0)
{
    if (!prep_loaded) loadData();
    if (parent == (Molecule*)NULL)
        throw MoleExcept(1050201,
                         "Atom::CONSTRUCTOR",
                         "No parent molecule given");
    if (res == (Residu*)NULL)
        throw MoleExcept(1050202,
                         "Atom::Atom",
                         "No residu given");
    try
    {
        links.reserve(6);
        atomlinked.reserve(6);

#ifdef GUI_DESIGN
        color[0] = 0;color[1] =0;color[2] =0;
        visible=true;
        showtype=ShowType::SMALLSPHERE;
#endif

        if (!atomicName.empty()) loadAtomicData();

    }
    catch (MoleExcept &e)
    {
        e.addTrace("Atom::Atom");
    }

}


Atom::Atom( const std::string& atomicName,
            Molecule *const parent,
            Residu *const res,
            const Coords &coords,
            const std::string& name,
            const std::string& mol2t)  throw(MoleExcept)
    :num(0),           fNum(0),                atomicNum(0) ,
      partial_charge(0),formal_charge(0),       BFactor(0) ,     data_loaded(false),
      inUse(true),inCycle(false),
      name(name),       atomicName(atomicName), mol2type(mol2t), bits(""),
      molecule(parent), residu(res),            fixpos(coords),   rotpos(0,0,0)

{
    if (!prep_loaded) loadData();
    if (parent == (Molecule*)NULL)
        throw MoleExcept(1050301,
                         "Atom::CONSTRUCTOR",
                         "No parent molecule given");
    if (res == (Residu*)NULL)
        throw MoleExcept(1050302,
                         "Atom::Atom",
                         "No residu given");
    try
    {
        links.reserve(6);
        atomlinked.reserve(6);

#ifdef GUI_DESIGN
        color[0] = 0;color[1] =0;color[2] =0;
        visible=true;
        showtype=ShowType::SMALLSPHERE;
#endif

        if (!atomicName.empty()) loadAtomicData();

    }
    catch (MoleExcept &e)
    {
        e.addTrace("Atom::Atom");

    }

}



Atom::~Atom()
{

    for (std::map<const Grid*, Box*>::iterator it= listBox.begin();it != listBox.end();++it)
    {

        Box *bx= (*it).second;

        bx->delAtom(*this);
    }
}






void  Atom::setUse          (const bool& newUse,
                             const bool& updateRes,
                             const bool& toBonds)
{
    inUse=newUse;
    if (toBonds)
    {
        for (ItBond itB = links.begin();itB != links.end(); itB++)
        {
            (*itB)->setUse(newUse,false);
        }
    }
    if (updateRes && residu != (Residu*)NULL) residu->checkUse();
}












////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// BondS ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



/*! \fn void  Atom::addBond     (const Bond  * const ed)
  * \brief Adding an Bond to the Atom
  * \param ed : Bond to add
  */
void  Atom::addBond    ( Bond  *const  ed)
{

    links.push_back(ed);
    atomlinked.push_back(&ed->getOtherAtom(this));
}





/*! \fn void  Atom::delBond    (const  Bond  *const ed )
  * \brief Delete an Bond to the Atom Bond list
  * \param ed : Bond to remove
  * Check whether the Bond is in the Bond list of the Atom and delete it if found
  * \warning It doesn't delete the Bond, just the link in the Atom.
  * \throw MoleExcept 1020101 when the given Bond is not part of this Atom
  */
void  Atom::delBond    (const  Bond  *const ed ) throw(MoleExcept)
{
    const ItBond ited = find(links.begin(),links.end(),ed);
    if (ited == links.end())
        throw MoleExcept(1050101,
                         "Atom::delBond",
                         "Given Bond is not part of this Atom");
    atomlinked.erase(atomlinked.begin()+
                     std::distance(links.begin(),ited));
    links.erase(ited);
}





/*! \fn bool  Atom::hasBondWith(Atom const &ve)
  * \brief Check if the given Atom shares an Bond with this Atom
  * \param ve : Atom to check
  * \return true when an Bond exists. False otherwise
  */
bool Atom::hasBondWith(const Atom &atom) const
{
    const size_t natms=atomlinked.size();
    for (size_t i=0;i < natms;i++) if (atomlinked.at(i) == &atom) return true;
    return false;
}







/*! \fn const Bond* Atom::getBondWith(Atom const &atom)
  * \brief Check if the given Atom shares an Bond with this
  * Atom, and returns the corresponding Bond
  * \param ve : Atom to check
  * \return Pointer of the Bond betwwen ve and this Atom,
  *  or (Bond*)NULL when no Bond exists.
  * \throw 1060201 - Bond::getOtherAtom
  */
const Bond* Atom::getBondWith(Atom const &atom) const throw(MoleExcept)
{
    try
    {
        const Bond *bond;
        for (ItCBond it = links.begin(); it != links.end(); it++)
        {
            bond = *it;
            if (&(bond->getOtherAtom(this)) == &atom) return bond;
        }
        return (Bond*)NULL;
    }catch (MoleExcept &e)
    {
        //1060201 - getOtherAtom
        e.addTrace("Atom::getBondWith");
        throw;
    }
}



Bond*  Atom::getBond(const size_t& pos)const  throw(MoleExcept)

{
    if (pos >= links.size())
        throw MoleExcept(000000,
                         "Atom::getBond",
                         "pos above the maximum number of bonds");
    return links.at(pos);
}



/*! \fn void Atom::cleanBond() throw(MoleExcept)
  * \brief Delete all Bonds of this Atom
  * Calls delBond() function of the parent Molecule to make a clean
  * deletion of all Bonds made with this Atom
  * \throw See Molecule::delBond(). If 1040201: Fatal error
  */
void Atom::cleanBond() throw(MoleExcept)
{

    try
    {
        molecule->delBond(links);
    }
    catch (MoleExcept &e)
    {
        e.addTrace("Atom::cleanBond()\nFrom : "+toString());
    }

    links.clear();

}




const Atom& Atom::getAtomLinked(const size_t& pos) const throw(MoleExcept)
{
    if (pos> atomlinked.size())
        throw MoleExcept(000000,
                         "Atom::getAtomLinked",
                         "pos is above the number of atom linked");
    return *atomlinked.at(pos);

}






////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// ATOMIC DATA /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/*! \fn void Atom::loadAtomicData() throw(MoleExcept)
 *  \brief Scan over AtomData struct a corresponding atomic Name
 *  Using the atomic name of the atom,
 * it will search over the AtomData struct the corresponding name to seek
 * for atomic number and to check whether the atom is a
 * biologically relevant atom
 *
 * \throw MoleExcept 1050201 - The given atom is not biologically relevant.
 * \throw MoleExcept 1050202 - Unrecognized atomic name

*/
void Atom::loadAtomicData()                         throw(MoleExcept)
{


    // FIRST SCAN - NO MODIFICATION OF atomicName :
    for (int Ivar1 =0; Ivar1<NB_ATMS; Ivar1++)
        if (AtomData[Ivar1].nom == atomicName)
        {
            atomicNum = AtomData[Ivar1].num;
            if (!AtomData[Ivar1].biol)
                throw MoleExcept(1050201,
                                 "Atom::loadAtomicData",
                                 getIdentifier()+
                                 " is not a biologically relevant atom");
            data_loaded=true;
            return;
        }



    // SECOND SCAN - MODIFICATION OF atomicName :
    // example : FE into Fe, MG into Mg
    string UpperLower="";

    for (size_t i=0; i<atomicName.length();i++)
    {
        if (i==0) UpperLower+=toupper(atomicName.at(i));
        else UpperLower+=tolower(atomicName.at(i));
    }

    for (int Ivar1 =0; Ivar1<NB_ATMS; Ivar1++)
        if (AtomData[Ivar1].nom == UpperLower)
        {
            atomicNum = AtomData[Ivar1].num;
            atomicName = UpperLower;
            if (!AtomData[Ivar1].biol)
                throw MoleExcept(1050201,
                                 "Atom::loadAtomicData",
                                 getIdentifier()+
                                 " is not a biologically relevant atom");
            data_loaded=true;
            return;
        }

    // LAST SCAN - SUBSTRING OF atomicName :
    // Use of the first letter => CAN CAUSE ISSUES

    for (int Ivar1 =0; Ivar1<NB_ATMS; Ivar1++)
        if (AtomData[Ivar1].nom == atomicName.substr(0,1))
        {
            atomicNum = AtomData[Ivar1].num;
            if (!AtomData[Ivar1].biol)
                throw MoleExcept(1050201,
                                 "Atom::loadAtomicData",
                                 getIdentifier()+
                                 " is not a biologically relevant atom");
            data_loaded=true;
            return;
        }
    throw MoleExcept(1050202,
                     "Atom::loadAtomicData",
                     "Unknown "+atomicName+
                     " for atomic name of atom : "+getIdentifier());

}







/** \fn void Atom::checkMOL2type() throw(MoleExcept)
  \brief Check the MOL2 type of the atom and update its properties
         according to MOL2 type and environment.
  \throw 1050301 - No MOL2 Type
  \throw 1050302 - MOL2 Type too long
  \throw 1050303 - MOL2 type given : Not recognized

  Check the MOL2_SPEC table to see if the type exists.
  Then it calls setAtomicNum() function to add atomic informations.
  It also clears atom physico-chemical properties. \n
  When the MOL2 type is aromatic, atom is set to aromatic.
  Also Hydrophobic if its a carbon and not linked to N.4, N.2, O.2 atoms
  by a double bond nor aromatic bond with O.co2.\n
  For pure hydrophobic, same as above.\n
  For metal, must be one of these values :
  Ca Cd Co Cu Fe Mg Mn Ni Zn CA CD CO CU FE MG MN NI ZN \n
  For acceptor, only when the formal charge is negative of
  if its either an oxygen, N.pl3, N.am or N.4\n
  For donor atom, only if its an oxygen, nitrogen or sulfur linked to an hydrogen.
  \note : called by setMOL2type
  */
void Atom::checkMOL2type() throw(MoleExcept)
{

    /* static const string Hyd_Db = " N.4 N.2 O.2 ";
    static const string A_Atm= " N.pl3 N.am N.4 ";
    static const string
=" Ca Cd Co Cu Fe Mg Mn Ni Zn CA CD CO CU FE MG MN NI ZN ";*/
    const size_t linkedatmSize=atomlinked.size();
    const size_t mol2typeLength=mol2type.length();
    bool wrong=false;
    if (mol2type.empty())
        throw MoleExcept(1050301,
                         "Atom::checkMOL2type",
                         "No MOL2 Type given");
    if (mol2type.length()>5)
        throw MoleExcept(1050302,
                         "Atom::checkMOL2type",
                         "MOL2 Type given is too long (max length:5)");

    if (mol2type=="Co")
    {
        mol2type="Co.oh";
    }

    // SCANNING MOL2 TABLE FOR THE GOOD MOL2TYPE:
    for (unsigned int Ivar1=0; Ivar1 < NB_MOL2;++Ivar1)
    {
        // Until we find the good one we continue:
        if (MOL2_SPEC[Ivar1].mol2type.find(mol2type)==string::npos ||
                MOL2_SPEC[Ivar1].mol2type.length() != mol2typeLength )continue;

        // Cleaning properties and setting atomicNum and Name
        props.clear();
        atomicNum=MOL2_SPEC[Ivar1].atomic_num;
        atomicName = AtomData[atomicNum-1].nom;

        // Aromatic : Must find Ar property
        if (MOL2_SPEC[Ivar1].property.find("Ar") != string::npos)
            props.setAromatic(true);

        // Hydrophobic : Must find H(ydrophobic) or Ar(omatic) Carbon
        if (MOL2_SPEC[Ivar1].property.find("H")  != string::npos
                || (MOL2_SPEC[Ivar1].property.find("Ar") != string::npos && atomicNum ==  6))
        {
            wrong=false;
            unsigned int NNpl3=0;
            for (size_t i=0; i<linkedatmSize; ++i) {

                const unsigned int &BType=links.at(i)->getBondType();
                const std::string &OthMOL2=atomlinked.at(i)->getMOL2Type();
                
                if (atomicNum == 6 ) {
                    if ((find(Hyd_Db.begin(),Hyd_Db.end(),OthMOL2) != Hyd_Db.end())
                    || (BType == BondType::AROMATIC && OthMOL2.compare("O.co2") == 0)) wrong=true;
                    
                    else if (BType == BondType::AROMATIC && OthMOL2.compare("N.pl3")==0)
                        NNpl3++;
                } 
                else {
                    if (find(Hyd_Db.begin(),Hyd_Db.end(),OthMOL2) != Hyd_Db.end()) {
                        wrong =true;
                    }
                }
            }
            if (!wrong && NNpl3 < 2)
                props.setHydrophobic(true);
        }
        // Apolar : must be np
        if (MOL2_SPEC[Ivar1].property.find("np")  !=string::npos)
            props.setApolar(true);

        // Cationic : must be PLUS
        if (MOL2_SPEC[Ivar1].property.find("PLUS")  !=string::npos)
            props.setCation(true);

        // Anionic : must be MOINS
        if (MOL2_SPEC[Ivar1].property.find("MOINS") !=string::npos)
            props.setAnion(true);

        // Dummy : must be NOP
        if (MOL2_SPEC[Ivar1].property.find("NOP")   !=string::npos)
            props.setDummy(true);

        // Metal : must be in the metallic above :
        if (find(Met_list.begin(),Met_list.end(),mol2type) != Met_list.end()) {
            props.setMetal(true);
        }

        // HBond acceptor : must be A(cceptor)
        // be in the above A_Atm list for Nitrogen, or be an oxygen or
        // have a formal charge below 0
        if (MOL2_SPEC[Ivar1].property.find("A")     !=string::npos)
        {
//            if (getParent().getMoleType()==MoleType::LIGAND) ;
            if(formal_charge<0){props.setAcceptor(true);}
            else if (formal_charge == 0 && atomicNum ==8 ) { props.setAcceptor(true);}
            else if(formal_charge==0
                    && (atomicNum ==  7  && find(A_Atm.begin(),A_Atm.end(),mol2type)== A_Atm.end())){
                //  && A_Atm.find(" "+mol2type+" ") == string::npos)))
                props.setAcceptor(true);
                for (size_t i=0; i<linkedatmSize;++i)
                    if (atomlinked.at(i)->getAtomicNum() == 1)
                    {props.setAcceptor(false);break;}
            }
        }

        // HBond donor : must be D(onor)
        // be in either oxygen, nitrogen or sulfur and have an hydrogen linked to
        if (MOL2_SPEC[Ivar1].property.find("D")      !=string::npos
                &&(atomicNum == 7 || atomicNum == 8 || atomicNum == 16 ))
        {


            for (size_t i=0; i<linkedatmSize;++i)
                if (atomlinked.at(i)->getAtomicNum() == 1)
                {props.setDonor(true);break;}
        }

        if (props.isAromatic())
        {
            for (size_t i=0; i<linkedatmSize;++i)
                if (atomicName =="N" && links.at(i)->getBondType() == BondType::AROMATIC)
                {props.setAcceptor(true);break;} // old weak acceptor
        }
        else
        {
            //            for (size_t i=0; i<linkedatmSize;++i)
            //               if (links.at(i)->getBondType() == BondType::DOUBLE
            //                 ||links.at(i)->getBondType() == BondType::TRIPLE
            //                       && atomicName !="C")
            //               {props.setAcceptor(true);break;}
            if ((atomicName=="S" || atomicName=="F" || atomicName=="Br" ||atomicName=="Cl" ||atomicName=="I") && formal_charge==0) props.setWeakAcceptor(true); // old weak acceptor
        }
        if (atomicName=="C")
        {
            bool wH=false;
            for (size_t i=0; i<linkedatmSize;++i)
                if (atomlinked.at(i)->getAtomicNum() == 1) {
                    wH=true;
                    break;
                }
            if ((props.isAromatic() || linkedatmSize==3 || linkedatmSize==2) && wH)
                props.setWeakDonor(true);
                
        }

        return;


    }
    throw MoleExcept(1050303,
                     "ATOM::Check_mol2type",
                     "MOL2 type given : "+mol2type+
                     " Not recognized for "+getIdentifier());

}








/*! \fn void Atom::setAtomicName(const string nAtm)
 *  \brief Load data according to the given atomic name
 *   \param nAtm : Atomic name to set the atom to
 *   When setting the atomic name for an atom,
 *   this function will search in AtomData table the corresponding id
 *   and update the radius, weight, atomic Num and van der Waals radius of the atom.
 * \throw MoleExcept 1050202 - No corresponding atom found
 * \throw MoleExcept 1050201 - The given atom is not biologically relevant.
 */
void Atom::setAtomicName(const string& nAtm)
{
    try
    {
        atomicName=nAtm;
        loadAtomicData();
    }
    catch (MoleExcept &e)
    {
        e.addTrace("Atom::setAtomicName");
        switch (e.getCode())
        {
        case 1050201: // Atom is not biologically relevant
            switch (ICMole::Moleaccess)
            {
            case Levels::STRICT: throw;break;
            case Levels::NOTICE: cerr << "NOTICE : "<< e.getData()<<endl;break;
            case Levels::WARNING: cerr<< "WARNING : "<< e.getData()<<endl;break;
            }break;
        case 1050202:// No corresponding atom found during loadAtomicData()
            throw;
        }
    }
}








void Atom::setAtomicNum (const unsigned int& atomicNumi)
{
    if (atomicNumi >= NB_ATMS
            || atomicNumi==0)
        throw MoleExcept(1050401,
                         "Atom::setAtomicNum",
                         "Given number does not correspond to an atomic number");
    atomicNum = atomicNumi;
    atomicName = AtomData[atomicNumi-1].nom;
}








/** \fn void Atom::setMOL2Type(const std::string&  MOL2)
     \brief Update the atom MOL2 type
     \param MOL2 : new MOL2 type to update to.
     Update the MOL2 type and call checkMOL2type() function.
     \throw 1050501 - If no MOL2 type given.
     \throw 1050302 - MOL2 Type too long
     \throw 1050303 - If unrecognized MOL2 type (from checkMol2Type
  */
void Atom::setMOL2Type  (const std::string&  MOL2) throw(MoleExcept)
{
    if (MOL2.empty())
        throw MoleExcept(1050501,
                         "Atom::setMOL2Type",
                         "No MOL2 Type given");
    mol2type=MOL2;
    try
    {
        checkMOL2type();
    }catch (MoleExcept &e)
    {
        e.addTrace("Atom::setMOL2Type");
        throw;
    }
}




const double& Atom::getWeight()
{
    if (!data_loaded) loadAtomicData();
    return AtomData[atomicNum-1].weight;

}

const double& Atom::getVdWRadius()
{
    if (!data_loaded) loadAtomicData();
    return AtomData[atomicNum-1].VDWradius;
}

const double& Atom::getRadius()
{
    if (!data_loaded) loadAtomicData();
    return AtomData[atomicNum-1].radius;
}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// BOOLEANS //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
bool Atom::isMetallic()   const
{ return (atomicNum == 12
          || atomicNum == 20
          || (atomicNum >= 25 && atomicNum <= 30)|| atomicNum == 48)
            ? true: false; }



bool Atom::isHalogen()    const
{ return (atomicNum ==  9
          || atomicNum == 17
          || atomicNum == 35
          || atomicNum == 53
          || atomicNum == 85)
            ? true:false;}


bool Atom::isNitrogen()   const {return (atomicNum ==  7)? true:false;}


bool Atom::isOxygen()     const {return (atomicNum ==  8)? true:false;}


bool Atom::isPhosphorus() const {return (atomicNum == 15)? true:false;}


bool Atom::isCarbon()     const {return (atomicNum ==  6)? true:false;}


bool Atom::isSulfur()     const {return (atomicNum == 16)? true:false;}


bool Atom::isHydrogen()   const {return (atomicNum == 1) ? true:false;}









const std::string& Atom::getResiduName() const
{
    if (residu == (Residu*)NULL)return "";
    return residu->getName();
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// MISCELLANEOUS ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/*! \fn std::string Atom::toString(const bool no_Bond) const
  * \brief Output verbose information about this Atom
  * \param no_Bond : Default TRUE. Also output linked Bond information
  */
const std::string Atom::toString(const bool &no_Bond, const bool &full_data) const
{

    std::ostringstream oss;
    oss << "----------- ATOM  "<< name<<"/"<<fNum<<"/"<<num;
    if (residu != (Residu*)NULL) oss << "/" << residu->getIdentifier();
    oss << " ----------"<<endl;
    oss << "BFactor"<<BFactor<<endl;
    oss<< "Charge : ";
    oss.width(5);oss.setf(ios::left); oss<<formal_charge;
    oss << " (Formal) / ";
    oss.width(5);oss<< partial_charge<<" (Partial)"<<endl;
    oss<< "Name   : ";
    oss.width(5); oss<<name;
    oss<< " (Name)   / ";
    oss.width(5); oss<< atomicName;
    oss<<" (Atomic)  / ";
    oss.width(5); oss<<mol2type;
    oss<< " (MOL2)"<<endl;
    oss<< "Coords : "<<fixpos.x<<";"<<fixpos.y<<";"<<fixpos.z<<" (Fixpos) / "
       <<rotpos.x<<";"<<rotpos.y<<";"<<rotpos.z<< " (Rotpos)"<<endl;
    if (!bits.empty())
        oss<< "Bits   : "<< bits <<endl;
    oss<< "Parent : "<<  molecule->getName()
       << "("<<molecule->getLongMoletype()<<")"<< molecule<<endl;


    if (atomicNum>=1 && full_data)
        oss << "Atomic Num : "        << atomicNum
            << "\t        Weight : "  << AtomData[atomicNum-1].weight
            <<"\tLong Name : "        << AtomData[atomicNum-1].longname<< endl;

    const std::string propstr=props.toString();
    if (propstr.length())oss << "Props  : " <<propstr<<endl;
    const std::string propFr=props.toFragType();
    if (propFr.length()) oss <<"FragType : "<< propFr<<endl;
    for (std::map<const Grid*, Box*>::const_iterator it=listBox.begin();it != listBox.end();++it)
    {
        cout << (*it).second->getId()<<"\t"<< (*it).second->fixpos.toString()<<"\t"<< (*it).second->fixpos.calcDist(fixpos)<<endl;
    }

    if (!no_Bond)
    {
        for (ItCBond it = links.begin() ; it != links.end(); ++it)
        {
            oss<< "  Bond : "<<(*it)->toString()<<endl;
        }
        for (size_t i=0; i< atomlinked.size();++i)
        {
            oss<<"Linked : "<< atomlinked.at(i)->getIdentifier()<<endl;
        }
    }
    oss<<endl;

    return oss.str();
}






/**
    \fn std::string Atom::getIdentifier()const
    \brief Return a small description of the atom, such as name, Id,
            and residu name, chain, and id.
    */
const std::string Atom::getIdentifier()const
{
    ostringstream oss;
    oss.width(5);oss.setf(ios::left); oss<< name;
    oss<<"|";
    oss.width(5); oss<< fNum;

    if (residu != (Residu*)NULL)
        oss << "|" << residu->getIdentifier();
    return oss.str();

}




/*! \fn void Atom::clone(Atom& Atom)
  * \brief Clone the data of the given Atom to this Atom
  * \param Atom : Atom to copy data from
  * \param with_linker : Default false. If true, also copy adress of related atom;
  * Clone Ids, names, bits, coordinates and all descriptors of the given Atom
  * to the current one. However, it does not copy  links.
  * Therefore no Bond will be created.
  */
void Atom::clone(const Atom& atom)
{
    num = atom.num;
    fNum = atom.fNum;
    atomicNum = atom.atomicNum;
    partial_charge = atom.partial_charge;
    formal_charge = atom.formal_charge;
    BFactor = atom.BFactor;
    name = atom.name;
    atomicName = atom.atomicName;
    mol2type = atom.mol2type;
    bits = atom.bits;
    fixpos = atom.fixpos;
    rotpos = atom.rotpos;
    props = atom.props;


}



double Atom::calcFixpos(const Atom& atom) const
{
    return atom.fixpos.calcDist(fixpos);
}

double Atom::calcFixpos(const Atom& atom, const double& threshold) const
{
    return atom.fixpos.calcDist(fixpos,threshold+0.1);
}


void Atom::addBox(const Grid * const grid,  Box* const box) throw(MoleExcept)
{
    if (grid == (Grid*)NULL)
        throw MoleExcept(1050701,"Atom::addBox","No grid given");
    if (box == (Box*)NULL)
        throw MoleExcept(1050703,"Atom::addBox","No Box Given");
    Grid* const gr2= const_cast<Grid* const>(grid);
    std::map<const Grid*, Box*>::iterator itp= listBox.find(gr2);
    if (itp != listBox.end())
    {
        if ((*itp).second == box)return;
        throw MoleExcept(1050703,
                         "Atom::addBox",
                         "A box from this grid is already associated to this atom");
    }

    listBox.insert(pair<Grid*,Box*>(gr2,box));
    box->addAtom(this);
}


Box* Atom::getBox(const Grid* const grid) const throw(MoleExcept)
{
    if (grid == (Grid*)NULL)
        throw MoleExcept(1050801,"Atom::getBox","No grid given");

    if (listBox.find(grid) == listBox.end())
        return (Box*)NULL;
    return listBox.at(grid);
}

void Atom::delBox(Box* const box)
{
    for (std::map<const Grid*, Box*>::iterator it= listBox.begin();it != listBox.end();++it)
    {
        if ((*it).second == box){
            listBox.erase(it);return;
        }
    }
}





void Atom::assignBits()
{
    const string DICT=" CB CD CD1 CD2 CE CE1 CE2 CE3 CG CG1 CG2 CH2 CZ CZ2 CZ3 HA3 HB HB1 HB2 HB3 HD11 HD12 HD13 HD23 HD3 HE3 HG11 HG12 HG13 HG2 HG21 HG22 HG23 HG3 HH2 HZ ND1 ND2 NE NE1 NE2 NH1 NH2 NZ OD1 OD2 OE1 OE2 OG OG1 OH SD SG ";
    const string DICT_ESSENTIAL= " HE HE21 HE22 HG1 HH HH11 HH12 HH21 HH22 HN1 HN2 HN3 HZ1 LPD1 LPD2 LPG1 LPG2 ";

    const string CAP_DICT = " HNCAP HOCAP OXT ";
    const string RName = residu->getName();
    if (name == "C"|| name=="CA" || name=="N"||name=="O")
    {
        if (residu->getResType()==ResType::STD_AA) bits="BACKBONE|DICT|DIRECT";
        else bits = "BACKBONE";
    }
    else if (name == "H")
    {
        if (RName!= "PHI") bits="BACKBONE|DICT|ESSENTIAL|DIRECT";
        else bits="BACKBONE";
    }
    else if (name == "HA")
    {
        if (RName != "PHI" &&RName!="HYP") bits="BACKBONE|DICT|DIRECT";
        else bits="BACKBONE";
    }
    else if (name =="HA2") bits="BACKBONE|DICT|DIRECT";
    else if (name == "HD1")
    {
        if (RName == "PHE" || RName == "TYR" || RName=="TRP") bits="DICT";
        else if (RName=="HIS" || RName=="HID")bits="DICT|ESSENTIAL";
    }
    else if (name =="HD2")
    {
        if (RName == "ASP") bits="DICT|ESSENTIAL";
        else bits="DICT";
    }
    else if (name == "HD21"|| name=="HD22")
    {
        if (RName == "ASN")bits="DICT|ESSENTIAL";
        else if (RName =="LEU")bits="DICT";
    }
    else if (name == "HE1"|| name=="HE2")
    {
        if (RName == "PHE"||RName=="MET"||RName=="TYR") bits="DICT";
        else if (name=="HE1")
        {
            if (RName == "HIS"||RName=="HID"||RName=="HIE")bits="DICT";
            else if (RName == "TRP") bits="DICT|ESSENTIAL";
        }
        else if (name=="HE2")
        {
            if (RName == "HIS"||RName=="HIE")bits="DICT|ESSENTIAL";
            else if (RName == "LYS") bits="DICT";

        }
    }
    else if ((name == "HZ2"|| name=="HZ3"))
    {
        if (RName =="LYS") bits ="DICT|ESSENTIAL";
        else if (RName=="TRP")bits="DICT";
    }
    else if (name=="HG")
    {
        if (RName == "SER"||RName=="CYS") bits="DICT|ESSENTIAL";
        else if (RName=="LEU") bits="DICT";
    }
    else if (CAP_DICT.find(" "+name+" ")!=string::npos)bits="CAP|DICT";
    else if (DICT_ESSENTIAL.find(" "+name+" ")!=string::npos)bits="DICT|ESSENTIAL";
    else if (DICT.find(" "+name+" ")!=string::npos)bits="DICT";




}


/**
\fn void Atom::setColor        (const double R,const double G, const double B)
\brief Set the color of the atom
\param R : Red component (value between 0 and 255)
\param G : Green component (value between 0 and 255)
\param B : Blue component (value between 0 and 255)
*/
void Atom::setColor        (const double &R,const double &G, const double &B)
{
    if (R > 255) color[0]=255; else if (R < 0) color[0]=0;else color[0]=R;
    if (G > 255) color[1]=255; else if (G < 0) color[1]=0;else color[1]=G;
    if (B > 255) color[2]=255; else if (B < 0) color[2]=0;else color[2]=B;

}

/**
\fn void Atom::setRcolor       (const double R)
\brief Set the red component color of the atom
\param R : Red component (value between 0 and 255)
*/
void Atom::setRcolor       (const double &R)      { if (R > 255) color[0]=255; else if (R < 0) color[0]=0;else color[0]=R;     }

/**
\fn void Atom::setGcolor       (const double G)
\brief Set the green component color of the atom
\param G : Green component (value between 0 and 255)
*/
void Atom::setGcolor       (const double &G)      { if (G > 255) color[1]=255; else if (G < 0) color[1]=0;else color[1]=G;      }

/**
\fn void Atom::setBcolor       (const double B)
\brief Set the blue component color of the atom
\param B : Blue component (value between 0 and 255)
*/
void Atom::setBcolor       (const double &B)      { if (B > 255) color[2]=255; else if (B < 0) color[2]=0;else color[2]=B;      }


//double getplp(Atom atom_prot) {
//    double dist_plp, plp_score ;
//    string mol2type;
//    if (mol2type == atom_prot.getMOL2Type())
//    {
////        dist_plp = fixpos.calcDist(atom_prot.fixpos);
//        if (dist_plp < 2.3){
////            plp_socre = F * (A-X) / A ;

//        } else if (dist_plp) {

//        }
//    }
//}

