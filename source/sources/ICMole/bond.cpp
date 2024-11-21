#include <iostream>
#include "headers/ICMole/bond.h"
#include "headers/ICMole/atom.h"
#include "headers/ICMole/residu.h"
using namespace std;
using namespace ICMole;


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////////////////////////////// CONSTRUCTORS /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*!
   *  \brief Standard constructor
   *  Can only be called by Molecule::addBond function
   *  \param V1 : First Atom of the Bond.
   *  \param V2 : Second Atom of the Bond
   *  \param atm : If the Bond represent an atom -> address of the atom
   *  \param num : Id of the Bond, as defined by the program. Default : 0 (depreciated 0 value)
   *  \param weight : For weighted Molecule, weight of the Bond. Default : 1
   */
Bond:: Bond(Atom  &nAtom1,
            Atom  &nAtom2,
            const unsigned int& num, const unsigned int &BType,
            Molecule& parent)
            :num(num),
             type(BType),
             atom1((nAtom1)),
             atom2(nAtom2),
             molecule(&parent),
             bits(""),
             inUse(true)
{
    atom1.links.push_back(this);atom1.atomlinked.push_back(&atom2);
    atom2.links.push_back(this);atom2.atomlinked.push_back(&atom1);
}

/*!
   *  \brief Destructor
   *
   *  Will delete Bond pointer from each Atom.
   */
Bond::~Bond() throw(MoleExcept)
{
    try
    {
        ItBond ited = find(atom1.links.begin(),atom1.links.end(),this);
        if (ited == atom1.links.end())
            throw MoleExcept(1060101,
                             "Bond::~Bond",
                             "Given Bond is not part of this Atom");


        atom1.links.erase(ited);
        atom1.atomlinked.erase(atom1.atomlinked.begin()+
                               std::distance(atom1.links.begin(),ited));

        ited = find(atom2.links.begin(),atom2.links.end(),this);
        if (ited == atom2.links.end())
            throw MoleExcept(1060101,
                             "Bond::~Bond",
                             "Given Bond is not part of this Atom");
        atom2.atomlinked.erase(atom2.atomlinked.begin()+
                               std::distance(atom2.links.begin(),ited));
        atom2.links.erase(ited);

    }
    catch (MoleExcept &e)
    {
        Moleaccess=Levels::FATAL;
        throw;
    }
}




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/////////////////////////////// MISCELLANEOUS /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*!
   *  \fn Atom&  getOtherAtom(Atom const &Atom) const
   *  \brief Given a Atom involved in the Bond, return the other Atom involved
   *  \param Atom : Atom to look at
   *  \return Atom pointer of the other Atom
   *  \throw MoleExcept 1030101 when the given Atom is not part of this Bond
   */
Atom&  Bond::getOtherAtom(Atom const &atom) const throw(MoleExcept)
{
         if (&atom == &atom1) return atom2;
    else if (&atom == &atom2) return atom1;
    throw MoleExcept(1060201,"Bond::getOtherAtom","Given Atom is not part of this Bond");
}


/*!
   *  \fn Atom& Bond::getOtherAtom(Atom *const Atom) const
   *  \brief Given a Atom involved in the Bond, return the other Atom involved
   *  \param Atom : Atom to look at
   *  \return Atom pointer
   *  \throw MoleExcept 1030102 when the given Atom is not part of this Bond
   */
Atom&  Bond::getOtherAtom(const Atom *const atom) const throw(MoleExcept)
{
         if (atom == &atom1) return atom2;
    else if (atom == &atom2) return atom1;
    throw MoleExcept(1060301,"Bond::getOtherAtom","Given Atom is not part of this Bond");
}

/*!
   * \fn Atom*  Bond::shareAtom(Bond const &ed) const
   * \brief Given another Bond, search the common Atom and return it
   * \param ed : Bond to look at
   * \return The Atom existing in the two Bonds or (Atom*)NULL
   */
Atom&   Bond::shareAtom(Bond const &bond) const throw(MoleExcept)
{
  if (&bond.getAtom1() == &atom1 || &bond.getAtom2()==&atom1) return atom1;
  if (&bond.getAtom1() == &atom2 || &bond.getAtom2()==&atom2) return atom2;
   throw MoleExcept(1060401,"Bond::shareAtom","No shared Atom");
}



/*!
   * \fn std::string Bond::toString() const
   *  \brief Returns a description of the Bond
   *
   */
const std::string Bond::toString(const Atom* atom) const
{
    ostringstream oss;


        oss << num << "/";
        oss.width(9);oss.setf(ios::left);
    switch (type)
    {
    case 101: oss << "SINGLE";   break;
    case 102: oss << "DOUBLE";   break;
    case 103: oss << "TRIPLE";   break;
    case 104: oss << "AROMATIC"; break;
    case 105: oss << "AMIDE";    break;
    case 106: oss << "ANY";      break;
    case 107: oss << "DUMMY";    break;
    case 108: oss << "UNDEFINED";break;
    case 109: oss << "LINKER";   break;
    case 110: oss << "SIDECHAIN";break;
    }
    oss << "/";
    if (atom != (Atom*)NULL)
    {
        if (atom == &atom1)
        {
            oss.setf(ios::left);oss<<atom1.getIdentifier()<<"/";
            oss.setf(ios::left);oss<<atom2.getIdentifier()<<"/";
        }
        else
        {
            oss.setf(ios::left);oss<<atom2.getIdentifier()<<"/";
            oss.setf(ios::left);oss<<atom1.getIdentifier()<<"/";
        }

    }
    else
    {
        oss.setf(ios::left);oss<<atom1.getIdentifier()<<"/";
        oss.setf(ios::left);oss<<atom2.getIdentifier()<<"/";
    }


    oss <<"BITS:"<<bits;
    if (!inUse) oss << " (Not Use)";
    return oss.str();
}

/*!
   * \fn void Bond::clone(Bond& bond)
   *  \brief Clone the Bond
   *
   *  Copy data of the given Bond into this one.
   * \warning It only copy Id, weight and label. It does not copy Atom pointer nor parent or bond related information.
   * \param ed : Bond to copy data from
   * \param with_parent : copy Bond links. It does not replace the initial Bond bond
   */
void Bond::clone(Bond const & bond)
{
    num=bond.num;
    type=bond.type;
}

/**
 * @brief Bond::setUse
 * @param newUse : True when must be kept for further processing
 * @param updateAtom : Update also both involved atom in the bond
 *
 * Say if the bond must be kept fur further processing.
 * will call Atom::setUse() to update its atoms when updateAtom is set to true.
 * The latter will also update residu use but not other bonds linked to the two
 * atoms.
 */
void  Bond::setUse(const bool& newUse, const bool& updateAtom)
{
    inUse=newUse;
    if (!updateAtom)return;
    atom1.setUse(newUse,true,false);
    atom2.setUse(newUse,true,false);
}


double Bond::getLength() const
{
 return atom1.calcFixpos(atom2);
}



void Bond::assignBits()
{
  bits="";
  string Name1, Name2;
  const string DICT=" CB CD CD1 CD2 CE CE1 CE2 CE3 CG CG1 CG2 CH2 CZ CZ2 HD1 HD2 HD21 HD22 HE HE1 HE2 HE21 HE22 HG HG1 HH HH11 HH12 HH21 HH22 HN1 HN2 HN3 HZ1 HZ2 HZ3 LPD1 LPD2 LPG1 LPG2 ";
  if (atom1.getName().compare(atom2.getName()) < 0)
    {
      Name1=atom1.getName(); Name2=atom2.getName();
    }
  else {Name1=atom2.getName(); Name1=atom1.getName();}

  if (Name1=="C")
    {
      if (Name2 == "CA" || Name2=="O")bits="BACKBONE|DICT";
      else if (Name2 == "N") bits="BACKBONE|DICT|INTERRES";
      else if (Name2== "OXT")bits="CAP|DICT";
    }
  else if (Name1 == "CA")
    {
      if (Name2 == "CB"||Name2=="HA3")bits="DICT";
      else if (Name2 == "HA"||Name2=="HA2")bits="BACKBONE|DICT";
      else if (Name2=="O") bits="BACKBONE|DICT|INTERRES";
      else if (Name2 == "N")
        {
          if (atom2.getResidu()->getName() =="AMN") bits="BACKBONE|DICT|INTERRES";
          else bits="BACKBONE|DICT";
        }
    }
  else if (Name1 == "H" && Name2 == "N")bits="BACKBONE|DICT";
  else if (Name1 == "HNCAP" && Name2 == "N")bits="CAP|DICT";
  else if (Name1 == "HOCAP" && Name2 == "OXT")bits="CAP|DICT";
  else if (Name1 == "N" && Name2 == "O") bits="BACKBONE|DICT|INTERRES";
  else if (DICT.find(" "+Name1+" ") != string::npos) bits="DICT";
}
