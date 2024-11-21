#include "headers/ICMole/cycle.h"
#include "headers/ICMole/molecule.h"
using namespace std;
using namespace ICMole;


Cycle::Cycle(Molecule *const parent):
    center("Du",parent,&parent->getAromaticRes(),"DuAr","Du"),
    aromatic(false),isMinimal(true),center_ok(false),level(0)
{

}

Cycle::~Cycle()
{

}



void Cycle::uniqueAtoms()
{
    sortAndUnique(atomlist);
}

/**
 * @brief Cycle::addAtom
 * @param atom : Atom to add
 *
 * Add this atom to the list of atom of this cycle.
 * This imply that the barycenter of the cycle is no longer up to date
 * so it will be calculated when calling to getFixpos(), getRotpos() or
 * getCenter() function
 *
 */
void Cycle::addAtom( Atom*const atom)
{
    if (atom == (Atom*)NULL)return;
    center_ok=false;
    atomlist.push_back(atom);
}


void Cycle::addAtoms( const AtomList& atmList)
{
    atomlist=atmList;
}

/**
 * @brief Cycle::delAtom
 * @param atom : atom to remove
 *
 * Remove the given atom from the atom list
 * This imply that the barycenter of the cycle is no longer up to date
 * so it will be calculated when calling to getFixpos(), getRotpos() or
 * getCenter() function
 */

void Cycle::delAtom(Atom *const atom)
{
    if (atom == (Atom*)NULL) return;
    ItAtom it=find(atomlist.begin(), atomlist.end(), atom);
    if (it == atomlist.end()) return;

    atomlist.erase(it);
    center_ok=false;
}




/**
 * @brief Cycle::getFixpos
 * @return The barycenter of the cycle
 *
 * Check whether an atom has been added or deleted from this cycle.
 * If not, the barycenter will be calculated. Otherwise it will return
 * the current barycenter.
 *
 * @warning Checks are only made on addition/deletion of an atom within the cycle
 * When an atom is moved, no check will be made and therefore the barycenter
 * can be false
 */
Coords& Cycle::getFixpos()
{
    if (!center_ok)
    {
        center.fixpos.clear();
        center.rotpos.clear();
        for (ItCAtom itA = atomlist.begin();itA != atomlist.end();itA++)

        {
            center.fixpos+=(*itA)->fixpos;
            center.rotpos+=(*itA)->rotpos;
        }
        if(atomlist.size() > 0) {
            center.fixpos/= (double)atomlist.size();
            center.rotpos/= (double)atomlist.size();
        }
        
        center_ok=true;
    }
    return center.fixpos;
}


/**
 * @brief Cycle::getRotpos
 * @return The barycenter of the cycle in the new axis
 *
 * Check whether an atom has been added or deleted from this cycle.
 * If not, the barycenter will be calculated. Otherwise it will return
 * the current barycenter.
 *
 * @warning Checks are only made on addition/deletion of an atom within the cycle
 * When an atom is moved, no check will be made and therefore the barycenter
 * can be false
 *
 */
Coords& Cycle::getRotpos()
{
    if (!center_ok)
    {
        center.fixpos.clear();
        center.rotpos.clear();
        for (ItCAtom itA = atomlist.begin();itA != atomlist.end();itA++)

        {
            center.fixpos+=(*itA)->fixpos;
            center.rotpos+=(*itA)->rotpos;
        }
        center.fixpos/= (double)atomlist.size();
        center.rotpos/= (double)atomlist.size();
        center_ok=true;
    }
    return center.rotpos;
}


/**
 * @brief Cycle::hasAtom
 * @param atom : Atom to check
 * @return Bool telling whether the given atom is within the cycle or not
 *
 * Check if the given atom is within the cycle or not. Return TRUE if so,
 * FALSE otherwise
 */
bool Cycle::hasAtom(const Atom &atom)const
{
    ItCAtom ita=find(atomlist.begin(),atomlist.end(), &atom);
    if (ita != atomlist.end()) return true; return false;
}



 void Cycle::setAromatic(const bool& ar)
 {
     aromatic=ar;
     for (ItAtom itA = atomlist.begin(); itA != atomlist.end();itA++)
     {
         (*itA)->props.setAromatic(ar);
     }
 }

/**
 * @brief Cycle::toString
 * @return A string describing the cycle.
 */
std::string Cycle::toString() const
{
    ostringstream oss;
    oss << "---- CYCLE ----- ";
    if (aromatic == true) oss << "AROMATIC";
    oss<<endl;
    for (ItCAtom itA = atomlist.begin();itA != atomlist.end();itA++)
        oss << " |-> "<<(*itA)->getIdentifier()
            << " : "<< (*itA)->props.toString() <<endl;

    oss <<endl;
    return oss.str();
}


void Cycle::calcVector()
{
    normal.clear();

    if (atomlist.size() < 2) {
        std::cerr << "Error: not enough atoms in the cycle to calculate vector" << std::endl;
        return;
    }

    for (ItAtom itA = atomlist.begin(); itA != atomlist.end();itA++)
      {
         Atom &atm = **itA;
        if (atm.getNumBond() < 2)
        {
            normal.push_back(Coords(0,0,0));
            continue;
        }
        const Atom &atm1 = atm.getAtomLinked(0);
        const Atom &atm2 = atm.getAtomLinked(1);
        normal.push_back(Coords(atm.fixpos.getNormal(atm1.fixpos,atm2.fixpos)));
      }
    const Atom&atm1 = *atomlist.at(0);
    const Atom&atm2 = *atomlist.at(1);
    norm_vector.setCoords(getFixpos().getNormal(atm1.fixpos,atm2.fixpos));

}
