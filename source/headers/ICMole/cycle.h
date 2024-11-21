#ifndef CYCLE_H
#define CYCLE_H


#include "atom.h"
namespace ICMole
{

  class Cycle
  {
    friend class Molecule;
  private:
    /*!< @brief atomlist : list of atoms within this cycle*/
    AtomList atomlist;
    /*!< \brief List of cycle included into this cycle */
    std::vector<Cycle> listCycle;
    /*!< @brief center : center of the atom. Saved in molecule::*/
    Atom    center;
    /*!< \brief true if the cycle is aromatic */
    bool  aromatic;
    /*!< \brief True when it's a minimal cycle */
    bool isMinimal;
    /*!< \brief True when barycenter no need to be updated*/
    bool center_ok;
    /*!< \brief How many sybcycle does it contains*/
    unsigned short level;
    /*!< \brief List of normal vector for all atoms involved in the cycle*/
    CoordList normal;
    /*!< \brief Normal vector of this cycle */
    Coords norm_vector;


  public:
   Cycle(Molecule *const parent);
   ~Cycle();

   bool operator<(const Cycle& rhs) const {

       if (rhs.atomlist.begin() < this->atomlist.begin() || rhs.atomlist.begin() == this->atomlist.begin() ){
//           std::cout << "inf" << std::endl;
           return true;
       }
//       std::cout << "sup" << std::endl;
       return false;
//       return rhs.atomlist.begin() < this->atomlist.begin() || (rhs.atomlist.begin() == this->atomlist.begin());
   }


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// ATOMS ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

   void addAtom( Atom*const atm);
   void addAtoms(const AtomList& atmList);
   void delAtom( Atom*const atm);

void uniqueAtoms();
   ItCAtom first()const { return atomlist.begin();}
   ItCAtom end()  const { return atomlist.end();}

   inline size_t getNumAtom() const {return atomlist.size();}

   bool hasAtom(const Atom &atm) const;
   Atom* getAtom(const size_t& pos)const {return atomlist.at(pos);}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// COORDS ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

   Atom& getCenter()  {return center;}
   Coords &getFixpos() ;
   Coords &getRotpos();


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// PROPERTIES //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



   inline bool isAromatic() const {return aromatic;}
    void setAromatic(const bool& ar);
   std::string toString() const;


   void calcVector() ;
   const Coords& getNormVector() const {return norm_vector;}
   const Coords& getNormVector(const size_t& pos) {
      if (pos >= normal.size()) {
         throw std::out_of_range("Invalid position in normal vector");
      }
      return normal.at(pos);
   }



  };

}

#endif // CYCLE_H
