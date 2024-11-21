#ifndef BOND_H
#define BOND_H

#include "global.h"
#include "physprop.h"
namespace ICMole
{


class Bond
{

    friend class Molecule;
private:
   /**
   * @brief num :  Id of the Bond, as defined by the program
   */
  unsigned int  num;

  /**
   * @brief Type of the bond. See BondType enumeration for the possibilities
   */
  unsigned int	type;


  /**
   * @brief First Atom involved in the Bond
   */
  Atom&         atom1;

  /**
   * @brief Second Atom involved in the Bond
   */
  Atom&         atom2;

  /**
   * @brief Parent of the Bond
   */
  Molecule*     molecule;

  /**
   * @brief Sybyl bits
   */
  std::string   bits;

  /**
   * @brief Tell whether this atom should be used for analysis or not
   */
  bool          inUse;



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// CONSTRUCTOR /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

      Bond(Bond const &);


      /*!< \brief Copy constructor private -> Unavailable */
      Bond& operator=(Bond const &);

      /*!< \brief Standard constructor. Can only be called by a molecule */
      Bond(Atom  &atom1,
           Atom  &atom2,
           const unsigned int& num,
           const unsigned int &BType,
           Molecule&  parent);


public:
      /*!< \brief Properties of the bond */
      PhysProp         props;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// DESCTRUCTOR /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
            ~Bond() throw(MoleExcept);









////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// IDs /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

///// Getters:

      /*!< \brief Return the Id of the Bond */
const unsigned int& getNum()      const {return num;    }


    /*!< \brief Return the type of the Bond */
const unsigned int& getBondType() const {return type; }


    /*!< \brief Return the parent Molecule */
          Molecule& getMolecule() const {return *molecule; }

          /*!< \brief Tell if this atom is used for processing */
      inline const         bool&  isUsed()      const       { return inUse;}

const std::string& getBits() const{ return bits;}


///// Setters

    /*!< \brief Set the Id of the Bond (as given by the Molecule) */
          void  setNum     (const unsigned int& n)  {   num = n; }


     /*!< \brief Set the weight of the Bond */
          void  setBondType(const unsigned int& n)  {type = n; }


    /*!< \brief Set the parent Molecule of the Bond */
          void  setMolecule ( Molecule  *const    parent)   { molecule = parent;}

    /*!< \brief Set Sybyl bond bits */
          void  setBits(const std::string &bit)    { bits=bit;}


          void  setUse(const bool& newUse, const bool& updateAtom=false);








////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// ATOMS ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

///// Getters :

    /*!< \brief Return the first Atom involved in the bond */
    Atom& getAtom1()    const {return atom1;}


    /*!< \brief Return the second Atom involved in the bond */
    Atom& getAtom2()    const {return atom2;}









////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// GENERICS ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/*!< \brief Return the other atom of the given one involved in this bond */
    Atom& getOtherAtom(Atom  const &atom) const throw(MoleExcept);


/*!< \brief Return the other atom of the given one involved in this bond */
    Atom& getOtherAtom(const Atom * const atom) const throw(MoleExcept);


/*!< \brief Return the atom that share this bond and the given on */
    Atom& shareAtom   (Bond    const &bond  ) const  throw(MoleExcept);


/*!< \brief Return a description of this bond : atoms involved and bond type */
    const std::string  toString(const Atom *atom=(Atom*)NULL) const;


/*!< \brief Copy this bond */
    void  clone(Bond const & bond);

     double getLength() const;

     void assignBits();
};


}
#endif // BOND_H
