#ifndef ATOM_H
#define ATOM_H

#include "global.h"
#include "coords.h"
#include "physprop.h"

namespace ICMole
{

class Atom
{
    friend class Bond;
    friend class Molecule;
    friend class MoleReader;
    friend class Cycle;

    /*!< \brief Static debugger */
#ifdef ICHEM_DEBUG
    static unsigned int debug;
#endif
private:
    static bool prep_loaded;

    static std::vector<std::string> Hyd_Db;

    static std::vector<std::string> Hyd_DbS;

    static std::vector<std::string> A_Atm;

    static std::vector<std::string> Met_list;

    //plp list
    static std::vector<std::string> plp_Alist;

    static std::vector<std::string> plp_Dlist;

    static std::vector<std::string> plp_Hlist;

    static std::vector<std::string> plp_Mlist;


    static void loadData();

    /*!< \brief  Number as assigned by the program */
    unsigned int        num;
    /*!< \brief Number as assigned by the input file */
    unsigned int        fNum;
    /*!< \brief Number as assigned by the interaction */
    unsigned int        tNum;
    /*!< \brief Atomic number of the atom : H => 1 ; He => 2 ... */
    unsigned int        atomicNum;
    /*!< \brief Partial charge of the atom */
    double              partial_charge;
    /*!< \brief Format charge of the atom */
    double              formal_charge;
    /*!< \brief BFactor of the atom (for PDB only) */
    double              BFactor;
    /*!< \brief True when loadAtomicData has been loaded and the atom checked */
    bool                data_loaded;
    /*!< \brief Tell whether this atom should be used for analysis or not */
    bool                inUse;
    /*!< \brief Tell whether this atom is involved in a cycle */
    bool                inCycle;
    /*!< \brief Name (as given in the file) */
    std::string         name;
    /*!< \brief Atomic name (standard Mendeleveiv notation : C, N, O) */
    std::string         atomicName;
    /*!< \brief SYBYL MOL2 atom type */
    std::string         mol2type;
    /*!< \brief SYBYL MOL2 bits*/
    std::string         bits;

    /*!< \brief  Parent Molecule*/
    Molecule*           molecule;
    /*!< \brief  Parent Residu */
    Residu*             residu;
    /*!< \brief  Parent box/ depending on the grid */
    std::map<const Grid*, Box*> listBox;
    /*!< \brief  List of Bonds */
    BondList            links;
    /*!< \brief List of linked atoms */
    AtomList            atomlinked;

    unsigned int        color[3];                       /*!< \brief RGB color : default coming from AtomData struct */
    ////////////////////////////////////////////////
    ////////////////////////////////////////////////
    ///////////////// CONSTRUCTORS /////////////////
    ////////////////////////////////////////////////
    ////////////////////////////////////////////////

    Atom();

    Atom( const std::string& atomicName,
           Molecule *const parent,
           Residu *const res,
           const std::string& name,
           const std::string& mol2t);

    Atom(const std::string& atomicName,
         Molecule *const parent,
         Residu *const res,
         const double& x,
         const double& y,
         const double& z,
         const std::string& name,
         const std::string& mol2t);

    Atom( const std::string& atomicName,
                Molecule *const parent,
                Residu *const res,
                const Coords &coords,
                const std::string& name,
                const std::string& mol2t);

     Atom( Atom const &);// No definition
     Atom& operator=(Atom const &);// No definition

    ~Atom();


#ifdef GUI_DESIGN
     /*!< \brief RGB color : default coming from AtomData struct */
     unsigned int        color[3];
     /*!< \brief GUI - Showing type - See enum GUIShowType */
     unsigned int        showtype;
     /*!< \brief GUI - True if visible */
     bool visible;
#endif

     ////////////////////////////////////////////////
     ////////////////////////////////////////////////
     /////////////// PRIVATE Functions //////////////
     ////////////////////////////////////////////////
     ////////////////////////////////////////////////

     void  addBond(Bond *const ed);
     void  delBond(const Bond *const ed);


     void  setParent(  Molecule  *const parent  );


     void setResidu(  Residu* const residu);

     void loadAtomicData();


public :

     Coords fixpos;
     Coords rotpos;
     PhysProp props;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// IDs /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

///// Getters:

/**
 * @brief Return the id of the atom, as given by the program.
 * @return atom ID
 */
inline const unsigned int& getNum()        const       { return num;      }


/**
 * @brief Return the parent molecule
 */
inline          Molecule& getParent()     const       { return *molecule;}


     /**
      * @brief Return the name of the atom, as given by the input file
      */
inline const std::string& getName()       const       { return name;     }


/*!< \brief Return the id of the Atom, as given by the file*/
inline const unsigned int& getFNum()      const       { return fNum;     }

/*!< \brief Return the id of the interaction coresponding to the point*/
inline const unsigned int& getTNum()      const       { return tNum;     }


    /*!< \brief Return the BFactor of the Atom, as given by the PDB file*/
inline const       double& getBFactor()   const       { return BFactor;  }


    /*!< \brief Return the atomic num of the atom */
inline const unsigned int& getAtomicNum() const       { return atomicNum;     }


    /*!< \brief Return the atomic name of the atom */
inline const  std::string& getAtomicName()const       { return atomicName;    }


     /*!< \brief Return the MOL2 type of the atom */
inline const  std::string&  getMOL2Type() const       { return mol2type;      }

    /*!< \brief Tell if this atom is used for processing */
inline const         bool&  isUsed()      const       { return inUse;}

inline const      double&  getPartialCharge() const   { return partial_charge;}
inline const      double&  getFormalCharge() const   { return formal_charge;}

inline const    std::string& getBits() const {return bits;}
       const double& getWeight() ;
       const double& getVdWRadius() ;
       const double& getRadius();

///// Setters:


    /*!< \brief Set the partial charge of the atom */
inline  void  setPartialCharge(const double& PCharge){partial_charge=PCharge;  }


    /*!< \brief Set the formal charge of the atom */
inline  void  setFormalCharge (const double& FCharge){ formal_charge=FCharge;  }


    /*!< \brief Set the BFactor of the atom */
inline  void  setBFactor      (const double& BF)     {       BFactor=BF;       }


    /*!< \brief Set the Id of the Atom */
inline  void  setNum          (const unsigned int &N){           num=N;        }


    /*!< \brief Set the Id of the Atom, as given by the file */
inline  void  setFNum         (const unsigned int& N){          fNum=N;        }

    /*!< \brief Set the interaction number id */
inline  void  setTNum         (const unsigned int& N){          tNum=N;        }


    /*!< \brief Set the Name of the Atom */
inline  void  setName         (const std::string& N) {         name=N;         }

        void  setUse          (const bool& newUse,
                               const bool &updateRes,
                               const bool& toBonds=true);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// ATOMIC DATA /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void checkMOL2type();
void setAtomicName(const std::string& nAtm);
void setAtomicNum (const unsigned int& atomicNum);
void setMOL2Type  (const std::string&  MOL2);





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// ENERGY ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


//double getplp(Coords atom_prot, int property_prot);



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// RESIDU ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

Residu* getResidu() const {return residu;}

const std::string getResiduName() const;
const std::string getChainName() const;






////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// ATOM LINKED /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


ItAtom            firstAtom() {return atomlinked.begin();}
ItAtom            lastAtom() {return atomlinked.end();}

ItCAtom            firstAtomC() {return atomlinked.begin();}
ItCAtom            lastAtomC() {return atomlinked.end();}

const Atom& getAtomLinked(const size_t& pos) const; 

size_t getNumBond() const {return atomlinked.size();}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// BONDS ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
bool  hasBondWith(const Atom &ve) const;
const Bond* getBondWith(Atom const &ve) const;
void  cleanBond();
Bond*  getBond(const size_t& pos) const;

ItBond            first() {return links.begin();}
ItBond            last() {return links.end();}



ItCBond            firstC() {return links.begin();}
ItCBond            lastC() {return links.end();}










////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// GENERICS //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

             void  clone           (const Atom& Atom);
const std::string  toString        (const bool &no_Bond=false,
                                    const bool &full_data=false) const;
const std::string  getIdentifier()const;

        void  assignBits();

double calcFixpos(const Atom& atom) const;
double calcFixpos(const Atom& atom, const double& threshold) const;








////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// PROPERTIES //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*!< \brief Determines if the atom is an hydrogen */
        bool  isHydrogen()       const;


/*!< \brief Determines if the atom is a metal
 *          (atomic number = 12, 20, 25, 26, 27, 28, 29, 30, 48) */
        bool  isMetallic()       const;


/*!< \brief Determines whethet the atom is a halogen */
        bool  isHalogen()        const;


/*!< \brief Determines if the atom is a nitrogen atom (atomic number = 7)*/
        bool  isNitrogen()       const;


/*!< \brief Determines if the atom is an oxygen atom (atomic number = 8)*/
        bool  isOxygen()         const;


/*!< \brief Determines if the atom is a phosphorus atom (atomic number = 15)*/
        bool  isPhosphorus()     const;


/*!< \brief Determines if the atom is a carbon atom (atomic number = 6)*/
        bool  isCarbon()         const;


/*!< \brief Determines if the atom is a sulfur atom (atomic number = 16)*/
        bool  isSulfur()         const;





        ////////////////////////////////////////////////
        ////////////////////////////////////////////////
        /////////////////// COLORS /////////////////////
        ////////////////////////////////////////////////
        ////////////////////////////////////////////////

        ////////////////////
        ///// SETTERS //////
        ////////////////////
void  setColor        (const double& Rx,const double& G, const double& B);//
void  setRcolor       (const double& Rx);//
void  setGcolor       (const double& G);//
void  setBcolor       (const double& B);//

        ////////////////////
        ///// GETTERS //////
        ////////////////////
    /*!< \brief Return the red color of the atom */
inline       double  getRcolor()        const       { return color[0];      }
    /*!< \brief Return the green color of the atom */
inline       double  getGcolor()        const       { return color[1];      }
    /*!< \brief Return the blue color of the atom */
inline       double  getBcolor()        const       { return color[2];      }



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// BOX /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

    void addBox(const Grid* const grid,  Box* const box);
    Box* getBox(const Grid * const grid) const;
    void delBox(Box* const box);




};

}


#endif // Atom_H
