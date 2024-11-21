#ifndef MOLECULE_H
#define MOLECULE_H


#include "residu.h"
#include "chain.h"
#include "coords.h"
#include "graphmatch.h"
namespace ICMole
{
struct templateData;


class Molecule
{
    friend class MoleReader;

private:

    static bool rules_loaded;

    Molecule(Molecule const &);           // No definition
    Molecule& operator=(Molecule const &);// No definition



protected:
    /*!< \brief List of Atoms */
    AtomList      Atoms;
    /*!< \brief List of Bonds */
    BondList      Bonds;
    /*!< \brief List of Residus */
    ResiduList    Residues;
    /*!< \brief List of Chains */
    ChainList    Chains;
    /*!< \brief List of Cycle */
    CycleList    Cycles;
    /*!< \brief Maximum Id for Atom*/
    unsigned int  maxNumAtom;
    /*!< \brief Maximum Id for Bond*/
    unsigned int  maxNumBond;
    /*!< \brief Maximum Id for Residu*/
    unsigned int  maxNumResidu;
    /*!< \brief Type of the molecule. See MoleType namespace */
    unsigned int  typemol;
    /*!< \brief Is the molecular fixpos is up to date or not */
    bool          bary_check;
    /*!< \brief Name of the molecule */
    std::string    name;


    Chain  dummyChain;
    Residu unknownRes;
    Residu aromaticRes;
    Residu cyclicRes;
    Complex *complex;
    Coords fixpos;
    Coords rotpos;



   std::vector<ResSet> ResSetList;
   std::vector<BondSet> BondSetList;
   std::vector<AtomSet> AtomSetList;
              AtomList ColorUpdate;
              AtomList ShowTypeUpdate;
   void addCycle(const AtomList& ListAt);

public:
   Coords extremH;
   Coords extremL;
   static unsigned int Rules[NB_MOLETYPE];
   static void loadRules();
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// CONSTRUCTORS ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
      Molecule(const unsigned int &NAtom  =50,
               const unsigned int &NBonds =100,
               const unsigned int MoleType=MoleType::PROTEIN);

    ////// Destructors :
    ~Molecule();










    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// COORDS ///////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    const Coords& getFixpos(const bool& withH=false);

    const Coords& getRotpos(const bool& withH=false);


      Coords getalpha(const int residu, const bool& withH=false);

      void getextrem(const bool& withH=false);

      Coords calcInertialMoment(const bool& withH=false,
                                const bool& onRotpos=false) ;

    void rotateMolecule(const std::vector<double>& matrix,
                        const Coords& transMobil,
                        const Coords& transRigid);



    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// ATOMS ////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    Atom& addAtom( const std::string& atomicName,
                   const std::string& name,
                   const std::string& mol2t,
                   Residu*            residu=(Residu*)NULL)
    throw(MoleExcept);//


    Atom& addAtom( const std::string& atomicName,
                   const double&      x,
                   const double&      y,
                   const double&      z,
                   const std::string& name="",
                   const std::string& mol2t="",
                   Residu*  res=(Residu*)NULL)
    throw(MoleExcept);//


    Atom& addAtom( const std::string& atomicName,
                   const Coords&     coords,
                   const std::string& name="",
                   const std::string& mol2t="",
                   Residu*  res=(Residu*)NULL)
    throw(MoleExcept);//



    void  addAtoms(const unsigned int &NtoAdd, Residu* const residu) throw(MoleExcept);


    void  delAtom( Atom *const Atom, const bool& with_maxnum=true)  throw(MoleExcept);


    void  delAtoms( AtomList& Atoms)                                throw(MoleExcept);


    void  renumAtom();


    const Atom& getAtom(const unsigned int &n,const bool &pos=true) const throw(MoleExcept);


    /*!< \brief Return the maximal Num for Atoms */
    const unsigned int&  getMaxNumAtom()   const  {return maxNumAtom;}

    /*!< \brief Return the number of Atom in this Molecule */
    unsigned int  numAtom()       const  {return (unsigned int)Atoms.size();}



    void reserveAtomSize(const unsigned int &NAtm);



    ////////////////////
    //// ITERATORS /////
    ////////////////////

    ItCAtom firstAtom()const {return Atoms.begin();}
    ItCAtom lastAtom() const {return Atoms.end();}







    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// BONDS ////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ///
    Bond&         addBond( Atom &Atom1, Atom &Atom2, const unsigned int &BondType) throw(MoleExcept);

    Bond&         addBond( Atom *const Atom1,
                           Atom *const Atom2,
                           const unsigned int &BondType) throw(MoleExcept);


    void          delBond(const Bond* const   Bond)          throw(MoleExcept);


    void          delBond(const   BondList& EdList)          throw(MoleExcept);


    void          delBond(const       Bond &  Bond)          throw(MoleExcept);


    void          renumBond();

    /*!< \brief Return the maximal Num for Bonds */
    const unsigned int&  getMaxNumBond()   const  {return maxNumBond;}

    /*!< \brief Return the number of Bonds in this Molecule */
    unsigned int  numBonds()        const  {return (unsigned int)Bonds.size();}

    const Bond& getBond(const size_t &n,const bool &pos=true) const throw(MoleExcept);



    void          reserveBondSize(const unsigned int &NBond);

    void cleanBond();

    ////////////////////
    //// ITERATORS /////
    ////////////////////

    ItCBond firstBond()const {return Bonds.begin();}
    ItCBond lastBond() const {return Bonds.end();}






    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////// RESIDUS ///////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////


    void          reserveResiduSize(const unsigned int &NResidu);

    Residu &      getUnknownRes()  {return unknownRes;}
    Residu &      getAromaticRes() {return aromaticRes;}

    size_t numResidus() const {return Residues.size();}
    const unsigned int& getMaxNumResidu()const  {return maxNumResidu;}
    void setResToAtom(Atom *const atm, Residu* const res) ;

    void moveResidu(Residu& residu, Molecule& molecule, const bool &delBonds=true) throw(MoleExcept);
    void moveResidu(ResiduList& residuList,
                    Molecule& molecule,
                    const bool &delBonds=true) throw(MoleExcept);
    size_t numChains() const {return Chains.size();}
    Chain& getChain(const size_t pos) throw(MoleExcept);

    ////////////////////
    ////// ADDERS //////
    ////////////////////
    Residu& addResidu(const std::string& chain,
                      const         int& resid,
                      const std::string& resname);

    void delResidu(Residu& residu);

    ////////////////////
    //// ITERATORS /////
    ////////////////////

    ItCRes firstResidu() const {return Residues.begin();}
    ItCRes lastResidu() const {return Residues.end();}
    Residu& getResidu(const size_t& pos) throw(MoleExcept);

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////// GENERICS //////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    std::string        toString(const bool& onlyUse=false)        const;
    std::string        getLongMoletype() const;
    const std::string& getName()         const {return name;}
    const unsigned int &getMoleType()     const {return typemol;}

    void setName(const std::string& new_name) {name=new_name;}

    void setMoleType(const unsigned int & type);

    void setComplex( Complex* cplx) {complex=cplx;}

    int SizeRes() { return Residues.size();}

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////// ALGORITHM //////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    /// Definitions in moletreat.cpp
    void createConnect(ResiduList& failedRes, AtomList& failedAtom,
                       const bool& cleanAll=false,
                       const std::string &verbose="")  throw(MoleExcept);

    void ringPerception() throw(MoleExcept);

    void createMoleGraph(Graph&,
                         const bool& onlyUsed=true, const bool& wHydrogen=true) const;


    void checkMOL2() throw(MoleExcept);

    void clear();

    void selChains(const std::vector<std::string> &allowedchains,
                   const bool &delete_all,
                   const std::string &verbose="") throw(MoleExcept);

    void matchTemplate(ResiduList &failedRes, AtomList &failedAtom, const std::string& verbose="",const bool& applyMOL2type=true) throw(MoleExcept);


    void cleanUnwanted(const std::string &verbose="");
    void selChains(const std::string& ST, const std::string &verbose="");
    void addHydrogen(const std::string &verbose="");
    void setUse(const bool& use);
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// CYCLES ///////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ItCCycle    firstCycle() const {return Cycles.begin();}
    ItCCycle    lastCycle()  const {return Cycles.end();}


    size_t numCycles()const {return Cycles.size();}

    size_t numUnkAtom()  const {return unknownRes.numAtom();}
    size_t numAromAtom() const {return aromaticRes.numAtom();}
    size_t numCyclAtom() const {return cyclicRes.numAtom();}
    Residu &getCycleRes() {return cyclicRes;}
    Residu &getAromRes() {return aromaticRes;}
    Residu &getUnkRes() {return unknownRes;}
    Cycle* getCycleFromCenter(Atom* const atom);
    void delCycle(Cycle*const cy) throw(MoleExcept);
    void renumResidu();
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////// CHAINS ///////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////

    ItCChain firstChain() const {return Chains.begin();}
    ItCChain lastChain() const {return Chains.end();}

    void setResiduTypes(const std::string&  verbose="") throw(MoleExcept);

    ////////////////////////////////////////////////
    ////////////////////////////////////////////////
    ////////////////////  SETS  ////////////////////
    ////////////////////////////////////////////////
    ////////////////////////////////////////////////
    void clearSets() {AtomSetList.clear();BondSetList.clear();;ResSetList.clear();;ColorUpdate.clear();}
    inline void addAtomShow(Atom* const atm){ShowTypeUpdate.push_back(atm);}
    inline void addAtomSet(AtomSet s){AtomSetList.push_back(s);}
    inline void addBondSet(BondSet s){BondSetList.push_back(s);}
    inline void addResSet (ResSet s){ResSetList.push_back(s);}
    inline unsigned int getNumSet() const {return (unsigned int)(AtomSetList.size()+BondSetList.size()+ResSetList.size());}
    inline ItAtom firstColAtom() {return ColorUpdate.begin();}
    inline ItAtom lastColAtom() {return ColorUpdate.end();}
    inline std::vector<AtomSet>::const_iterator firstAtomSet() const {return AtomSetList.begin();}
    inline std::vector<AtomSet>::const_iterator lastAtomSet()  const {return AtomSetList.end();}
    inline std::vector<BondSet>::iterator firstBondSet() {return BondSetList.begin();}
    inline std::vector<BondSet>::iterator lastBondSet() {return BondSetList.end();}
    inline std::vector<ResSet>::iterator firstResSet() {return ResSetList.begin();}
    inline std::vector<ResSet>::iterator lastResSet() {return ResSetList.end();}
    inline void addAtomColor(Atom* const atm){if (atm != (Atom*)NULL) ColorUpdate.push_back(atm);}

    /////////////////////////////////////////////////
    /////////////////////////////////////////////////
    //////////////////// ENERGIE ////////////////////
    /////////////////////////////////////////////////
    /////////////////////////////////////////////////

    void outplp() const;




};



}
#endif // MOLECULE_H
