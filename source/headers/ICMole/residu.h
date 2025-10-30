#ifndef RESIDU_H
#define RESIDU_H
#include "graph.h"

namespace ICMole
{


struct HetData
{

    std::string HETreplace;
    bool hasTemplate;
    double MW;
    unsigned int ResT;
};

  class Residu
  {
    friend class MoleReader;
    friend class Molecule;
    friend class Chain;
private:
    static std::map<std::string,HetData> HETClass;
    static bool rules_loaded;
    static bool loadHETClass() ;

    /*!< \brief List of atoms included in this residu */
    AtomList      atoms;
    /*!< \brief String containing Residu name/Residu ID/Chain Name*/
    std::string   identifier;
    /*!< \brief Name of the residu */
    std::string   name;
    /*!< \brief Sybyl Bits of the residu */
    std::string   bits;
    /*!< \brief ID of this residu as given by the toolkit */
    int           num;
    /*!< \brief ID of this residu as given by the input file */
    int           fileNum;
    /*!< \brief Type of the residu. See ResType namespace for the list */
    unsigned int  resType;
    /*!< \brief List of inter substructure bond */
    ResiduList   InterRes;
    /*!< \brief Root atom of the residu i.e. alpha carbon  */
    Atom*         rootAtom;
    /*!< \brief Molecule where this residu is included in  */
    Molecule*     molecule;
    /*!< \brief Chain related to the residu */
    Chain*        chain;
    /*!< \brief Data associated to this residu */
    HetData*      hetdata;
    /*!< \brief Tell whether this atom should be used for analysis or not */
    bool          inUse;







////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// CONSTRUCTORS /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

    Residu( Molecule     &molecule,
           const  int         &Num,
           const  int         &FNum,
           const std::string  &Name,
           Chain              &chain,
           const unsigned int &ResType);

    Residu(Molecule *const    molecule,
           const  int         &Num,
           const  int         &FNum,
           const std::string  &Name,
           Chain              &chain,
           const unsigned int &ResType);



public:

    static unsigned int Rules[NB_MOLETYPE][NB_RESTYPE];
    static void loadRules(const bool &force=false);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// ATOMS ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

    /*!< \brief Add an atom to this residu */
    void addAtom(Atom& atom);


    /*!< \brief Remove an atom to this residu */
    void delAtom(const Atom &atom);


    /*!< \brief Return the number of atoms in this residu */
     size_t numAtom() const {return atoms.size();}

    ItAtom firstAtom() {return atoms.begin();}
    ItAtom lastAtom() {return atoms.end();}
    ItCAtom firstCAtom() const {return atoms.begin();}
    ItCAtom lastCAtom() const {return atoms.end();    }




 void setFNum(const int& newNum){fileNum=newNum;}
 void setName(const std::string &newName){name=newName;}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// CHAIN ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

    /*!< \brief Return the chain from which this residu is included */
      Chain* getChain() const {return chain;}


    /*!< \brief Return the name of the residu */
      const std::string& getChainName() const;










////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// GENERICS ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

      /*!< \brief Return the name, id and chain name of this residu*/
   const std::string& getIdentifier() const {return identifier;}


    /*!< \brief Return the residu number as given by IChem */
   const int &        getNum()        const {return num;}


   /*!< \brief Return the residu number as given by the input file */
   const int &        getFNum()       const {return fileNum;}


   /*!< \brief Return the name of the residu */
   const std::string& getName()       const {return name;}

    Atom& getAtom(const size_t&) const;

   /**
    * @brief addInterResidu : add a connected residu of this residu
    */
   void addInterResidu( Residu* const);

   size_t getNInterResidu() const {return InterRes.size();}

   ItCRes firstRes() const {return InterRes.begin();}
   ItCRes lastRes() const {return InterRes.end();}
   Residu& getInterRes(const size_t& pos) const {return *InterRes.at(pos);}

   const unsigned int& getResType() const {return resType;}

   /*!< \brief Return a string describing this residu and its atoms */
   const std::string  toString()      const;


   /*!< \brief A string describing the residu type (Amino Acid/Cofactor/Ligand ...)*/
   const std::string  getLongName()   const;

   Atom* getRootAtom() const {return rootAtom;}

   const bool& isUsed() const {return inUse;}

   void checkUse();

   void  setUse          (const bool& newUse,
                          const bool& updateAtom=true,
                          const bool& toBonds=true);

   void setNum(const unsigned int& number){num=number;}

   void checkResidu(bool wSetType=false) ;


   void setHetData(HetData* const newHetData){hetdata=newHetData;}
   void setResiduType(const unsigned int& newResType ){resType=newResType;}

   void applyResiduType(const std::string verbose="") ;

   Molecule * getParent() const {return molecule;}

   void getGraph(Graph &Gr) const;


   void getDistMatrix(std::vector<unsigned short>& matrix,
                      std::map<const Atom*,unsigned short>& order,
                      const bool& wHydrogen=false) const;

   const double& getWeight();
  };
}











#endif // RESIDU_H
