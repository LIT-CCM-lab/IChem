#ifndef COMPLEX_H
#define COMPLEX_H
#include "grid.h"

#define N_GRID_CPLX 4

namespace ICMole
{



/** \class Complex
   \brief Handle multiple molecule, such as protein/ligand complexes
   \version 3.0

 */
class Complex
{

    friend class Molecule;
    /*!< \brief List of all molecules */
    MoleList allmolecules;
    /*!< \brief Molecule grouped by MoleType */
    MoleList  molecules[NB_MOLETYPE];
    /*!< \brief Last molecule stored by MoleType */
    Molecule* lastmolecule[NB_MOLETYPE];
    /*!< \brief Atom List of all molecules */
    AtomList  atomlist;
    /*!< \brief Bond List of all molecules */
    BondList bondlist;
    /*!< \brief Residu list of all molecules */
    ResiduList residulist;
    /*!< \brief Chain List of all molecules */
    ChainList chainlist;
    /*!< \brief Grid points that represent the complex */
    Grid  grid[N_GRID_CPLX];
    /*!< \brief Tell whether grid are up to date */
    bool uptoGrid;
    /*!< \brief Length of each grid box */
    double boxStep[N_GRID_CPLX];


public:

    Complex();
    ~Complex();
    void addMolecule(Molecule *const molecule) ;
    size_t getNumMolecule(const unsigned int &MoleType) const ;
    size_t getNumMolecule() const;
    Molecule * getMole(const size_t& MoleType )const {return lastmolecule[MoleType];}
    void  clear() ;
    void  deleteMole(Molecule *const molecule) ;
    void removeMole(Molecule *const molecule) ;
    bool isMoleIn(Molecule* const molecule) const;
    bool isMoleIn(const Molecule&  molecule) const;

    void loadData(const bool &only_used=false);

    void updateMoleType(Molecule* const molecule,const unsigned int & type);

    void addAtom  ( Atom*   const atom  ) ;
    void addBond  ( Bond*   const bond  ) ;
    void addResidu( Residu* const residu) ;
    void addChain (Chain * const chain ) ;

    void deleteAtom  (const Atom*   const atom  ) ;
    void deleteBond  (const Bond*   const bond  ) ;
    void deleteResidu(const Residu* const residu) ;
    void deleteChain (const Chain*  const chain ) ;

    ItCAtom firstAtom() const {return atomlist.begin();}
    ItCAtom lastAtom()  const {return atomlist.end();}
    ItCBond firstBond() const {return bondlist.begin();}
    ItCBond lastBond()  const {return bondlist.end();}
    size_t getNumBond() const { return bondlist.size();}
    size_t getNumAtom() const { return atomlist.size();}
    size_t getNumResidu() const { return residulist.size();}
    ItCChain firstChain() const {return chainlist.begin();}
    ItCChain lastChain()  const {return chainlist.end();}

    ItCMole firstMole() const {return allmolecules.begin();}
    ItCMole lastMole() const {return allmolecules.end();}
    ItCMole firstMole(const size_t& moletype) const {return molecules[moletype].begin();}
    ItCMole lastMole(const size_t& moletype) const {return molecules[moletype].end();}

    ItCRes firstResidu() const {return residulist.begin();}
    ItCRes lastResidu() const {return residulist.end();}
    Grid& genGrid(const double& boxStep, const bool &wHydrogen=false) ;
    void clearGrid();
    int getMaxGrid() const;
    Grid& getGrid(const unsigned short& n) {return grid[n];}

    size_t numResidu() const {return residulist.size();}
    Molecule& getMoleAtPos(const size_t& pos) const {return * allmolecules.at(pos);}

    void splitMoleculeInToComplex(Molecule &molecule,const std::string &verbose="") ;
    std::string toString() const;
    void selectWater(const std::string &verbose="")  ;
    void selectLigand(const std::string& PDB_DIR, const ResiduList& failedRes,const AtomList& failedAtm, const std::string &verbose="", const bool &wUnDrugg=false);
    template <class IT,class BA> double calcMoleSmallDist(const IT& ref, const BA& comp, const double& best_dist, const bool& wHydrogen=false) const;
template<class BA> bool calcSmallDist(const Molecule& ref, const BA& comp, const double& best_dist, const bool& wHydrogen=false) ;
    void moveMoleToComplex(Molecule * const molecule, Complex &complex);
};
}

#endif // COMPLEX_H
