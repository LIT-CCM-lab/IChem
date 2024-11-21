#ifndef GRID_H
#define GRID_H
#include <set>
#include "global.h"
#include "coords.h"

#define R 72.9
#define PI 3.14159265
#define NB_PROJ 120
#define U 0.5
//#define ICHEM_DEBUG

namespace ICMole
{
class Grid
{

    friend class VolSite;

friend class Complex;

private:
    /*!< \brief Number of cube per range on the i axis. Default : 0 */
    int           numPtRangeI;

    /*!< \brief Number of cube per range on the j axis. Default : 0 */
    int           numPtRangeJ;

     /*!< \brief Number of cube per range on the k axis. Default : 0 */
    int           numPtRangeK;

    /*!< \brief Total number of cubes into the main cube. Default : 0 */
    int           maxCubeNum;

    int           adjLevel;

    /**
     * @brief status : 0- No Calculation made
     *  1 - Grid generated
     *  2 - Color cube generated and up to date
     */
    short        status;

    /*!< \brief Take hydrogen into account */
    bool          wHydrogen;

    /*!< \brief True for VolSite Calculation, false otherwise */
    bool         wVolSite;

    /*!< \brief Length of a cube. Default : 1.5 Angstroems */
    double        boxStep;

    /*!< \brief Rotation matrix to switch from x,y,z to i,j,k axis */
    double       rotMatrix[9];


    double      projLength;

    /*!< \brief i axis vector */
    Coords       vect_i;

    /*!< \brief j axis vector */
    Coords       vect_j;

    /*!< \brief k axis vector */
    Coords       vect_k;

    /*!< \brief Starter position of the main cube */
    Coords       X8;

    /*!< \brief Center of the grid */
    Coords       center;

    /*!< \brief Size of the grid */
    Coords       box_size;

    /*!< \brief All boxes */
    BoxList      AllBoxes;

    void prepareMatrix(const Coords& center,
                             const Coords& pt1,
                             const Coords& pt2);
    void  genBoxes(const double& dist1,
                   const double& dist2,
                   const double& dist3,
                   const double& dist4,
                   const double& dist5,
                   const double& dist6 );

public :
    ////////////////////////////////////////////////
    ////////////////////////////////////////////////
    ///////////////// CONSTRUCTORS /////////////////
    ////////////////////////////////////////////////
    ////////////////////////////////////////////////
    Grid();
    ~Grid();


    ////////////////////////////////////////////////
    ////////////////////////////////////////////////
    /////////////// GETTERS/SETTERS ////////////////
    ////////////////////////////////////////////////
    ////////////////////////////////////////////////
    ////////////////////
    ///// SETTERS //////
    ////////////////////
    inline   void setProjLength(const double      &L)   {projLength=L;   }
    ////////////////////
    ///// GETTERS //////
    ////////////////////
    inline    const    int&  getNumPtRangeI() const {return numPtRangeI;}			                                                                      /*!< Nomber of cube per range on the i axis */
    inline    const    int&  getNumPtRangeJ() const {return numPtRangeJ;}
    inline    const    int&  getNumPtRangeK() const {return numPtRangeK;}
    inline    const    int&  getMaxCubeNum()  const {return maxCubeNum;}
    inline    const double&  getBoxStep()     const {return boxStep;}
    inline    const Coords&  getcenter()     const {return center;}

    ////////////////////////////////////////////////
    ////////////////////////////////////////////////
    //////////////// MISCELLANEOUS /////////////////
    ////////////////////////////////////////////////
    ////////////////////////////////////////////////
    ///////////////////////
    // GRID CONSTRUCTION //
    ///////////////////////
    void  createComplexBoxes  (Complex& complex,
                               const double& box_Step)
    throw(MoleExcept);
    void  createMolecularBoxes(Molecule& molecule,
                               const double& box_Step)
    throw(MoleExcept);
    void  createComplexBoxes  (Molecule&  protein,
                               Molecule & ligand,
                               const double & box_Step,
                               const Coords& box_size)
    throw(MoleExcept);
    void  createComplexBoxes  (Complex&  protein,
                               Molecule & ligand,
                               const double & box_Step,
                               const Coords& box_size)
    throw(MoleExcept);
    void createParentGrid(const Grid& grid,
                          const double&box_Step)
    throw(MoleExcept);
    ///////////////////////
    //////// UTILS ////////
    ///////////////////////
    void rotating(Coords &ini) const;
    void rotateMole(Molecule & mole,

                    const bool& find_box=true, const bool &withResSelRules=true)const;
    void rotateAtoms(AtomList &atomlist,const bool& find_box=true)const;
    ///////////////////////
    /////// EXPORTS ///////
    ///////////////////////
    void printInFile(const std::string& f,
                     const BoxList& liste,
                     const bool& wBorder=false, const bool& onlyUse=false) const;
    ///////////////////////
    /////// BOX PROCESS ///
    ///////////////////////

    inline void setWHydrogen    (const bool& H){wHydrogen=H;}

    void setRule(const unsigned int& MType,
                 const unsigned int& RType,
                 const unsigned int& RVal);
    void setIniValues();

    ItBox  first(){return AllBoxes.begin();}
    ItBox  end  (){return AllBoxes.end();}
    Box* getBox(const size_t& N)const {return AllBoxes.at(N);}
    Box* getBox(int i, int j, int k) const {return AllBoxes.at(i-1+(j)*numPtRangeI+(k)*numPtRangeI*numPtRangeJ);}

    BoxList getAllBox() { return AllBoxes;}

    void getAdjacency(std::multimap<unsigned int, Box *> &liste, const Box * const box,
                      const  int& min_level,
                      const  int& max_level);
    void getAdjacency(std::vector<Box*>& liste, const Box* const box,
                      const  double &dist, const bool& spherical=false,
                      const bool& borderonly=false);
    void getFullAdjacency(std::vector<Box*>& liste, const Box* const box,
                      const  double &dist) const;

    void colorCube(Molecule & mole,
                   const bool& withRotation=true,
                   const bool& withResSelRules=true,
                   const bool& withBorder=false)throw(MoleExcept);
    Box* getBox(const Coords& pos)const;
    void getAdjacentAtoms(std::vector<Atom*>& liste,const Box&  box, const double &dist,const bool &wHydrogen=false) const;
    void calcBSA(Molecule &protein, Molecule& ligand, double&bsa, double &volume_ligand);
    void calcBSA(Complex& complex, Molecule& ligand, double&bsa, double &volume_ligand);
    void getAdjacentAtoms(std::vector<Atom*>& liste,const Atom&  atom, const double &dist,const bool &wHydrogen=false) const;
void clear();
};

}

#endif // GRID_H
