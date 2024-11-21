#ifndef MOLEEXCEPT_H
#define MOLEEXCEPT_H

#include <string>
#include <sstream>

/*
 1010101 - Coords::operator/= (const Coords & coords) : x, y or z value within coords is set to 0
 1010201 - Coords::operator/= : a is equal to 0
 1010301 - Coords::operator/ :  a is equal to 0
 1020101 - Vertex::delEdge : Given edge is not part of this vertex
 1030101 - Edge::~Edge : Given edge is not part of this vertex
 1030201 - Edge::getOtherVertex : Given vertex is not part of this edge
 1030301 - Edge::getOtherVertex : Given vertex is not part of this edge | Vertex::getEdgeWith | Vertex::toString (not following)
 1030401 - Edge::shareVertex : No shared vertex
 1040101 - Graph::addVertex : Bad allocation exception - exiting.
 1040102 - Graph::addVertexs : Bad allocation exception - exiting.
 1040301 - Graph::delVertex : No vertex given
 1040302 - Graph::delVertex : This graph doesn't contains the following vertex
 1040401 - Graph::delVertexs : No vertex given
 1040402 - Graph::delVertexs : This graph doesn't contains the following vertex
 1040501 - Graph::getVertex :  given value above the number of vertexs within the graph
 1040502 - Graph::getVertex : No vertex found with the given id
 1040601 - Graph::addEdge : Vertexs given for the creation of the edge are the same
 1040602 - Graph::addEdge : First Vertex is not part of this graph
 1040603 - Graph::addEdge : Second Vertex is not part of this graph
 1040604 - Graph::addEdge : Bad allocation
 1040701 - Graph::delEdge : No edge given
 1040702 - Graph::delEdge : This graph doesn't contains the following edge
 1040801 - Graph::delEdge : No edge given
 1040802 - Graph::delEdge : This graph doesn't contains the following edge
 1040901 - Graph::delEdge : This graph doesn't contains the following edge
 1050101 - Atom::delBond : Given Bond is not part of this Atom
 1050201 - Atom::loadAtomicData : Atom is not a biologically relevant atom
 1050202 - Atom::loadAtomicData : Unknown atomic name of atom
 1050301 - Atom::checkMOL2type : No MOL2 Type given | Atom::setMOL2Type
 1050302 - Atom::checkMOL2type : MOL2 Type given is too long (max length:5) | Atom::setMOL2Type
 1050303 - Atom::Check_mol2type : MOL2 type given : Not recognized
 1050401 - Atom::setAtomicNum : Given number does not correspond to an atomic number
 1050501 - Atom::checkMOL2type : No MOL2 Type given
 1060101 - Bond::~Bond : Given Bond is not part of this Atom
 1060201 - Bond::getOtherAtom : Given Atom is not part of this Bond | Atom::getBondWith
 1060301 - Bond::getOtherAtom : Given Atom is not part of this Bond
 1060401 - Bond::shareAtom : No shared Atom
 1090102 - Chain::Chain : No chain name given
 1090201 - Chain::Chain : No parent molecule given
 1090202 - Chain::Chain : No chain name given
 1090301 - Chain::addResidu : No residu given
 1090401 - Chain::addResidu : No residu given
 1090402 - Chain::addResidu : Residu not found in Chain
 1100101 - Complex::getNumMolecule : Given MoleType is above the maximum
 1100201 - Complex::delMole : Given molecule is not part of this complex
 1110101 - GLOBAL::loadHETClass : No ICHEM_LIB DEFINED
 1100102 - GLOBAL::loadHETClass : Cannot open HETLIST file
 1120101 - Grid::createComplexBoxes : Box step is null or negative. MUST BE POSITIVE VALUE
 1120102 - Grid::createComplexBoxes : Molecule in complex has no atoms
 1120103 - Grid::createComplexBoxes : No molecule in complex
 1120104 - Grid::createComplexBoxes : Number of small boxes less than 5. Must increase box_size or reduce box_step
 1120201 - Grid::createComplexBoxes : When no protein given or no atoms in it
 1120202 - Grid::createComplexBoxes : When no ligand given or no atoms in it
 1120203 - Grid::createComplexBoxes : Box step is below 0 or equals 0
 1120204 - Grid::createComplexBoxes : Not enough boxes on one axis. Minimal value is 5.
 1120301 - Grid::createMolecularBoxes : When no molecule given or no atoms in it
 1120302 - Grid::createMolecularBoxes : Box step is below 0 or equals 0
 1120303 - Grid::createMolecularBoxes : Selection rules for this molecule is set to 'ignor
 1120303 - Grid::createMolecularBoxes : Not enough boxes on one axis. Minimal value is 5.
 1120401 - Grid::genBoxes : Bad allocation exception - exiting

*/


namespace ICMole{
/**
 * \class MoleExcept
 * \brief Handle molecular exceptions
 * \author DESAPHY Jeremy
 * \version 4.0
 * \date 22 March 2013
 *
 * For the integer value, 7 digits are required : <br/>
 *  <ul><li> First digit - NAMESPACE </li>
 *      <li> Second and Thrid digits - CLASS </li>
 *      <li> Fourth and Fifth digit - FUNCTION</li>
 *      <li> Sixth & Seventh digits - ERROR TYPE</li>
 *</ul>
 *
 * ICMole: 1xxxxxx
 *   Coords   : 101xxxx
 *   Vertex   : 102xxxx
 *   Edge     : 103xxxx
 *   Graph    : 104xxxx
 *   Atom     : 105xxxx
 *   Bond     : 106xxxx
 *   Molecule : 107xxxx
 *   Residu   : 108xxxx
 *   Chain    : 109xxxx
 *   Complex  : 110xxxx
 *   Global   : 111xxxx
 *   Grid     : 112xxxx
 *   GraphMatch:113xxxx
 * ICPars: 2xxxxxx
 *   MoleReader : 201xxxx
 * ICCalcs: 3xxxxxx
 *   VolSite  : 301xxxx
 *   Interactions: 302xxxx
 *   PDBConvert : 303xxxx
 * ICTools: 9xxxxxx
 *  switch :    901xxxx
 *  genKey :    902xxxx
 *  BSACalc:    903xxxx
 *  ints   :    904xxxx
 *  grim   :    905xxxx
 *  IFP    :    907xxxx
 *  IFPscreen : 908xxxx
 */


class MoleExcept
{
private:
    unsigned int cause;  /*!< Integer specifying the issue*/
     std::string source; /*!< Where the issue comes from*/
     std::string data;   /*!< What kind of issue it is */
     std::string trace;
public:
     /*!
        *  \brief Constructor
        *
        *  Constructor for exception handling
        *  \param code : Integer specifying the issue - 5 digits.
        *  \param source : Where the issue comes from (class and function)
        *  \param data :  What kind of issue it is
        */
     MoleExcept(const unsigned int code, const std::string source="", const std::string data=""):cause(code),source(source),data(data){ trace+="  ||--> "+source; trace +=" "; trace+="\n";}
    /*!
       *  \brief Standard Constructor
       *
       */
    MoleExcept(const MoleExcept &e) : cause(e.cause),source(e.source),data(e.data) {}
    /*!
       *  \brief Return the integer specifying the issue
       *
       */
    inline unsigned int getCode()           {return cause;}
    inline std::string getSource()          {return source;}
    inline std::string getData()            {return data;}
    inline void setCode(const unsigned int c){cause=c;}
    inline void addTrace(std::string s){trace+="  ||--> "+s; trace +=" "; trace+="\n";}
    inline std::string getTrace()const {return trace;}
    inline std::string what()const {
        std::ostringstream N;
        N<<"  Code : "<<cause<<std::endl
         <<"Source : "<<source<<std::endl
         <<"  Data : "<<data<<std::endl<<std::endl
         <<" Trace : "<<std::endl<< trace<<std::endl; return N.str();}
};
}
#endif // MOLEEXCEPT_H
