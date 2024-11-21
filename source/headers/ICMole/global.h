#ifndef GLOBAL_H
#define GLOBAL_H

//#define ICHEM_DEBUG 1

#include <vector>
#include <map>
#include <fstream>
#include <algorithm>
#include "moleexcept.h"
#include "moledata.h"


namespace Levels {static const unsigned int NONE=0,NOTICE=1, WARNING=2, STRICT=3,FATAL=4;}
namespace ICMole
{

static unsigned int Moleaccess=Levels::NOTICE;

class Coords;
class Vertex;
class Edge;
class Graph;
class Atom;
class Bond;
class Molecule;
class Residu;
class Chain;
class Complex;
class Cycle;
class Grid;
class Box;
class Grim;
class Fingerprint;
class Similarity;
class PharmProp;

typedef std::vector<Edge*>    EdgeList;
typedef std::vector<Vertex*>  VertexList;
typedef std::vector<Coords>   CoordList;
typedef std::vector<Atom*>    AtomList;
typedef std::vector<Bond*>    BondList;
typedef std::vector<Residu*>  ResiduList;
typedef std::vector<Chain*>   ChainList;
typedef std::vector<Molecule*> MoleList;
typedef std::vector<Cycle*>    CycleList;
typedef std::vector<Box*>    BoxList;
typedef std::vector<PharmProp*> PharmList;


typedef BondList::iterator              ItBond;
typedef BondList::const_iterator        ItCBond;
typedef AtomList::iterator              ItAtom;
typedef AtomList::const_iterator        ItCAtom;
typedef EdgeList::iterator              ItEdge;
typedef EdgeList::const_iterator        ItCEdge;
typedef VertexList::iterator            ItVert;
typedef VertexList::const_iterator      ItCVert;
typedef ResiduList::iterator            ItRes;
typedef ResiduList::const_iterator      ItCRes;
typedef ChainList::iterator             ItChain;
typedef ChainList::const_iterator       ItCChain;
typedef MoleList::iterator             ItMole;
typedef MoleList::const_iterator       ItCMole;
typedef CycleList::iterator             ItCycle;
typedef CycleList::const_iterator       ItCCycle;
typedef BoxList::iterator             ItBox;
typedef BoxList::const_iterator       ItCBox;
typedef PharmList::const_iterator     ItcPharm;


extern std::string get_IChem_LIB_Path(const std::string alternate="");
extern void tokenStr(const std::string& str,
                     std::vector<std::string>& tokens,
                     const std::string& delimiters);
extern void tokenize(const std::string& str,
                      std::vector<unsigned int>& tokens,
                      const std::string& delimiters);

template<class T> extern void sortAndUnique(std::vector<T*> & n)
{
    sort(n.begin(),n.end());
//    unique(n.begin(),n.end());
    n.resize(std::distance(n.begin(),std::unique(n.begin(), n.end())));
}


template<class T> extern void sortAndUnique(std::vector<T> & n)
{
   sort(n.begin(),n.end());
//   unique(n.begin(),n.end());
   n.resize(std::distance(n.begin(),std::unique(n.begin(), n.end())));

}


template<class T> extern size_t getPos(std::vector<T> &n, T seekv)
{
    return std::distance(n.begin(),find(n.begin(),n.end(),seekv));
}

template <typename T>void SortUniqueVector(std::vector<T>& toCheck)
{
    sort(toCheck.begin(), toCheck.end());
    unique (toCheck.begin(), toCheck.end());
    toCheck.resize(std::distance(toCheck.begin(),std::unique (toCheck.begin(), toCheck.end())));
}

//bool cycleinf(Cycle *i, Cycle *j) {
//    if (i->atomlist.begin() < j->atomlist.begin() )
//        return true;
//    return false;
//}


//template <typename T>void SortUniqueCycle(std::vector<T>& toCheck)
//{
//    sort(toCheck.begin(), toCheck.end(),cycleinf);

//    toCheck.resize(std::distance(toCheck.begin(),std::unique (toCheck.begin(), toCheck.end())));
//}


double getplp(double plpA, double plpB, double plpC, double plpD, double plpE, double distance_plp) ;

typedef struct
{
    AtomList list;
    std::string name;
    std::string type;
    std::string subtype;
} AtomSet;
typedef struct
{
    BondList list;
    std::string name;
    std::string type;
    std::string subtype;
} BondSet;
typedef struct
{
    ResiduList list;
    std::string name;
    std::string type;
    std::string subtype;
} ResSet;



}
#endif // GLOBAL_H


