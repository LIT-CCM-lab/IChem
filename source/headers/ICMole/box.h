
#ifndef BOX_H
#define BOX_H
#include "global.h"
#include "coords.h"
#include "physprop.h"
namespace ICMole
{

class Box
{
    friend class Grid;
private:
    unsigned int id;                              /*!< \brief Id of the box */
    bool Margin;
    bool used;
    bool pha_Used;
    int clust;
    AtomList atomlist;

    std::multimap<unsigned int,Box*> adjacents;
    Box(const double &x, const double &y, const double &z);
    Box(Coords& ,double length);

public:
    Coords fixpos;
    Coords rotpos;
    PhysProp Props;
    ~Box();

    void setUse(const bool& use){used=use;  }
    const bool& getUse() const { return used;}

    inline       ItCAtom firstAtom  ()  const                               { return atomlist.begin();}                         /*!< \brief Return the first atom of this box */
    inline       ItCAtom lastAtom   ()  const                               { return atomlist.end();  }                         /*!< \brief Return the last atom of this box */


    inline         void addAtom     (Atom* const atom)                       { atomlist.push_back(atom);} /*!< \brief Add an atom which is inside this box */
    inline         void addAtom     (Atom& atom)                       { atomlist.push_back(&atom);} /*!< \brief Add an atom which is inside this box */
    void delAtom     (Atom& atm);
    inline         void setId       (const unsigned int _id)           { id=_id;}                                          /*!< \brief Set the id of the box*/
    inline         void setClust      (int _id)           { clust=_id;}                                          /*!< \brief Set the cluster id of the box*/
    inline unsigned int getClust       ()                           const { return clust;}                                       /*!< \brief Return the cluster id of the box */
    inline         void clearAtoms  ()                                 { atomlist.clear();}                                /*!< \brief Clear the list of atoms inside the box */
    inline unsigned int numAtoms    ()                           const { return atomlist.size(); }                         /*!< \brief Return the number of atoms inside this box */
    inline unsigned int getId       ()                           const { return id;}                                       /*!< \brief Return the id of the box */
    inline         void setMargin   (const bool isM)                   {Margin=isM;}
    inline         void setphaUsed  (const bool isU)                   {pha_Used=isU;}
    bool            isMargin()  const {return Margin;}
    bool            isMainBox() const {return !Margin;}
    bool            isUsed() const { return pha_Used;} /*!< \brief Return the pharma use status */
    void addAdjacent(const unsigned int& adj,  Box* const bx){adjacents.insert(std::pair<unsigned int,Box*>(adj,bx));}
    bool isAdjacent  (const Box * const _box) const;

};
}

#endif // BOX_H
