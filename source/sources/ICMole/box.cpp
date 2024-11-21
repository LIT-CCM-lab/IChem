#include "headers/ICMole/box.h"
#include "headers/ICMole/atom.h"
using namespace ICMole;
using namespace std;



/**
  * \brief Box constructor
  *
  * \param _x : x coordinate of the box
  * \param _y : y coordinate of the box
  * \param _z : z coordinate of the box
*/

Box::Box(const double &x,
         const double &y,
         const double &z):Margin(false),used(true),pha_Used(false),atomlist(0),fixpos(x,y,z) {}


Box::~Box()
{
    const size_t listsize=atomlist.size();
    for (size_t i=0;i <listsize;++i)
    {
        Atom& atom = *atomlist.at(i);
        atom.delBox(this);
    }
}





void Box::delAtom     (Atom& atm)
{
    ItAtom it = find(atomlist.begin(),atomlist.end(),&atm);
    if (it != atomlist.end())
    {
        atomlist.erase(it);

    }
}


/**
  \fn bool Box::isAdjacent  (Box * const _box)
  \param _box : Box to look at
  Return true when the given box is adjacent to this box. Return false otherwise.
  */
bool Box::isAdjacent  (const Box * const _box) const
{
    if (_box == (Box*)NULL)return false;
    const double x(_box->rotpos.x-rotpos.x);
    const double y=(_box->rotpos.y-rotpos.y);
    const double z=(_box->rotpos.z-rotpos.z);
    if (fabs(x)<=1 && fabs(y)<=1 && fabs(z)<=1)return true;
    //const Coords diff(_box->rotpos-rotpos);
    //if (diff.x >= -1 && diff.x <= 1 && diff.y >= -1 && diff.y <= 1 && diff.z >= -1 && diff.z <= 1)return true;
    return false;
    /*if (find(Adjacent.begin(),Adjacent.end(),_box) == Adjacent.end()) return false;
    return true;*/
}
