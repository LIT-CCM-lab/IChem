#ifndef PHARMPROP_H
#define PHARMPROP_H

#include <string>

#include "global.h"
#include "coords.h"

namespace ICMole
{
/*! \class PharmProp
  * \brief Representation of physico-chemical properties in a pharmacophore
*/



class PharmProp
{
private :
    unsigned int TypePha;
    Coords c_base;
    Coords c_1A;
    Coords c_3A;
    int tol;


public:
    PharmProp( const unsigned int PharType,const Coords base,const Coords c1,const Coords c3,int tol2);

    Coords getBase() const {return c_base;}
    Coords get_1A() const { return c_1A;}
    Coords get_3A() const { return c_3A;}
    int get_tol() const {return tol;}
    inline unsigned int prop() const {return TypePha;}
    PharmProp& operator= (const PharmProp& _prop);
    const std::string toString() const; /*!< \brief Return a string with all properties */
    void settol(int tole){tol=tole;}

};
}
#endif // PHARMPROP_H
