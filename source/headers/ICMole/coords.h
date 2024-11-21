#ifndef COORDS_H
#define COORDS_H
#include <math.h>
#include "moleexcept.h"
#include <iostream>

/**
 * \file coords.h
 * \brief Representation of 3D coordinates
 * \author DESAPHY Jeremy
 * \version 4.0
 * \date March 22th 2013
 *
 * Allow to handle 3D coordinates of any point: atoms, molecules ...
 *
 */


/*! \namespace ICMole
 *
 * Anything related to the handling of molecular representation is in this namespace
 */
namespace ICMole
{


/*! \class Coords
 * \brief Class handling 3D coordinates
 */
class Coords
{    


public:
     double x, y, z; /*!< Represent the position in each axis */
    ///////////////
    // CONSTRUCTORS
    ///////////////
    Coords();
    Coords(const double px, const double py, const double pz);
    Coords(Coords const &,const bool& norm=false);
    ~Coords();
    /////////////////////////////
    /////::Accessors/Mutators
    /////////////////////////////

    /*!
       *  \brief Reset coordinates to 0
       *
       */
    inline void clear() { x=0; y=0; z=0;}

    void getCoords() {
      std::cout << "CoordX: " << x << std::endl;
      std::cout << "CoordY: " << y << std::endl;
      std::cout << "CoordZ: " << z << std::endl;
    };


    /**
       \fn void setCoords(const double x_, const double y_, const double z_)
         \brief Setting new values

         Set the coordinates of x,y,z  values<br/>
         How to use : Coords c1(2,3,4); c1.Set_coords(1,2,3);
         \param x_ : The new x value
         \param y_ : The new y value
         \param z_ : The new z value
       */
    inline void setCoords(const double &px, const double &py, const double &pz) { x=px; y=py; z=pz;}

    /**
        \fn  void setCoords(const Coords &coords)
         \brief Setting coordinates

         Set the coordinates of the parameter object into the called object<br/>
         How to use : Coords c1(2,3,4); Coords c2(1,2,3); c2.Set_coords(c1);<br/>
         Same as Coords& operator=  (const Coords& _coords);
         \param c : The Coords object to set with
       */
    inline void setCoords(const Coords &coords) {x=coords.x; y=coords.y;z=coords.z;}
    void swap(Coords &coords);

    ////////////
    // OPERATORS
    /////////////
    void  operator-= (const Coords &coords);
    void  operator+= (const Coords &coords);
    void  operator*= (const Coords &coords);
    void  operator/= (const Coords &coords) throw(MoleExcept);
    bool  operator== (const Coords& coords) const;
    bool  operator!= (const Coords& coords) const;

    Coords& operator=  (const Coords& coords);

    Coords operator-   (const Coords &coords) const { return Coords(x-coords.x,y-coords.y,z-coords.z);}
    Coords operator+   (const Coords &coords) const { return Coords(x+coords.x,y+coords.y,z+coords.z);}
    Coords operator*   (const Coords &coords) const { return Coords(x*coords.x,y*coords.y,z*coords.z);}

void normalize();
    /////////////////////////////
    /////::Rotation
    /////////////////////////////
    void  rotateX (const double& rad);
    void  rotateY (const double& rad);
    void  rotateZ (const double& rad);
    void  rotate(const double &rad_x, const double &rad_y, const double &rad_z);

    /////////////////////////////
    /////::Calculators
    /////////////////////////////
    Coords  getNormal(const Coords &B, const Coords &C) const;

    double  calcDist (const Coords &_coords) const;
    double  calcDist (const Coords &_coords, const double &max) const;
    double  calcAngle(const Coords &_coords1, const Coords &_coords2) const;

    /////////////////////////////
    /////::Translation
    /////////////////////////////
    inline void  translateX(const double& _x) {x+=_x;}                               /*!< Translate the x coordinate */
    inline void  translateY(const double& _y) {y+=_y;}                               /*!< Translate the y coordinate */
    inline void  translateZ(const double& _z) {z+=_z;}                               /*!< Translate the z coordinate */
    inline void  translate(const double &_x, const double &_y, const double &_z) {x+=_x;y+=_y;z+=_z;}



    /////////////////////////////
    /////Operator templates :
    /////////////////////////////
    /*!
       *  \brief Multiplication affectation operator
       *
       *  This operator allows you to multiplicate the coordinates by a scalar, either int, float or double and affect directly in it<br/>
       *  How to use : int t = 4; Coords c2(1,2,3); c2*=t;
       *  \param c : The object to multiplicate by
       *  \param a : The scalar
       */
    template < typename T > friend  void operator*=( Coords &c, const T &a)
    {
    c.setCoords(c.x*a, c.y*a, c.z*a);
    }
    /*!
      *  \brief Multiplication affectation operator
      *
      *  This operator allows you to multiplicate the coordinates by a scalar, either int, float or double and affect directly in it<br/>
      *  How to use : int t = 4; Coords c2(1,2,3); c2*=t;
      *  \param c : The object to multiplicate by
      *  \param a : The scalar
      *  \throw MoleExcept : when the scalar equals 0
      */
    template < typename T > friend  void operator/=( Coords &c, const T &a) throw(MoleExcept)
    {
        if (a ==0 ) throw MoleExcept(1010201,"Coords::operator/=","a is equal to 0");

    c.setCoords(c.x/a, c.y/a, c.z/a);
    }
    /*!
       *  \brief Substraction affectation operator
       *
       *  This operator allows you to substract the coordinates by a scalar, either int, float or double and affect directly in it<br/>
       *  How to use : int t = 4; Coords c2(1,2,3); c2-=t;
       *  \param c : The object to substract by
       *  \param a : The scalar
       */
    template < typename T > friend  void operator-=( Coords &c, const T &a)
    {
        c.x-=a; c.y-=a; c.z-=a;
    }
    /*!
       *  \brief Addition affectation operator
       *
       *  This operator allows you to add the coordinates by a scalar, either int, float or double and affect directly in it<br/>
       *  How to use : int t = 4; Coords c2(1,2,3); c2+=t;
       *  \param c : The object to add by
       *  \param a : The scalar
       */
    template < typename T > friend  void operator+=( Coords &c, const T &a)
    {
        c.x+=a; c.y+=a; c.z+=a;
    }
    /*!
       *  \brief Multiplication operator
       *
       *  This operator allows you to multiplicate the coordinates by a scalar, either int, float or double<br/>
       *  How to use : int t = 4; Coords c2(1,2,3); Coords c(2,3,4);  c2=c*t;
       *  \param c : The object to multiplicate by
       *  \param a : The scalar
       */
    template < typename T > friend  Coords operator*(const Coords &c, const T &a)
    {
        return Coords(c.x*a, c.y*a, c.z*a);
    }
    /*!
       *  \brief Addition operator
       *
       *  This operator allows you to add the coordinates by a scalar, either int, float or double<br/>
       *  How to use : int t = 4; Coords c2(1,2,3); Coords c(2,3,4);  c2=c+t;
       *  \param c : The object to add by
       *  \param a : The scalar
       */
    template < typename T > friend  Coords operator+(const Coords &c, const T &a)
    {
        return Coords(c.x+a,c.y+a,c.z+a);
    }
    /*!
       *  \brief Division operator
       *
       *  This operator allows you to divide the coordinates by a scalar, either int, float or double<br/>
       *  How to use : int t = 4; Coords c2(1,2,3); Coords c(2,3,4);  c2=c/t;
       *  \param c : The object to divide by
       *  \param a : The scalar
       */
    template < typename T > friend  Coords operator/(const Coords &c, const T &a) throw(MoleExcept)
    {
        if (a==0) throw MoleExcept(1010301,"Coords::operator/","a is equal to 0");
        return Coords(c.x/a, c.y/a,c.z/a);

    }
    /*!
       *  \brief Substract operator
       *
       *  This operator allows you to substract the coordinates by a scalar, either int, float or double<br/>
       *  How to use : int t = 4; Coords c2(1,2,3); Coords c(2,3,4);  c2=c-t;
       *  \param c : The object to substract by
       *  \param a : The scalar
       */
    template < typename T > friend  Coords operator-(const Coords &c, const T &a)
    {
        return Coords(c.x-a,c.y-a,c.z-a);
    //Coords result = c; result.setCoords(result.x-a, result.y-a, result.z-a); return result;
    }

    std::string toString() const;

};
}
#endif // COORDS_H
