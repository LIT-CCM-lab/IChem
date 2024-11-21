#include "headers/ICMole/coords.h"


using namespace ICMole;



////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// CONSTRUCTOR || DESTRUCTORS //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////


/*!
     *  \brief Constructor
     *
     *  Coords Constructor class.
     *  Set coordinates to 0.
     *  Not recommended
     */
Coords::Coords():x(0),y(0),z(0) {}
/*!
     *  \brief Constructor
     *
     *  Coords Constructor class.
     *  \param px : x axis value
     *  \param py : y axis value
     *  \param pz : z axis value
     */
Coords::Coords(const double px, const double py, const double pz):x(px),y(py),z(pz){}

/*!
     *  \brief Copy constructor
     *
     *  Coords copy constructor
     *  \param coords : Other coordinates object
     */
Coords::Coords(Coords const &coords, const bool &normalize)
    {
        x = coords.x;
        y=coords.y;
        z=coords.z;
        if (normalize)
        {
        const double norm = sqrt(x*x+y*y+z*z);
        x/=norm; y/=norm; z/=norm;
        }
    }


Coords::   ~Coords(){}


/*!
   *  \brief Swapping content
   *
   *  This function allows you to switch values between two Coords object.
   *  How to use : Coords c1(2,3,4); Coords c2(1,2,3); c1.Swap(c2);
   *   c2 will have for coordinates (2,3,4) and c1 (1,2,3)
   *  \param coords : The object to swap with
   */
void Coords::swap(Coords &pcoords)
{
    const Coords t=Coords(pcoords);
    pcoords=*this;
    x=t.x;y=t.y;z=t.z;
}


void Coords::normalize()
{
    const double norm = sqrt(x*x+y*y+z*z);
    x/=norm; y/=norm; z/=norm;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// OPERATORS //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////


/*!
   *  \brief Substract affectation operator
   *
   *  This operator allows you to substract the coordinates of another Coords object
   *  How to use : Coords c1(2,3,4); Coords c2(1,2,3); c1-=c2;
   *  \param coords : The object to substract to the reference
   */

void Coords:: operator-= (const Coords &coords)  { x-= coords.x; y-= coords.y; z-= coords.z; }//CHECK
/*!
   *  \brief Multiplication affectation operator
   *
   *  This operator allows you to multiply the coordinates of another Coords object
   *  How to use : Coords c1(2,3,4); Coords c2(1,2,3); c1*=c2;
   *  \param coords : The object to multiply to the reference
   */
void Coords:: operator*= (const Coords & coords)  { x*= coords.x;  y*= coords.y;  z*= coords.z;  }//CHECK
/*!
   *  \fn void Coords:: operator/= (const Coords & coords)
   *\brief Division affectation operator
   *
   *  This operator allows you to divide the coordinates of another Coords object
   *  How to use : Coords c1(2,3,4); Coords c2(1,2,3); c1/=c2;
   *  \param coords : The object to divide to the reference
   *  \throw MoleExcept 1010101 - When x,y or z is equal to 0
   */
void Coords:: operator/= (const Coords & coords)  throw(MoleExcept){
  if (coords.x ==0 || coords.y == 0 || coords.z == 0) throw MoleExcept(1010101,"Coords::operator/=","x,y or z is equal to 0");
  x/= coords.x;  y/= coords.y;  z/= coords.z;  }//CHECK
/*!
   *  \brief Addition affectation operator
   *
   *  This operator allows you to add the coordinates of another Coords object
   *  How to use : Coords c1(2,3,4); Coords c2(1,2,3); c1+=c2;
   *  \param coords : The object to add to the reference
   */
void Coords:: operator+= (const Coords & coords)  { x+= coords.x;  y+= coords.y;  z+= coords.z;  }//CHECK

/*!
   *  \brief Comparison operator
   *
   *  This operator allows you to compare two Coords object.
   *  Return 1 if the two objects have the same coordinates, 0 otherwise
   *  How to use : Coords c1(2,3,4); Coords c2(1,2,3); if(c1==c2) cout << "THE SAME"<<endl;
   *  \param v : The Coords object to compare with
   */
bool Coords:: operator== (const Coords & v) const {  return (x == v.x && y == v.y && z == v.z);  }//CHECK
/*!
   *  \brief Opposite Comparison operator
   *
   *  This operator allows you to compare two Coords object.
   *  Return 0 if the two objects have the same coordinates, 1 otherwise
   *  How to use : Coords c1(2,3,4); Coords c2(1,2,3); if(c1!=c2) cout << "NOT THE SAME"<<endl;
   *  \param v : The Coords object to compare with
   */
bool Coords:: operator!= (const Coords & v) const {  return !(x == v.x && y == v.y && z == v.z); }//CHECK
/*!
   *  \brief Affectation operator
   *
   *  This operator allows you to affect the value of one Coords object into another.
   *  Return a new Coords object. No need to delete it.
   *  How to use : Coords c1(2,3,4); Coords c2(1,2,3); c1=c2;
   *  \param _coords : The Coords object to set to
   */
Coords& Coords::operator=  (const Coords& _coords)
{
    x = _coords.x; y=_coords.y; z=_coords.z;
    return *this;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// ROTATION ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/*!
   *  \brief Rotation around x axis
   *
   *  This operator allows you to rotate the point (so all the 3D coordinates) around the x axis. <br/>
   *  How to use : Coords c1(2,3,4); Coords c2(1,2,3); c1.Rot_x(0.2);
   * <table style="border:none;text-align:center;"><tr><td>[</td><td> 1</td><td>0 </td><td> 0</td><td>]</td></tr>
   * <tr><td>  [ </td><td>0  </td><td>cos(rad)</td><td> -sin(rad)</td><td>]</td></tr>
   * <tr><td>   [</td><td> 0 </td><td> sin(rad) </td><td> cos(rad)</td><td>]</td></tr></table>
   *  \param rad : rotating of rad radian
   */

void Coords::rotateX(const double& rad) { const double _x= x; const double _y= cos(rad)*y-sin(rad)*z; const double _z=  sin(rad)*y+cos(rad)*z; x=_x;y=_y;z=_z; }

/*!
   *  \brief Rotation around y axis
   *
   *  This operator allows you to rotate the point (so all the 3D coordinates) around the y axis. <br/>
   *  How to use : Coords c1(2,3,4); Coords c2(1,2,3); c1.Rot_y(0.3);
   * <table style="border:none;text-align:center;"><tr><td>[</td><td> cos(rad)</td><td>0 </td><td> sin(rad)</td><td>]</td></tr>
   * <tr><td>   [ </td><td>0  </td><td>1</td><td> 0</td><td>]</td></tr>
   * <tr><td>   [</td><td> -sin(rad) </td><td> 0 </td><td> cos(rad)</td><td>]</td></tr></table>
   *  \param rad : rotating of rad radian
   */
void Coords::rotateY(const double& rad) {const double _y= y; const double _x= cos(rad)*x+sin(rad)*z; const double _z= -sin(rad)*x+cos(rad)*z; x=_x;y=_y;z=_z; }

/*!
   *  \brief Rotation around z axis
   *
   *  This operator allows you to rotate the point (so all the 3D coordinates) around the z axis. <br/>
   *  How to use : Coords c1(2,3,4); Coords c2(1,2,3); c1.Rot_z(0.3);
   * <table style="border:none;text-align:center;"><tr><td>[</td><td> cos(rad)</td><td>-sin(rad) </td><td> 0</td><td>]</td></tr>
   * <tr><td>   [ </td><td>sin(rad) </td><td>cos(rad)</td><td> 0</td><td>]</td></tr>
   * <tr><td>   [</td><td> 0 </td><td> 0 </td><td> 1</td><td>]</td></tr></table>
   *  \param rad : rotating of rad radian
   */
void Coords::rotateZ(const double& rad) { const double _z= z; const double _y=sin(rad)*x+cos(rad)*y; const double _x=  cos(rad)*x-sin(rad)*y; x=_x;y=_y;z=_z; }


/*!
   *  \brief Rotation around all axis axis
   *
   *  This operator allows you to rotate the point (so all the 3D coordinates) around all axis. <br/>
   *  Will call Rot_x(), Rot_y(), and Rot_z()
   *  \param rad_x : rotation of rad_x radian around x axis
   *  \param rad_y : rotation of rad_y radian around y axis
   *  \param rad_z : rotation of rad_z radian around z axis
   */
void Coords::rotate(const double &rad_x, const double &rad_y, const double &rad_z) { rotateX(rad_x);rotateY(rad_y) ; rotateZ(rad_z);}




///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// CALCULATORS //////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////



/*!
   *  \brief Calculate the angle between the caller object and two other Coords objects
   *
   *  Calculate angle between 3 Coords objects. The caller object is the center of the angle
   *  \param _coords1 : First Coords object
   *  \param _coords2 : Second Coords object
   *  \return angle : in radian
   */
double Coords::calcAngle(const Coords &_coords1, const Coords &_coords2) const
{
  const double A1=sqrt((this->x- _coords1.x)*(this->x- _coords1.x)
                      +(this->y- _coords1.y)*(this->y- _coords1.y)
                      +(this->z- _coords1.z)*(this->z- _coords1.z))*
                  sqrt((this->x- _coords2.x)*(this->x- _coords2.x)
                      +(this->y- _coords2.y)*(this->y- _coords2.y)
                      +(this->z- _coords2.z)*(this->z- _coords2.z));
    return      acos(((_coords1.x-this->x)*(_coords2.x-this->x)+
                      (_coords1.y-this->y)*(_coords2.y-this->y)+
                      (_coords1.z-this->z)*(_coords2.z-this->z))/A1);
}




/*!
   *  \brief Calculate the normal vector between the caller object and two other Coords objects
   *
   *  Calculate the normal between 2 Coords objects and the caller. The caller object is the center of the normal vector
   *  \param B : First Coords object
   *  \param C : Second Coords object
   *  \param res : The normal vector is return in the parameter
   */

Coords  Coords::getNormal(const Coords &B, const Coords& C) const
{
            const Coords V1(B.x-x,B.y-y,B.z-z);
            const Coords V2(C.x-x,C.y-y,C.z-z);
                Coords res( (V1.y*V2.z-V1.z*V2.y),
                            (V1.z*V2.x-V1.x*V2.z),
                            (V1.x*V2.y-V1.y*V2.x));
            const double norm = sqrt(res.x*res.x+ res.y*res.y+res.z*res.z);
            res.setCoords(res.x/norm+x,res.y/norm+y,res.z/norm+z);
            return res;
}


/*!
       *  \brief Calculate the distance with another Coords object
       *  Calculate the distance between two coords objects.
       *  \param _coords : The other Coords object
       *  \return  the distance between the two objects
       */
double  Coords::calcDist (const Coords &_coords) const {
  const Coords Xd(x-_coords.x,y-_coords.y,z-_coords.z); 
  return sqrt(Xd.x*Xd.x + Xd.y*Xd.y + Xd.z*Xd.z);
}

// double Coords::calcDist(const Coords &_coords) const {

//     if (this == nullptr) {

//         std::cerr << "Error: 'this' is a null pointer!" << std::endl;

//         return -1;

//     }

//     // Assuming _coords is passed by reference, lets check its validity (printing vlaues)

//     std::cout << "Coords being compared - This: (" << this->x << ", " << this->y << ", " << this->z

//               << "), Other: (" << _coords.x << ", " << _coords.y << ", " << _coords.z << ")" << std::endl;

//     const Coords Xd(x - _coords.x, y - _coords.y, z - _coords.z);

//     return sqrt(Xd.x * Xd.x + Xd.y * Xd.y + Xd.z * Xd.z);

// }


double Coords::calcDist(const Coords &_coords, const double& max) const {
   double dist = 0;
   double diff = this->x - _coords.x;

   if (diff > max || diff < -max) 
      return max;

   dist = diff*diff;
   diff = this->y - _coords.y;

   if (diff > max || diff < -max) 
      return max;

   dist += diff*diff;
   diff = this->z- _coords.z;

   if (diff > max || diff < -max)
      return max;

   dist += diff*diff;

   return sqrt(dist);
}




std::string Coords::toString() const
{
    std::ostringstream oss;
    oss << "(x= " << x << " ; y=" << y << " ; z=" << z <<")";
    return oss.str();
}
