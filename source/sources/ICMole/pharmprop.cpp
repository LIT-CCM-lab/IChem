#include "headers/ICMole/pharmprop.h"
#include "headers/ICMole/global.h"


using namespace ICMole;
using namespace std;


PharmProp::PharmProp(const unsigned int PharType, const Coords base, const Coords c1, const Coords c3, int tol2)
{
    TypePha=PharType;
    c_base=base;
    c_1A=c1;
    c_3A=c3;
    tol=tol2;
}

