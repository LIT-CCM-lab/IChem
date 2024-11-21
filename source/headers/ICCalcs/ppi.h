#ifndef PPI_H
#define PPI_H

#include "headers/ICMole/global.h"
#include "headers/ICMole/coords.h"
#include "headers/ICCalcs/interaction.h"
#include "headers/ICMole/complex.h"
#include "headers/ICCalcs/volsite.h"


namespace ICMole{

class PPI
{

public:

    struct Interaction
    {
        Atom *atomA;
        Atom *atomB;
        double dist;

    };
    struct Interface{
        std::string chain1;
        std::string chain2;
        int size;
        Interface( std::string chainA, std::string chainB, int sizeT):chain1(chainA),chain2(chainB),size(sizeT) {}
    };

    struct InterDescriptor{

//% n'est plus utilise dans la nouvelle version du code.
//        InterDescriptor();
        InterDescriptor();

        ChainList vecChain;
        Molecule Inter;
        Interactions Ints;

        InterResults results;

        double area;
        Coords barycentre;
        Coords posmin;
        Coords posmax;
        inline   void  calcArea() { area =0; } //% formule de calcul de l'air de contact

        inline Chain* getchain(int pos) { return vecChain[pos];}
        inline Molecule& getInterMol() {return Inter; }
        inline Interactions& getInterClass() { return Ints; }
        inline InterResults& getResult() {return results; }
        void calculcoords();
        void calculcoords(int i);
        void calculcoordsintermol();
        void calculbari();
        void calculbari(int i);
        //    inline   void  setPosMin(const double x,const double y,const double z){posmin(x,y,z);}
        //    inline   void  setPosMax(const double x,const double y,const double z){posmax(x,y,z);}

    };


    Complex &cp_orig;
    PPI(Complex &cp,std::string name);
    ////////////////////////////////////////////////
    ////////////////////////////////////////////////
    ////////////////// Function ////////////////////
    ////////////////////////////////////////////////
    ////////////////////////////////////////////////


    void mapinitialisation();
    void mapdetection();
    void makeinterface();
    void calcenfouiss(bool outputall);
    void writeinterface(PPI::InterDescriptor& inter);
    void detectcavity(bool outputall);
    void detectPPi(std::vector<std::string> &vecChain, int nchaint, bool outputall = false);
    void defPropsppi(Molecule& protein);
    bool appartient( std::string name, const std::string& liste);

private:

    std::vector<PPI::Interface*> intersize;
    double interfacetype[45];
    std::vector<ICMole::Interactions*> interfaces;
    std::vector<ICMole::InterResults*> scores;
    std::string PDB_name;
    std::map<std::string,std::map<std::string,std::vector<PPI::Interaction> > >mapinteraction;


};

}
#endif // PPI_H
