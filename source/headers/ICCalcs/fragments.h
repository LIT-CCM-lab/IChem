#ifndef FRAGMENTS_H
#define FRAGMENTS_H

#include "headers/ICCalcs/interaction.h"
namespace ICMole
{



class FragGenerator
{
    struct Fragment
    {
       AtomList atoms;
       BondList anchor;
       bool wCycle;
       std::vector<unsigned int> anchorID;
       unsigned short id;
       InterResults inters;

    };
private:
    Fragment* listFragments;
    Complex& complex;
    unsigned int NCL;
    std::map<unsigned int,unsigned int> NCLUS;
    BondList anchorBond;
public:
    FragGenerator(Complex &complex);
    void calcInteractions();
    ~FragGenerator();

    void fragmentLigand();
    void saveFragments(const std::string& name="");







};
}

#endif // FRAGMENTS_H
