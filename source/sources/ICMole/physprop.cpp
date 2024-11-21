#include "headers/ICMole/physprop.h"
#include "headers/ICMole/global.h"
using namespace ICMole;
using namespace std;

PhysProp::PhysProp():
    aromatic(false),    acceptor(false),
    donor(false),      cation(false),
    anion(false), hydrophobic(false),
    metal(false),       dummy(false),
    ring(false),substituent(false),    scaffold(false),
    linker(false),weakDonor(false),weakAcceptor(false),inProtein(false),
    inLigand(false),closeProtein(false),
    inCavity(false),   outCavity(false),noInfo(false)

{}



PhysProp& PhysProp::operator=  (const PhysProp& _prop)
{
    aromatic     = _prop.aromatic;
    acceptor     = _prop.acceptor;
    donor        = _prop.donor;
    cation       = _prop.cation;
    anion        = _prop.anion;
    metal        = _prop.metal;
    dummy        = _prop.dummy;
    ring         = _prop.ring;
    hydrophobic  = _prop.hydrophobic;
    inProtein    = _prop.inProtein;
    closeProtein = _prop.closeProtein;
    inCavity     = _prop.inCavity;
    outCavity    = _prop.outCavity;
    substituent  = _prop.substituent;
    scaffold     = _prop.scaffold;
    linker       = _prop.linker;
    noInfo       = _prop.noInfo;
    weakDonor    = _prop.weakDonor;
    weakAcceptor = _prop.weakAcceptor;
    inLigand     = _prop.inLigand;


    return *this;
}

bool PhysProp::operator==  (const PhysProp& _prop)
{
  if (aromatic   == _prop.aromatic &&
    acceptor     == _prop.acceptor &&
    donor        == _prop.donor &&
    cation       == _prop.cation &&
    anion        == _prop.anion &&
    metal        == _prop.metal &&
    dummy        == _prop.dummy &&
    ring         == _prop.ring &&
    hydrophobic  == _prop.hydrophobic &&
    inProtein    == _prop.inProtein &&
    closeProtein == _prop.closeProtein &&
    inCavity     == _prop.inCavity &&
    outCavity    == _prop.outCavity &&
    substituent  == _prop.substituent &&
    scaffold     == _prop.scaffold &&
    linker       == _prop.linker &&
    noInfo       == _prop.noInfo &&
    weakDonor    == _prop.weakDonor &&
    weakAcceptor == _prop.weakAcceptor &&
    inLigand     == _prop.inLigand )
  {

 return true;
  }else
      return false;
}

void PhysProp::clear()
{
    aromatic     = false;
    apolar       = false;
    acceptor     = false;
    donor        = false;
    cation       = false;
    anion        = false;
    metal        = false;
    dummy        = false;
    ring         = false;
    hydrophobic  = false;
    inProtein    = false;
    closeProtein = false;
    inCavity     = false;
    outCavity    = false;
    substituent  = false;
    scaffold     = false;
    linker       = false;
    noInfo       = false;
    weakDonor    = false;
    weakAcceptor = false;
    inLigand     = false;


}

const std::string PhysProp::toCavString() const
{
    ostringstream oss;
    if (noInfo) oss << "No Info|";
    if (inProtein)oss << "Protein|";
    if (closeProtein)oss <<"Close|";
    if (inCavity)oss <<"Cavity|";
    if (outCavity)oss <<"Out Cavity|";
    return oss.str();
}

const std::string PhysProp::toString() const
{
    ostringstream oss;
    if (hydrophobic) oss << "Hydrophobic|";
    if (apolar) oss <<"Apolar|";
    if (aromatic) oss << "Aromatic|";
    if (acceptor) oss << "Acceptor|";
    if (donor) oss << "Donor|";
    if (cation) oss << "Cation|";
    if (anion) oss << "Anion|";
    if (metal) oss << "Metal|";
    if (metalA) oss << "MetalAcceptor|";
    if (ring) oss << "In Ring|";
//    if (weakAcceptor) oss << "Weak Acceptor|";
//    if (weakDonor) oss << "Weak Donor|";
    return oss.str().substr(0,oss.str().length()-1);

}
const std::string PhysProp::toFragType() const
{
    if (substituent) return "Substituent";
    if (scaffold)    return "Scaffold";
    if (linker)      return "Linker";
    return "";
}



