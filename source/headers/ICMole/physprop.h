#ifndef PHYSPROP_H
#define PHYSPROP_H

#include <string>

namespace ICMole
{
/*! \class PhysProp
  * \brief Representation of physico-chemical properties
*/
class PhysProp
{
private:
    bool aromatic;
    bool acceptor;
    bool donor;
    bool cation;
    bool anion;
    bool hydrophobic;
    bool apolar;
    bool metalA;
    bool metal;
    bool dummy;
    bool ring;
    bool substituent;
    bool scaffold;
    bool linker;
    bool weakDonor;
    bool weakAcceptor;
    bool inProtein;    /*!< \brief Has an atom center within the box (former isProt) */
    bool inLigand;     /*!< \brief Has a ligand atom center within the box (former isSit)*/
    bool closeProtein; /*!< \brief Is this within vdW radius of a protein atom ? (former Out)*/
    bool inCavity;     /*!< \brief Is this within the site ? */
    bool outCavity;
    bool noInfo;


public:
    PhysProp();
    inline const bool& isAromatic()        const { return aromatic;    }
    inline const bool& isAcceptor()        const { return acceptor;    }
    inline const bool& isDonor()           const { return donor;       }
    inline const bool& isCation()          const { return cation;      }
    inline const bool& isAnion()           const { return anion;       }
    inline const bool& isHydrophobic()     const { return hydrophobic; }
    inline const bool& isApolar()          const { return apolar;      }
    inline const bool& isMetalA()          const { return metalA;      }
    inline const bool& isMetal()           const { return metal;       }
    inline const bool& isScaffold()        const { return scaffold;    }
    inline const bool& isInRing()          const { return ring;        }
    inline const bool& isNoInfo()          const { return noInfo;      }
    inline const bool& isweakDonor()       const { return weakDonor;   }
    inline const bool& isweakAcceptor()    const { return weakAcceptor;}
    inline const bool& isInProtein()       const { return inProtein;   }
    inline const bool& isInLigand()        const { return inLigand;    }
    inline const bool& isCloseProtein()    const { return closeProtein;}
    inline const bool& isInCavity()        const { return inCavity;    }
    inline const bool& isOutCavity()       const { return outCavity;   }
    inline const bool& isSubstituent()     const { return substituent; }
    inline const bool& isLinker()          const { return linker;      }
    inline const bool& isDummy()           const { return dummy;       }


    inline void setAcceptor     (const bool& prop) {    acceptor=prop;  }
    inline void setDonor        (const bool& prop) {       donor=prop;  }
    inline void setCation       (const bool& prop) {      cation=prop;  }
    inline void setAnion        (const bool& prop) {       anion=prop;  }
    inline void setHydrophobic  (const bool& prop) { hydrophobic=prop;  }
    inline void setApolar       (const bool& prop) {      apolar=prop;  }
    inline void setMetal        (const bool& prop) {       metal=prop;  }
    inline void setMetalA       (const bool& prop) {      metalA=prop;  }
    inline void setDummy        (const bool& prop) {       dummy=prop;  }
    inline void setRing         (const bool& prop) {        ring=prop;  }
    inline void setAromatic     (const bool& prop) {    aromatic=prop;  }
    inline void setSubstituent  (const bool& prop) { substituent=prop;  }
    inline void setScaffold     (const bool& prop) {    scaffold=prop;  }
    inline void setLinker       (const bool& prop) {      linker=prop;  }
    inline void setWeakDonor    (const bool& prop) {   weakDonor=prop;  }
    inline void setWeakAcceptor (const bool& prop) {weakAcceptor=prop;  }

    inline void setInLigand     (const bool& prop) {    inLigand=prop;  }
    inline void setInProtein    (const bool& prop) {   inProtein=prop;  }
    inline void setCloseProtein (const bool& prop) {closeProtein=prop;  }
    inline void setInCavity     (const bool& prop) {    inCavity=prop;  }
    inline void setOutCavity    (const bool& prop) {   outCavity=prop;  }
    inline void setNoInfo       (const bool& prop) {      noInfo=prop;  }

    PhysProp& operator=  (const PhysProp& _prop);/*!< \brief Copy constructor */
    bool operator== (const PhysProp& _prop);
    const std::string toString() const; /*!< \brief Return a string with all properties which are set to True */
    const std::string toFragType() const;
    const std::string toCavString() const;
    void clear();/*!< \brief Set all properties to false */
};

}

#endif // PHYSPROP_H
