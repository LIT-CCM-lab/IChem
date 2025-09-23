#include <iostream>
#include "headers/ICCalcs/interaction.h"
#include "headers/ICMole/complex.h"
#include "headers/ICMole/box.h"
#include "headers/ICMole/cycle.h"
#include "headers/ICCalcs/volsite.h"
#include "headers/ICPars/molereader.h"

// #include "headers/ICCalcs/neighborsearch.hpp"


//#define ICHEM_DEBUG
using namespace std;
using namespace ICMole;



double *Interactions::path_dist;
unsigned int Interactions::path_dist_size;
unsigned int Interactions::size_triplets;
bool Interactions::Load_triplet=false;
int Interactions::vect_list[1000][1000];


Interactions::Interactions(ICMole::Complex &cp) throw(ICMole::MoleExcept): complex(cp),grid(complex.genGrid(4.5)) {
    if (cp.getMole(MoleType::PROTEIN)== (Molecule*)NULL)
        throw MoleExcept(3020101, "Interactions::Interactions", "No protein found in complex");

    for(short i = 0; i < NB_INTTYPE; ++i) {
        wInterType[i] = true;
    }
} 

void Interactions::addInteraction(InterResults& interResult, Atom& atomP, Atom& atomL, double dist, int& NInter, double* angle, unsigned int interactionType) const {

    double ang = angle ? *angle : -100000.0; // <optional> C++17 library

    interResult.listInters.emplace_back(NInter++, &atomP, &atomL, (atomL.fixpos + atomP.fixpos) / 2, interactionType, dist, ang);
    
    switch (interactionType) {
        case InterType::HBOND_PROT:
            interResult.N++;
            break;
        case InterType::HBOND_LIG:
            interResult.O++;
            break;
        case InterType::METAL:
            interResult.Zn++;
            break;
        case InterType::HYDROPHOBIC:
            interResult.CA++; 
            break;
        case InterType::ARFACEFACE:
        case InterType::AREDGEFACE:
            interResult.CZ++; 
            break;
        case InterType::IONIC_PROT:
            interResult.NZ++; 
            break;
        case InterType::IONIC_LIG:
            interResult.OD1++; 
            break;
        default: break;
    }
}

// Special case: nitrogen linked to sulfonamide
void Interactions::checkMetalNitrogenSulfonamideCase(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const {
    
    if (atomL.isNitrogen() && dist < params.Dist_Metal) {
        bool linkedToSulfur = false;

        for (size_t i = 0; i < atomL.getNumBond(); ++i) {

            const Atom& linkedAtom = atomL.getAtomLinked(i);
            if (linkedAtom.isSulfur()) {
                linkedToSulfur = true;
                break;
            }
        }
        if (linkedToSulfur) {
            addInteraction(interResult, atomP, atomL, dist, NInter, nullptr, InterType::METAL);
        }
    }
}


void Interactions::checkMetalInteractions(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const {
    
    if(!atomL.props.isAcceptor())
        return;

    if (!wInterType[InterType::METAL] || !atomP.props.isMetal())
        return;

    // Normal case
    if (dist <= params.Dist_Metal && dist >= params.dist_Metal) {
        addInteraction(interResult, atomP, atomL, dist, NInter, nullptr, InterType::METAL);
        return;
    }

    if(  dist < 3.4  && ( atomL.props.isMetal() && (atomP.props.isAcceptor() || atomP.props.isAnion())
        || (atomL.props.isAnion() || atomL.props.isAcceptor() ) && atomP.props.isMetal()
        || (atomP.isMetallic() && atomL.getAtomicName() == "N"))) {

        atomL.props.setMetalA(true);
    }
}


void Interactions::checkHydrogenBondLigandAcceptor(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const  {
    
    if (wInterType[InterType::HBOND_PROT] && atomP.props.isDonor() && dist <= params.Dist_H && dist >= params.dist_H) {
        
        for (size_t i = 0; i < atomP.getNumBond(); ++i) {
            const Atom& atomP2 = atomP.getAtomLinked(i);
            double angle = atomP2.fixpos.calcAngle(atomP.fixpos, atomL.fixpos);

            if (!atomP2.isHydrogen()) 
                continue;
                
            if (angle <= params.Angl_H - params.AngT_H || angle >= params.Angl_H + params.AngT_H) 
                continue;

            addInteraction(interResult, atomP, atomL, dist, NInter, &angle, InterType::HBOND_PROT);
        }
    }
}

void Interactions::checkHydrogenBondLigandDonor(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const {
    
    if (wInterType[InterType::HBOND_LIG] && atomP.props.isAcceptor() && dist <= params.Dist_H && dist >= params.dist_H) {
        for (size_t i=0; i < atomL.getNumBond(); ++i) {

            const Atom &atomL2 = atomL.getAtomLinked(i);
            if (!atomL2.isHydrogen())     
                continue;

            double angle = atomL2.fixpos.calcAngle(atomP.fixpos, atomL.fixpos);

            if (angle <= params.Angl_H - params.AngT_H || angle >= params.Angl_H + params.AngT_H)  
                continue;
            
            addInteraction(interResult, atomP, atomL, dist, NInter, &angle, InterType::HBOND_LIG);
        }
    }
}

void Interactions::checkWeakHydrogenBondLigandAcceptor(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const {
    
    if (wInterType[InterType::WHBOND_PROT] && atomP.props.isweakDonor() && dist <= params.Dist_WHBond) {
        for (size_t i=0; i < atomP.getNumBond(); ++i) {
            
            const Atom &atomP2 = atomP.getAtomLinked(i);
            double angle = atomP2.fixpos.calcAngle(atomP.fixpos,  atomL.fixpos);

            if (!atomP2.isHydrogen())     
                continue;

            if (angle <= params.Angl_WHBond - params.AngT_WHBond || angle >= params.Angl_WHBond + params.AngT_WHBond)  
                continue;

            addInteraction(interResult, atomP, atomL, dist, NInter, &angle, InterType::WHBOND_PROT);
        }
    }
}

void Interactions::checkWeakHydrogenBondLigandWeakAcceptor(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const {
    
    if (atomL.props.isweakAcceptor() && wInterType[InterType::WHBOND_PROT] && dist < params.Dist_WHBond && (atomP.props.isweakDonor() || atomP.props.isDonor())) {
        
        for (size_t i=0; i < atomP.getNumBond(); ++i) {
            const Atom &atomP2 = atomP.getAtomLinked(i);
            double angle = atomP2.fixpos.calcAngle(atomP.fixpos,  atomL.fixpos);

            if (!atomP2.isHydrogen())     
                continue;

            if (angle <= params.Angl_WHBond - params.AngT_WHBond || angle >= params.Angl_WHBond + params.AngT_WHBond)
                continue;

            addInteraction(interResult, atomP, atomL, dist, NInter, &angle, InterType::WHBOND_PROT);
        }
    }
}

void Interactions::checkWeakHydrogenBondLigandDonorProteinWeakAcceptor(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const {

    if (wInterType[InterType::WHBOND_LIG] && atomP.props.isweakAcceptor() && dist <= params.Dist_WHBond) {
        for (size_t i=0; i< atomL.getNumBond(); ++i) {
            const Atom &atomL2 = atomL.getAtomLinked(i);
            double angle = atomL2.fixpos.calcAngle(atomP.fixpos,  atomL.fixpos);

            if (!atomL2.isHydrogen())     
                continue;

            if (angle <= params.Angl_WHBond - params.AngT_WHBond || angle >= params.Angl_WHBond + params.AngT_WHBond) 
                continue;
            
            addInteraction(interResult, atomP, atomL, dist, NInter, &angle, InterType::WHBOND_LIG);
        }
    }
}


void Interactions::checkWeakHydrogenBondLigandWeakDonorProteinAcceptor(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const {
    
    if (!atomL.props.isweakDonor())
        return;
    if (!wInterType[InterType::WHBOND_LIG])
        return;
    if (!(atomP.props.isAcceptor() || atomP.props.isweakAcceptor()))
        return;
    if (dist > params.Dist_WHBond) 
        return;

    for (size_t i = 0; i < atomL.getNumBond(); ++i) {
        const Atom& atomL2 = atomL.getAtomLinked(i);
        double angle = atomL2.fixpos.calcAngle(atomP.fixpos, atomL.fixpos);

        if (!atomL2.isHydrogen()) 
            continue;

        if (angle <= params.Angl_WHBond - params.AngT_WHBond || 
            angle >= params.Angl_WHBond + params.AngT_WHBond) {
            continue;
        }

        addInteraction(interResult, atomP, atomL, dist, NInter, &angle, InterType::WHBOND_LIG);
    }
}


void Interactions::checkIonicProteinInteractions(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const {

    if (!wInterType[InterType::IONIC_PROT])
        return;
    if (!atomL.props.isAnion())
        return;
    if (!atomP.props.isCation() || atomP.props.isMetal())
        return;
    if (dist > params.Dist_Ionic || dist < params.dist_Ionic)
        return;
    
    addInteraction(interResult, atomP, atomL, dist, NInter, nullptr, InterType::IONIC_PROT);
}

void Interactions::checkIonicLigandInteractions(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const {
    
    if (!wInterType[InterType::IONIC_LIG])
        return;
    if (!atomP.props.isAnion())
        return;
    if (dist > params.Dist_Ionic || dist < params.dist_Ionic)
        return;

    addInteraction(interResult, atomP, atomL, dist, NInter, nullptr, InterType::IONIC_LIG);
}

void Interactions::checkPiCationInteraction(Atom& atomL, Atom& atomP, double dist, InterResults& interResult, int& NInter) const {
    if (!wInterType[InterType::PICATION]) 
        return;
    if (atomP.getName() != "DuCy") 
        return;
    if (dist >= params.Dist_PiCation) 
        return;

    Cycle* cycleP = atomP.getParent().getCycleFromCenter(&atomP);
    if (cycleP == nullptr) 
        return;

    cycleP->calcVector();

    double angle = cycleP->getCenter().fixpos.calcAngle(cycleP->getNormVector(), atomL.fixpos);
    if (angle <= params.Angl_PiCation - params.AngT_PiCation || angle >= params.Angl_PiCation + params.AngT_PiCation)
        return;

    addInteraction(interResult, cycleP->getCenter(), atomL, dist, NInter, &angle, InterType::PICATION);
}


void Interactions::detectInteractions(Molecule& ligand, InterResults& interResult, bool wMerge, bool oldh) const {
    
    double dist;
    int NInter = 0;
    double min_allowed_dist = 1.5;
    double max_allowed_dist = std::max({ params.Dist_H, params.Dist_Hyd, params.Dist_Ionic, params.Dist_Metal, params.Dist_Arom, params.Dist_PiCation });

    
    // Build KD-tree on protein atoms
    NeighborSearch neighborSearch;
    std::vector<std::array<double,3>> proteinPoints;
    std::vector<Atom*> proteinAtoms;

    
    for (auto itM = complex.firstMole(); itM != complex.lastMole(); ++itM) {
        Molecule* mol = *itM;

        for(auto itC = mol->firstCycle(); itC != mol->lastCycle(); ++itC) {
            (*itC)->getFixpos();
        }
    }
    
    // Helper to add an atom
    auto push = [&](Atom* a) { 
        proteinPoints.push_back({ a->fixpos.x, a->fixpos.y, a->fixpos.z });
        proteinAtoms.push_back(a);
    };

    for (ItCAtom it = complex.firstAtom(); it != complex.lastAtom(); ++it) {
        Atom* a = *it;
        if (!a->isHydrogen() && a->isUsed() && a->getParent().getMoleType() == MoleType::PROTEIN)
            push(a);
    }

    // Push manually the dummy atoms
    for (auto m = complex.firstMole(MoleType::PROTEIN); m != complex.lastMole(MoleType::PROTEIN); ++m) {
        
        Molecule* prot = *m;
        for (auto itC = prot->firstCycle(); itC != prot->lastCycle(); ++itC) {
            Cycle* cy = *itC;
            if (!cy->isAromatic()) 
                continue;   

            cy->getFixpos();                       
            Atom* ctr = &cy->getCenter();          
            push(ctr);
        }
    }

    neighborSearch.build(proteinPoints);

    // Prepare ligand points    
    std::vector<std::array<double,3>> ligandPoints;
    std::vector<Atom*> ligandAtoms;

    for (ItCAtom itA = ligand.firstAtom(); itA != ligand.lastAtom(); ++itA) {
        Atom* a = *itA;
        if (a->isUsed() && !a->isHydrogen() && a->getName() != "DuCy") {
            ligandPoints.push_back({ a->fixpos.x, a->fixpos.y, a->fixpos.z });
            ligandAtoms.push_back(a);
        }
    }

    // Query all ligand points at once
    auto pairs = neighborSearch.query(ligandPoints, max_allowed_dist);

    // Sort by ligand_idx so that contacts of the same ligand atom
    // are consecutive so we can clear hydlist at each change
    std::sort(pairs.begin(), pairs.end(), [](const NeighborSearch::Contact& a, const NeighborSearch::Contact& b) {
        return a.ligand_idx < b.ligand_idx;
    });

    // Helper to clear hydrophobic list
    auto flushHydList = [&](std::map<Residu*,resbest>& hydlist) {
        for (auto& [res, rb] : hydlist) {
            if (!rb.atmP || !rb.atmL) 
                continue;
            if (rb.atmL->getResidu()->getIdentifier() == rb.atmP->getResidu()->getIdentifier())
                continue;

            addInteraction(interResult, *rb.atmP, *rb.atmL, rb.dist, NInter, nullptr, InterType::HYDROPHOBIC);
        }
        hydlist.clear();
    };

    std::map<Residu*,resbest> hydlist;
    int currentLigand = -1;


    for (const auto& pair : pairs) {
        // We can the new ligand atom we flush previous hydlist
        if (pair.ligand_idx != currentLigand) {
            flushHydList(hydlist);
            currentLigand = pair.ligand_idx;
        }

        Atom& atomL = *ligandAtoms [pair.ligand_idx];
        Atom& atomP = *proteinAtoms[pair.protein_idx];
        dist = std::sqrt(pair.distance_squared);

        if (&atomP.getParent() == &ligand || !atomP.isUsed() || atomP.isHydrogen())
            continue;
        if (dist > max_allowed_dist)
            continue;

        // All interactions check
        checkMetalNitrogenSulfonamideCase(atomL, atomP, dist, interResult, NInter);
        checkMetalInteractions(atomL, atomP, dist, interResult, NInter);

        if (atomL.props.isAcceptor()) {
            checkHydrogenBondLigandAcceptor(atomL, atomP, dist, interResult, NInter);
            checkWeakHydrogenBondLigandAcceptor(atomL, atomP, dist, interResult, NInter);
        }
        checkWeakHydrogenBondLigandWeakAcceptor(atomL, atomP, dist, interResult, NInter);

        if (atomL.props.isDonor()) {
            checkHydrogenBondLigandDonor(atomL, atomP, dist, interResult, NInter);
            checkWeakHydrogenBondLigandDonorProteinWeakAcceptor(atomL, atomP, dist, interResult, NInter);
        }
        checkWeakHydrogenBondLigandWeakDonorProteinAcceptor(atomL, atomP, dist, interResult, NInter);

        checkIonicProteinInteractions(atomL, atomP, dist, interResult, NInter);

        if (atomL.props.isCation()) {
            checkIonicLigandInteractions(atomL, atomP, dist, interResult, NInter);
            checkPiCationInteraction(atomL, atomP, dist, interResult, NInter);
        }
        
        processHydrophobicInteraction(atomL, atomP, dist, oldh, interResult, hydlist, neighborSearch, proteinAtoms);
    }

    // flush the last ligand atom’s hydrophobics
    flushHydList(hydlist);

    if (wInterType[InterType::HYDROPHOBIC]) {
        auto* protein = complex.getMole(MoleType::PROTEIN);
        if (protein) {
            checkAromaticHydrophobicInteractions(ligand, *protein, interResult, NInter, params.dist_Hyd, params.Dist_Hyd);
        }
    }

    processAromaticInteractions(ligand, neighborSearch, proteinAtoms, max_allowed_dist, interResult, wInterType, NInter, min_allowed_dist);

    if (wMerge) {
        mergeInteractions(interResult);
    }
}


void Interactions::checkAromaticHydrophobicInteractions(Molecule& ligand, Molecule& protein, InterResults& interResult, int& NInter, double dist_H, double Dist_H) const {
    
    for (ItCCycle ligandCycle = ligand.firstCycle(); ligandCycle != ligand.lastCycle(); ++ligandCycle) {
        Cycle& ligCycle = **ligandCycle;
        if (!ligCycle.isAromatic()) 
            continue;

        for (ItCCycle proteinCycle = protein.firstCycle(); proteinCycle != protein.lastCycle(); ++proteinCycle) {
            Cycle& protCycle = **proteinCycle;
            if (!protCycle.isAromatic())
                continue;

            for (size_t i = 0; i < ligCycle.getNumAtom(); ++i) {
                Atom& atomL = *ligCycle.getAtom(i);
                if (!atomL.props.isHydrophobic()) 
                    continue;

                for (size_t j = 0; j < protCycle.getNumAtom(); ++j) {
                    Atom& atomP = *protCycle.getAtom(j);
                    if (!atomP.props.isHydrophobic()) 
                        continue;

                    double distance = atomL.calcFixpos(atomP);
                    if (distance >= dist_H && distance <= Dist_H) {
                        addInteraction(interResult, atomP, atomL, distance, NInter, nullptr, InterType::HYDROPHOBIC);
                    }
                }
            }
        }
    }
}


void Interactions::processAromaticInteractions(Molecule& ligand, NeighborSearch& neighborSearch, const std::vector<Atom*>& proteinAtoms, double max_allowed_dist, InterResults& interResult, const bool* wInterType, int& NInter, double min_allowed_dist) const {
    
    using Contact = NeighborSearch::Contact;

    for (ItCCycle itLC = ligand.firstCycle(); itLC != ligand.lastCycle(); ++itLC) {
        Cycle& ligCycle = **itLC;
        if (!ligCycle.isAromatic())
            continue;

        Atom& ligCtr = ligCycle.getCenter();
        ligCycle.calcVector();

        // KD-tree query around ligand ring center
        std::array<double,3> q{{ligCtr.fixpos.x, ligCtr.fixpos.y, ligCtr.fixpos.z}};
        double searchRadius = std::max(1.5 * max_allowed_dist, 10.0);
        std::vector<Contact> hits = neighborSearch.query({q}, searchRadius);

        // We collect unique protein rings
        std::vector<Cycle*> protCycles; 

        // We examine every hits
        for (const Contact& h : hits) {

            Atom& prot = *proteinAtoms[h.protein_idx];
            double dist = std::sqrt(h.distance_squared);

            // Not the same molecule
            if (&prot.getParent() == &ligand)
                continue;

            // PI-cation check: to refactor into its own function
            if (wInterType[InterType::PICATION] && prot.props.isCation() && dist <= params.Dist_PiCation && dist >= min_allowed_dist) {
                double ang = ligCtr.fixpos.calcAngle(ligCycle.getNormVector(), prot.fixpos);

                bool ok = (ang >= params.Angl_PiCation - params.AngT_PiCation && ang <= params.Angl_PiCation + params.AngT_PiCation) || (ang >= params.Angl_PiCation - params.AngT_PiCation - M_PI && ang <= params.Angl_PiCation + params.AngT_PiCation - M_PI);

                if (ok)
                    addInteraction(interResult, prot, ligCtr, dist, NInter, &ang, InterType::PICATION);
            }

            // Aromatic center vector protCycles (DuAr)
            if (prot.getName() == "DuAr" && dist >= min_allowed_dist && dist <= 10.0) {
                if (Cycle* c = prot.getParent().getCycleFromCenter(&prot))
                    protCycles.push_back(c);
            }
        }

        // Filter to unique protein cycles
        std::sort(protCycles.begin(), protCycles.end());
        protCycles.erase(std::unique(protCycles.begin(), protCycles.end()), protCycles.end());

        // Check for interactions with ligand/protein ring pair
        for (Cycle* cp : protCycles) {
            
            Cycle& cycleP = *cp;
            double centerDist = cycleP.getFixpos().calcDist(ligCycle.getFixpos());
            if(centerDist > max_allowed_dist)
                continue;

            // Face-face / Edge-face
            cycleP.calcVector();

            bool face = false; 
            int edge = 0;
            bool abort = false;

            for (size_t i = 0; i < ligCycle.getNumAtom() && !abort; ++i) {
                
                Atom& aL = *ligCycle.getAtom(i);
                if (!aL.isCarbon())
                    continue;
                
                const Coords& nL = ligCycle.getNormVector(i);

                for (size_t j = 0; j < cycleP.getNumAtom(); ++j) {
                    
                    Atom& aP = *cycleP.getAtom(j);
                    if (!aP.isCarbon()) 
                        continue;

                    double d = aL.calcFixpos(aP);
                    
                    if (d > params.Dist_Arom * 3) { 
                        abort = true; 
                        break; 
                    }

                    const Coords& nP = cycleP.getNormVector(j);
                    double ang = acos((aL.fixpos.x - nL.x) * (aP.fixpos.x - nP.x) +
                                      (aL.fixpos.y - nL.y) * (aP.fixpos.y - nP.y) +
                                      (aL.fixpos.z - nL.z) * (aP.fixpos.z - nP.z));

                    if (d < params.Dist_Arom && d > params.dist_Arom) {
                        // face-face
                        if ((ang > params.Angl_AromFF - params.AngT_AromFF &&
                             ang < params.Angl_AromFF + params.AngT_AromFF) ||
                            (ang > M_PI + params.Angl_AromFF - params.AngT_AromFF &&
                             ang < M_PI + params.Angl_AromFF + params.AngT_AromFF) ||
                            (ang > params.Angl_AromFF - params.AngT_AromFF - M_PI &&
                             ang < params.Angl_AromFF + params.AngT_AromFF - M_PI))
                            face = true;

                        // edge-face
                        if ((ang > params.Angl_AromEF - params.AngT_AromEF &&
                             ang < params.Angl_AromEF + params.AngT_AromEF) ||
                            (ang > M_PI + params.Angl_AromEF - params.AngT_AromEF &&
                             ang < M_PI + params.Angl_AromEF + params.AngT_AromEF) ||
                            (ang > params.Angl_AromEF - params.AngT_AromEF - M_PI &&
                             ang < params.Angl_AromEF + params.AngT_AromEF - M_PI))
                            ++edge;
                    }
                }
            }

            if (face && wInterType[InterType::ARFACEFACE])
                addInteraction(interResult, cycleP.getCenter(), ligCycle.getCenter(), centerDist, NInter, nullptr, InterType::ARFACEFACE);
            else if (edge > 5 && wInterType[InterType::AREDGEFACE])
                addInteraction(interResult, cycleP.getCenter(), ligCycle.getCenter(), centerDist, NInter, nullptr, InterType::AREDGEFACE);
        }
    }
}


void Interactions::processHydrophobicInteraction(Atom& atomL, Atom& atomP, double dist, bool oldh, InterResults& interResult, std::map<Residu*, ICMole::resbest>& hydlist, NeighborSearch& neighborSearch, const std::vector<Atom*>& proteinAtoms) const {
    
    if (!wInterType[InterType::HYDROPHOBIC])
        return;
    if (!atomL.props.isHydrophobic() || !atomP.props.isHydrophobic())
        return;
    if (atomL.props.isAromatic() && atomP.props.isAromatic())
        return;
    if (dist < params.dist_Hyd || dist > params.Dist_Hyd)
        return;

    
    if (!oldh) {
        std::array<double,3> q {{ atomP.fixpos.x, atomP.fixpos.y, atomP.fixpos.z }};
        std::vector<std::array<double,3>> query {q};

        const double density_radius = 4.5;
        auto neighbours = neighborSearch.query(query, density_radius);

        int nbatm = 0, nbhyd = 0;
        for (const auto& hit : neighbours) {
            Atom* nb = proteinAtoms[ hit.protein_idx ];
            if (nb->isHydrogen())
                continue;

            ++nbatm;
            if (nb->props.isHydrophobic()) 
                ++nbhyd;
        }

        if (nbatm == 0 || (100.0 * nbhyd / nbatm) <= 50.0)
            return;
    }

    Residu* res = atomP.getResidu();
    auto it = hydlist.find(res);

    if (it == hydlist.end()) {
        ICMole::resbest r;
        r.atmP = &atomP;
        r.atmL = &atomL;
        r.dist = dist;
        r.id   = 0;
        hydlist.emplace(res, r);
    }

    else if (dist < it->second.dist) {
        it->second.atmP = &atomP;
        it->second.atmL = &atomL;
        it->second.dist = dist;
    }
}


/*
    Merge close proximity interactions within a molecule for optimization purposes.
    It looks for hydrophobic interactions between atoms in the molecule and combine them if they are close
*/
void Interactions::mergeInteractions(InterResults& interResult) const
{
    int NInter;
    Coords new_center;
    Atom *new_atmP = nullptr, *new_atmL = nullptr;
    bool modif = true ; double dist;
    
    while (modif)
    {
        //        cout << interResult.listInters.max_size() << endl;
        modif = false;
        // tmp vector for new interpoints
        std::vector<InterPoint> new_interpoints;
        for (size_t inter_i=0; inter_i < interResult.listInters.size();++inter_i)
        {
            InterPoint &intPi = interResult.listInters.at(inter_i);
            //            cout << interResult.listInters.size()  << " work on : " << intPi.point << endl;

            if (intPi.interaction != InterType::HYDROPHOBIC
                    ||  intPi.merged_to != -1) continue;

            const size_t sizeList=interResult.listInters.size();

            for (size_t inter_j=inter_i+1; inter_j < sizeList;++inter_j)
            {
                InterPoint &intPj = interResult.listInters.at(inter_j);

                if (intPj.interaction != InterType::HYDROPHOBIC
                        ||  intPj.merged_to != -1 || intPi.point==intPj.point) continue;
                dist = intPi.center.calcDist(intPj.center,1.1);
                if ( dist > 1 ) continue;
                //                if (interResult.listInters.size() == 2048) {
                //                    cout << "maximum size"  << endl;
                //                }

                modif=true;
                NInter = (int)interResult.listInters.size();
                new_center=(intPj.center+intPi.center)/2;

                if (new_center.calcDist(intPi.Prot_Ref->fixpos) <
                        new_center.calcDist(intPj.Prot_Ref->fixpos) )
                {
                    new_atmP=intPi.Prot_Ref;
                }else new_atmP=intPj.Prot_Ref;

                if (new_center.calcDist(intPi.Lig_Ref->fixpos) < new_center.calcDist(intPj.Lig_Ref->fixpos ))
                {
                    new_atmL=intPi.Lig_Ref;
                }else new_atmL=intPj.Lig_Ref;
                InterPoint IntP(NInter,
                                new_atmP,
                                new_atmL,
                                new_center,
                                InterType::HYDROPHOBIC,(intPi.dist<intPj.dist)? intPi.dist:intPj.dist);

                new_interpoints.push_back(IntP);
                if (intPi.merged_to != -1)
                {
                    for (size_t inter_k=0; inter_k < interResult.listInters.size();++inter_k)
                    {
                        InterPoint &intPk = interResult.listInters.at(inter_k);
                        if (intPk.merged_to == intPi.merged_to) intPk.merged_to=NInter;

                    }
                }
                intPi.merged_to=NInter;
                if (intPj.merged_to != -1)
                {
                    for (size_t inter_k=0; inter_k < interResult.listInters.size();++inter_k)
                    {
                        InterPoint &intPk = interResult.listInters.at(inter_k);
                        if (intPk.merged_to == intPj.merged_to) intPk.merged_to=NInter;

                    }
                }
                intPj.merged_to=NInter;
                NInter++;

            } //END inter_j
        } //END inter_i
        interResult.listInters.insert(interResult.listInters.end(), new_interpoints.begin(), new_interpoints.end());
    } //END WHILE
}



    //    double bestdist;InterPoint *toMove;
    //    for (size_t inter_i=0; inter_i < interResult.listInters.size();++inter_i)
    //    {
    //        InterPoint &intPi = interResult.listInters.at(inter_i);
    //        bestdist=100;
    //        if (intPi.interaction != InterType::IONIC_PROT
    //                &&intPi.interaction != InterType::IONIC_LIG) continue;

    //        for (size_t inter_j=inter_i+1; inter_j < interResult.listInters.size();++inter_j)
    //        {
    //            InterPoint &intPj = interResult.listInters.at(inter_j);
    //            if ((intPi.interaction == InterType::IONIC_PROT
    //                 && intPj.interaction != InterType::HBOND_PROT)
    //                    ||(intPi.interaction == InterType::IONIC_LIG
    //                       && intPj.interaction != InterType::HBOND_LIG) )continue;

    //            dist = intPi.center.calcDist(intPj.center,1.1);
    //            if (dist < bestdist) {bestdist = dist;toMove=&intPj;}


    //        }
    //        if (bestdist <1)
    //        {

    //#ifdef ICHEM_DEBUG
    //                cout << "#########################################"<<endl;
    //                if (intPi.interaction == InterType::IONIC_PROT)
    //                    cout << "MOVING IONIC PROT CATION TO H-BOND PROT DONOR POSITION"<<endl;
    //                else
    //                    cout << "MOVING IONIC LIG CATION TO H-BOND LIG DONOR POSITION"<<endl;
    //                cout << "Ionic pos: "
    //                     <<intPi.Lig_Ref->getIdentifier()
    //                    << "  "
    //                    << intPi.Prot_Ref->getIdentifier()
    //                    <<  intPi.center.toString()<<endl
    //                     << "HBond pos: "
    //                     <<toMove->Lig_Ref->getIdentifier()
    //                    << "  "
    //                    << toMove->Prot_Ref->getIdentifier()
    //                    << toMove->center.toString()<<endl
    //                       ;
    //#endif
    //            intPi.center = toMove->center;

    //        }
    //    }

// }

void Interactions::mergeSpeInts(InterResults& interResult)const
{
    int NInter;
    Coords new_center;
    Atom *new_atmP, *new_atmL;
    bool modif=true;double dist;

    while (modif)
    {
        modif =false;
        for (size_t inter_i=0; inter_i < interResult.listInters.size();++inter_i)
        {
            InterPoint &intPi = interResult.listInters.at(inter_i);
            if (intPi.interaction != InterType::HYDROPHOBIC
                    ||  intPi.merged_to != -1)continue;

            const size_t sizeList=interResult.listInters.size();
            for (size_t inter_j=inter_i+1; inter_j < sizeList;++inter_j)
            {
                InterPoint &intPj = interResult.listInters.at(inter_j);
                if (intPj.interaction != InterType::HYDROPHOBIC
                        ||  intPj.merged_to != -1 || intPi.point==intPj.point)continue;
                dist =intPi.center.calcDist(intPj.center,1.1);
                if (dist > 0.1 )continue;


                modif=true;
                NInter = (int)interResult.listInters.size();
                new_center=(intPj.center+intPi.center)/2;

                if (new_center.calcDist(intPi.Prot_Ref->fixpos) <
                        new_center.calcDist(intPj.Prot_Ref->fixpos) )
                {
                    new_atmP=intPi.Prot_Ref;
                }else new_atmP=intPj.Prot_Ref;

                if (new_center.calcDist(intPi.Lig_Ref->fixpos) < new_center.calcDist(intPj.Lig_Ref->fixpos ))
                {
                    new_atmL=intPi.Lig_Ref;
                }else new_atmL=intPj.Lig_Ref;
                InterPoint IntP(NInter,
                                new_atmP,
                                new_atmL,
                                new_center,
                                InterType::HYDROPHOBIC,(intPi.dist<intPj.dist)? intPi.dist:intPj.dist);

                interResult.listInters.push_back(IntP);
                if (intPi.merged_to != -1)
                {
                    for (size_t inter_k=0; inter_k < interResult.listInters.size();++inter_k)
                    {
                        InterPoint &intPk = interResult.listInters.at(inter_k);
                        if (intPk.merged_to == intPi.merged_to) intPk.merged_to=NInter;

                    }
                }
                intPi.merged_to=NInter;
                if (intPj.merged_to != -1)
                {
                    for (size_t inter_k=0; inter_k < interResult.listInters.size();++inter_k)
                    {
                        InterPoint &intPk = interResult.listInters.at(inter_k);
                        if (intPk.merged_to == intPj.merged_to) intPk.merged_to=NInter;

                    }
                }
                intPj.merged_to=NInter;
                NInter++;

            }//END inter_j
        }//END inter_i
    }//END WHILE

    double bestdist;InterPoint *toMove;
    for (size_t inter_i=0; inter_i < interResult.listInters.size();++inter_i)
    {
        InterPoint &intPi = interResult.listInters.at(inter_i);
        bestdist=100;
        if (intPi.interaction != InterType::IONIC_PROT
                &&intPi.interaction != InterType::IONIC_LIG) continue;

        for (size_t inter_j=inter_i+1; inter_j < interResult.listInters.size();++inter_j)
        {
            InterPoint &intPj = interResult.listInters.at(inter_j);
            if ((intPi.interaction == InterType::IONIC_PROT
                 && intPj.interaction != InterType::HBOND_PROT)
                    ||(intPi.interaction == InterType::IONIC_LIG
                       && intPj.interaction != InterType::HBOND_LIG) )continue;

            dist = intPi.center.calcDist(intPj.center,1.1);
            if (dist < bestdist) {bestdist = dist;toMove=&intPj;}


        }
        if (bestdist <1)
        {

#ifdef ICHEM_DEBUG
            cout << "#########################################"<<endl;
            if (intPi.interaction == InterType::IONIC_PROT)
                cout << "MOVING IONIC PROT CATION TO H-BOND PROT DONOR POSITION"<<endl;
            else
                cout << "MOVING IONIC LIG CATION TO H-BOND LIG DONOR POSITION"<<endl;
            cout << "Ionic pos: "
                 <<intPi.Lig_Ref->getIdentifier()
                << "  "
                << intPi.Prot_Ref->getIdentifier()
                <<  intPi.center.toString()<<endl
                 << "HBond pos: "
                 <<toMove->Lig_Ref->getIdentifier()
                << "  "
                << toMove->Prot_Ref->getIdentifier()
                << toMove->center.toString()<<endl
                   ;
#endif
            intPi.center = toMove->center;

        }
    }

}


std::string Interactions::toString(const InterResults& interResult) const {
    ostringstream oss;

    // Stores ligand cycle atoms per interaction index
    std::map<int, std::vector<std::string>> ligandCycleAtoms; 

    for (std::vector<InterPoint>::const_iterator it = interResult.listInters.begin();
         it != interResult.listInters.end(); it++)
    {
        const InterPoint& interpt = *it;
        if (interpt.merged_to != -1) continue;

        switch (interpt.interaction)
        {
        case InterType::HBOND_LIG:   oss << "HBond_LIG         "; break;
        case InterType::HBOND_PROT:  oss << "HBond_PROT        "; break;
        case InterType::IONIC_LIG:   oss << "Ionic_LIG         "; break;
        case InterType::IONIC_PROT:  oss << "Ionic_PROT        "; break;
        case InterType::METAL:       oss << "Metal_Acceptor    "; break;
        case InterType::HYDROPHOBIC: oss << "Hydrophobic       "; break;
        case InterType::ARFACEFACE:  oss << "Aromatic_Face/Face"; break;
        case InterType::AREDGEFACE:  oss << "Aromatic_Edge/Face"; break;
        case InterType::PICATION:    oss << "Pi/Cation         "; break;
        case InterType::WHBOND_LIG:  oss << "Weak_HBond_LIG    "; break;
        case InterType::WHBOND_PROT: oss << "Weak_HBond_PROT   "; break;
        }

        if (interpt.interaction != InterType::AREDGEFACE
            && interpt.interaction != InterType::ARFACEFACE
            && interpt.interaction != InterType::PICATION)
        {
            oss << "\t"; oss.setf(ios::left);
            oss << "|" << interpt.Prot_Ref->getIdentifier();
            oss << "\t|" << interpt.Lig_Ref->getIdentifier();
            oss << "\t|" << interpt.point;
            oss << "\t|"; oss.width(7); oss.setf(ios::right); oss << interpt.dist;
            oss << "\t|";
            if (interpt.angle != -100000)
                oss << (interpt.angle * 180 / M_PI);
            else
                oss << "/";
        }
        else if (interpt.interaction == InterType::PICATION)
        {
            if (interpt.Prot_Ref->getName() == "DuAr")
            {
                Cycle* cyc = interpt.Prot_Ref->getParent().getCycleFromCenter(interpt.Prot_Ref);
                oss << "\t"; oss.setf(ios::left);
                oss << "|"; oss.width(5); oss.setf(ios::left); oss << cyc->getCenter().getName();
                oss << "|"; oss.width(5); oss.setf(ios::left); oss << "NULL";
                oss << "|"; oss << cyc->getAtom(0)->getResidu()->getIdentifier();
                oss << "\t|" << interpt.Lig_Ref->getIdentifier();
                oss << "\t|" << interpt.point;
                oss << "\t|"; oss.width(7); oss.setf(ios::right); oss << interpt.dist;
                oss << "\t|";
                if (interpt.angle != -100000)
                    oss << (interpt.angle * 180 / M_PI);
                else
                    oss << "/";
            }
            else
            {
                oss << "\t"; oss.setf(ios::left);
                oss << "|" << interpt.Prot_Ref->getIdentifier();
                Cycle* cyc = interpt.Lig_Ref->getParent().getCycleFromCenter(interpt.Lig_Ref);
                oss << "\t";
                oss << "|"; oss.setf(ios::left, ios::adjustfield); oss.width(5); oss << cyc->getCenter().getName();
                oss << "|"; oss.setf(ios::left); oss.width(5); oss << "0";
                oss << "|"; oss << "UNK1001-XX";
                oss << "\t|" << interpt.point;
                oss << "\t|"; oss.width(7); oss.setf(ios::right); oss << interpt.dist;
                oss << "\t|";
                oss << (interpt.angle * 180 / M_PI);
            }
        }
        else
        {
            Cycle* cyc = interpt.Prot_Ref->getParent().getCycleFromCenter(interpt.Prot_Ref);
            oss << "\t";
            oss << "|"; oss.setf(ios::left, ios::adjustfield); oss.width(5); oss << cyc->getCenter().getName();
            oss << "|"; oss.setf(ios::left); oss.width(5); oss << "NULL";
            oss << "|"; oss << cyc->getAtom(0)->getResidu()->getIdentifier();

            cyc = interpt.Lig_Ref->getParent().getCycleFromCenter(interpt.Lig_Ref);
            oss << "\t";
            oss << "|"; oss.setf(ios::left, ios::adjustfield); oss.width(5); oss << cyc->getCenter().getName();
            oss << "|"; oss.setf(ios::left); oss.width(5); oss << "NULL";
            oss << "|"; oss << cyc->getAtom(0)->getResidu()->getIdentifier();

            oss << "\t|" << interpt.point;
            oss << "\t|"; oss.width(7); oss.setf(ios::right); oss << interpt.dist;
            oss << "\t|/";
        }

        oss << endl;

        // Store ligand cycle atoms for ARFACEFACE and AREDGEFACE interactions mapped by interaction index
        if (interpt.interaction == InterType::ARFACEFACE || interpt.interaction == InterType::AREDGEFACE)
        {
            Cycle* cyc = interpt.Lig_Ref->getParent().getCycleFromCenter(interpt.Lig_Ref);
            if (cyc != nullptr)
            {
                for (ItCAtom itA = cyc->first(); itA != cyc->end(); ++itA)
                {
                    ligandCycleAtoms[interpt.point].push_back((*itA)->getName());
                }
            }
        }
    }

    // Ligand cycle atoms printed, grouped by interaction index
    if (!ligandCycleAtoms.empty())
    {
        oss << "\n==== Ligand atoms involved in the aromatic cycle: ====\n";
        for (const auto& entry : ligandCycleAtoms)
        {
            oss << "Interaction " << entry.first << ": ";
            for (const auto& atom : entry.second)
            {
                oss << atom << " ";
            }
            oss << endl;
        }
    }
    return oss.str();
}


void Interactions::interToMOL2(  InterResults& interResult,
                                 const bool&Ligand,
                                 const bool&Protein,
                                 const bool&Center,
                                 const bool&Merge)const
{


    if (interResult.listInters.empty())return;



    bool exists=false;
    Atom *atml;

    const std::string AtmNames[NB_INTTYPE]= {"","N"   ,"O"  ,"NZ" ,"OD1","CA","Zn","CZ","CZ","","", "Zn"};
    const std::string mol2Names[NB_INTTYPE]={"","N.am","O.2","N.4","O.co2"  ,"C.3","Zn","C.ar","C.ar","","","O.co2"};
    const std::string AtomicName[NB_INTTYPE]={"","N"   ,"O"  ,"N" ,"O","C","Zn","C","C","","","O"};
    unsigned int NAtm=1;
    Molecule& mole=interResult.Ints;
    mole.clear();
    if (Merge) mole.setMoleType(MoleType::INT_MERG);
    else if (Ligand && !Protein && !Center) mole.setMoleType(MoleType::INT_LIG);
    else if (!Ligand && Protein && !Center) mole.setMoleType(MoleType::INT_PROT);
    else if (!Ligand && !Protein && Center) mole.setMoleType(MoleType::INT_CENT);
    else mole.setMoleType(MoleType::INT_MERG);


    if (Ligand || Merge)
    {



        Residu *listRes[NB_INTTYPE]={(Residu*)NULL,

                                     &mole.addResidu("Z",1,"SEL"),
                                     &mole.addResidu("Z",2,"ALL"),
                                     &mole.addResidu("Z",3,"LYL"),
                                     &mole.addResidu("Z",4,"ASL"),
                                     &mole.addResidu("Z",5,"GLL"),
                                     &mole.addResidu("Z",7,"ZIL"),
                                     &mole.addResidu("Z",6,"PHL"),
//                                     &mole.addResidu("Z",6,"PHL"),
                                     (Residu*)NULL,(Residu*)NULL,(Residu*)NULL};
        listRes[InterType::AREDGEFACE]=listRes[InterType::ARFACEFACE];

        for (size_t i=0; i< interResult.listInters.size();++i)
        {
            InterPoint &interP=interResult.listInters.at(i);
            if (interP.merged_to != -1)continue;
            exists=false;
            for (ItCAtom itC = mole.firstAtom();itC != mole.lastAtom();itC++)
            {
                Atom &atm = **itC;
                if (atm.fixpos.calcDist(interP.Lig_Ref->fixpos) > 0.005)continue;
                if (atm.getName()== AtmNames[interP.interaction]){ exists=true;break;}
            }
            if (exists)continue;
            if (listRes[interP.interaction] == (Residu*)NULL)continue;
            atml=&mole.addAtom(AtomicName[interP.interaction],
                    interP.Lig_Ref->fixpos,
                    AtmNames[interP.interaction],
                    mol2Names[interP.interaction],
                    listRes[interP.interaction]);
            atml->setFNum(NAtm);++NAtm;
            atml->setTNum(i);
            atml->setPartialCharge(interP.Lig_Ref->getPartialCharge());

        }

    }
    if (Protein || Merge)
    {

        Residu *listRes[NB_INTTYPE]={(Residu*)NULL,

                                     &mole.addResidu("Z",1,"SEP"),
                                     &mole.addResidu("Z",2,"ALP"),
                                     &mole.addResidu("Z",3,"LYP"),
                                     &mole.addResidu("Z",4,"ASP"),
                                     &mole.addResidu("Z",5,"GLP"),
                                     &mole.addResidu("Z",7,"ZIP"),
                                     &mole.addResidu("Z",6,"PHP"),
                                     //                                     &mole.addResidu("Z",6,"PHP"),
                                     (Residu*)NULL,(Residu*)NULL,(Residu*)NULL};
        listRes[InterType::AREDGEFACE]=listRes[InterType::ARFACEFACE];

        for (size_t i=0; i< interResult.listInters.size();++i)
        {
            InterPoint &interP=interResult.listInters.at(i);
            if (interP.merged_to != -1)continue;
            exists=false;
            for (ItCAtom itC = mole.firstAtom();itC != mole.lastAtom();itC++)
            {
                Atom &atm = **itC;
                if (atm.fixpos.calcDist(interP.Prot_Ref->fixpos) > 0.005)continue;
                if (atm.getName()== AtmNames[interP.interaction]){ exists=true;break;}
            }
            if (exists)continue;
            if (listRes[interP.interaction] == (Residu*)NULL)continue;
            atml=&mole.addAtom(AtomicName[interP.interaction],
                    interP.Prot_Ref->fixpos,
                    AtmNames[interP.interaction],
                    mol2Names[interP.interaction],
                    listRes[interP.interaction]);
            atml->setFNum(NAtm);++NAtm;
            atml->setTNum(i);

        }

    }
    if (Center || Merge)
    {


        Residu *listRes[NB_INTTYPE]={(Residu*)NULL,

                                     &mole.addResidu("Z",1,"SEC"),
                                     &mole.addResidu("Z",2,"ALC"),
                                     &mole.addResidu("Z",3,"LYC"),
                                     &mole.addResidu("Z",4,"ASC"),
                                     &mole.addResidu("Z",5,"GLC"),
                                     &mole.addResidu("Z",7,"ZIC"),
                                     &mole.addResidu("Z",6,"PHC"),
                                     //                                     &mole.addResidu("Z",6,"PHC"),
                                     (Residu*)NULL,(Residu*)NULL,(Residu*)NULL};
        listRes[InterType::AREDGEFACE]=listRes[InterType::ARFACEFACE];

        for (size_t i=0; i< interResult.listInters.size();++i)
        {
            exists = false;
            InterPoint &interP=interResult.listInters.at(i);
            if (interP.merged_to != -1)continue;
            if (listRes[interP.interaction] == (Residu*)NULL)continue;
            for (ItCAtom itC = mole.firstAtom();itC != mole.lastAtom();itC++)
            {
                Atom &atm = **itC;
                if (atm.fixpos.calcDist(interP.center) > 0.005)continue;
                if (atm.getName()== AtmNames[interP.interaction]){
                    exists=true;
//                    cout << "superposé Center ";
                    atm.setTNum(interP.point);
//                    for (int flo:atm.getTNum())
//                        cout << flo << " ";
//                    cout << endl;
                    break;

                }
            }
            if (exists)continue;
            atml=&mole.addAtom(AtomicName[interP.interaction],
                    interP.center,
                    AtmNames[interP.interaction],
                    mol2Names[interP.interaction],
                    listRes[interP.interaction]);
            atml->setFNum(NAtm);++NAtm;
//            atml->setTNum(i);
            atml->setTNum(interP.point);
        }
    }




}

void Interactions::calcenfouiss(string mol_name) // calcul l'enfouissement des interfaces
{

    cout << "calcul enfouiss"<< endl;
    VolSite volsite(complex);
    volsite.proj(55-20,false);
    MoleReader mread2;

    string PDB_name =     complex.getMole(0)->getName();
    string enffile = "../../ints.enfs"; //*$

    ofstream offenf;
    if (fopen(enffile.c_str(),"r")==NULL)
    {
        offenf.open(enffile.c_str(),ios::out|ios::app);
    }
    else
    {
        offenf.open(enffile.c_str(),ios::out|ios::app);
    }

    //        cout << A << "-" <<  B << " avec "  << intersize[i]->size << endl;
    //        cout << PDB_name+"_ints_"+A+"_"+B+".mol2" << endl;
    Molecule test;
    mread2.loadNewFile(mol_name+".mol2");
    mread2.loadNextMolecule(test,MoleType::LIGAND);


    //        cout << "ma proteine " << test.toString() << endl;
    //        return;

    vector<double> FGPCount;
    vector<vector<double> > FGPBuried;
    for (size_t j=0;j<8;++j)
    {
        FGPCount.push_back(0);
        if (j==0)continue;
        vector<double> tmp;
        for (size_t j=0;j<10;++j)tmp.push_back(0); FGPBuried.push_back(tmp);
    }
    for (size_t iAtm=0; iAtm < test.numAtom();++iAtm)
    {
        const Atom& atm= test.getAtom(iAtm);
        Box* box = volsite.getGrid().getBox(atm.fixpos);
        //cout << atm.getName()<<"\t"<<box->getId()<< " " << atm.fixpos.calcDist(box->fixpos)<< " " << volsite.getProjValue(*box)<<endl;
        double pos= floor((volsite.getProjValue(*box)-30)/10);
        if (offenf.is_open())
        {
            offenf  << atm.getName()<< " " << volsite.getProjValue(*box)  << endl;
        }

        //            cout << volsite.getProjValue(*box) << "\t col : " << pos+1 << endl;
        //pas de distinction protein-lignadd
        //passer tous les <5 a <7
        if (atm.getName()=="CA")      {FGPBuried[0][pos]+=1;FGPCount[1]++;FGPCount[0]++;}
        else if (atm.getName()=="O")  {FGPBuried[1][pos]+=1;FGPCount[2]++;FGPCount[0]++;}
        else if (atm.getName()=="OD1"){FGPBuried[2][pos]+=1;FGPCount[3]++;FGPCount[0]++;}
        else if (atm.getName()=="N")  {FGPBuried[1][pos]+=1;FGPCount[2]++;FGPCount[0]++;}
        else if (atm.getName()=="NZ") {FGPBuried[2][pos]+=1;FGPCount[3]++;FGPCount[0]++;}
        else if (atm.getName()=="CZ") {FGPBuried[3][pos]+=1;FGPCount[4]++;FGPCount[0]++;}
        else if (atm.getName()=="Zn") {FGPBuried[4][pos]+=1;FGPCount[5]++;FGPCount[0]++;}
    }
    for (size_t i=0; i<4;++i)
    {

        for (size_t j=0;j<10;++j)
        {
            if (FGPBuried[i][j]==0)continue;
            //                    cout << i << ":" << j << "valeur " << FGPBuried[i][j]<< endl;
            FGPBuried[i][j]=FGPBuried[i][j];

        }
        FGPCount[i+1]=FGPCount[i+1]/FGPCount[0]*100.0;

    }

    //    cout <<complex.ge<<"\t";
    ofstream ofscore;
    ofscore.open("../../descriptor.sre", ios::out|ios::app);
    for (size_t i=0;i<=4;++i) {
        cout << FGPCount[i]<< " ";
        ofscore << FGPCount[i]<< " ";
    }
    for (size_t i=0; i<4;++i) // 4 sans metal 5 avec metal
    {

        for (size_t j=0;j<10;++j)
        {

            cout << FGPBuried[i][j]<< " ";
            ofscore <<FGPBuried[i][j]<< " ";
        }

    }


    cout <<"#"<< PDB_name << endl;

    ofscore <<"#"<< PDB_name<< endl;
    ofscore.close();
    offenf.close();

    return;
}


void Interactions::loadTripletsComb(const bool full)
{
    ostringstream oss;
    ifstream Finput;
    oss << ICMole::get_IChem_LIB_Path()<<"/datas/PH4_list_6_7"<< ((full)? "F":"M")<<".csv";
    Finput.open(oss.str().c_str(),ios::in);
    if (!Finput.is_open()) {cerr << "ERROR WHILE OPENING " << oss.str()<<endl;return;}


    const unsigned int _pos=7;
    const unsigned int dist_size=6;
    map<char,int> Prop_To_Int;
    string ligne;
    Prop_To_Int.insert(pair<char,int>('A',InterType::HBOND_PROT));
    Prop_To_Int.insert(pair<char,int>('B',InterType::HBOND_LIG));
    Prop_To_Int.insert(pair<char,int>('C',InterType::IONIC_PROT));
    Prop_To_Int.insert(pair<char,int>('D',InterType::IONIC_LIG));
    Prop_To_Int.insert(pair<char,int>('E',InterType::AREDGEFACE));
    Prop_To_Int.insert(pair<char,int>('F',InterType::HYDROPHOBIC));
    Prop_To_Int.insert(pair<char,int>('G',InterType::METAL));

//    HBOND_PROT=1,  /*!< HBond Protein donor        */
//    HBOND_LIG=2,  /*!< HBond Ligand donor         */
//    IONIC_PROT=3,  /*!< Ionic Protein cation       */
//    IONIC_LIG=4,  /*!< Ionic Ligand cation        */
//    HYDROPHOBIC=5,  /*!< Hydrophobic                */
//    METAL=6,  /*!< Metal acceptor interaction */
//    ARFACEFACE=7,  /*!< Aromatic face to face      */
//    AREDGEFACE=8,  /*!< Aromatic edge to face      */
//    PICATION=9,  /*!< PI cation                  */
//    WHBOND_LIG=10,
//    WHBOND_PROT=11,

    path_dist = new double[7];
    path_dist[0]=0;
    path_dist[1]=4;
    path_dist[2]=6;
    path_dist[3]=9;
    path_dist[4]=13;
    path_dist[5]=17;
    path_dist[6]=32;

    path_dist_size=6;




    //  if (debug){
    //      cout << "DEBUG ******************************************************"<<endl
    //           << "DEBUG ********** LOADING TRIPLETS COMBINAISONS *************"<<endl
    //           << "DEBUG ******************************************************"<<endl;
    //      cout << oss.str() <<endl;
    //      cout << "DEBUG Number of dist_size : " << dist_size <<endl
    //           << "DEBUG Number of interaction groups : " << _pos <<endl;
    //    }

    for (unsigned int i=0; i<1000; i++)
        for (unsigned int j=0; j<1000; j++) vect_list[i][j]=-1;


    std::getline( Finput, ligne ); size_triplets=atoi(ligne.c_str());

    //  if (debug)    cout << "Number of possible triplets : "<< (size_triplets+1) <<endl;
    int ni=0,ni2=0,pos=0;

    while (!Finput.eof())
    {
        std::getline( Finput, ligne );
        if (ligne.length() < 8) continue; // Idk why this line is here. But without it, "core dumped" issues occurred
        ni = atoi(ligne.substr(0,3).c_str());
        oss.str("");
        oss << Prop_To_Int[ligne[4]]<<Prop_To_Int[ligne[5]]<<Prop_To_Int[ligne[6]];
        ni2 = atoi(oss.str().c_str());
        pos = atoi(ligne.substr(8).c_str());
        vect_list[ni][ni2]=pos;
        //      if (debug)cout << ni << " " << ni2 << " " << pos <<endl;
        if (ligne[4] == 'E')
        {
            oss.str("");
            oss << InterType::ARFACEFACE<<Prop_To_Int[ligne[5]]<<Prop_To_Int[ligne[6]];
            ni2 = atoi(oss.str().c_str());
            vect_list[ni][ni2]=pos;
            //          if (debug)cout << ni << " " << ni2 << " " << pos <<endl;
        }
        if (ligne[5] == 'E')
        {
            oss.str("");
            oss << Prop_To_Int[ligne[4]]<<InterType::ARFACEFACE<<Prop_To_Int[ligne[6]];
            ni2 = atoi(oss.str().c_str());
            vect_list[ni][ni2]=pos;
            //          if (debug)       cout << ni << " " << ni2 << " " << pos <<endl;
        }
        if (ligne[6] == 'E')
        {
            oss.str("");
            oss << Prop_To_Int[ligne[4]]<<Prop_To_Int[ligne[5]]<<InterType::ARFACEFACE;
            ni2 = atoi(oss.str().c_str());
            vect_list[ni][ni2]=pos;
            //          if (debug)     cout << ni << " " << ni2 << " " << pos <<endl;
        }
    }
    Finput.close();
    Interactions::Load_triplet=true;


}

Fingerprint& Interactions::generateTriplets(InterResults& interResult,const bool full)
{




    if (!Load_triplet) loadTripletsComb(full);

    Fingerprint *fgp = new Fingerprint(size_triplets+1,true);
    double dist1, dist2, dist3;
    double dist11, dist12, dist13;
    const unsigned int size = interResult.listInters.size();
    if (size < 3) { return *fgp;}
    double TableDist[size][size];
    unsigned int IV1=0, IV2=0,IV3=0;

    // Compute the pairwise IPA distances
    for (std::vector<InterPoint>::const_iterator pt1=interResult.listInters.begin();pt1 != interResult.listInters.end();pt1++){
        //  for (ItPInter pt1 = interResult.listInters.begin(); pt1 != interResult.listInters.end(); pt1++){
        IV2=IV1;
        for (std::vector<InterPoint>::const_iterator pt2=pt1 ; pt2 != interResult.listInters.end(); pt2++)	{

            dist1 = getCenter(interResult,pt1->point).calcDist(getCenter(interResult,pt2->point));
            TableDist[IV1][IV2]=dist1;
            TableDist[IV2][IV1]=dist1;
            IV2++;

        }
        IV1++;
    }

    //cout << "#############  BEGIN "<<endl;
    int num=0,num1=0;
    IV1=0; IV2=0;IV3=0;
    int intdist1,intdist2,intdist3;
    unsigned int intP1, intP2,intP3;
    unsigned int intP11, intP12,intP13;
    int nd=0, np=0;
    for (std::vector<InterPoint>::const_iterator pt1=interResult.listInters.begin();pt1 != interResult.listInters.end();pt1++){
        IV2=IV1+1;intP1 = pt1->interaction;
        // TODO: Correct the bug with pt2 = pt1 + 1
        for (std::vector<InterPoint>::const_iterator pt2=pt1 ; pt2 != interResult.listInters.end(); pt2++)	{
            IV3=IV2+1;intP2 = pt2->interaction;
            for (std::vector<InterPoint>::const_iterator pt3 = pt2+1; pt3 != interResult.listInters.end(); pt3++) {
                //  if (!(pt1->interaction == InterType::HYDROPHOBIC|| pt2->interaction == InterType::HYDROPHOBIC || pt3->interaction == InterType::HYDROPHOBIC))continue;
                //  if (pt1->interaction != InterType::HYDROPHOBIC && pt2->interaction != InterType::HYDROPHOBIC && pt3->interaction != InterType::HYDROPHOBIC)continue;
                // cout << "ST"<<endl;
                num1++;
                intP3 = pt3->interaction;
                dist1 = TableDist[IV1][IV2];
                dist2 = TableDist[IV2][IV3];
                dist3 = TableDist[IV3][IV1];


                //                cout << dist1<< " " << dist2 << " " << dist3 << " " << pt1->interaction << " " << pt2->interaction << " " <<pt3->interaction<<endl;
                if (dist1 > dist2+dist3){cerr << "TRIANGULAR INEGALITY NOT VERIFIED !! " <<endl;continue;}
                if (dist2 > dist1+dist3){cerr << "TRIANGULAR INEGALITY NOT VERIFIED !! " <<endl;continue;}
                if (dist3 > dist1+dist2){cerr << "TRIANGULAR INEGALITY NOT VERIFIED !! " <<endl;continue;}


                dist11=dist1;
                intP11=intP1;

                if (dist2 > dist1) {
                    dist11 = dist2;
                    intP11=intP2;
                    dist12=dist1;
                    intP12=intP1;
                } else {
                    dist12=dist2;
                    intP12=intP2;
                }

                if (dist3 > dist11){
                    dist13=dist12;
                    intP13=intP12;
                    dist12=dist11;
                    intP12=intP11;
                    dist11=dist3;
                    intP11=intP3;
                } else if (dist3 > dist12) {
                    dist13=dist12;
                    intP13=intP12;
                    dist12=dist3;
                    intP12=intP3;
                } else {
                    dist13=dist3;
                    intP13=intP3;
                }

                if (dist13 == 0) dist13 = 0.001;
                intdist1 = -1; intdist2=-1;intdist3=-1;
                for (unsigned int pos = 0 ; pos < path_dist_size; pos++) if (dist11 >= path_dist[pos] && dist11 < path_dist[pos+1])   intdist1=(pos+1);
                for (unsigned int pos = 0 ; pos < path_dist_size; pos++) if (dist12 >= path_dist[pos] && dist12 < path_dist[pos+1])   intdist2=(pos+1);
                for (unsigned int pos = 0 ; pos < path_dist_size; pos++) if (dist13 >= path_dist[pos] && dist13 < path_dist[pos+1])   intdist3=(pos+1);

                if (intdist1 == -1 || intdist2 == -1 || intdist3 == -1) continue;
                nd = intdist1*100+intdist2*10+intdist3;
                np = intP11*100+intP12*10+intP13;
                //                if (debug)   cout << nd << " " << np << " " << vect_list[nd][np]<<endl;

               if (vect_list[nd][np] != -1) {
                    fgp->pushOne(vect_list[nd][np]);
                }


                num++;
                IV3++;
            }
            IV2++;
        }
        IV1++;
    }

    //    if (debug) cout << num << " " << num1<<endl;// <<" " << fgp->__ToString()<<endl;

    return *fgp;
}


// Map of residue pointers used to build the fingerprint IFP, we iterate over residues in NumToRes and turn detected interactions into bits // our problem for NA sodium is that the size of NumToRes = 0
void Interactions::genIFP(InterResults& interResult, const unsigned int& fgpType) const {

    static const int intToPos[5][NB_INTTYPE]= {
        {-1,3,4,5,6, 0,-1, 1, 2,-1,-1,-1,-1,-1},
        {-1,0,1,2,3,-1,4,-1,-1,-1,-1,-1,-1,-1},
        {-1,3,4,5,6, 0,8, 1, 2, 7, -1, -1,-1,-1},
        {-1,0,1,2,3,-1,7,-1,-1, 6, 4, 5,-1,-1},
        {-1,-1,-1,-1,-1,-1,0,-1,-1,-1,-1,-1,-1}
    };

    static const short length[5]={7,5,9,8,1};
    //     U ,H,H,I,I, H,M , A, A
    //     N ,B,B,O,O, Y,E , R, R
    //     D ,P,L,P,L, D,T , F, E
    //     E , , , , ,  ,  , F, F
    //     F
    
    map<int,Residu*> NumtoRes;
    for (ItCRes itR = complex.firstResidu();itR != complex.lastResidu();++itR) 
    {
        Residu *res =*itR;

        

        if (!res->isUsed() || res->getParent()->getMoleType()== MoleType::LIGAND) 
            continue;
        
        if (Residu::Rules[res->getParent()->getMoleType()][res->getResType()] == MoleType::UNDEFINED || Molecule::Rules[res->getParent()->getMoleType()] == MoleType::UNDEFINED)
            continue;
        
        NumtoRes.insert(pair<int,Residu*>(res->getNum(),res));
    }


    multimap<Residu*,InterPoint*> listRes;

    for (size_t i=0; i< interResult.listInters.size();++i)
    {
        InterPoint& interP = interResult.listInters.at(i);
        if (interP.merged_to != -1)continue;
        Residu *res =interP.Prot_Ref->getResidu();
        if (interP.interaction==InterType::AREDGEFACE || interP.interaction==InterType::ARFACEFACE)
        {
            Cycle* cyc=interP.Prot_Ref->getParent().getCycleFromCenter(interP.Prot_Ref);
            res=cyc->getAtom(0)->getResidu();

        }
        listRes.insert(pair<Residu*,InterPoint*>(res,&interP));

    }

    interResult.IFP=Fingerprint(NumtoRes.size()*length[fgpType]);
    int NRes=0;
    ostringstream oss;
    for (map<int,Residu*>::iterator it =NumtoRes.begin(); it != NumtoRes.end();++it ) {

        std::pair <std::multimap<Residu*,InterPoint*>::iterator, std::multimap<Residu*,InterPoint*>::iterator> ret;
        ret = listRes.equal_range((*it).second);

        for (std::multimap<Residu*,InterPoint*>::iterator it2=ret.first; it2!=ret.second; ++it2) {
            if (intToPos[fgpType][(*it2).second->interaction] == -1)
                continue;

            interResult.IFP.bitOn(NRes*length[fgpType]+intToPos[fgpType][(*it2).second->interaction]);
        }
        interResult.IFPString+="|";
        oss.str("");
        if ((*it).second->getResType() == ResType::STD_AA)
        {
            bool found=false;
            for (size_t pData=0;pData <NBAA;++pData)
            {

                if (AAcid[pData].name==(*it).second->getName())
                {
                    oss << AAcid[pData].code << (*it).second->getFNum();
                    found=true;
                    break;
                }
            }
            if (!found){
                oss << (*it).second->getName() << (*it).second->getFNum();
            }
        }
        else {      oss <<(*it).second->getName()<<(*it).second->getFNum();
        }
        //        cout << oss.str() << endl;
        interResult.IFPString+=(*it).second->getChainName();
        for (int ni=0; ni < length[fgpType]-(int)oss.str().length()-2;ni++)
            interResult.IFPString += " ";
        interResult.IFPString += oss.str();


        NRes++;
    }
}







