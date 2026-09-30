// Unit tests of ICMole::Molecule editing (test/dataset ligands)

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "headers/ICMole/atom.h"
#include "headers/ICMole/cycle.h"
#include "headers/ICMole/molecule.h"
#include "headers/ICPars/molereader.h"

using namespace ICMole;

namespace
{

int failures = 0;

#define CHECK(...)                                                                     \
    do {                                                                               \
        if (!(__VA_ARGS__)) {                                                          \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #__VA_ARGS__ \
                      << std::endl;                                                    \
            ++failures;                                                                \
        }                                                                              \
    } while (0)

std::size_t ringsWith(const Molecule& molecule, const Atom& atom)
{
    std::size_t n = 0;
    for (ItCCycle it = molecule.firstCycle(); it != molecule.lastCycle(); ++it)
        if ((*it)->hasAtom(atom)) ++n;
    return n;
}

// Every atom of a ring pseudo-residue must be the centre of a remaining ring
bool centresValid(Molecule& molecule, Residu& pseudo)
{
    for (ItCAtom itA = pseudo.firstCAtom(); itA != pseudo.lastCAtom(); ++itA)
    {
        bool found = false;
        for (ItCCycle it = molecule.firstCycle(); it != molecule.lastCycle(); ++it)
            found = found || &(*it)->getCenter() == *itA;
        if (!found) return false;
    }
    return true;
}

// Deleting an atom shared by fused rings removes all of them, and their centres
void testDelAtomInFusedRings(const std::string& path)
{
    MoleReader reader;
    reader.loadNewFile(path);
    Molecule molecule(100, 100, MoleType::LIGAND);
    reader.loadNextMolecule(molecule, MoleType::LIGAND);
    molecule.checkMOL2();
    molecule.ringPerception();

    Atom* fused = nullptr;
    for (ItCAtom it = molecule.firstAtom(); it != molecule.lastAtom() && fused == nullptr; ++it)
        if (ringsWith(molecule, **it) >= 2) fused = *it;
    CHECK(fused != nullptr);
    if (fused == nullptr) return;

    const std::size_t before = molecule.numCycles();
    const std::size_t involved = ringsWith(molecule, *fused);
    molecule.delAtom(fused);

    CHECK(molecule.numCycles() == before - involved);
    for (ItCCycle it = molecule.firstCycle(); it != molecule.lastCycle(); ++it)
        for (ItCAtom itA = (*it)->first(); itA != (*it)->end(); ++itA) CHECK(*itA != fused);
    CHECK(centresValid(molecule, molecule.getAromRes()));
    CHECK(centresValid(molecule, molecule.getCycleRes()));
}

}

int main()
{
    // Legacy reading needs the IChem data (HETLIST) of this source tree
    setenv("ICHEM_LIB", ICHEM_SOURCE_DIR, 1);

    try
    {
        // Ligands with fused rings
        testDelAtomInFusedRings(std::string(ICHEM_SOURCE_DIR) + "/test/dataset/fraglit/1HNN_1_2_ligand.mol2");
        testDelAtomInFusedRings(std::string(ICHEM_SOURCE_DIR) + "/test/dataset/fraglit/1H9Z_1_7_ligand.mol2");
    }
    catch (MoleExcept& e)
    {
        std::cerr << "Unexpected exception: " << e.getData() << std::endl;
        ++failures;
    }

    if (failures != 0)
    {
        std::cerr << "[FAILED] " << failures << " molecule check(s)" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "[OK] molecule tests passed" << std::endl;
    return EXIT_SUCCESS;
}
