#include <iostream>
#include <headers/ICTools/switch.h>

using namespace std;
using namespace ICMole;


void IChemSwitch::helpAtomProps() const {
    cout << " AtomProps - Print atom properties of ligands\n"
              << "    Usage: IChem AtomProps protein.mol2 ligand1.mol2 [ligand2.mol2 ...]\n"
              << "    Options:\n"
              << "      protein.mol2:       The protein structure file\n"
              << "      ligand.mol2        One or more ligand files\n"
              << "    Description:\n"
              << "      This command prints the atom properties (ID, MOL2 type, properties, and charge)\n"
              << "      for each ligand molecule in the provided files.\n"
              << "\n";
}


void IChemSwitch::printAtomInfo(const ICMole::Atom& atom) const {
    cout << "ATOM \t" << atom.getIdentifier() << "\n"
              << "\t type: " << atom.getMOL2Type() << "\n"
              << "\t props : " << atom.props.toString() << "\n"
              << "\t charge: " << atom.getFormalCharge() << "\n\n";
}

/* 
    - The command must have this following argument structure: Protein ligand
    - Ligand can be multimol2 in one file, but it can also be many files (each either single or multiple mol2)
    - The command can take many arguments: ligand1, ligand2, ligand3, ..., ligandn
    - The result must display all of these properties: Identifier - type - props - charge
 */


void IChemSwitch::AtomProps() const {
    if (Input_Values.size() < 2) {
        throw MoleExcept(9020101, "IChemSwitch::AtomProps", "Must provide at least one protein file and one ligand file");
    }

    // Protein input
    const std::string& fProtein = Input_Values.at(0);

    Complex icomplex;
    MoleReader iread;

    try {
        // Load the protein file
        iread.loadNewFile(fProtein);
        iread.loadInComplex(icomplex, MoleType::PROTEIN);
        Molecule* protein = icomplex.getMole(MoleType::PROTEIN);

        if (protein == nullptr) {
            throw MoleExcept(9020102, "IChemSwitch::AtomProps", "No protein found in " + fProtein);
        }

        // Display protein properties
        cout << "\t\t\tProtein: " << protein->getName() << " (from file: " << fProtein << ")" << "\n\n";
        for(ItCAtom itPA = protein->firstAtom(); itPA != protein->lastAtom(); ++itPA ) {
            printAtomInfo(**itPA);
        }

        // Process each ligand file
        for (size_t i = 1; i < Input_Values.size(); ++i) {
            const std::string& fLigand = Input_Values.at(i); // Current ligand

            // Load the ligand file
            iread.loadNewFile(fLigand);

            // Looping through all molecules in the file
            while (!iread.isEOF()) {
                Molecule ligand;
                iread.loadNextMolecule(ligand, MoleType::LIGAND);

                cout << "Ligand: " << ligand.getName() << " from file: " << fLigand << "\n\n";

                // Print atom properties for the current ligand
                for (ItCAtom itLA = ligand.firstAtom(); itLA != ligand.lastAtom(); ++itLA) {
                    printAtomInfo(**itLA);
                }
            }
        }
    } catch (MoleExcept& e) {
        // throw runtime_error("Processing error \n");
        cerr << "Error: " << e.getCode() << " " << e.getSource() << " " << e.getData() << "\n";
    }
}