#include <iostream>
#include <headers/ICTools/switch.h>

using namespace std;
using namespace ICMole;


void IChemSwitch::helpAtomProps() const {
    cout << " AtomProps - Print atom properties of ligands\n"
              << "    Usage: IChem AtomProps protein.mol2 ligand1.mol2 [ligand2.mol2 ...]\n"
              << "    Description:\n"
              << "      This command prints the atom properties (ID, MOL2 type, properties, and charge)\n"
              << "      for each mol2 file provided as argument\n"
              << "\n";
}


void IChemSwitch::printAtomInfo(const ICMole::Atom& atom) const {
    cout << "ATOM \t" << atom.getIdentifier() << "\n"
              << "\t type: " << atom.getMOL2Type() << "\n"
              << "\t props : " << atom.props.toString() << "\n"
              << "\t charge: " << atom.getFormalCharge() << "\n\n";
}

/*
    * Main function of AtomProps command, allowing the loading of multimol2 file (protein, ligand)
    * Able to take as many arguments as possible  
*/
void IChemSwitch::AtomProps() const {
    
    Complex icomplex;
    MoleReader iread;
    int ligandCount = 1;

    try {
        
        // Process each ligand file
        for (size_t i = 0; i < Input_Values.size(); ++i) {
            const std::string& fLigand = Input_Values.at(i); // Current ligand

            // Load the ligand file
            iread.loadNewFile(fLigand);

            // Looping through all molecules in the file
            while (!iread.isEOF()) {
                Molecule ligand;
                iread.loadNextMolecule(ligand, MoleType::LIGAND);

                cout << "Ligand " << ligandCount++ << " with name " << ligand.getName() << " from file: " << fLigand << "\n\n";
                // cout << "Ligand " << ligandCount++ << " from file: " << fLigand << "\n\n";
                // Print atom properties for the current ligand
                for (ItCAtom itLA = ligand.firstAtom(); itLA != ligand.lastAtom(); ++itLA) {
                    printAtomInfo(**itLA);
                }
            }
        }
    } catch (MoleExcept& e) {
        cerr << "AtomProps processing files failed: " << e.getCode() << " " << e.getSource() << " " << e.getData() << "\n";
    }
}
