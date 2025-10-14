#include <iostream>
#include <headers/ICTools/switch.h>

using namespace std;
using namespace ICMole;


void IChemSwitch::helpAtomProps() const {
        cout  << "NAME\n"
              << "      IChem AtomProps - Display atom-level properties from MOL2 files\n\n\n"
              << "SYNOPSIS:\n"
              << "      IChem AtomProps <mol2_file1> [<mol2_file2> ... <mol2_filen]\n\n\n"
              << "DESCRIPTION:\n"
              << "      Parses one or more MOL2 files and displays the atom-level properties for each structure\n\n"
              << "  Supported input:\n"
              << "      - Any valid MOL2 file\n\n"
              << "  Displayed properties per atom:\n"
              << "      - Atom ID\n"
              << "      - MOL2 atom type\n"
              << "      - Chemical property (e.g. apolar, hydrophobic, acceptor)\n"
              << "      - Charge\n\n\n"
              << "EXAMPLES:\n"
              << "      IChem AtomProps ligand.mol2\n"
              << "      IChem AtomProps protein.mol2 ligand.mol2 \n"
              << "\n";
}


void IChemSwitch::printAtomInfo(const ICMole::Atom& atom) const {
    cout << "ATOM \t" << atom.getIdentifier() << "\n"
              << "\t type : " << atom.getMOL2Type() << "\n"
              << "\t props : " << atom.props.toString() << "\n"
              << "\t charge: " << atom.getPartialCharge() << "\n\n";
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
        
        // Process each molecule file
        for (size_t i = 0; i < Input_Values.size(); ++i) {
            const std::string& fLigand = Input_Values.at(i); // Current molecule file

            // Load the molecule file
            iread.loadNewFile(fLigand);

            // Looping through all molecules
            while (!iread.isEOF()) {
                Molecule ligand;
                iread.loadNextMolecule(ligand, MoleType::LIGAND);

                cout << "Molecule " << ligandCount++ << " with name " << ligand.getName() << " from file: " << fLigand << "\n\n";
                // Print atom properties for the current molecule
                for (ItCAtom itLA = ligand.firstAtom(); itLA != ligand.lastAtom(); ++itLA) {
                    printAtomInfo(**itLA);
                }
            }
        }
    } catch (MoleExcept& e) {
        cerr << "AtomProps processing files failed: " << e.getCode() << " " << e.getSource() << " " << e.getData() << "\n";
    }
}
