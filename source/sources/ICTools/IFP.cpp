#include "headers/ICTools/switch.h"
#include "headers/ICTools/ifp_module.hpp"

using namespace std;

void IChemSwitch::helpIFP() const {
    cout
        << "IFP - Interaction FingerPrint\n"
        << "Usage: IChem [options] IFP protein ligand\n"
        << "       IChem [options] IFP protein ligand ligand_ref\n"
        << "       IChem [options] IFP protein1 ligand1 protein2 ligand2\n"
        << "\n"
        << "Description:\n"
        << "  Compute interaction fingerprints between a protein and a ligand\n"
        << "  Compare a reference ligand and detect various interaction types\n"
        << "\n"
        << "Options:\n"
        << "  -name N (LIG)  Name of the fingerprint (default: ligand name)\n"
        << "  --polar        Detect and output only polar interactions\n"
        << "  --metal        Detect and output only metal interactions\n"
        << "  --extended     Include in the fingerprint:\n"
        << "                    * Metal/Acceptor interactions\n"
        << "                    * Weak Hydrogen bonds\n"
        << "                    * Pi-Cation interactions\n"
        << "  --ligD         Print all possible ligand interactions\n"
        << "  --bitstringOFF Do not output the bitstring representation\n"
        << "\n"
        << "Thresholds and angles:\n"
        << "  -D_Hb    N (3.5)  Hbond length (Å)\n"
        << "  -D_Hyd   N (4.5)  Hydrophobic length (Å)\n"
        << "  -D_Io    N (4.0)  Ionic length (Å)\n"
        << "  -D_Me    N (2.8)  Metal/Acceptor length (Å)\n"
        << "  -D_Ar    N (5.0)  Aromatic interaction length (Å)\n"
        << "  -D_Pic   N (5.0)  Pi-Cation interaction length (Å)\n"
        << "  -a_H     N (Pi)   Hbond angle (rad)\n"
        << "  -at_H    N (Pi/3) Hbond tolerance (rad)\n"
        << "  -a_ArFF  N (Pi)   Aromatic Face-to-Face angle (rad)\n"
        << "  -at_ArFF N (Pi/6) Aromatic FF tolerance (rad)\n"
        << "  -a_ArEF  N (Pi/2) Aromatic Edge-to-Face angle (rad)\n"
        << "  -at_ArEF N (Pi/3) Aromatic EF tolerance (rad)\n"
        << "  -a_Pic   N (Pi)   Pi-Cation angle (rad)\n"
        << "  -at_Pic  N (Pi/6) Pi-Cation tolerance (rad)\n"
        << "\n"
        << "Example:\n"
        << "  IChem IFP protein.mol2 ligand.mol2\n"
        << "  IChem IFP protein.mol2 ligand.mol2 ligand_ref.mol2\n"
        << "\n"
        << "---------------------------------------------------------------------------\n\n";
}

void IChemSwitch::IFP() const {
    IFPModule::runIFP(Input_Values, Opt_Values);
}
