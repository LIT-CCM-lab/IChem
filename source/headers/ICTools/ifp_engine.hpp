#ifndef IFP_ENGINE_HPP
#define IFP_ENGINE_HPP

#include <string>
#include <vector>

#include "headers/ICTools/ifp_module.hpp"
#include "headers/ICCalcs/interaction.h"
#include "ichemo.h"

namespace IFPInternal {

// Reuse the same OptionMap as the CLI module
using OptionMap = IFPModule::OptionMap;

using ICMole::Interactions;
using ICMole::InterResults;
using ICMole::Molecule;
using ICMole::Fingerprint;

// Threshold overrides
struct InteractionOverrides {
    // Max distances
    bool has_D_Hb   = false; double D_Hb  = 0.0;
    bool has_D_Hyd  = false; double D_Hyd = 0.0;
    bool has_D_Io   = false; double D_Io  = 0.0;
    bool has_D_Me   = false; double D_Me  = 0.0;
    bool has_D_Ar   = false; double D_Ar  = 0.0;
    bool has_D_Pic  = false; double D_Pic = 0.0;
    bool has_D_WHb  = false; double D_WHb = 0.0;

    // Min distances
    bool has_d_Hb   = false; double d_Hb  = 0.0;
    bool has_d_Hyd  = false; double d_Hyd = 0.0;
    bool has_d_Io   = false; double d_Io  = 0.0;
    bool has_d_Me   = false; double d_Me  = 0.0;
    bool has_d_Ar   = false; double d_Ar  = 0.0;
    bool has_d_Pic  = false; double d_Pic = 0.0;
    bool has_d_WHb  = false; double d_WHb = 0.0;

    // Angles
    bool has_a_H    = false; double a_H    = 0.0;
    bool has_at_H   = false; double at_H   = 0.0;
    bool has_a_ArFF = false; double a_ArFF = 0.0;
    bool has_at_ArFF= false; double at_ArFF= 0.0;
    bool has_a_ArEF = false; double a_ArEF = 0.0;
    bool has_at_ArEF= false; double at_ArEF= 0.0;
    bool has_a_Pic  = false; double a_Pic  = 0.0;
    bool has_at_Pic = false; double at_Pic = 0.0;
};

// Options used by both CLI and API
struct IFPOptions {
    std::string fingerprintName;

    bool includeSolvent = true;
    bool includeCofactor = true;
    bool oldHydrophobic  = true;

    // 11 bit mask
    unsigned bitMask = 0;

    InteractionOverrides overrides;
};

// One ligand IFP result
struct IFPEntry {
    std::string name;
    Fingerprint fp;
    std::string fpString;
};

// Parse CLI options into IFPOptions
IFPOptions parseIFPOptions(const OptionMap& optionsValues);

// Apply overrides to an Interactions object
void applyOverrides(Interactions& interactions, const InteractionOverrides& overrideInteractionParams);

// Apply solvent/cofactor residue rules
void configureResidueRules(const IFPOptions& options);

// Compute interactions + IFP for a single ligand
void computeIFPForLigand(Interactions& interactions, Molecule& ligand, const IFPOptions& options, InterResults& output);

// Utility to compute all IFPs from (protein, ligandFile)
// Used by both CLI 4 arg mode and the API
std::vector<IFPEntry> computeIFPsFromFiles(const std::string& proteinFile, const std::string& ligandFile, const IFPOptions&  options);

} // namespace IFPInternal

#endif // IFP_ENGINE_HPP
