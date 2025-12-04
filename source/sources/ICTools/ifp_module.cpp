#include "headers/ICTools/ifp_module.hpp"
#include "headers/ICCalcs/interaction.h"
#include "headers/ICMole/similarity.h"
#include "ichemo.h"

// #include <map>
// #include <vector>
// #include <string>
// #include <sstream>

using namespace std;
using namespace ICMole;

namespace {

using OptionMap = IFPModule::OptionMap;

// inline unsigned chooseIFPType(bool polarOnly, bool extended, bool metalOnly) {
//     // 0 - standard
//     // 1 - polar only
//     // 2 - extended
//     // 3 - extended polar
//     // 4 - metal only
//     if (!polarOnly && !extended && !metalOnly) return 0;
//     if ( polarOnly && !extended)               return 1;
//     if (!polarOnly &&  extended)               return 2;
//     if ( metalOnly)                            return 4;
//     return 3;
// }

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

inline void applyOverrides(Interactions& interactions, const InteractionOverrides& overrideInteractionParams) {
    
    // Max distances
    if (overrideInteractionParams.has_D_Hb)   interactions.setDist_H(overrideInteractionParams.D_Hb);
    if (overrideInteractionParams.has_D_Hyd)  interactions.setDist_Hyd(overrideInteractionParams.D_Hyd);
    if (overrideInteractionParams.has_D_Io)   interactions.setDist_Ionic(overrideInteractionParams.D_Io);
    if (overrideInteractionParams.has_D_Me)   interactions.setDist_Metal(overrideInteractionParams.D_Me);
    if (overrideInteractionParams.has_D_Ar)   interactions.setDist_Arom(overrideInteractionParams.D_Ar);
    if (overrideInteractionParams.has_D_WHb)  interactions.setDist_WHBond(overrideInteractionParams.D_WHb);
    if (overrideInteractionParams.has_D_Pic)  interactions.setDist_PICation(overrideInteractionParams.D_Pic);

    // Min distances
    if (overrideInteractionParams.has_d_Hb)   interactions.setMinDist_H(overrideInteractionParams.d_Hb);
    if (overrideInteractionParams.has_d_Hyd)  interactions.setMinDist_Hyd(overrideInteractionParams.d_Hyd);
    if (overrideInteractionParams.has_d_Io)   interactions.setMinDist_Ionic(overrideInteractionParams.d_Io);
    if (overrideInteractionParams.has_d_Me)   interactions.setMinDist_Metal(overrideInteractionParams.d_Me);
    if (overrideInteractionParams.has_d_Ar)   interactions.setMinDist_Arom(overrideInteractionParams.d_Ar);
    if (overrideInteractionParams.has_d_WHb)  interactions.setMinDist_WHBond(overrideInteractionParams.d_WHb);
    if (overrideInteractionParams.has_d_Pic)  interactions.setMinDist_PICation(overrideInteractionParams.d_Pic);

    // Angles
    if (overrideInteractionParams.has_a_H)      interactions.setAngl_H(overrideInteractionParams.a_H);
    if (overrideInteractionParams.has_at_H)     interactions.setAngl_Tol_H(overrideInteractionParams.at_H);
    if (overrideInteractionParams.has_a_ArFF)   interactions.setAngl_AromFF(overrideInteractionParams.a_ArFF);
    if (overrideInteractionParams.has_at_ArFF)  interactions.setAngl_Tol_AromFF(overrideInteractionParams.at_ArFF);
    if (overrideInteractionParams.has_a_ArEF)   interactions.setAngl_AromEF(overrideInteractionParams.a_ArEF);
    if (overrideInteractionParams.has_at_ArEF)  interactions.setAngl_Tol_AromEF(overrideInteractionParams.at_ArEF);
    if (overrideInteractionParams.has_a_Pic)    interactions.setAngl_PICation(overrideInteractionParams.a_Pic);
    if (overrideInteractionParams.has_at_Pic)   interactions.setAngl_Tol_PICation(overrideInteractionParams.at_Pic);
}

struct IFPOptions {
    std::string fingerprintName;

    bool includeSolvent = true;
    bool includeCofactor = true;
    bool oldHydrophobic  = true;
    bool outputBitstring = true;

    // 11 bit mask controlling which interaction types are encoded:
    // bit 0 : hydrophobic
    // bit 1 : aromatic face-to-face
    // bit 2 : aromatic edge-to-face
    // bit 3 : H-bond protein
    // bit 4 : H-bond ligand
    // bit 5 : ionic protein
    // bit 6 : ionic ligand
    // bit 7 : pi-cation
    // bit 8 : metal
    // bit 9 : weak H-bond protein
    // bit 10: weak H-bond ligand
    unsigned bitMask = 0;

    InteractionOverrides overrides;
};

struct IFPEntry {
    std::string name;
    Fingerprint fp;
    std::string fpString;
};


// At most 3 digits after the decimal point
double parseIFPNumeric(const std::string& value, const std::string& optName)
{
    if (value.empty()) {
        std::ostringstream oss;
        oss << "Option " << optName << " requires a numeric value";
        throw MoleExcept(9020101, "IChem::IFP", oss.str());
    }

    std::size_t i = 0;
    const std::size_t n = value.size();

    // Optional leading sign
    if (value[i] == '+' || value[i] == '-') {
        ++i;
        if (i == n) {
            std::ostringstream oss;
            oss << "Option " << optName << " has an invalid numeric value: " << value;
            throw MoleExcept(9020101, "IChem::IFP", oss.str());
        }
    }

    bool hasDigitsBeforeDot = false;
    while (i < n && std::isdigit(static_cast<unsigned char>(value[i]))) {
        hasDigitsBeforeDot = true;
        ++i;
    }

    // We require at least one digit before the decimal point
    if (!hasDigitsBeforeDot) {
        std::ostringstream oss;
        oss << "Option " << optName
            << " has an invalid numeric value (expected digits before decimal point): "
            << value;
        throw MoleExcept(9020101, "IChem::IFP", oss.str());
    }

    bool hasDot = false;
    int  decCount = 0;

    if (i < n && value[i] == '.') {
        hasDot = true;
        ++i;

        // At least one digit after the dot, and at most 3
        while (i < n && std::isdigit(static_cast<unsigned char>(value[i])) && decCount < 3) {
            ++decCount;
            ++i;
        }

        if (decCount == 0) {
            std::ostringstream oss;
            oss << "Option " << optName
                << " has an invalid numeric value (expected digits after decimal point): "
                << value;
            throw MoleExcept(9020101, "IChem::IFP", oss.str());
        }

        if (i < n && std::isdigit(static_cast<unsigned char>(value[i]))) {
            // There are more than 3 decimal digits
            std::ostringstream oss;
            oss << "Option " << optName
                << " accepts at most 3 digits after the decimal point: " << value;
            throw MoleExcept(9020101, "IChem::IFP", oss.str());
        }
    }

    // No extra garbage characters allowed
    if (i != n) {
        std::ostringstream oss;
        oss << "Option " << optName
            << " has an invalid numeric value: " << value;
        throw MoleExcept(9020101, "IChem::IFP", oss.str());
    }

    return std::stod(value);
}


// Parse Opt_Values into an IFPOptions struct
IFPOptions parseIFPOptions(const OptionMap& optionsValues) {
    IFPOptions options;

    // For the profile, one is mandatory
    enum class Profile {
        None,   // no profile picked yet
        Basic,
        All,
        WeakH,
        PiCat,
        Metal,
        Old
    };

    Profile profile = Profile::None;

    auto selectProfile = [&](Profile p, const std::string& optName) {
        if (profile != Profile::None) {
            throw MoleExcept(9020101,"IChem::IFP","Options --all, --weakh, --picat, --metal, --basic and --old are mutually exclusive");
        }
        profile = p;
    };

    for (const auto& [key, vals] : optionsValues) {
        const bool hasValue = !vals.empty();
        const std::string value = hasValue ? vals.front() : std::string{};

        // Options
        if (key == "-name") {
            if (!hasValue || value.empty()) {
                throw MoleExcept(9020101,"IChem::IFP","Option -name requires a non-empty value");
            }
            options.fingerprintName = value;
        }
        else if (key == "--solvent") {
            options.includeSolvent  = false;
        }
        else if (key == "--cofactor") {
            options.includeCofactor = false;
        }
        else if (key == "--newH") {
            options.oldHydrophobic  = false;
        }

        // Profiles (exactly one, can't have two at the same time)
        else if (key == "--all") {
            selectProfile(Profile::All, key);
        }
        else if (key == "--weakh") {
            selectProfile(Profile::WeakH, key);
        }
        else if (key == "--picat") {
            selectProfile(Profile::PiCat, key);
        }
        else if (key == "--metal") {
            selectProfile(Profile::Metal, key);
        }
        else if (key == "--basic") {
            selectProfile(Profile::Basic, key);
        }
        else if (key == "--old") {
            selectProfile(Profile::Old, key);
        }

        // Distances max
        else if (key == "-D_Hb") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -D_Hb requires a numeric value");
            }
            options.overrides.has_D_Hb = true;
            options.overrides.D_Hb = parseIFPNumeric(value, key);
        }
        else if (key == "-D_Hyd") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -D_Hyd requires a numeric value");
            }
            options.overrides.has_D_Hyd = true;
            options.overrides.D_Hyd = parseIFPNumeric(value, key);
        }
        else if (key == "-D_Io") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -D_Io requires a numeric value");
            }
            options.overrides.has_D_Io = true;
            options.overrides.D_Io = parseIFPNumeric(value, key);
        }
        else if (key == "-D_Me") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -D_Me requires a numeric value");
            }
            options.overrides.has_D_Me = true;
            options.overrides.D_Me = parseIFPNumeric(value, key);
        }
        else if (key == "-D_Ar") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -D_Ar requires a numeric value");
            }
            options.overrides.has_D_Ar = true;
            options.overrides.D_Ar = parseIFPNumeric(value, key);
        }
        else if (key == "-D_Pic") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -D_Pic requires a numeric value");
            }
            options.overrides.has_D_Pic = true;
            options.overrides.D_Pic = parseIFPNumeric(value, key);
        }

        else if (key == "-D_WHb") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -D_WHb requires a numeric value");
            }
            options.overrides.has_D_WHb = true;
            options.overrides.D_WHb = parseIFPNumeric(value, key);
        }

        // Distances min
        else if (key == "-d_Hb") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -d_Hb requires a numeric value");
            }
            options.overrides.has_d_Hb = true;
            options.overrides.d_Hb = parseIFPNumeric(value, key);
        }
        else if (key == "-d_Hyd") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -d_Hyd requires a numeric value");
            }
            options.overrides.has_d_Hyd = true;
            options.overrides.d_Hyd = parseIFPNumeric(value, key);
        }
        else if (key == "-d_Io") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -d_Io requires a numeric value");
            }
            options.overrides.has_d_Io = true;
            options.overrides.d_Io = parseIFPNumeric(value, key);
        }
        else if (key == "-d_Me") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -d_Me requires a numeric value");
            }
            options.overrides.has_d_Me = true;
            options.overrides.d_Me = parseIFPNumeric(value, key);
        }
        else if (key == "-d_Ar") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -d_Ar requires a numeric value");
            }
            options.overrides.has_d_Ar = true;
            options.overrides.d_Ar = parseIFPNumeric(value, key);
        }

        else if (key == "-d_Pic") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -d_Pic requires a numeric value");
            }
            options.overrides.has_d_Pic = true;
            options.overrides.d_Pic = parseIFPNumeric(value, key);
        }

        else if (key == "-d_WHb") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -d_WHb requires a numeric value");
            }
            options.overrides.has_d_WHb = true;
            options.overrides.d_WHb = parseIFPNumeric(value, key);
        }

        // Angles
        else if (key == "-a_H") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -a_H requires a numeric value");
            }
            options.overrides.has_a_H = true;
            options.overrides.a_H = parseIFPNumeric(value, key);
        }
        else if (key == "-at_H") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -at_H requires a numeric value");
            }
            options.overrides.has_at_H = true;
            options.overrides.at_H = parseIFPNumeric(value, key);
        }
        else if (key == "-a_ArFF") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -a_ArFF requires a numeric value");
            }
            options.overrides.has_a_ArFF = true;
            options.overrides.a_ArFF = parseIFPNumeric(value, key);
        }
        else if (key == "-at_ArFF") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -at_ArFF requires a numeric value");
            }
            options.overrides.has_at_ArFF = true;
            options.overrides.at_ArFF = parseIFPNumeric(value, key);
        }
        else if (key == "-a_ArEF") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -a_ArEF requires a numeric value");
            }
            options.overrides.has_a_ArEF = true;
            options.overrides.a_ArEF = parseIFPNumeric(value, key);
        }
        else if (key == "-at_ArEF") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -at_ArEF requires a numeric value");
            }
            options.overrides.has_at_ArEF = true;
            options.overrides.at_ArEF = parseIFPNumeric(value, key);
        }
        else if (key == "-a_Pic") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -a_Pic requires a numeric value");
            }
            options.overrides.has_a_Pic = true;
            options.overrides.a_Pic = parseIFPNumeric(value, key);
        }
        else if (key == "-at_Pic") {
            if (!hasValue) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -at_Pic requires a numeric value");
            }
            options.overrides.has_at_Pic = true;
            options.overrides.at_Pic = parseIFPNumeric(value, key);
        }

        // Unknown option => return exception
        else {
            std::ostringstream oss;
            oss << "Unknown IFP option: " << key
                << ". Allowed options include: "
                << "--all, --basic, --weakh, --picat, --metal, --old, "
                << "--solvent, --cofactor, --newH, "
                << "-name, -D_*, -d_*, -a_*, -at_*";
            throw MoleExcept(9020101, "IChem::IFP", oss.str());
        }
    }

    // Build the 11 bit mask according to the chosen profile
    constexpr unsigned BIT_HYDRO   = 1u << 0;
    constexpr unsigned BIT_AR_FF   = 1u << 1;
    constexpr unsigned BIT_AR_EF   = 1u << 2;
    constexpr unsigned BIT_HB_P    = 1u << 3;
    constexpr unsigned BIT_HB_L    = 1u << 4;
    constexpr unsigned BIT_ION_P   = 1u << 5;
    constexpr unsigned BIT_ION_L   = 1u << 6;
    constexpr unsigned BIT_PIC     = 1u << 7;
    constexpr unsigned BIT_METAL   = 1u << 8;
    constexpr unsigned BIT_WH_P    = 1u << 9;
    constexpr unsigned BIT_WH_L    = 1u << 10;

    constexpr unsigned MASK_BASE   = BIT_HYDRO | BIT_AR_FF | BIT_AR_EF
                                   | BIT_HB_P  | BIT_HB_L
                                   | BIT_ION_P | BIT_ION_L;

    constexpr unsigned MASK_WEAKH  = BIT_WH_P | BIT_WH_L;
    constexpr unsigned MASK_PIC    = BIT_PIC;
    constexpr unsigned MASK_METAL  = BIT_METAL;
    constexpr unsigned MASK_ALL    = MASK_BASE | MASK_PIC | MASK_METAL | MASK_WEAKH;

    // old 7 bit layout
    constexpr unsigned FLAG_OLD_LAYOUT = 1u << 31;

    switch (profile) {
        case Profile::None:
            // If no profile: forbidden: the user MUST pick one interaction profile
            throw MoleExcept(9020101,"IChem::IFP","You must specify one of --all, --basic, --weakh, --picat, --metal or --old");
        case Profile::Basic:
            // 11-bit layout, only base bits allowed
            options.bitMask = MASK_BASE;
            break;
        case Profile::All:
            options.bitMask = MASK_ALL;
            break;
        case Profile::WeakH:
            options.bitMask = MASK_WEAKH;
            break;
        case Profile::PiCat:
            options.bitMask = MASK_PIC;
            break;
        case Profile::Metal:
            options.bitMask = MASK_METAL;
            break;
        case Profile::Old:
            // Same base interactions, but mark that we want the old 7 bit layout
            options.bitMask = MASK_BASE | FLAG_OLD_LAYOUT;
            break;
    }

        // Check consistency of min/max pairs when both are provided => min can't be > to max
    auto checkInterval = [&](bool hasMin, double minVal,
                             bool hasMax, double maxVal,
                             const char* label) {
        if (hasMin && hasMax && minVal > maxVal) {
            std::ostringstream oss;
            oss << "For " << label
                << " minimal distance (" << minVal
                << ") is greater than maximal distance (" << maxVal << ")";
            throw MoleExcept(9020101, "IChem::IFP", oss.str());
        }
    };

    checkInterval(options.overrides.has_d_Hb,  options.overrides.d_Hb,
                  options.overrides.has_D_Hb,  options.overrides.D_Hb,
                  "H-bond (-d_Hb / -D_Hb)");

    checkInterval(options.overrides.has_d_Hyd, options.overrides.d_Hyd,
                  options.overrides.has_D_Hyd, options.overrides.D_Hyd,
                  "Hydrophobic (-d_Hyd / -D_Hyd)");

    checkInterval(options.overrides.has_d_Io,  options.overrides.d_Io,
                  options.overrides.has_D_Io,  options.overrides.D_Io,
                  "Ionic (-d_Io / -D_Io)");

    checkInterval(options.overrides.has_d_Me,  options.overrides.d_Me,
                  options.overrides.has_D_Me,  options.overrides.D_Me,
                  "Metal (-d_Me / -D_Me)");

    checkInterval(options.overrides.has_d_Ar,  options.overrides.d_Ar,
                  options.overrides.has_D_Ar,  options.overrides.D_Ar,
                  "Aromatic (-d_Ar / -D_Ar)");

    checkInterval(options.overrides.has_d_Pic,  options.overrides.d_Pic,
                  options.overrides.has_D_Pic,  options.overrides.D_Pic,
                  "Aromatic (-d_Pic / -D_Pic)");

    checkInterval(options.overrides.has_d_WHb,  options.overrides.d_WHb,
                  options.overrides.has_D_WHb,  options.overrides.D_WHb,
                  "Weak H-bond (-d_WHb / -D_WHb)"); 
                             

    return options;
}

// Apply residue rules once
inline void configureResidueRules(const IFPOptions& options) {
    if (options.includeSolvent) {
        Residu::Rules[MoleType::PROTEIN][ResType::WATER] = MoleType::PROTEIN;
    }
    if (options.includeCofactor) {
        Residu::Rules[MoleType::PROTEIN][ResType::COFACTOR] = MoleType::PROTEIN;
    }
}

void computeIFPForLigand(Interactions& interactions, Molecule& ligand, const IFPOptions& options, InterResults& output) {
   
    // RingDetection not run yet
    if (ligand.firstCycle() == ligand.lastCycle()) {
        ligand.ringPerception();
    }

    interactions.detectInteractions(ligand, output, true, options.oldHydrophobic);
    interactions.genIFP(output, options.bitMask);
}

// 3 arguments mode (protein, ligand, ligand_ref)
// We compute IFPs for all reference ligands first, then all query ligands
void processReferenceMode(const std::string& ligandFile, const std::string& refLigandFile, Interactions& interactions, const IFPOptions& options)
{
    std::vector<IFPEntry> docked; // ligandFile
    std::vector<IFPEntry> refs;   // refLigandFile

    MoleReader reader;

    // docked ligands
    reader.loadNewFile(ligandFile);
    while (!reader.isEOF()) {
        Molecule ligand;
        reader.loadNextMolecule(ligand, MoleType::LIGAND);
        ligand.checkMOL2();

        InterResults result;
        computeIFPForLigand(interactions, ligand, options, result);

        docked.push_back(IFPEntry{
            ligand.getName(),
            result.IFP,
            result.IFPString
        });
    }

    // reference ligands
    reader.loadNewFile(refLigandFile);
    while (!reader.isEOF()) {
        Molecule ligand;
        reader.loadNextMolecule(ligand, MoleType::LIGAND);
        ligand.checkMOL2();

        InterResults result;
        computeIFPForLigand(interactions, ligand, options, result);

        refs.push_back(IFPEntry{
            ligand.getName(),
            result.IFP,
            result.IFPString
        });
    }

    // Print all IFPs
    if (options.outputBitstring) {
        // docked ligands
        for (const auto& d : docked) {
            cout << d.name << '\t' << d.fpString     << '\n'
                 << d.name << '\t' << d.fp.toString() << '\n';
        }
        // all references
        for (const auto& r : refs) {
            cout << r.name << '\t' << r.fpString     << '\n'
                 << r.name << '\t' << r.fp.toString() << '\n';
        }
    }

    // Similarities: each docked vs ref
    Similarity sims(false);
    for (const auto& d : docked) {
        sims.setRef(d.fp);
        for (const auto& r : refs) {
            sims.setComp(r.fp);
            cout << d.name << '\t'
                 << r.name << '\t'
                 << sims.Tanimoto() << '\n';
        }
    }
}


// 2 arguments mode: protein + ligand file (single or multi ligand)
void processLigandMode(const std::string& ligandFile, Complex& complex, Interactions& interactions, const IFPOptions& options) {
    
    MoleReader reader;
    reader.loadNewFile(ligandFile);
    const std::size_t numLigands = reader.getNumMolecules();

    if (numLigands == 1) {
        // Single ligand, loaded in the complex
        reader.loadInComplex(complex, MoleType::LIGAND);

        if (complex.getMole(MoleType::LIGAND) == nullptr)
            throw MoleExcept(9020102, "IChem::IFP", "No ligand found in " + ligandFile);

        Molecule& ligand = *complex.getMole(MoleType::LIGAND);
        ligand.checkMOL2();

        InterResults interRes;
        computeIFPForLigand(interactions, ligand, options, interRes);

        cout << interactions.toString(interRes) << '\n';
        if (options.outputBitstring) {
            cout << interRes.IFPString << '\n' // Residue list
                 << interRes.IFP.toString() << '\n'; // IFP
        }
    }
    else {
        // Multi ligand file, not loaded in the complex
        reader.loadNewFile(ligandFile);

        std::vector<IFPEntry> entries;
        entries.reserve(numLigands);

        while (!reader.isEOF()) {
            Molecule ligand;
            reader.loadNextMolecule(ligand, MoleType::LIGAND);
            ligand.checkMOL2();

            InterResults res;
            computeIFPForLigand(interactions, ligand, options, res);

            entries.push_back(IFPEntry{
                ligand.getName(),
                res.IFP,
                res.IFPString
            });
        }

        if (options.outputBitstring) {
            for (const auto& entry : entries) {
                cout << entry.name << '\t' << entry.fpString << '\n'
                     << entry.name << '\t' << entry.fp.toString() << '\n';
            }
        }
    }
}

// Used for the 4 argument mode
std::vector<IFPEntry> computeIFPsFromFiles(const std::string& proteinFile, const std::string& ligandFile, const IFPOptions&  options) {
    
    std::vector<IFPEntry> entries;

    // Protein
    Complex complex;
    MoleReader proteinReader;

    proteinReader.loadNewFile(proteinFile);
    proteinReader.get_format_file();
    proteinReader.loadInComplex(complex, MoleType::PROTEIN);

    if (complex.getMole(MoleType::PROTEIN) == nullptr)
        throw MoleExcept(9020102, "IChem::IFP", "No protein found in " + proteinFile);
    
    // Interactions object for this protein
    Interactions interactions(complex);
    applyOverrides(interactions, options.overrides);

    // 1 ou many ligands
    MoleReader ligandReader;
    ligandReader.loadNewFile(ligandFile);

    while (!ligandReader.isEOF()) {
        
        Molecule ligand;
        ligandReader.loadNextMolecule(ligand, MoleType::LIGAND);
        ligand.checkMOL2();

        InterResults result;
        computeIFPForLigand(interactions, ligand, options, result);

        IFPEntry entry;
        entry.name = ligand.getName();
        entry.fp = result.IFP;
        entry.fpString = result.IFPString;

        entries.push_back(std::move(entry));
    }

    if (entries.empty()) {
        throw MoleExcept(9020102, "IChem::IFP", "No ligand found in " + ligandFile);
    }

    return entries;
}


void processFourArgMode(const std::string& protein1File, const std::string& ligand1File, const std::string& protein2File, const std::string& ligand2File, const IFPOptions&  options) {
    
    // Convention:
    // (protein1, ligand1) = DOCKED (query)
    // (protein2, ligand2) = REF (reference)

    std::vector<IFPEntry> dockedEntries = computeIFPsFromFiles(protein1File, ligand1File, options);
    std::vector<IFPEntry> refEntries = computeIFPsFromFiles(protein2File, ligand2File, options);

    const std::size_t numDocked = dockedEntries.size();
    const std::size_t numRef = refEntries.size();

    if (options.outputBitstring) {
        // DOCKED block
        for (const auto& entry : dockedEntries) {
            std::cout << entry.name << '\t' << entry.fpString      << '\n'
                      << entry.name << '\t' << entry.fp.toString() << '\n';
        }

        // REF block
        for (const auto& entry : refEntries) {
            std::cout << entry.name << '\t' << entry.fpString      << '\n'
                      << entry.name << '\t' << entry.fp.toString() << '\n';
        }
    }

    // Similarity
    Similarity sims(false);

    for (std::size_t i = 0; i < numDocked; ++i) {
        const auto& docked = dockedEntries[i];

        for (std::size_t j = 0; j < numRef; ++j) {
            const auto& ref = refEntries[j];

            // Fingerprint lengths must match
            if (docked.fpString.size() != ref.fpString.size())
                throw MoleExcept( 9020105, "IChem::IFP", "Cannot compare fingerprints with different sizes in 4-argument mode");

            sims.setRef(ref.fp);
            sims.setComp(docked.fp);

            std::cout << docked.name << '\t' << ref.name    << '\t' << sims.Tanimoto() << '\n';
        }
    }
}



} // anonymous namespace

// unsigned chooseIFPType_forTests(bool polarOnly, bool extended, bool metalOnly) {
//     return chooseIFPType(polarOnly, extended, metalOnly);
// }

namespace IFPModule {

void runIFP(const std::vector<std::string> &Input_Values, const OptionMap &Opt_Values) {
    
    const std::size_t inputSize = Input_Values.size();

    if (inputSize < 2 || inputSize > 4)
        throw MoleExcept(9020101, "IChem::IFP", "Number of parameters must be 2, 3 or 4");

    const IFPOptions options = parseIFPOptions(Opt_Values);

    try {
        // Global rules, once per call
        Residu::loadRules();
        Molecule::loadRules();
        configureResidueRules(options);

        // 4 argument mode
        if (inputSize == 4) {
            const std::string& protein1File = Input_Values[0];
            const std::string& ligand1File  = Input_Values[1];
            const std::string& protein2File = Input_Values[2];
            const std::string& ligand2File  = Input_Values[3];

            processFourArgMode(protein1File, ligand1File, protein2File, ligand2File, options);
            return;
        }

        // 2 and 3 argument modes
        const std::string& proteinFile = Input_Values.at(0);
        const std::string& ligandFile = Input_Values.at(1);
        const bool withReference = (inputSize == 3);

        Complex complex;
        MoleReader proteinReader;
        proteinReader.loadNewFile(proteinFile);
        proteinReader.get_format_file();
        proteinReader.loadInComplex(complex, MoleType::PROTEIN);

        if (complex.getMole(MoleType::PROTEIN) == nullptr)
            throw MoleExcept(9020102, "IChem::IFP", "No protein found in " + proteinFile);

        Interactions interactions(complex);
        applyOverrides(interactions, options.overrides);

        if (withReference) {
            const std::string& refLigandFile = Input_Values.at(2);
            processReferenceMode(ligandFile, refLigandFile, interactions, options);
        } else {
            processLigandMode(ligandFile, complex, interactions, options);
        }
    }
    catch (MoleExcept& exception) {
        exception.addTrace("IChem::IFP");
        throw;
    }
}

} // namespace IFPModule