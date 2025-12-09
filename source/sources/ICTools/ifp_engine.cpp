#include "headers/ICTools/ifp_engine.hpp"

#include "headers/ICMole/similarity.h"
#include <sstream>
#include <cctype>

using namespace std;
using namespace ICMole;

namespace IFPInternal {

// Numeric parser
static double parseIFPNumeric(const std::string& value, const std::string& optName) {
    
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

    if (!hasDigitsBeforeDot) {
        std::ostringstream oss;
        oss << "Option " << optName << " has an invalid numeric value (expected digits before decimal point): " << value;
        throw MoleExcept(9020101, "IChem::IFP", oss.str());
    }

    int  decCount = 0;

    if (i < n && value[i] == '.') {
        ++i;

        while (i < n && std::isdigit(static_cast<unsigned char>(value[i])) && decCount < 3) {
            ++decCount;
            ++i;
        }

        if (decCount == 0) {
            std::ostringstream oss;
            oss << "Option " << optName << " has an invalid numeric value (expected digits after decimal point): " << value;
            throw MoleExcept(9020101, "IChem::IFP", oss.str());
        }

        if (i < n && std::isdigit(static_cast<unsigned char>(value[i]))) {
            std::ostringstream oss;
            oss << "Option " << optName << " accepts at most 3 digits after the decimal point: " << value;
            throw MoleExcept(9020101, "IChem::IFP", oss.str());
        }
    }

    if (i != n) {
        std::ostringstream oss;
        oss << "Option " << optName << " has an invalid numeric value: " << value;
        throw MoleExcept(9020101, "IChem::IFP", oss.str());
    }

    return std::stod(value);
}

// Overrides for our options 
void applyOverrides(Interactions& interactions, const InteractionOverrides& overrideInteractionParams)
{
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

// Parse options
IFPOptions parseIFPOptions(const OptionMap& optionsValues) {
    
    IFPOptions options;

    enum class Profile {
        None,
        Basic,
        All,
        WeakH,
        PiCat,
        Metal,
        Old
    };

    Profile profile = Profile::None;

    auto selectProfile = [&](Profile p) {
        if (profile != Profile::None) {
            throw MoleExcept(9020101, "IChem::IFP", "Options --all, --weakh, --picat, --metal, --basic and --old are mutually exclusive");
        }
        profile = p;
    };

    for (const auto& keyvalue : optionsValues) {
        const std::string& key = keyvalue.first;
        const auto& vals = keyvalue.second;
        const bool hasValue = !vals.empty();
        const std::string value = hasValue ? vals.front() : std::string{};

        // Generic options
        if (key == "-name") {
            if (!hasValue || value.empty()) {
                throw MoleExcept(9020101, "IChem::IFP", "Option -name requires a non-empty value");
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

        // Profiles
        else if (key == "--all") {
            selectProfile(Profile::All);
        }
        else if (key == "--weakh") {
            selectProfile(Profile::WeakH);
        }
        else if (key == "--picat") {
            selectProfile(Profile::PiCat);
        }
        else if (key == "--metal") {
            selectProfile(Profile::Metal);
        }
        else if (key == "--basic") {
            selectProfile(Profile::Basic);
        }
        else if (key == "--old") {
            selectProfile(Profile::Old);
        }

        // Max distances
        else if (key == "-D_Hb") {
            options.overrides.has_D_Hb = true;
            options.overrides.D_Hb = parseIFPNumeric(value, key);
        }
        else if (key == "-D_Hyd") {
            options.overrides.has_D_Hyd = true;
            options.overrides.D_Hyd = parseIFPNumeric(value, key);
        }
        else if (key == "-D_Io") {
            options.overrides.has_D_Io = true;
            options.overrides.D_Io = parseIFPNumeric(value, key);
        }
        else if (key == "-D_Me") {
            options.overrides.has_D_Me = true;
            options.overrides.D_Me = parseIFPNumeric(value, key);
        }
        else if (key == "-D_Ar") {
            options.overrides.has_D_Ar = true;
            options.overrides.D_Ar = parseIFPNumeric(value, key);
        }
        else if (key == "-D_Pic") {
            options.overrides.has_D_Pic = true;
            options.overrides.D_Pic = parseIFPNumeric(value, key);
        }
        else if (key == "-D_WHb") {
            options.overrides.has_D_WHb = true;
            options.overrides.D_WHb = parseIFPNumeric(value, key);
        }

        // Min distances
        else if (key == "-d_Hb") {
            options.overrides.has_d_Hb = true;
            options.overrides.d_Hb = parseIFPNumeric(value, key);
        }
        else if (key == "-d_Hyd") {
            options.overrides.has_d_Hyd = true;
            options.overrides.d_Hyd = parseIFPNumeric(value, key);
        }
        else if (key == "-d_Io") {
            options.overrides.has_d_Io = true;
            options.overrides.d_Io = parseIFPNumeric(value, key);
        }
        else if (key == "-d_Me") {
            options.overrides.has_d_Me = true;
            options.overrides.d_Me = parseIFPNumeric(value, key);
        }
        else if (key == "-d_Ar") {
            options.overrides.has_d_Ar = true;
            options.overrides.d_Ar = parseIFPNumeric(value, key);
        }
        else if (key == "-d_Pic") {
            options.overrides.has_d_Pic = true;
            options.overrides.d_Pic = parseIFPNumeric(value, key);
        }
        else if (key == "-d_WHb") {
            options.overrides.has_d_WHb = true;
            options.overrides.d_WHb = parseIFPNumeric(value, key);
        }

        // Angles
        else if (key == "-a_H") {
            options.overrides.has_a_H = true;
            options.overrides.a_H = parseIFPNumeric(value, key);
        }
        else if (key == "-at_H") {
            options.overrides.has_at_H = true;
            options.overrides.at_H = parseIFPNumeric(value, key);
        }
        else if (key == "-a_ArFF") {
            options.overrides.has_a_ArFF = true;
            options.overrides.a_ArFF = parseIFPNumeric(value, key);
        }
        else if (key == "-at_ArFF") {
            options.overrides.has_at_ArFF = true;
            options.overrides.at_ArFF = parseIFPNumeric(value, key);
        }
        else if (key == "-a_ArEF") {
            options.overrides.has_a_ArEF = true;
            options.overrides.a_ArEF = parseIFPNumeric(value, key);
        }
        else if (key == "-at_ArEF") {
            options.overrides.has_at_ArEF = true;
            options.overrides.at_ArEF = parseIFPNumeric(value, key);
        }
        else if (key == "-a_Pic") {
            options.overrides.has_a_Pic = true;
            options.overrides.a_Pic = parseIFPNumeric(value, key);
        }
        else if (key == "-at_Pic") {
            options.overrides.has_at_Pic = true;
            options.overrides.at_Pic = parseIFPNumeric(value, key);
        }

        // Unknown
        else {
            std::ostringstream oss;
            oss << "Unknown IFP option: " << key
                << ". Allowed options include: "
                << "--all, --basic, --weakh, --picat, --metal, --old, "
                << "--solvent, --cofactor, --newH,"
                << "-name, -D_*, -d_*, -a_*, -at_*";
            throw MoleExcept(9020101, "IChem::IFP", oss.str());
        }
    }

    // Bitmask
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

    constexpr unsigned FLAG_OLD_LAYOUT = 1u << 31;

    switch (profile) {
        case Profile::None:
            throw MoleExcept(9020101, "IChem::IFP", "You must specify one of --all, --basic, --weakh, --picat, --metal or --old");
        case Profile::Basic:
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
            options.bitMask = MASK_BASE | FLAG_OLD_LAYOUT;
            break;
    }

    // min/max consistency checks
    auto checkInterval = [&](bool hasMin, double minVal,
                             bool hasMax, double maxVal,
                             const char* label)
    {
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

    checkInterval(options.overrides.has_d_Pic, options.overrides.d_Pic,
                  options.overrides.has_D_Pic, options.overrides.D_Pic,
                  "Aromatic (-d_Pic / -D_Pic)");

    checkInterval(options.overrides.has_d_WHb,  options.overrides.d_WHb,
                  options.overrides.has_D_WHb,  options.overrides.D_WHb,
                  "Weak H-bond (-d_WHb / -D_WHb)");

    return options;
}

// Residue rules
void configureResidueRules(const IFPOptions& options)
{
    if (options.includeSolvent) {
        Residu::Rules[MoleType::PROTEIN][ResType::WATER] = MoleType::PROTEIN;
    }
    if (options.includeCofactor) {
        Residu::Rules[MoleType::PROTEIN][ResType::COFACTOR] = MoleType::PROTEIN;
    }
}

// Compute interactions + IFP for a single ligand
void computeIFPForLigand(Interactions& interactions, Molecule& ligand, const IFPOptions& options, InterResults& output) {
    if (ligand.firstCycle() == ligand.lastCycle()) {
        ligand.ringPerception();
    }

    interactions.detectInteractions(ligand, output, true, options.oldHydrophobic);
    interactions.genIFP(output, options.bitMask);
}

// (protein, ligandFile) -> IFPEntry list, used by 4 argument mode 
std::vector<IFPEntry> computeIFPsFromFiles(const std::string& proteinFile, const std::string& ligandFile, const IFPOptions& options) {
    
    std::vector<IFPEntry> entries;

    Complex complex;
    MoleReader proteinReader;

    proteinReader.loadNewFile(proteinFile);
    proteinReader.get_format_file();
    proteinReader.loadInComplex(complex, MoleType::PROTEIN);

    if (complex.getMole(MoleType::PROTEIN) == nullptr) {
        throw MoleExcept(9020102, "IChem::IFP", "No protein found in " + proteinFile);
    }

    Interactions interactions(complex);
    applyOverrides(interactions, options.overrides);

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

} // namespace IFPInternal
