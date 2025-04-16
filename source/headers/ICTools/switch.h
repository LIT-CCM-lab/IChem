#ifndef SWITCH_H
#define SWITCH_H
#include "ichemo.h"

#define ICHEM_VERSION "5.3.0"
#define ICHEM_RELEASE "2024-12-09"
#define ICHEM_TOOLS   14

namespace IChemTools
{
    const unsigned int GENKEY  =0;
    const unsigned int REALIGN =1;
    const unsigned int BSACalc =2;
    const unsigned int IFP =3;
    const unsigned int ints=4;
    const unsigned int Grim =5;
    const unsigned int volsite =6;
    const unsigned int pdbconv=7;
    const unsigned int patch=8;
    const unsigned int utils=9;
    const unsigned int sims=10;
    const unsigned int scoring=11;
    const unsigned int atomProps = 12;
//    const unsigned int cavLig=13;
}
class IChemSwitch;

class IChemSwitch
{
    public:

        typedef void (IChemSwitch::*pf1 ) () const;

    private:

        static pf1 listFunc[ICHEM_TOOLS];
        static pf1 listHelp[ICHEM_TOOLS];
        static bool verbose;                  /*!< Verbose mode used for standard output */
        static bool licence_loaded;           /*!< Boolean to check whether the licence is loaded */
        static bool allowed[ICHEM_TOOLS];     /*!< Boolean array. True when the tool is allowed. False when not */
        const static std::string moduluskeys;              /*!< For RSA Encryption */
        const static std::string privatexp;                /*!< For RSA Encryption */
        const static std::string publicexp;                /*!< For RSA Encryption */
        const static std::string possInputs;               /*!< List of possible arguments and options, separated by spaces before AND after */
        const static std::vector<std::string> possPrgs;                 /*!< List of possible tools separated by spaces before AND after */
        std::map<std::string,std::vector<std::string> > Opt_Values;               /*!< List of Options given by the user (before the tool name) */
        std::vector<std::string> Input_Values;             /*!< List of Input parameters given by the user (after the tool name */

        void calcGrim(ICMole::Grim& grim)const;
        void AtomProps() const;
        void helpAtomProps() const;
        void printAtomInfo(const ICMole::Atom& atom) const;
        
        void helpGenKey()    const;
        void helpRealign()   const;
        void helpBSAcalc()   const;
        void helpIFP()       const;
        void helpIFPScreen() const;
        void helpints()      const;
        void helpgrim()      const;
        void helpVolSite()   const;
        void helpPDBConvert() const;
        void helpPatch() const;
        void helpUtils() const;
        void helpFGPS() const;
        void helpScoring() const;
        void helpdetectPPI() const;
        void helpPharma() const;
        void helpconvertCav() const;
        void helpTest() const;
//                      void helpGrScreen() const;

        void genKey()    const  throw(ICMole::MoleExcept);
        void realign()   const throw(ICMole::MoleExcept);
        void BSAcalc()   const throw(ICMole::MoleExcept);
        void IFP()       const throw(ICMole::MoleExcept);
        void IFPScreen() const throw(ICMole::MoleExcept);
        void ints() const throw(ICMole::MoleExcept);
        void grim()      const throw(ICMole::MoleExcept);
        void volsite()   const throw(ICMole::MoleExcept);
        void PDBConversion() const throw(ICMole::MoleExcept);
        void patch()   const throw(ICMole::MoleExcept);
        void utils() const throw(ICMole::MoleExcept);
        void runFGPS() const throw(ICMole::MoleExcept);
        void scoring()   const throw(ICMole::MoleExcept);
        void test() const throw(ICMole::MoleExcept);
        void pharma() const throw(ICMole::MoleExcept);
        void convertCav() const throw(ICMole::MoleExcept);
//                      void GrScreen() const throw(ICMole::MoleExcept);

        void checkLicence()    throw(ICMole::MoleExcept);
        void yourlicence();

        static const std::vector<std::string> FetchPossInput() ;
    public:

        IChemSwitch(const int &argc, const std::vector<std::string> &argv) throw(ICMole::MoleExcept);
        IChemSwitch(const int& argc, char* argv[]) throw(ICMole::MoleExcept);
        void help();
};

#endif // SWITCH_H
