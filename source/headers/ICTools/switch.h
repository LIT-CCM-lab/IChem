#ifndef SWITCH_H
#define SWITCH_H
#include "ichemo.h"

#define ICHEM_VERSION "5.3.11"
#define ICHEM_RELEASE "2026-09-25"
#define ICHEM_TOOLS   12

namespace IChemTools
{
    const unsigned int REALIGN =0;
    const unsigned int BSACalc =1;
    const unsigned int IFP =2;
    const unsigned int ints=3;
    const unsigned int Grim =4;
    const unsigned int volsite =5;
    const unsigned int pdbconv=6;
    const unsigned int patch=7;
    const unsigned int utils=8;
    const unsigned int sims=9;
    const unsigned int scoring=10;
    const unsigned int atomProps = 11;
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
        const static std::string possInputs;               /*!< List of possible arguments and options, separated by spaces before AND after */
        const static std::vector<std::string> possPrgs;                 /*!< List of possible tools separated by spaces before AND after */
        std::map<std::string,std::vector<std::string> > Opt_Values;               /*!< List of Options given by the user (before the tool name) */
        std::vector<std::string> Input_Values;             /*!< List of Input parameters given by the user (after the tool name */

        void calcGrim(ICMole::Grim& grim)const;
        void AtomProps() const;
        void helpAtomProps() const;
        void printAtomInfo(const ICMole::Atom& atom) const;
        
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

        void realign()   const ;
        void BSAcalc()   const ;
        void IFP()       const ;
        void IFPScreen() const ;
        void ints() const ;
        void grim()      const ;
        void volsite()   const ;
        void PDBConversion() const ;
        void patch()   const ;
        void utils() const ;
        void runFGPS() const ;
        void scoring()   const ;
        void test() const ;
        void pharma() const ;
        void convertCav() const ;
//                      void GrScreen() const ;

        static const std::vector<std::string> FetchPossInput() ;
    public:

        IChemSwitch(const int &argc, const std::vector<std::string> &argv) ;
        IChemSwitch(const int& argc, char* argv[]) ;
        void help();
};

#endif // SWITCH_H
