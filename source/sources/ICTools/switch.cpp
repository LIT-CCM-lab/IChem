#include "headers/ICTools/switch.h"
#include <sstream>

using namespace std;
using namespace ICMole;



//IChemSwitch::pf1 te=&IChemSwitch::test;

bool IChemSwitch::verbose       =false;
const std::string IChemSwitch::possInputs    =
        " realign BSAcalc IFP grim ints fgps volsite pdbconv patch utils sims scoring \
        -type -name -logf -d_Hb -d_Hyd -d_Io -d_Me -d_Ar -a_H -at_H -a_ArFF -at_ArFF -a_ArEF -at_ArEF --enf --noMerge --d\
        -boxStep -name \
        -step -boxS -b -n -nPTS -l -p -pha -phaC -phaP --weight --tol --outpoutsite --cofactor --dna --solvent --ligand  --hydrogen --svm --desc -drog \
        --wMOL2 --wUnDrug --HARMSIZE  --MSEMET --CSECYS --MOVHET --ALTATM --NUMATM --UPDMAS --process \
        --BONDSE --CLNUNW --MOL2TY --SPLITM -SelChain --noLig --SELWAT --SELLIG --TOMOL2 --RESTYP \
        -mutant -bornes -c -c1 -c2 -n \
        --wInts --small --binary -metric";

        const vector<std::string> IChemSwitch::possPrgs=IChemSwitch::FetchPossInput();

const vector<string> IChemSwitch::FetchPossInput()
{
    vector<string> RX;
    RX.push_back("realign");
    RX.push_back("BSAcalc");
    RX.push_back("IFP");
    RX.push_back("ints");
    RX.push_back("grim");
    RX.push_back("volsite");
    RX.push_back("pdbconv");
    RX.push_back("patch");
    RX.push_back("utils");
    RX.push_back("sims");
    RX.push_back("scoring");
    RX.push_back("AtomProps");
//    RX.push_back("cavLig");
    return RX;
}


IChemSwitch::pf1 IChemSwitch::listFunc[ICHEM_TOOLS]={&IChemSwitch::realign,
                                                     &IChemSwitch::BSAcalc,
                                                     &IChemSwitch::IFP,
                                                     &IChemSwitch::ints,
                                                     &IChemSwitch::grim,
                                                     &IChemSwitch::volsite,
                                                     &IChemSwitch::PDBConversion,
                                                     &IChemSwitch::patch,
                                                     &IChemSwitch::utils,
                                                     &IChemSwitch::runFGPS,
                                                     &IChemSwitch::scoring,
                                                     &IChemSwitch::AtomProps
                                                };
//                                                     &IChemSwitch::convertCav};
IChemSwitch::pf1 IChemSwitch::listHelp[ICHEM_TOOLS]={&IChemSwitch::helpRealign,
                                                     &IChemSwitch::helpBSAcalc,
                                                     &IChemSwitch::helpIFP,
                                                     &IChemSwitch::helpints,
                                                     &IChemSwitch::helpgrim,
                                                     &IChemSwitch::helpVolSite,
                                                     &IChemSwitch::helpPDBConvert,
                                                     &IChemSwitch::helpPatch,
                                                     &IChemSwitch::helpUtils,
                                                     &IChemSwitch::helpFGPS,
                                                     &IChemSwitch::helpScoring,
                                                     &IChemSwitch::helpAtomProps
                                                     };
//                                                     &IChemSwitch::helpconvertCav};


IChemSwitch::IChemSwitch(const int& argc, char *argv[]) 
{
    if (argc == 1) {help();return;}

    // Global aromaticity model: reset on every invocation (-F batch mode reuses
    // the process), --oldAro switches back to the legacy bond-count rules.
    Molecule::setAromaticityMode(AromaticityMode::SP2_PLANAR);

    int TaskType=-1;
    std::string tmpStr,tmpStr2;bool opts=true;
    vector<string> options;
    for (int i=1; i < argc; i++)
    {

        tmpStr = argv[i];
        if (!opts) {Input_Values.push_back(tmpStr); continue;}
        // Boolean option :
        if (tmpStr.substr(0,2)=="--")
        {
            if (tmpStr =="--verb") {verbose=true;continue;}
            if (tmpStr =="--oldAro")
            {
                Molecule::setAromaticityMode(AromaticityMode::BOND_COUNT);
                continue;
            }
            options.clear();
            options.push_back("TRUE");
            Opt_Values.insert(pair<string,vector<string> >(tmpStr,options));
        }
        else if (tmpStr.substr(0,1)=="-")
        {

            if (i+1 ==argc)
                throw MoleExcept(9010102,
                                 "IChem::CONSTRUCTOR",
                                 "No value given for "+tmpStr);
            tmpStr2 = argv[i+1];

            if (possInputs.find(" "+tmpStr2+" ") != string::npos)
                throw MoleExcept(9010103,
                                 "IChem::CONSTRUCTOR",
                                 "No value given for "+tmpStr);
            else if (tmpStr == "-debug")
            {

                vector<string> tokens;
                tokenStr(tmpStr2,tokens,"-");
                continue;
            }
            else
            {
                options.clear();
                while(tmpStr2.substr(0,1) != "-" &&
                      i+2!=argc &&
                      find(possPrgs.begin(),possPrgs.end(), tmpStr2) == possPrgs.end())
                {
                    options.push_back(tmpStr2);
                    i++;
                    tmpStr2=argv[i+1];
                }
                i--;
                Opt_Values.insert(pair<string,vector<string> >(tmpStr,options));
            }
            i++;
        }
        else
        {
            const vector<string>::const_iterator pos =
                    find(possPrgs.begin(),possPrgs.end(), tmpStr);
            if (pos != possPrgs.end())
            {
                TaskType = std::distance(possPrgs.begin(),pos);
                opts=false;
            }
            else {
                std::ostringstream oss;
                oss << "Unknown option or misplaced token before tool name: '" << tmpStr << "'";
                throw MoleExcept(9010106, "IChem::CONSTRUCTOR", oss.str());
            }
        }


    }

    if (verbose)
    {
        for (std::map<std::string,vector<std::string> >::iterator
             it = Opt_Values.begin();
             it != Opt_Values.end();
             it++)
        {
            cout << "OPTION : "<< (*it).first ;
            for (vector<string>::iterator
                 it2= (*it).second.begin();
                 it2 != (*it).second.end();
                 it2++) cout <<"\t"<< (*it2)<<endl;
        }
        cout << "TASKTYP : "<< TaskType<<endl;
        for (std::vector<string>::iterator
             it = Input_Values.begin();
             it != Input_Values.end();
             it++)
        {
            cout << "INPUT : "<< (*it)<<endl;
        }
    }

    if (TaskType == -1)
    {
        if (!Input_Values.empty()|| !Opt_Values.empty()){ throw MoleExcept(9010104,
                                                                           "IChem::CONSTRUCTOR",
                                                                           "No running type given (tool missing)");
        }
        else
        {help();return;}
    }
    try{
        (this->*listFunc[TaskType])();

    }
    catch (MoleExcept &e){throw;}
    catch (exception &e){throw;}



}
IChemSwitch::IChemSwitch(const int &argc,const std::vector<std::string>& argv) 
{
    if (argc == 1) {help();return;}

    // Global aromaticity model: reset on every invocation (-F batch mode reuses
    // the process), --oldAro switches back to the legacy bond-count rules.
    Molecule::setAromaticityMode(AromaticityMode::SP2_PLANAR);

    int TaskType=-1;
    std::string tmpStr,tmpStr2;bool opts=true;
    vector<string> options;
    for (int i=0; i < argc; i++)
    {

        tmpStr = argv[i];
        if (!opts) {Input_Values.push_back(tmpStr); continue;}
        // Boolean option :
        if (tmpStr.substr(0,2)=="--")
        {
            if (tmpStr =="--verb") {verbose=true;continue;}
            if (tmpStr =="--oldAro")
            {
                Molecule::setAromaticityMode(AromaticityMode::BOND_COUNT);
                continue;
            }
            options.clear();
            options.push_back("TRUE");
            Opt_Values.insert(pair<string,vector<string> >(tmpStr,options));
        }
        else if (tmpStr.substr(0,1)=="-")
        {

            if (i+1 == argc)
                throw MoleExcept(9010102,
                                 "IChem::CONSTRUCTOR",
                                 "No value given for "+tmpStr);
            tmpStr2 = argv[i+1];

            if (possInputs.find(" "+tmpStr2+" ") != string::npos)
                throw MoleExcept(9010103,
                                 "IChem::CONSTRUCTOR",
                                 "No value given for "+tmpStr);
            else if (tmpStr == "-debug")
            {

                vector<string> tokens;
                tokenStr(tmpStr2,tokens,"-");
                continue;
            }
            else
            {
                options.clear();
                while(tmpStr2.substr(0,1) != "-" &&
                      i+1!=argc &&
                      find(possPrgs.begin(),possPrgs.end(), tmpStr2) == possPrgs.end())
                {
                    options.push_back(tmpStr2);
                    i++;
                    tmpStr2=argv[i+1];
                }
                i--;
                Opt_Values.insert(pair<string,vector<string> >(tmpStr,options));
            }
            i++;
        }
        else
        {
            const vector<string>::const_iterator pos =
                    find(possPrgs.begin(),possPrgs.end(), tmpStr);
            if (pos != possPrgs.end())
            {
                TaskType = std::distance(possPrgs.begin(),pos); 
                opts=false;
            }
            else {
                std::ostringstream oss;
                oss << "Unknown option or misplaced token before tool name: '" << tmpStr << "'";
                throw MoleExcept(9010106, "IChem::CONSTRUCTOR", oss.str());
            }
        }


    }

    if (verbose)
    {
        for (std::map<std::string,vector<std::string> >::iterator
             it = Opt_Values.begin();
             it != Opt_Values.end();
             it++)
        {
            cout << "OPTION : "<< (*it).first ;
            for (vector<string>::iterator
                 it2= (*it).second.begin();
                 it2 != (*it).second.end();
                 it2++) cout <<"\t"<< (*it2)<<endl;
        }
        cout << "TASKTYP : "<< TaskType<<endl;
        for (std::vector<string>::iterator
             it = Input_Values.begin();
             it != Input_Values.end();
             it++)
        {
            cout << "INPUT : "<< (*it)<<endl;
        }
    }

    if (TaskType == -1)
    {
        if (!Input_Values.empty()|| !Opt_Values.empty()){ throw MoleExcept(9010104,
                                                                           "IChem::CONSTRUCTOR",
                                                                           "No running type given");
        }
        else
        {help();return;}
    }
    try{
        (this->*listFunc[TaskType])();

    }
    catch (MoleExcept &e){throw;}
    catch (exception &e){throw;}



}


void IChemSwitch::help()
{
    cout << "IChem Version : " << ICHEM_VERSION <<endl;
    cout << "Date : " << ICHEM_RELEASE << endl<<endl;
    cout << "Global options (before the tool name):" << endl
         << "  --verb     Verbose output" << endl
         << "  --oldAro   Legacy ring aromaticity perception (ring double/aromatic" << endl
         << "             bond count). Default: all ring atoms sp2 AND ring planar" << endl
         << endl;
    for (size_t i=0; i< ICHEM_TOOLS;i++)
    {

        (this->*listHelp[i])();

    }
}
