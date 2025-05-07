#include "headers/ICTools/switch.h"

using namespace std;
using namespace ICMole;



//IChemSwitch::pf1 te=&IChemSwitch::test;

bool IChemSwitch::verbose       =false;
bool IChemSwitch::licence_loaded=false;
// bool IChemSwitch::allowed[ICHEM_TOOLS]= {false,false,false,false,false,false,false,false,false,false,false,false,false,false};
bool IChemSwitch::allowed[ICHEM_TOOLS]= {true,true,true,true,true,true,true,true,true,true,true,true,true};
const std::string IChemSwitch::moduluskeys   ="1722935102953";
const std::string IChemSwitch::privatexp     ="817040704029";
const std::string IChemSwitch::publicexp     ="65537";
const std::string IChemSwitch::possInputs    =
        " realign genKey BSAcalc IFP grim ints fgps volsite pdbconv patch utils licence sims scoring \
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
    RX.push_back("genKey");
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


IChemSwitch::pf1 IChemSwitch::listFunc[ICHEM_TOOLS]={&IChemSwitch::genKey,
                                                     &IChemSwitch::realign,
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
IChemSwitch::pf1 IChemSwitch::listHelp[ICHEM_TOOLS]={&IChemSwitch::helpGenKey,
                                                     &IChemSwitch::helpRealign,
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


IChemSwitch::IChemSwitch(const int& argc, char *argv[]) throw(ICMole::MoleExcept)
{
    // checkLicence();

    if (argc == 1) {help();return;}



    int TaskType=-1;
    std::string tmpStr,tmpStr2;bool opts=true;
    vector<string> options;
    for (int i=1; i < argc; i++)
    {

        tmpStr = argv[i];
        if (tmpStr == "license") { yourlicence(); return;}
        if (!opts) {Input_Values.push_back(tmpStr); continue;}
        // Boolean option :
        if (tmpStr.substr(0,2)=="--")
        {
            if (tmpStr =="--verb") {verbose=true;continue;}
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
                TaskType = std::distance(possPrgs.begin(),pos); opts=false;
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
        if (!allowed[TaskType])
            throw MoleExcept(9010105,
                             "IChemSwitch::CONSTRUCTOR",
                             "Your licence is not valid for this tool");
        (this->*listFunc[TaskType])();

    }
    catch (MoleExcept &e){throw;}
    catch (exception &e){throw;}



}
IChemSwitch::IChemSwitch(const int &argc,const std::vector<std::string>& argv) throw(ICMole::MoleExcept)
{
    // checkLicence();

    if (argc == 1) {help();return;}



    int TaskType=-1;
    std::string tmpStr,tmpStr2;bool opts=true;
    vector<string> options;
    for (int i=0; i < argc; i++)
    {

        tmpStr = argv[i];
        if (tmpStr == "license") { yourlicence(); return;}
        if (!opts) {Input_Values.push_back(tmpStr); continue;}
        // Boolean option :
        if (tmpStr.substr(0,2)=="--")
        {
            if (tmpStr =="--verb") {verbose=true;continue;}
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
                TaskType = std::distance(possPrgs.begin(),pos); opts=false;
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
        if (!allowed[TaskType])
            throw MoleExcept(9010105,
                             "IChemSwitch::CONSTRUCTOR",
                             "Your licence is not valid for this tool");
        (this->*listFunc[TaskType])();

    }
    catch (MoleExcept &e){throw;}
    catch (exception &e){throw;}



}


/**
 * @brief IChemSwitch::checkLicence
 * Load ICHEM_LIC licence file to check if user can use tools
 *
 * @throw 9010301 : No ICHEM_LIC found
 * @throw 9010302 : Cannot open licence file
 * @throw 9010303 : Your licence is not valid
 * @throw 9010304 : Your licence has expired
 * @throw 9010305 : Your licence has expired
 * @throw 9010306 : Your licence has expired
 */
void IChemSwitch::checkLicence()    throw(ICMole::MoleExcept)
{

    if (licence_loaded) return;

    // Fetching ICHEM_LIC Value from bash :
    char * pPath;
    std::string ligne;
    std::vector<std::string> Licence_Tokens;

    const string var_env = "ICHEM_LIC";
    pPath = getenv (var_env.c_str());
    if (pPath==NULL)
        throw MoleExcept(9010301,
                         "IChem::checkLicence",
                         "No ICHEM_LIC found");
    ifstream ifs(pPath,ios::in);        if (!ifs.is_open())
        throw MoleExcept(9010302,
                         "IChem::checkLicence",
                         "Cannot open licence file ");

    //Reading line of comments:
    getline(ifs,ligne);
    getline(ifs,ligne);
    // Reading the good line :
    getline(ifs,ligne);
    ifs.close();


    // Decrypting the licence :
    const std::string RSAdecrypted = ec::RSA::Decrypt(ligne, moduluskeys, privatexp);
    // Example of outputs : 1|1|1|1|1|1|1|1|1|1|1|2020|12|31|

    // Splitting line into a string array regarding |
    tokenStr(RSAdecrypted,Licence_Tokens,"|");

    // Checking size :
    if (Licence_Tokens.size() < 3+ICHEM_TOOLS)
        throw MoleExcept(9010303,
                         "IChemSwitch::genKey",
                         "Your licence is not valid");


    // Updating IChem allowed tools :
    for (size_t i=0; i< ICHEM_TOOLS;i++)
    {

        allowed[i]=(Licence_Tokens.at(i) == "1")?true:false;

    }

    // Checking date :
    const int year =atoi(Licence_Tokens.at(ICHEM_TOOLS).c_str());
    const int month=atoi(Licence_Tokens.at(ICHEM_TOOLS+1).c_str());
    const int day  =atoi(Licence_Tokens.at(ICHEM_TOOLS+2).c_str());
    licence_loaded=true;

    // Getting current date :
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );


    // Comparing date :
    if (year < now->tm_year+1900)
        throw MoleExcept(9010304,
                         "IChem::checkLicence",
                         "Your licence has expired on "+
                         Licence_Tokens.at(ICHEM_TOOLS)  +" "+
                         Licence_Tokens.at(ICHEM_TOOLS+1)+" "+
                         Licence_Tokens.at(ICHEM_TOOLS+2));
    else if (year > now->tm_year+1900) return;


    if (month < now->tm_mon+1)
        throw MoleExcept(9010305,
                         "IChem::checkLicence",
                         "Your licence has expired on "+
                         Licence_Tokens.at(ICHEM_TOOLS)  +" "+
                         Licence_Tokens.at(ICHEM_TOOLS+1)+" "+
                         Licence_Tokens.at(ICHEM_TOOLS+2));
    else if (month > now->tm_mon+1)return;


    if (day < now->tm_mday)
        throw MoleExcept(9010306,
                         "IChem::checkLicence",
                         "Your licence has expired on "+
                         Licence_Tokens.at(ICHEM_TOOLS)  +" "+
                         Licence_Tokens.at(ICHEM_TOOLS+1)+" "+
                         Licence_Tokens.at(ICHEM_TOOLS+2));
    cerr << "Warning: Your licence expires in "
         <<Licence_Tokens.at(ICHEM_TOOLS)<<" "
        <<Licence_Tokens.at(ICHEM_TOOLS+1)<<" "
       <<Licence_Tokens.at(ICHEM_TOOLS+2)<<endl;

}


void IChemSwitch::help()
{
    cout << "IChem Version : " << ICHEM_VERSION <<endl;
    cout << "Date : " << ICHEM_RELEASE << endl<<endl;
    cout << "./IChem license\nTo show your license\n############################"<<endl << endl;
    for (size_t i=0; i< ICHEM_TOOLS;i++)
    {

        if (!allowed[i]) continue;
        (this->*listHelp[i])();

    }
}
