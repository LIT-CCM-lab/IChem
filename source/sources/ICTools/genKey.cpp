#include "headers/ICTools/switch.h"

using namespace std;
using namespace ICMole;


/**
 * @brief IChemSwitch::helpGenKey
 * Display help for generate licence file:
 */
void IChemSwitch::helpGenKey() const {
    cout
        << "genKey - License key creation" << endl
        << "Usage: IChem genKey year month day allowed_tools" << endl
        << endl
        << "Description:" << endl
        << "  Generate a license key defining the expiry date and enabled tools" << endl
        << endl
        << "Arguments:" << endl
        << "  year          Expiring year  (4 digits, e.g. 2014)" << endl
        << "  month         Expiring month (2 digits, e.g. 09)" << endl
        << "  day           Expiring day   (2 digits, e.g. 23)" << endl
        << "  allowed_tools List of 0/1 flags defining allowed tools, in this order:" << endl
        << "                1. Licence Key generator" << endl
        << "                2. Realigning molecules" << endl
        << "                3. BSA Calculation" << endl
        << "                4. IFP generator" << endl
        << "                5. Interaction detection" << endl
        << "                6. Graph Interaction Matching" << endl
        << "                7. VolSite" << endl
        << "                8. PDB to MOL2 conversion" << endl
        << "                9. Patching MOL2 conversion" << endl
        << "                10. Utils" << endl
        << "                11. sims" << endl
        << "                12. scoring" << endl
        << "                13. Detection of PPi" << endl
        << endl
        << "Example:" << endl
        << "  IChem genKey 2025 10 20 1 1 1 1 1 1 1 1 1 1 1 1 1" << endl
        << endl
        << "---------------------------------------------------------------------------" << endl
        << endl;
}

void IChemSwitch::genKey()const  
{

    if (Input_Values.size() != 3+ICHEM_TOOLS)
    {
        cout << "GOT : "<<Input_Values.size()<<endl;
        cout << "EXPECTED : "<< (3+ICHEM_TOOLS)<<endl;
        throw MoleExcept(9020101,
                         "IChemSwitch::GenKey",
                         "Not enough entry");
    }


    ostringstream licence_decoded;
    int val;
    for (size_t i=3; i<3+ICHEM_TOOLS;i++ )
    {
        val = atoi(Input_Values.at(i).c_str());
        if (val != 0 && val != 1)
            throw MoleExcept(9020102,
                             "IChemSwitch::GenKey",
                             "Invalid Value. Must be 0 or 1");
        switch (i)
        {
        case 3: cout<< "      Licence key generator : ";break;
        case 4: cout<< "      Molecular realignment : ";break;
        case 5: cout<< "            BSA Calculation : ";break;
        case 6: cout<< "              IFP generator : ";break;
        case 7: cout<< "      Interaction detection : ";break;
        case 8: cout<< " Graph Interaction Matching : ";break;
        case 9: cout<< "                    VolSite : ";break;
        case 10:cout<< "     PDB to MOL2 conversion : ";break;
        case 11:cout<< "   Patching MOL2 conversion : ";break;
        case 12:cout<< "                      Utils : ";break;
        case 13:cout<< "                       Sims : ";break;
        case 14:cout<< "                    Scoring : ";break;
        case 15:cout<< "  Detection, analyse of PPi : ";break;


        }
        cout << val<<endl;
        licence_decoded<< val<<"|";
    }
    const int year = atoi(Input_Values.at(0).c_str());
    cout << "            Expiration YEAR : "<< year<<endl;
    if (year < 2000 || year > 3000)
        throw MoleExcept(9020103,
                         "IChemSwitch::GenKey",
                         "Invalid Value for the year");
    licence_decoded<< year<<"|";

    const int month = atoi(Input_Values.at(1).c_str());
    cout << "           Expiration MONTH : "<< month<<endl;
    if (month ==0 || month > 12)
        throw MoleExcept(9020104,
                         "IChemSwitch::GenKey",
                         "Invalid Value for the month. Must be between 1 and 12");
    licence_decoded<< month<<"|";

    const int day = atoi(Input_Values.at(2).c_str());
    cout << "             Expiration DAY : "<< day<<endl;
    if (day < 1 && day >31)
        throw MoleExcept(9020105,
                         "IChemSwitch::GenKey",
                         "Invalid Value. Must be between 1 or 31");
    licence_decoded<< day<<"|";
    cout << licence_decoded.str()<<endl;

    // RSA encryption :
    std::string RSAencrypted = ec::RSA::Encrypt(licence_decoded.str(),
                                                moduluskeys,
                                                publicexp);


    cout << endl<<endl<<endl;
    cout << "Copy the following into the licence file :"<<endl;
    cout << "# Licence Key IChem : "<<endl;
    cout << "# Expiration date (D/M/Y) "<< day<< "/"<<month<<"/"<<year<<endl;
    cout << RSAencrypted<<endl;
}






/**
 * \brief Show to the user the content of it's licence
 *
 * Load ICHEM_LIC path and it's license content. Then decrypt it and analyse the results.
 */
void IChemSwitch::yourlicence()
{

    char * pPath;
    const string var_env = "ICHEM_LIC";
    std::vector<std::string> Licence_Tokens;


    // Fetching ICHEM_LIC Value from bash :
    pPath = getenv (var_env.c_str());   if (pPath==NULL) throw MoleExcept(9010301, "IChem::checkLicence", "No ICHEM_LIC found");
    ifstream ifs(pPath,ios::in);        if (!ifs.is_open())throw MoleExcept(9010302, "IChem::checkLicence", "Cannot open licence file ");
    std::string ligne;
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
    if (Licence_Tokens.size() < 3+ICHEM_TOOLS) throw MoleExcept(9010303,"IChemSwitch::genKey","Your licence is not valid");

    // the -1 is for licence key generator that should not be shown to users
    for (size_t i=0; i< ICHEM_TOOLS;i++)
    {
        if (Licence_Tokens.at(i)=="1"){
            switch (i)
            {
            case 0: cout<< "      Licence key generator : ";break;
            case 1: cout<< "      Molecular realignment : ";break;
            case 2: cout<< "            BSA Calculation : ";break;
            case 3: cout<< "              IFP generator : ";break;
            case 4: cout<< "      Interaction detection : ";break;
            case 5: cout<< " Graph Interaction Matching : ";break;
            case 6: cout<< "                    VolSite : ";break;
            case 7: cout<< "     PDB to MOL2 conversion : ";break;
            case 8: cout<< "   Patching MOL2 conversion : ";break;
            case 9: cout<< "                      Utils : ";break;
            case 10:cout<< "                       Sims : ";break;
            case 11:cout<< "                    Scoring : ";break;
            case 12:cout<< " Detection, analysis of PPi : ";break;

            }
            cout<< ((Licence_Tokens.at(i)=="1")?"YES":"NO")<<endl;
        }
    }
    const int year =atoi(Licence_Tokens.at(ICHEM_TOOLS).c_str());
    const int month=atoi(Licence_Tokens.at(ICHEM_TOOLS+1).c_str());
    const int day  =atoi(Licence_Tokens.at(ICHEM_TOOLS+2).c_str());
    cout << "Your licence expires on (Y/M/D): " << year << "/" << month<< "/"<<day<<endl;
}
