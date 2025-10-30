#include <iostream>
#include "headers/ICPars/molereader.h"
#include "headers/ICMole/molecule.h"
#include "headers/ICMole/atom.h"
#include "headers/ICMole/coords.h"

using namespace std;
using namespace ICMole;

//#define ICHEM_DEBUG


/**
 * @brief Explode a line from a PDB File into an associative array
 * @param str : string to explode
 * @param tokens : associative array that will be filled with the string data
 *
 * Theses are the list of tokens outputted:
 * <li> Type : line header </li>
 * <li> Num : Id of the atom </li>
 * <li> Name : Name of the atom </li>
 * <li> Res : Residu name</li>
 * <li> X : X coordinates</li>
 * <li> Y : Y coordinates</li>
 * <li> Z : Z coordinates</li>
 * <li> Occup : Occupancy factor </li>
 * <li> Heat : Heat factor </li>
 * <li> AtmName : Atomic name</li>
 *
 */
void MoleReader::tokenize( const string& str, map<string,string>& tokens)
{
     const string delimiters = ", \n, ;, (, )";
    const struct
    {
        unsigned int beg;
        unsigned int len;
        string Name;
    } Pos[10]={{0,6,"Type"},
               {6,5,"Num"},
               {12,4,"Name"},
               {17,3,"Res"},
               {30,8,"X"},
               {38,8,"Y"},
               {46,8,"Z"},
               {54,6,"Occup"},
               {60,6,"Heat"},
               {76,2,"AtmName"}};

    //function to set colons with a valid content, there is a probleme when colons overlap
   /* string col;
    string::size_type pos;
    string::size_type lastPos;*/

    const size_t length = str.length();

    for (size_t ppos=0; ppos < 10; ppos++)
    {
        if (Pos[ppos].beg+Pos[ppos].len > length) break;
        const std::string col = str.substr(Pos[ppos].beg, Pos[ppos].len);
        const string::size_type lastPos = col.find_first_not_of(delimiters, 0);
        const string::size_type pos = col.find_first_of(delimiters, lastPos);
        if (pos-lastPos <= 0) continue;
        tokens.insert(pair<string,string>(Pos[ppos].Name,col.substr(lastPos, pos-lastPos)));
    }
    std::string col = str.substr(21,1);
    if (col.compare(" ") != 0)
    {
        tokens.insert(pair<string,string>("Chain",col));
        col = str.substr(22,5);
        tokens.insert(pair<string,string>("ResId",col));//col.substr(lastPos, pos-lastPos)
    }
    else
    {
        col = str.substr(21,6);
        const string::size_type  lastPos = col.find_first_not_of(delimiters, 0);
        const string::size_type  pos = col.find_first_of(delimiters, lastPos);
        tokens.insert(pair<string,string>("ResId",col.substr(lastPos, pos-lastPos)));
    }
#ifdef ICHEM_DEBUG
        for(map<string,string>::iterator
            it = tokens.begin();
            it != tokens.end() ;
            it++)
            cout << (*it).first<<":"<<(*it).second<<" || "; cout <<endl;
#endif

}


void MoleReader::loadPDB(Molecule& molecule,const unsigned int& TypeMol) 
{
    if (!finput.is_open() && !lineinfile)
        throw MoleExcept(2020801,
                         "MoleReader::loadPDB",
                         "No file opened");

    const string delimiters(", \n, ;, (, )");
    const string whitespaces (" \t\f\v\n\r");


    string            col2;
    string            buffer;
    string::size_type pos;
    string::size_type lastPos;

    int Res_Id=0;
    unsigned int atm1,atm2,atm3,atm4;
    Residu* res;
    ostringstream oss;
    map<string,string>col;

    std::map<std::string,Residu*>::iterator ITTI;
    std::map<unsigned int,Atom*> AtmL;
    mapRes.clear();
    molecule.setMoleType(TypeMol);
    try
    {

       while(!EOFile)
        {
            getLine(buffer);
            if (EOFile||buffer == "END") break;

            if (buffer.substr(0,6) == "HEADER")
            {
                lastPos =buffer.find_last_not_of(whitespaces);
                buffer.erase(lastPos+1,buffer.length()-lastPos-1);
                buffer.erase(0,7);
                lastPos =buffer.find_first_not_of(whitespaces);
                buffer.erase(0,lastPos);
                molecule.setName(buffer);
            }
            if (buffer.substr(0,4) == "ATOM"
             || buffer.substr(0,6).find("HETATM") != string::npos)
            {
                try{
                col.clear();
                tokenize(buffer, col);
                if (col.size()<8) continue;

                Res_Id = atoi(col["ResId"].c_str());

                res=(Residu*)NULL;
                if ( Res_Id >=0 && col["Res"].length() > 0)
                {
                    oss.str("");
                    oss << Res_Id<< col["Res"]<<col["Chain"];

                    ITTI = mapRes.find(oss.str());
                    if (ITTI == mapRes.end()){

                        res= &molecule.addResidu(col["Chain"],Res_Id, col["Res"]);
                        mapRes.insert(pair<string,Residu*>(oss.str(),res));
                    }
                    else res = (*ITTI).second;
                }
                if (res== (Residu*)NULL) res = &molecule.getUnknownRes();


                if (col["AtmName"] == "FE") col["AtmName"]="Fe";
                else if (col["AtmName"] == "MG") col["AtmName"]="Mg";
                else if (col["AtmName"] == "MN") col["AtmName"]="Mn";
                else if (col["AtmName"] == "ZN") col["AtmName"]="Zn";
                else if (col["AtmName"] == "MO") col["AtmName"]="Mo";
                else if (col["AtmName"] == "CA") col["AtmName"]="Ca";
                else if (col["AtmName"] == "CU") col["AtmName"]="Cu";
                else if (col["AtmName"] == "AL") col["AtmName"]="Al";
                else if (col["AtmName"] == "BR") col["AtmName"]="Br";
                else if (col["AtmName"] == "CL") col["AtmName"]="Cl";
                else if (col["AtmName"] == "NA") col["AtmName"]="Na";// PATCH ERR::3

                Coords coords = Coords(atof(col["X"].c_str()), atof(col["Y"].c_str()), atof(col["Z"].c_str()));

                Atom &atm = molecule.addAtom(
                        (((col["AtmName"]).length() > 0) ? col["AtmName"] : col["Name"]  ),
                        coords,
                        (col["Name"].length()>0) ? col["Name"]:col["AtmName"],
                        "", res);

//                Atom &atm =
//                        molecule.addAtom((((col["AtmName"]).length() > 0) ? col["AtmName"] : col["Name"]  ),
//                               atof(col["X"].c_str()),
//                               atof(col["Y"].c_str()),
//                               atof(col["Z"].c_str()),
//                               (col["Name"].length()>0) ? col["Name"]:col["AtmName"],
//                                "",
//                               res);
                atm.setFNum(atoi(col["Num"].c_str()));

                atm.setBFactor(atof(col["Heat"].c_str()));

                AtmL.insert(pair<unsigned int,Atom*>(atm.getFNum(),&atm));
                }catch (MoleExcept &e){
                    cerr << e.getCode()<<"\t"<< e.getData()<<endl;
                }
            }
            else if (buffer.substr(0,6) == "CONECT")
            {
try{
                buffer = buffer.substr(0,buffer.find_last_not_of(whitespaces)+2);

                const string::size_type sized = buffer.length();
#ifdef ICHEM_DEBUG
                cout << "CONNECT : "<< buffer.length()<<endl;
#endif
                atm2=0; atm3=0; atm4=0;
                col2 = buffer.substr(6,5);
#ifdef ICHEM_DEBUG
                cout << " ; Atom1 = " << col2 ;
#endif
                lastPos = col2.find_first_not_of(delimiters, 0);
                pos = col2.find_first_of(delimiters, lastPos);
                atm1=atoi(col2.substr(lastPos, pos-lastPos).c_str());
                col2 = buffer.substr(11,5);
#ifdef ICHEM_DEBUG
                cout << " ; Atom2 = " << col2 ;
#endif
                lastPos = col2.find_first_not_of(delimiters, 0);
                pos = col2.find_first_of(delimiters, lastPos);
                atm2=atoi(col2.substr(lastPos, pos-lastPos).c_str());
                if (sized >= 21)
                {
                    col2 = buffer.substr(16,5);
#ifdef ICHEM_DEBUG
                    cout << " ; Atom3 = " << col2 ;
#endif
                    lastPos = col2.find_first_not_of(delimiters, 0);
                    pos = col2.find_first_of(delimiters, lastPos);
                    atm3=atoi(col2.substr(lastPos, pos-lastPos).c_str());
                }
                if (sized >= 26)
                {
                    col2 = buffer.substr(21,5);
#ifdef ICHEM_DEBUG
                    cout << " ; Atom4 = " << col2 ;
#endif
                    lastPos = col2.find_first_not_of(delimiters, 0);
                    pos = col2.find_first_of(delimiters, lastPos);
                    atm4=atoi(col2.substr(lastPos, pos-lastPos).c_str());
                }

                if (AtmL.find(atm1) == AtmL.end()) continue;

                if (atm2 != 0 && atm2< atm1) {
                    if (AtmL.find(atm2) == AtmL.end()) continue;
#ifdef ICHEM_DEBUG
                    if (debug)cout << "Create bond : "
                                   << AtmL.at(atm1)->getIdentifier()
                                   <<" - "
                                  << AtmL.at(atm2)->getIdentifier() <<endl;
#endif
                    molecule.addBond(AtmL.at(atm1),
                                     AtmL.at(atm2),
                                     BondType::SINGLE);}

                if (atm3 != 0 && atm3< atm1) {
                    if (AtmL.find(atm3) == AtmL.end()) continue;
#ifdef ICHEM_DEBUG
                    if (debug)cout << "Create bond : "
                                   << AtmL.at(atm1)->getIdentifier()
                                   <<" - "
                                  << AtmL.at(atm2)->getIdentifier() <<endl;
#endif
                    molecule.addBond(AtmL.at(atm1),
                                     AtmL.at(atm3),
                                     BondType::SINGLE);}
                if (atm4 != 0 && atm4< atm1) {
                    if (AtmL.find(atm4) == AtmL.end()) continue;
#ifdef ICHEM_DEBUG
                    cout << "Create bond : "
                                   << AtmL.at(atm1)->getIdentifier()
                                   <<" - "
                                  << AtmL.at(atm2)->getIdentifier() <<endl;
#endif
                    molecule.addBond(AtmL.at(atm1)
                                     ,AtmL.at(atm4),
                                     BondType::SINGLE);
                }
#ifdef ICHEM_DEBUG
                cout <<endl;
#endif
                }catch (MoleExcept &e){}
            }
        }

    }
    catch ( MoleExcept &e)
    {
        cerr << e.getCode()<<endl;
    }


}

