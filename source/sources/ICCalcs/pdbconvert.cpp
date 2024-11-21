#include <iostream>
#include <math.h>
#include "headers/ICCalcs/pdbconvert.h"
#include "headers/ICPars/molereader.h"
using namespace ICMole;
using namespace std;

PDBConvert::PDBConvert(const string &name)
{
    PDBName=name;
    master=-1;

}


bool PDBConvert::for_scPDB=false;


void findAndReplace( std::string& ligne, const string& from, const string& to)
{
  const size_t from_size= from.size();
  size_t input_size = ligne.length();
  for (size_t n_pos = 0; n_pos<input_size-from_size; n_pos++)
    {
      if (ligne.substr(n_pos,from_size).compare(from) == 0) ligne.replace(n_pos,from_size,to);
    }

}


/**
 * @brief Scan a PDB File and save selected lines
 * @param Path : Path of the Protein.pdb file
 *
 * Open the PDB File and saves line positions of ATOM, HETATM,TER,REMARK,
 * CONECT,SEQRES,MODEL,ENDMDL, HEADER,MASTER,END, ORIGX ,SCALE,MTRIX block
 *
 * Select only the first model in NMR pdb files
 *
 */
void PDBConvert::loadPDB(const std::string& Path)
throw(ICMole::MoleExcept)
{
    unsigned short nmodel=0;
    filelines.clear();
    int fileline=0;
    std::string ligne,header;
    std::vector<unsigned int> listatm;   // Creating another list for atoms to handle models

// STEP 1 - LOADING FILE :
    if (Path.length()==0)
        throw MoleExcept(3030101,
                         "PDBConvert::Constructor",
                         "No file given");

    filestream.open(Path.c_str());
    if (!filestream.is_open())
        throw MoleExcept(3030102,
                         "PDBConvert::Constructor",
                         "Cannot open "+Path);



// STEP 2 - Reading file and handling model
    while(!filestream.eof())
    {
   // Getting line :
        std::getline(filestream,ligne);
    // Fetching header :
        header=ligne.substr(0,6);
        if ((    header.compare("ATOM  ") == 0
               ||header.compare("HETATM") == 0
               ||header.compare("TER   ")  ==0) && nmodel<=1) listatm.push_back(fileline);
        else if (header.compare("REMARK") == 0)  remark.push_back(fileline);
        else if (header.compare("CONECT") == 0) connect.push_back(fileline);
        else if (header.compare("SEQRES") == 0)  seqres.push_back(fileline);
        else if (header.compare("MODEL ") == 0) { nmodel++; }
        else if (header.compare("ENDMDL") == 0) { if (nmodel==1) atoms=listatm;}
        // 1 line per entry :
        else if (header.compare("HEADER") == 0) header_line = fileline;
        else if (header.compare("MASTER") == 0) master = fileline;
        else if (header.compare("END   ") == 0) break;
        else{
            header = ligne.substr(0,5);
            if (header.compare("ORIGX") == 0
                ||header.compare("SCALE") == 0
                ||header.compare("MTRIX") == 0) { coords.push_back(fileline);}
          }
        filelines.push_back(ligne);
    fileline++;
    }// END WHILE

// STEP 3 - Getting PDB ID :
    header =filelines.at(header_line);
    PDBName = header.substr(header.find_last_not_of(" ")-3,4);

// STEP 4 - In case of not NMR file => get the whole atom list
    if (nmodel == 0 && listatm.size() != 0) atoms= listatm;

filestream.close();


    if (!for_scPDB) return;
    cout << PDBName<<"|READER|HEADER\t##### READING FILE #####\n";
    cout << PDBName<<"|READER|STATIS\t"<< remark.size() << "\t REMARK lines\n"<<
            PDBName<<"|READER|STATIS\t"<<  atoms.size() << "\t ATOM|HETATM|TER lines\n"<<
            PDBName<<"|READER|STATIS\t"<<connect.size() << "\t CONECT lines\n"<<
            PDBName<<"|READER|STATIS\t"<< seqres.size() << "\t SEQRES lines\n"<<
            PDBName<<"|READER|STATIS\t"<< coords.size() << "\t ORIGX|SCALE|MTRIX lines"<<endl;



}


/**
 * @brief Harmonize length of each line
 * @param length : string length that all lines must fit to
 *
 *
 *  Check whether each line as a length at least equal to length parameter (default 80)
 */
void PDBConvert::harmonizeSizeLine(const size_t& length)
{
    const size_t nlines = filelines.size();
   size_t strlen=0;
  unsigned short nmodif=0;

  if (for_scPDB)
      cout << PDBName <<"|HARMON|START\t\n"
           << PDBName <<"|HARMON|HEADER\t###### Harmonizing size line\n" ;
  for (size_t nline = 0;
              nline < nlines;
            ++nline)  // Scanning every line
    {

      if (filelines.at(nline).length() >= 80) continue; // Checking length of string
      string &ligne=filelines.at(nline);
      strlen = ligne.length();

      nmodif++;

      if (for_scPDB) cout << PDBName <<"|HARMON|MODIFS\t"<< nline<<" INI\t" <<strlen;
      for (size_t ni = strlen; ni < length; ni++)   // Adding spaces to get the good length
        {
          ligne+= " ";
        }
      if (for_scPDB)cout << PDBName <<"|HARMON|MODIFS\tEND\t" <<ligne.length()<<"\n";
    }
  if (for_scPDB) cout << PDBName <<"|HARMON|STATIS\t"<< nmodif <<"\tline modified\n"<< PDBName<<"|HARMON|END\t"<<endl;
}


/**
  * @brief Changing Selenomethionine to methonine
  *
  * Check SEQRES, ATOM|HETATM|TER block and change MSE to MET.
  * Change in SEQRES block " MSE" => " MET"
  * Change in HETATM block :
  *   -> Header from "HETATM" to "ATOM  "
  *   -> Residu name from "MSE" to "MET"
  *   -> Atom name for "SE" from "SE" to "SD" and atom type to "S"
*/
void PDBConvert::changeMSEtoMET()
{
  const string MSE_seqres = "MSE ";
  const string MET_seqres = "MET ";
  unsigned int nmodif=0;

  const size_t nAtomsSize=atoms.size();

  if (for_scPDB) cout << PDBName <<"|MSEMET|START\t\n"<<PDBName<<"|MSEMET|HEADER\t##### Changing MSE to MET "<<endl ;

  // STEP 1 - Scanning ATOM LINES :
  for (size_t pos=0; pos < nAtomsSize;pos++)
  {
  // Setting LIGNE:
    std::string &ligne=filelines.at(atoms.at(pos));
  // NOT interesting line => continue
    if (ligne.substr(17,3).compare("MSE") != 0) continue;
    nmodif++;
  // Modify  line's header :
    if (ligne.substr(0,6).compare("HETATM")==0) ligne.replace(0,6,"ATOM  ");
  // Modifying "MSE" to "MET"
    ligne.replace(17,3,"MET");
  // Modifying atom name and type
     if (ligne.substr(12,4).compare("SE  ")!=0) continue;
     ligne.replace(12,4," SD ");
     ligne.replace(76,2," S");

     cout <<PDBName <<"|MSEMET|MODIFS\tNEW\t"<<atoms.at(pos)<<"\t"<< ligne<<"\n";

  }
  if (nmodif == 0)
  {
      if (for_scPDB){cout << PDBName<<"|MSEMET|STATIS\t0\tNo modified line\n";
                   cout << PDBName<<"|MSEMET|END\t"<<endl;}
      return ;
  }


  const size_t nSeqRes=seqres.size();
  for (size_t pos=0; pos < nSeqRes;pos++)
  {
  // Setting LIGNE:
    std::string &ligne=filelines.at(atoms.at(pos));

    if (ligne.find(MSE_seqres)==string::npos)continue;
    if (for_scPDB) { cout <<PDBName <<"|MSEMET|MODIFS\tOLD\t"<<atoms.at(pos) << "\t"<< ligne<<"\n";}
    nmodif++;
    findAndReplace(ligne,MSE_seqres,MET_seqres);
    if (for_scPDB) { cout <<PDBName <<"|MSEMET|MODIFS\tNEW\t"<<atoms.at(pos) << "\t"<< ligne<<"\n";}
  }


  if (for_scPDB) cout << PDBName<<"|MSEMET|STATIS\t"<< nmodif<<"\tLine modified\n"<<PDBName<<"|MSEMET|END\t"<<endl;

}







/**
 * @brief  Changing SelenoCystein to cystein
  *
  * Check SEQRES, ATOM|HETATM|TER block and change CSE to CYS.
  * Change in SEQRES block " CSE" => " CYS"
  * Change in HETATM block :
  *   -> Header from "HETATM" to "ATOM  "
  *   -> Residu name from "CSE" to "CYS"
  *   -> Atom name for "SG" from "SG" to "SE" and atom type to "S"
 */
void PDBConvert::changeCSEtoCYS()
{
  const string CSE_seqres = " CSE";
  const string CYS_seqres = " CYS";
  unsigned int nmodif=0;
const size_t nAtomsSize=atoms.size();
  if (for_scPDB)
      cout << PDBName <<"|CSECYS|START\t"<<endl
           << PDBName <<"|CSECYS|HEADER\t##### Changing CSE to CYS \n" ;
  // STEP 1 - Scanning ATOM LINES :
  for (size_t pos=0; pos < nAtomsSize;pos++)
  {
  // Setting LIGNE:
    std::string &ligne=filelines.at(atoms.at(pos));
  // NOT interesting line => continue
    if (ligne.substr(17,3).compare("CSE") != 0) continue;
    nmodif++;
  // Modify  line's header :
    if (ligne.substr(0,6).compare("HETATM")==0) ligne.replace(0,6,"ATOM  ");
  // Modifying "MSE" to "MET"
    ligne.replace(17,3,"CYS");
  // Modifying atom name and type
     if (ligne.substr(12,4).compare("SE  ")!=0) continue;
     ligne.replace(12,4," SG ");
     ligne.replace(76,2," S");

    if (for_scPDB) cout <<PDBName <<"|CSECYS|MODIFS\tNEW\t"
                        <<atoms.at(pos)<<"\t"<< ligne<<"\n";

  }
  if (nmodif == 0)
  {
      if (for_scPDB){cout << PDBName<<"|CSECYS|STATIS\tNo modified line\n"
                          << PDBName<<"|CSECYS|END\t"<<endl;}
      return ;
  }


  const size_t nSeqRes=seqres.size();
  for (size_t pos=0; pos < nSeqRes;pos++)
  {
  // Setting LIGNE:
    std::string &ligne=filelines.at(atoms.at(pos));

    if (ligne.find(CSE_seqres)==string::npos)continue;
    if (for_scPDB) { cout <<PDBName <<"|CSECYS|MODIFS\tOLD\t"
                          <<atoms.at(pos) << "\t"<< ligne<<endl;}
    nmodif++;
    findAndReplace(ligne,CSE_seqres,CYS_seqres);
    if (for_scPDB) { cout <<PDBName <<"|CSECYS|MODIFS\tNEW\t"
                          <<atoms.at(pos) << "\t"<< ligne<<endl;}
  }



  if (for_scPDB) { cout <<PDBName <<"|CSECYS|STATIS\t"
                       << nmodif <<"\tline modified \n"
                       <<PDBName<<"|CSECYS|END\t"<<endl;}
}











/**
 * @brief Move all HETATM to the end of ATOM block
 *   Warning : Should be run AFTER running changeCSEtoCYS() and changeMSEtoMET() functions.
  * Otherwise there will be a gap in residues.
  * This function will scan every line of ATOM|HETATM|TER blocks and split in two ATOM and HETATM.
  * hetatoms are grouped by chains.
  * For TER lines, it will depends on whether the corresponding residu is an het residu or not.
  * At last, all ATOM lines stays and HETATM are put at the end of atom list, ordered by chain
 */
void PDBConvert::moveHETATMtoend()
{
  listLine             atms;          // new line list for ATOMS
  map<string, listLine >           hets;          // new line list for HET. 1 list for each chain
  map<string, listLine >::iterator Chain_hets;    //
  string curr_res="";

  // curr_res_het tells whether the current residu is indeed a HETATM
  bool curr_res_het=false;
  const size_t nAtomsSize=atoms.size();


  if (for_scPDB) {
      cout <<PDBName <<"|HETEND|START\t\n"
          <<PDBName<<"|HETEND|HEADER\t##### Moving HETATM to end of file \n" ;}



  for (size_t pos=0; pos < nAtomsSize;pos++)
  {
  // Getting LIGNE:
    const unsigned int &Lpos = atoms.at(pos);
    std::string &ligne=filelines.at(Lpos);

 // Getting residu name for this atom :
    const std::string atomRes = ligne.substr(17,3);

    // Case either HETATM or TER (with the same current residu name)
    if (ligne.substr(0,6).compare("HETATM")==0
      ||(ligne.substr(0,3).compare("TER")
         && atomRes.compare(curr_res)==0
         && curr_res_het))
    {
        curr_res_het=true;
        if (for_scPDB) {
            cout <<PDBName <<"|HETEND|MODIFS\t"<< Lpos<<"\t"<<ligne <<"\n";
        }

        // Getting chain name:
        const std::string atomChain = ligne.substr(21,1);
        // Searching chain in the list of chains
        Chain_hets = hets.find(atomChain);

        // Not in the list => adding it
        if (Chain_hets==hets.end())
        {
            listLine Newchain;
            Newchain.push_back(Lpos);
            // Include the given HETATM line into this chain :
            hets.insert(pair<string,listLine>(atomChain,Newchain));
        }
        // Otherwise Include the given HETATM line into this chain :
        else (*Chain_hets).second.push_back(Lpos);
        curr_res=atomRes;

    }
    // And if it's just an atom we put it in atms list
    else{
        atms.push_back(Lpos);
        curr_res_het=false;}
     curr_res=atomRes;

  }

    // STEP 2 : Updating the atoms list
    atoms.clear();// cleaning list
    atoms.insert( atoms.end(), atms.begin(), atms.end() );// adding the new atms (which contains only ATOM and related TER lines)

    // adding the Hetatm (which contains only HETATOM and related TER lines)
    for (Chain_hets = hets.begin();
         Chain_hets != hets.end();
         Chain_hets++)
      {
        if (for_scPDB)
            cout <<PDBName <<"|HETEND|STATIS\t"
                << (*Chain_hets).second.size()<< "\t"<<(*Chain_hets).first << "\t atom moved"<<endl;
        atoms.insert( atoms.end(),
                      (*Chain_hets).second.begin(),
                      (*Chain_hets).second.end() );
      }
    if (for_scPDB) cout << PDBName<<"|HETEND|END\t"<<endl;


}








void PDBConvert::selAltAtm()
{
  if (for_scPDB){
      cout <<PDBName<<"|ALTPOS|START\t\n"
           <<PDBName<<"|ALTPOS|HEADER\t##### Selecting alternative atoms : "<<endl;
    }
  string                    curr_altloc  ="",
                            tmp_name     ="",
                            tmp_reschain ="",
                            tmp_altloc   ="",
                            curr_reschain="";

  // The two maps below have for key the occupancy name "A","B","C","D"
  // and for value, the occupancy sum for occup and the atom list line for atom_occ
  std::map<string,double>   occups;
  std::map<string,listLine> atom_occ;
  listLine                  new_atoms,
                            LT;
  double                    best_occup   =0,
                            sum_occup    =0;

  unsigned int              NLine        =0,
                            NRes         =0;
  const size_t              nAtomsSize   =atoms.size();

  //Scanning atom lines for alternative atoms
  for (size_t pos=0; pos < nAtomsSize;pos++)
  {
  // Getting LIGNE:
    const unsigned int &Lpos = atoms.at(pos);
    std::string &ligne=filelines.at(Lpos);
  // Getting alternative location
    tmp_altloc   = ligne.substr(16,1);
  // No alternative position => Adding to atom and continue
    if (tmp_altloc == " ") {new_atoms.push_back(Lpos); continue;}


////////// STEP1: Setting alternative data information:

  // Getting atom name:
    tmp_name     = ligne.substr(12,4);
  // Getting atom chain:
    tmp_reschain = ligne.substr(17,10);

  // Clearing occupancy maps :
    occups.clear();
    atom_occ.clear();

  // Adding occupancy value of the first alternative atom:
    occups.insert(pair<string,double>(tmp_altloc,atof(ligne.substr(54,6).c_str())));

 // Cleaning the listline vector and adding the line of the first alternative atom:
    LT.clear();
    LT.push_back(Lpos);
 // Since an insert within map copy the given data (i.e. the vector)
 // We don't need to case about former occupancy states:
    atom_occ.insert(pair<string,listLine>(tmp_altloc,LT));

 // Setting up current information:
    curr_altloc  =tmp_altloc;
    curr_reschain=tmp_reschain;


////////// STEP2: Getting all alternative atom:
    // So we read the next lines to get all the alternate positions.
    // It end when atome name and chain are different
    while(1)
      {

        pos++;
        if (pos == nAtomsSize) break;
        const unsigned int &LposA = atoms.at(pos);
        const std::string& ligneA = filelines.at(LposA);

        tmp_name     = ligneA.substr(12,4);      // Getting atom name
        tmp_reschain = ligneA.substr(17,10);     // Getting atom chain
        tmp_altloc   = ligneA.substr(16,1);      // Getting alternative location

        if (ligneA.substr(0,3)=="TER" && curr_altloc=="") curr_altloc=tmp_altloc;

        // We check that we are currently in an alternate position line
        // Otherwise we go back to the last line and stop the while
        if ((((tmp_altloc=="A" && curr_altloc != "A") ||
             (tmp_altloc=="1"&& curr_altloc != "1")) ||
                                     tmp_altloc==" " ||
                         curr_reschain != tmp_reschain)) {--pos;break;}


        // When we switch from altloc A to altloc B ...etc
        if (curr_altloc != tmp_altloc)
          {
            // We create a new list of line position
            // That will keep all line for this altenative position
            listLine LT;
            // We add the related line
            LT.push_back(LposA);
            // we associate into atom_occ the altloc value to the list of line
            atom_occ.insert(pair<string,listLine>(tmp_altloc,LT));
            //  we associate into occups the altloc value and the occupency value
            occups.insert(pair<string,double>(tmp_altloc,atof(ligneA.substr(54,6).c_str())));
            //  now we set the current altloc value to the current altloc
            curr_altloc=tmp_altloc;
          }
        else // Since there is no change in the altloc :
          {
            // We sum the corresponding occupency
            occups.at(tmp_altloc)+=atof(ligneA.substr(54,6).c_str());
            // And add the line to the list for this altloc
            atom_occ.at(tmp_altloc).push_back(LposA);
          }
      }// END WHILE

////////// STEP3 : Selection
    // Now we have all alternative atoms
    // We search for the best one, i.e. the one with the best occupancy

     best_occup =0;
     sum_occup =0;
     tmp_altloc="";// Contains the name of the best alternative position
     // Here we sum all occupancy
     // And select the alternative position set based on the best occupancy
     for (map<string,double>::iterator
          it = occups.begin();
          it != occups.end();
          ++it)
       {
         // the second condition is for alternate position
         // with occupancy at 0 (ex:NZ ALYS E 676 in 1N62)
         if ((*it).second > best_occup
         || (tmp_altloc=="" && (*it).second==0))
         {best_occup=(*it).second;
          tmp_altloc=(*it).first;}

         sum_occup+=(*it).second;

       }
     for (map<string,listLine>::iterator
          it = atom_occ.begin();
          it != atom_occ.end();
        ++it)
       {
         // Occupancy sum must be the same value has the number of alternative atoms
         // Otherwise there is an issue
         if ((fabs(sum_occup-(double)(*it).second.size()) > 1))
             cerr<< PDBName<<"|ALTLOC|WARNI\tOccup sum does not equal to the number of atoms for "
                 << curr_reschain << " altloc "<< (*it).first
                 << " :: SUM:" << sum_occup
                 << " ATOMS:"<< (*it).second.size()<<endl;

      }

     listLine& sellist = atom_occ.at(tmp_altloc);
     const size_t nAtmAP=sellist.size();
     for (size_t posAP = 0; posAP < nAtmAP; ++posAP)
     {
         const unsigned int &LposA = sellist.at(posAP);
         std::string& ligneA = filelines.at(LposA);
         ligneA.replace(16,1," ");
         new_atoms.push_back(LposA);
         if (for_scPDB) cout << PDBName<<"|ALTPOS|MODIFS\tIN \t" << ligneA<<"\n";

     }
     if (for_scPDB)
         for (map<string,listLine>::iterator
              it = atom_occ.begin();
              it != atom_occ.end();
              it++)
           {
             if ((*it).first == tmp_altloc)continue;
             listLine& sellist = (*it).second;
             const size_t nAtmAP=sellist.size();
             for (size_t posAP = 0; posAP < nAtmAP; ++posAP)
             {
                 const unsigned int &LposA = sellist.at(posAP);
                 std::string& ligneA = filelines.at(LposA);
                   cout << PDBName<<"|ALTPOS|MODIFS\tOUT \t" << ligneA<<"\n";

             }
           }




  }// END for size_t pos



  atoms.clear();
  atoms.insert( atoms.end(), new_atoms.begin(), new_atoms.end() );
if (for_scPDB)
    cout <<PDBName<<"|ALTPOS|STATIS\t"<<NLine<<"\t Number of alternative atoms \n"
     <<PDBName<<"|ALTPOS|STATIS\t"<<NRes<<"Number of alternative residues \n"
     <<PDBName<<"|ALTPOS|END\t"<<endl;
}





void PDBConvert::renumAtms()
{
  if (for_scPDB){ cout << PDBName<<"|ATMNUM|START\t"<<endl<<PDBName<<"|ATMNUM|HEADER\t##### Renum atoms : "<<endl;
    }
  map<int,int> Old_to_New;
  unsigned int nA=1, oldnA;
  ostringstream oss;

  const size_t              nAtomsSize   =atoms.size();


  //Scanning atom lines
  for (size_t pos=0; pos < nAtomsSize;pos++)
  {
  // Getting LIGNE:
    const unsigned int &Lpos = atoms.at(pos);
    std::string &ligne=filelines.at(Lpos);

  // Take the old value
      oldnA = atoi(ligne.substr(6,5).c_str());
  // Storing in a table the link between old and new value
      Old_to_New.insert(pair<int,int>(oldnA,nA));


      oss.str("");
      if (nA < 10000)  oss<<" ";
      if (nA <  1000)  oss<<" ";
      if (nA <   100)  oss<<" ";
      if (nA <    10)  oss<<" ";
      oss << nA;
      ligne.replace(6,5,oss.str());// Update the atom line
      nA++;
  }


  int new_vals[5];
  vector<string> vals;
  // Updating the connect line
  for (size_t pos=0; pos < connect.size();pos++)
  {
  // Getting LIGNE:

    const unsigned int &Lpos = connect.at(pos);
    std::string &ligne=filelines.at(Lpos);
    // A connect is defined as follow :
    //        COLUMNS   DATA TYPE   FIELD   DEFINITION
    //        1-6       Record      name    "CONECT"
    //        7-11      Integer     serial  Atom serial number
    //        12-16     Integer     serial  Serial number of bonded atom
    //        17-21     Integer     serial  Serial number of bonded atom
    //        22-26     Integer     serial  Serial number of bonded atom
    //        27-31     Integer     serial  Serial number of bonded atom
    // STEP 1-So we first check if the first atom serial number still exist
    //        (it can be deleted with discardAtm function).=> delete line
    // STEP 2-We check all the next values. We delete values of the deleted atoms and updating others
    // STEP 3-If all next values are deleted, we delete the line
    // STEP 4-Recreating the CONECT line
    // STEP 5-Adding space to have 80 characters



    string new_line = ligne;
    vals.clear();
    vals.push_back(ligne.substr(0,6));      //0
    vals.push_back(ligne.substr(6,5));      //1
    vals.push_back(ligne.substr(11,5));     //2
    vals.push_back(ligne.substr(16,5));     //3
    vals.push_back(ligne.substr(21,5));     //4
    vals.push_back(ligne.substr(26,5));     //5


    //STEP1
    if (Old_to_New.find(atoi(vals.at(1).c_str())) == Old_to_New.end())
    {
        if(for_scPDB)cout << PDBName<<"|RENA\t atm del:"<<new_line<<endl;

       connect.erase(connect.begin()+pos);
       --pos;
        continue;
    }
    else new_vals[0] = Old_to_New.at(atoi(vals.at(1).c_str()));

    new_vals[1]=-100;
    new_vals[2]=-100;
    new_vals[3]=-100;
    new_vals[4]=-100;


    //STEP2
    for (unsigned int nT = vals.size()-1; nT >0; nT--)
      {
        if (Old_to_New.find(atoi(vals.at(nT).c_str())) == Old_to_New.end()) vals.erase(vals.begin()+nT);
        else new_vals[nT-1]=Old_to_New.at(atoi(vals.at(nT).c_str()));
      }
    // if (for_scPDB) cout << vals.size()<<":"<< new_vals[0]<<"\t"<< new_vals[1]<<"\t"<< new_vals[2]<<"\t"<< new_vals[3]<<"\t"<< new_vals[4]<<"\t"<<endl;
    //STEP3
    if (vals.size() <= 2){
        if(for_scPDB)cout << PDBName<<"|ATMNUM|MODIFS\t del connect"<<endl;
        connect.erase(connect.begin()+pos);
        --pos;
        continue;
    }

    //STEP4
    oss.str("");
    oss<< "CONECT";
    for (unsigned int nT = 0; nT <= 4; nT++)
      {
        if (new_vals[nT] == -100) continue;
        if (new_vals[nT] < 10000)  oss<<" ";
        if (new_vals[nT] <  1000)  oss<<" ";
        if (new_vals[nT] <   100)  oss<<" ";
        if (new_vals[nT] <    10)  oss<<" ";
        oss << new_vals[nT];

      }

    //STEP5
    while (oss.str().length() < 80) oss <<" ";
    if(for_scPDB)
      cout << PDBName<<"|ATMNUM|MODIFS\tOUT\t" <<ligne<<endl
           << PDBName<<"|ATMNUM|MODIFS\tIN \t" << oss.str()<<endl;
    ligne = oss.str();

  }


  if (for_scPDB)cout <<PDBName<<"|ATMNUM|END\t"<<endl;
}


void PDBConvert::toComplex(Complex &complex) const
{
    vector<string> curr_lines;
    curr_lines.push_back(filelines.at(header_line));
    for (size_t pos=0; pos < atoms.size();++pos) curr_lines.push_back(filelines.at(atoms.at(pos)));
    for (size_t pos=0; pos < connect.size();++pos) curr_lines.push_back(filelines.at(connect.at(pos)));
    if (master != -1)curr_lines.push_back(filelines.at(master));
    curr_lines.push_back("END");
//for (size_t pos=0; pos< curr_lines.size(); ++pos) cout << curr_lines.at(pos)<<endl;
    MoleReader mr(curr_lines,FileFormat::PDB);

    mr.loadInComplex(complex,MoleType::PROTEIN);

}






