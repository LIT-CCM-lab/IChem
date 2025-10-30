#include <stdio.h>
#include "headers/ICTools/switch.h"
#include "headers/ICCalcs/volsite.h"
#include "headers/ICPars/pharmwriter.h"

using namespace std;
using namespace ICMole;


void IChemSwitch::helpVolSite() const {
    cout
        << "volsite - Cavity detection and druggability prediction" << endl
        << "Usage: IChem [options] volsite prot [lig]" << endl
        << endl
        << "Description:" << endl
        << "  Detect cavities on the protein surface and estimate their druggability" << endl
        << "  If a ligand is provided, VolSite restricts detection to its binding site" << endl
        << "  Optionally outputs descriptors, SVM features, or cavity-based pharmacophores" << endl
        << endl
        << "Options:" << endl
        << "  -step N (1.5)   Grid spacing (Å)" << endl
        << "  -boxS N (20)    Main box edge length (Å)" << endl
        << "  -b N (55)       Minimal buriedness threshold" << endl
        << "  -name N         Prefix for output cavity names" << endl
        << "  -n N (5)        Min neighbors for buried boxes" << endl
        << "  -nPTS N (70)    Min number of cubes to define a cavity" << endl
        << "  -NPTS N (400)   Max number of cubes to define a cavity" << endl
        << "  --dna           Include DNA in the protein" << endl
        << "  --cofactor      Include cofactors" << endl
        << "  --solvent       Include solvent molecules" << endl
        << "  --hydrogen      Include hydrogens in detection" << endl
        << "  --desc          Output cavity descriptors (descriptor.txt)" << endl
        << "  --svm           Output SVM features for druggability model" << endl
        << "  --pharm         Generate cavity-based pharmacophores (.chm, .mol2, .pml, .plp)" << endl
        << "  --outExclu      Add exclusion spheres to pharmacophore output" << endl
        << "  -drog N         Set observed druggability value" << endl
        << endl
        << "Example:" << endl
        << "  IChem volsite protein.mol2              # Detect all cavities" << endl
        << "  IChem volsite protein.mol2 ligand.mol2  # Detect ligand-binding site cavities" << endl
        << "  IChem --pharm volsite protein.mol2 ligand.mol2  # Generate pharmacophore" << endl
        << endl
        << "---------------------------------------------------------------------------" << endl
        << endl;
}


// Too many options commented, to be revised
// Meanwhile, I will have my version with no comments
// void IChemSwitch::helpVolSite() const {
//     cout << "Volsite help" <<endl
//             << "VOLSITE - Cavity detection in a mol2 file"                                  <<endl
//          << "    volsite prot lig                                                    (1)"<<endl
//          << "    volsite prot                                                        (2)"<<endl
//          << endl
//          << "   [General options]"<<endl
//          << "      -step  N (1.5)    Edge length of each box                 (Angstroem)"<<endl
//          << "      -boxS  N (20)     Edge length of the main box             (Angstroem)"<<endl
//          << "      -b     N (55)     Minimal threshold for buriedness "                  <<endl
//          << "      -name  N          PDB Name for output cavity names "                  <<endl
//          << "      -n     N  (5)     Minimal neighbours for buried cavity boxes"                  <<endl
//          << "      -nPTS  N (70)     Minimal number of cubes to consider it a cavity "   <<endl
//          << "      -NPTS  N (400)    Maximal number of cubes to consider it a cavity "   <<endl
//          << "      --dna             Consider DNA as part of the protein"                <<endl
//          << "      --cofactor        Consider cofactor as part of the protein "          <<endl
//          << "      --solvent         Consider solvent as part of the protein "           <<endl
//             //         << "      --ligand          Keep all ligand cubes as part of cavities (non fonctionnel)" <<endl
//             //<< "      --outputsite      Output protein site around cavities"<<endl
//          << "      --hydrogen        Consider hydrogens "                                <<endl
//          << "      --desc            Write a descriptor file name descriptor.txt"        <<endl
//          << "      --svm             Build a svm property file "                       <<endl
//          << "      -drog  N          Observed druggability "                             <<endl
//             //         << endl
// //         << "      --pharm           Generate a pharmacophore (.chm) from cavity "  <<endl
// //         << "      --outExclu        Output exclusion sphere in pharmacophore file" <<endl
//             //         << "             [1]        Closest cavity cubes "                            <<endl
//             //         << "             [2]        All cavity cubes "                                <<endl
//             //         << "             [3]        All cavity cubes with buriedness > 80"            <<endl
//             //         << "             [4]        Centered cubes with buriedness > 80"              <<endl
//             //         << "             [5]        Optimized pharmacophore"              <<endl

//             //         << "             [3]          Barycenter of cube shared property (-b optimized)"    <<endl
//             //         << "             [4]          Closest cavity cube to protein     (-b optimized)"    <<endl
//             //         << "             [5]          Closest cavity cube to barycenter  (-b optimized)"   <<endl
//             //         << "             [6]          3 with less hydrophobic features   (-b optimized)"   <<endl
//             //         << "             [7]          3 (hydrophobics features must be link to at least 2 residues)   (-b optimized)"   <<endl
//             //         << "   [Pharmacophore options]"<<endl
//             //         << "      -phaP  N (20)     Lowest percentage of standard cavity keep to generate pharmacophore" <<endl
//             //         << "      -phaC  N (50)     Number of cavity cube keep to generate pharmacophore "     <<endl
//             //         << "      --tol             Tolerance varying with numbers of neighbors"     <<endl
//             //         << "      --weight  N (50)  Weight varying with features percentage"     <<endl
//          << "###########################################################################"<<endl
//          <<endl;

// }

void IChemSwitch::volsite() const 
{
    Complex cplx;
    MoleReader mread;
    const unsigned int InputSize = (unsigned int)Input_Values.size();
    if (InputSize ==0) {helpVolSite(); return;}

    if (InputSize != 1 && InputSize != 2 && InputSize != 3) {helpVolSite();throw MoleExcept(9010501,"IChem::runVolSite","Not enough parameters");}
    bool w_lig = false;
    bool w_chm = false;
    bool pharma_end = false;
    if (InputSize >= 2) w_lig=true;
    if (InputSize == 3) w_chm=true;
    const std::string RProt = Input_Values.at(0);
    mread.loadNewFile(RProt);
    mread.loadInComplex(cplx,MoleType::PROTEIN);
    if (cplx.getNumMolecule(MoleType::PROTEIN)==0){helpVolSite(); throw MoleExcept(9010502,"IChem::runVolSite","No protein found in "+ RProt);}
    if (InputSize < 2 ) cout << "Run VolSite on : "<< RProt << endl;


    if (w_lig)
    {
        const std::string RLig  = Input_Values.at(1);
        mread.loadNewFile(RLig);  mread.loadInComplex(cplx,MoleType::LIGAND);
        if (cplx.getNumMolecule(MoleType::LIGAND)==0) {helpVolSite(); throw MoleExcept(9010503,"IChem::runVolSite","No ligand found in "+ RLig);}
        cout << "Ligand file : "<< RLig << endl;
    }

    if (w_chm){ //conversion de chm en mol2 à la volée
        const std::string RChm  = Input_Values.at(2);
        mread.loadNewFile(RChm);  mread.loadInComplex(cplx,MoleType::PHARMACOPHORE);
        if (cplx.getNumMolecule(MoleType::PHARMACOPHORE)==0) {helpVolSite(); throw MoleExcept(9010503,"IChem::runVolSite","No pharma found in "+ RChm);}
        cout << "Chm file : "<< RChm << endl;
    }

    Coords box_size(20,20,20);
    string name="";
    double  box_step=1.5;
    double pharma_step=1.5;
    int projs=55,nPTS=70,NPTS=400,drog=0,lissage=20,proj_len=8;
    unsigned int minNeighbors=5;
    int PHA_mod=1,PHA_P=20,PHA_C=50;
    bool cofactor=false,ligand=false,solvent=false,hydrogen=false,svm=false,outsite=false,nucleic =false,pharmacophore=false,outputexclu=false,toler=false,weight=false;
    bool Bweight = false;

    /* FOR Franck TRIALS */
    for (std::map<std::string,std::vector<std::string> >::const_iterator
         it = Opt_Values.begin();
         it != Opt_Values.end();
         it++)
    {
        const std::string &opt_name = (*it).first;
        const vector<string> & opt_val = (*it).second;
        const std::string& value = opt_val.at(0);
        //  -rn -cn --values -sim -outInt -match -max -size --all_cliques -score
        if (opt_name.compare("-step")                 == 0) { box_step=atof(value.c_str());                           }
        else if (opt_name.compare("-boxS")            == 0) { double nI=atof(value.c_str()); box_size.setCoords(nI,nI,nI);}
        else if (opt_name.compare("-b")               == 0) { projs=atof(value.c_str());                                 }
        else if (opt_name.compare("-name")            == 0) { name=value;                                                }
        else if (opt_name.compare("--dna")            == 0) { nucleic = true;                                            }
        else if (opt_name.compare("--cofactor")       == 0) { cofactor=true;                                             }
        else if (opt_name.compare("--ligand")         == 0) { ligand=true;                                               }
        else if (opt_name.compare("--solvent")        == 0) { solvent=true;                                              }
        else if (opt_name.compare("--hydrogen")       == 0) { hydrogen=true;                                             }
        else if (opt_name.compare("-nPTS")            == 0) { nPTS=atof(value.c_str());                                  }
        else if (opt_name.compare("-NPTS")            == 0) { NPTS=atof(value.c_str());                                  }
        else if (opt_name.compare("--desc")           == 0) { svm =true;                                                 }
        else if (opt_name.compare("--svm")            == 0) { svm=true;                                                  }
        else if (opt_name.compare("--outputsite")     == 0) { outsite =true;                                             }
        else if (opt_name.compare("-drog")            == 0) { drog =atof(value.c_str());                                 }
        else if (opt_name.compare("-n")               == 0) { minNeighbors=atof(value.c_str());                          }
        else if (opt_name.compare("-l")               == 0) { lissage=atof(value.c_str());                               }
        else if (opt_name.compare("-p")               == 0) { proj_len=atof(value.c_str());                              }
        else if (opt_name.compare("--pharm")          == 0) { pharmacophore=true;PHA_mod=34;               }
        else if (opt_name.compare("--outExclu")       == 0) { outputexclu=true;                                          }
        else if (opt_name.compare("-pha")             == 0) { pharmacophore=true;PHA_mod=atof(value.c_str());            }
        else if (opt_name.compare("-phaP")            == 0) { PHA_P=atof(value.c_str());                                 }
        else if (opt_name.compare("-phaC")            == 0) { PHA_C=atof(value.c_str());                                 }
        else if (opt_name.compare("--tol")            == 0) { toler=true;                                                }
        else if (opt_name.compare("--weight")         == 0) { weight=true;                                               }

        else throw MoleExcept(9010504,"IChem::runVolSite","Unrecognized optionl : "+opt_name);
    }

    //% faire la verbose
//        cout << "Assigning properties to protein atoms"<<endl;
//        if (PHA_mod == 1) { PHA_mod = 1;}
//        else if (PHA_mod == 2){ PHA_mod = 8;}
//        else if (PHA_mod == 3){ PHA_mod = 18;}
//        else if (PHA_mod == 4){ PHA_mod = 28;}
//        else if (PHA_mod == 5){ PHA_mod = 33;}


    Grid grid;
    Grid grid_pha;
    if (cofactor){ cout << "Keeping cofactor as protein atoms"<<endl; Residu::Rules[MoleType::PROTEIN][ResType::COFACTOR]=MoleType::PROTEIN;}
    if (nucleic) { cout << "Keeping Nucleic as protein atoms" <<endl; Residu::Rules[MoleType::PROTEIN][ResType::NUCLEIC]=MoleType::PROTEIN;}
    if (solvent) {
        cout << "Keeping solvent as protein atoms "<<endl; Residu::Rules[MoleType::PROTEIN][ResType::WATER]=MoleType::PROTEIN;Residu::Rules[MoleType::LIGAND][ResType::WATER]=MoleType::LIGAND;}
    if (hydrogen){ cout << "Keeping hydrogens as protein atoms"<<endl;grid.setWHydrogen(true);}
    if (pharmacophore && w_lig) {
        cout << "Generate pharmacophore from cavity" << endl;
        box_step =1;
        if ( box_step != 1.5 ){
            pharma_step = box_step;
            box_step = 1.5;
        }
        if (PHA_mod >= 8 ){
            PHA_mod -=8;


            VolSite volsite_pha(*cplx.getMole(MoleType::PROTEIN),*cplx.getMole(MoleType::LIGAND),grid_pha,pharma_step,box_size,proj_len,pharmacophore);

            if (nPTS != 70) volsite_pha.setMinSize(nPTS);
            if (NPTS != 400)volsite_pha.setMaxSize(NPTS);
            if (drog !=0) volsite_pha.setDruggabilty(drog);
            volsite_pha.setName(cplx.getMole(MoleType::PROTEIN)->getName());
            volsite_pha.proj(projs-lissage,false);

            volsite_pha.treatments(minNeighbors,2.4);
            //            cout << "treatmnts here"<<endl;
            volsite_pha.setCavityKP(PHA_P);
            volsite_pha.setCavityKC(PHA_C);
            volsite_pha.setTolerance(toler);
            volsite_pha.setWeigth(weight);
            // creation d'un cubetoPHA directement ici
            // Pharma .mol2


            // Generation de l'ensemble des poitns utilisable
            // Un seule bloc continu, plusieurs départ possible si la cavité n'est pas centré.
            volsite_pha.clusterPHA((16+box_size.x)/pharma_step/2,(16+box_size.y)/pharma_step/2,(16+box_size.z)/pharma_step/2);
            if (!volsite_pha.getPhaSuccess()) {
                volsite_pha.clusterPHA((6+16+box_size.x)/pharma_step/2,(16+box_size.y)/pharma_step/2,(16+box_size.z)/pharma_step/2);
            }
            if (!volsite_pha.getPhaSuccess()) {
                volsite_pha.clusterPHA((16-6+box_size.x)/pharma_step/2,(16+box_size.y)/pharma_step/2,(16+box_size.z)/pharma_step/2);
            }
            if (!volsite_pha.getPhaSuccess()) {
                volsite_pha.clusterPHA((16+box_size.x)/pharma_step/2,(16+6+box_size.y)/pharma_step/2,(16+box_size.z)/pharma_step/2);
            }
            if (!volsite_pha.getPhaSuccess()) {
                volsite_pha.clusterPHA((16+box_size.x)/pharma_step/2,(16-6+box_size.y)/pharma_step/2,(16+box_size.z)/pharma_step/2);
            }
            if (!volsite_pha.getPhaSuccess()) {
                volsite_pha.clusterPHA((16+box_size.x)/pharma_step/2,(16+box_size.y)/pharma_step/2,(16+6+box_size.z)/pharma_step/2);
            }
            if (!volsite_pha.getPhaSuccess()) {
                volsite_pha.clusterPHA((16+box_size.x)/pharma_step/2,(16+box_size.y)/pharma_step/2,(16-6+box_size.z)/pharma_step/2);
            }
            if (!volsite_pha.getPhaSuccess()) {
                volsite_pha.clusterPHA((12+16+box_size.x)/pharma_step/2,(16+box_size.y)/pharma_step/2,(16+box_size.z)/pharma_step/2);
            }
            if (!volsite_pha.getPhaSuccess()) {
                volsite_pha.clusterPHA((16-12+box_size.x)/pharma_step/2,(16+box_size.y)/pharma_step/2,(16+box_size.z)/pharma_step/2);
            }
            if (!volsite_pha.getPhaSuccess()) {
                volsite_pha.clusterPHA((16+box_size.x)/pharma_step/2,(16+12+box_size.y)/pharma_step/2,(16+box_size.z)/pharma_step/2);
            }
            if (!volsite_pha.getPhaSuccess()) {
                volsite_pha.clusterPHA((16+box_size.x)/pharma_step/2,(16-12+box_size.y)/pharma_step/2,(16+box_size.z)/pharma_step/2);
            }
            if (!volsite_pha.getPhaSuccess()) {
                volsite_pha.clusterPHA((16+box_size.x)/pharma_step/2,(16+box_size.y)/pharma_step/2,(16+12+box_size.z)/pharma_step/2);
            }
            if (!volsite_pha.getPhaSuccess()) {
                volsite_pha.clusterPHA((16+box_size.x)/pharma_step/2,(16+box_size.y)/pharma_step/2,(16-12+box_size.z)/pharma_step/2);
            }

            if (!volsite_pha.getPhaSuccess()) {
                cout << "Cavity detection failed" << endl;
            }

            if (PHA_mod == 13) {
                cout << "kmeans begin" <<endl;
                volsite_pha.KmeansPHA(cplx.getMole(MoleType::PROTEIN)->numAtom()/20, false);
                cout << "kmeand end" << endl;
            }
            MoleWriter mw;
            Molecule mole_pha;
            PharmList pharmaco;
            ostringstream ossname;
            ossname.str("");
            string size="ALL";
            if (name.length() != 0) ossname << name << "_Pharmacophore";
            else ossname << "Pharmacophore";

            ossname << ".mol2";

            Coords const &toto = cplx.getMole(MoleType::LIGAND)->getFixpos();


            // Calcul du pharmacophore a partir des points de la cavité
            //STD pha_mod = 26
            cplx.getMole(MoleType::PROTEIN)->ringPerception();
            volsite_pha.cubeToPHA(mole_pha,pharmaco,(*cplx.getMole((MoleType::LIGAND))) , ossname.str(),size,PHA_mod,true,outputexclu);//PHA_mod);
            //            cout << mole_pha.tocd String() <<endl;
            mw.newFile(ossname.str(),FileFormat::MOL2); mw.writeMOL2(&mole_pha);


            // Ecris le pharmacophore en chm
            ossname.str("");
            if (name.length() != 0) ossname << name << "_Pharmacophore";
            else ossname << "Pharmacophore";
            ossname<<".chm";
            PharmWriter pw(ossname.str(),FileFormat::CHM);
            string name_str = ossname.str();
            pw.writeCHM(pharmaco,name_str.c_str(),Bweight); //% Recupéré Bweightcelui de volsit
            ossname.str("");
            ossname << "Pharmacophore.pml";
            string name_pml = ossname.str();
            pw.writePML(pharmaco,name_pml.c_str(),Bweight);
            ossname.str("");ossname << "Pharmacophore";
            //            ossname << "2";   ossname<<".chm";
//            name_str = ossname.str();
//            pw.writeCHMss(pharmaco,name_str,Bweight); //% Recupéré Bweightcelui de volsite


            // Ecris le pharmacophore en mol2 et en subdivision
            //            MoleWriter mwss;
            //            mwss.writeMOL2ss(&mole_pha, PHA_mod+8);



            //            MoleWriter lipha;
            //            ostringstream ossnamelipha;
            //            if (name.length() != 0) ossnamelipha << name << "_mod";
            //            else ossnamelipha << "ligand_pharma";
            //            ossnamelipha << ".mol2";
            //            lipha.newFile(ossnamelipha.str(),FileFormat::MOL2); lipha.writePharmaLig(cplx.getMole(MoleType::LIGAND));



            pharmacophore = false;
            pharma_end = true;

            //            cout << cplx.getMole(MoleType::PHARMACOPHORE)->toString() << endl;
//            cout << "Cavity center : " << toto.toString() << endl;
            if (w_chm && PHA_mod == 34 && pharma_step == 1){

                // selection de k parmis n
                //                const int n =10; const int k = 2;
                //                int comb[20] = {0};
                //                int i = 0;
                //                while (i >= 0) {
                //                    if (comb[i] < n + i - k + 1) { comb[i]++;
                //                        if (i == k - 1) {
                //                            for (int j = 0; j < k; j++) {
                //                                printf("%d ", comb[j]);
                //                            }
                //                            printf("\n");
                //                        } else { comb[++i] = comb[i - 1]; }
                //                    } else i--; }

                //                cout << "i'm here" << endl;
                mw.newFile("chmmodified.mol2",FileFormat::MOL2); mw.writeMOL2(cplx.getMole(MoleType::PHARMACOPHORE));


                ofstream otfs;
                otfs.open("../distanceCHM.dsc",ios::out|ios::app);
                ofstream ofs;
                ofs.open("../descripteurCHM.dsc",ios::out|ios::app);
                Molecule &pha_desc = *cplx.getMole(MoleType::PHARMACOPHORE);
                for (ItCAtom itAtm= pha_desc.firstAtom(); itAtm != pha_desc.lastAtom(); ++itAtm)
                {
                    Atom &atm = **itAtm;
                    otfs << atm.getIdentifier() << "\t" <<atm.fixpos.calcDist(toto) << atm.getFormalCharge() << endl;
                    for (ItCAtom itAtm2= itAtm+1; itAtm2 != pha_desc.lastAtom(); ++itAtm2)
                    {
                        Atom &atm2 = **itAtm2;
                        //                    cout << atm.getName() << "-" <<atm2.getName() << "   ";
                        //                        cout << atm.getName() << " avec " << atm2.getName() << " " << atm.getFormalCharge() * atm2.getFormalCharge() << " d: "<< atm.fixpos.calcDist(atm2.fixpos) << endl;
                        ofs  << atm.getFormalCharge() * atm2.getFormalCharge() << " ";
                        ofs  << atm.fixpos.calcDist(atm2.fixpos) << endl;
                    }
                }
            }
        }
    }


    //    VolSite volsite(*cplx.getMole(MoleType::PROTEIN),*cplx.getMole(MoleType::LIGAND),grid,box_step);


    VolSite volsite(*cplx.getMole(MoleType::PROTEIN),grid,box_step,proj_len);

    if (ligand) volsite.setproteinOnly(false);

    if (svm)     { cout << "Build the file descriptor.txt"<<endl;}

    //volsite.setWithLigand(true);
    if (nPTS != 70) volsite.setMinSize(nPTS);
    if (drog !=0) volsite.setDruggabilty(drog);
    volsite.setName(cplx.getMole(MoleType::PROTEIN)->getName());

    volsite.proj(projs-lissage,false);

    volsite.treatments(minNeighbors,2.5);

    volsite.clusterCAV(projs);

    if (w_lig && !pharma_end)
    {
        double dist;
        volsite.output_cavs_lig(*cplx.getMole(MoleType::LIGAND),name,svm,dist,pharmacophore,PHA_mod);
        if (box_step != 1.5 || projs != 55)
        {
            cout << "WARNING - Druggability value only valid with default parameter"<<endl;
        }
        cout.precision(3);
        cout << "Protein name : " << ((name.length()) ? name : "Unknown")<<endl
             << "Druggability : " << dist << endl<< "Prediction : ";
        if (dist >0 ) cout << "Is Druggable"<<endl;
        else cout << "Not druggable"<<endl;
    }
    else
    {
        if (pharma_end == false ){
            volsite.output_cavs_all(name,svm,pharmacophore,PHA_mod,Input_Values);
            cout << "End" << endl;
        }
    }

}


