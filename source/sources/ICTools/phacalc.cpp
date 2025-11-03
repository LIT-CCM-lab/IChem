#include <stdio.h>
#include "headers/ICTools/switch.h"
#include "headers/ICCalcs/pharmaco.h"

using namespace std;
using namespace ICMole;
void IChemSwitch::helpPharma() const
{
    cout << "Pharmacophore - Generation of pharmacophore descriptor"                    <<endl
         << "    pharma pharmacophore.chm                                                    (1)"<<endl
         << endl
         << "   [General options]"<<endl
         << "   --div generate small pharmacophore "<<endl
         << "###########################################################################"<<endl
         <<endl;

}

void IChemSwitch::pharma() const 
{
    Complex cplx;
    MoleReader mread;
    const unsigned int InputSize = (unsigned int)Input_Values.size();
    if (InputSize ==0) {helpVolSite(); return;}

    if (InputSize != 1 ) {helpVolSite();throw MoleExcept(9010502,"IChem::runPharma","Too much parameters");}
    const std::string RChm  = Input_Values.at(2);
    mread.loadNewFile(RChm);  mread.loadInComplex(cplx,MoleType::PHARMACOPHORE);
    if (cplx.getNumMolecule(MoleType::PHARMACOPHORE)==0) {helpVolSite(); throw MoleExcept(9010503,"IChem::runPharma","No pharmacophore found in "+ RChm);}
    cout << "Run on: "<< RChm << endl;

    bool generate_div = false;

    for (std::map<std::string,std::vector<std::string> >::const_iterator
         it = Opt_Values.begin();
         it != Opt_Values.end();
         it++)
    {
        const std::string &opt_name = (*it).first;
        const vector<string> & opt_val = (*it).second;
        const std::string& value = opt_val.at(0);
        //  -rn -cn --values -sim -outInt -match -max -size --all_cliques -score
        if (opt_name.compare("--div")            == 0) { generate_div=true;                           }

        else throw MoleExcept(9010504,"IChem::runPharma","Unrecognized optionl : "+opt_name);
    }


    MoleWriter mw;
    mw.newFile("Chm2mol2.mol2",FileFormat::MOL2); mw.writeMOL2(cplx.getMole(MoleType::PHARMACOPHORE));
    ofstream ofs;
    ofs.open("descriptorCHM.dsc",ios::out|ios::app);
    Molecule &pha_desc = *cplx.getMole(MoleType::PHARMACOPHORE);
    int nb_atm = 0;
    for (ItCAtom itAtm= pha_desc.firstAtom(); itAtm != pha_desc.lastAtom(); ++itAtm)
    {
        Atom &atm = **itAtm;
        nb_atm++;
        for (ItCAtom itAtm2= itAtm+1; itAtm2 != pha_desc.lastAtom(); ++itAtm2)
        {
            Atom &atm2 = **itAtm2;
            //                    cout << atm.getName() << "-" <<atm2.getName() << "   ";
            cout << atm.getName() << " vs " << atm2.getName() << " " << atm.getFormalCharge() * atm2.getFormalCharge() << " dist: "<< atm.fixpos.calcDist(atm2.fixpos) << endl;
            ofs  << atm.getFormalCharge() * atm2.getFormalCharge() << " ";
            ofs  << atm.fixpos.calcDist(atm2.fixpos) << endl;
        }
    }
    if (generate_div){
        if ( nb_atm < 8 ){
            int natm=0;
            //            for (natm=0; natom < nb_atm; natm++ ){

        }
    }else {

    }
}


