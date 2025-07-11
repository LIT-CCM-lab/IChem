#include "headers/ICTools/switch.h"
#include "headers/ICCalcs/interaction.h"
#include "headers/ICPars/molereader.h"
#include "headers/ICCalcs/volsite.h"
#include "headers/ICMole/box.h"
using namespace std;
using namespace ICMole;

void IChemSwitch::helpScoring()     const
{
    cout
//            << " SCORING FUNCTION GENERATOR : "                                                  <<endl
//            << "    scoring prot lig name "                                                      <<endl

//            << "###########################################################################"<<endl
            << endl;
}

void IChemSwitch::scoring()      const throw(ICMole::MoleExcept)
{
    const unsigned int InputSize = (unsigned int)Input_Values.size();
    if (InputSize!=3){helpScoring();throw MoleExcept(9010301,"IChem::scoring","Not enough parameters");}
    const std::string protein = Input_Values.at(0);
    const std::string ligand  = Input_Values.at(1);
    const std::string name  = Input_Values.at(2);
    Complex cplx;
    MoleReader mread;


    try{
        cplx.clear();
        if (verbose) cout <<"# READING MOLECULES"<<endl;
        mread.loadNewFile(protein); mread.loadInComplex(cplx,MoleType::PROTEIN);
        mread.loadNewFile(ligand);  mread.loadInComplex(cplx,MoleType::LIGAND);

        if (cplx.getNumMolecule(MoleType::LIGAND) ==0 || cplx.getMole(MoleType::LIGAND)->numAtom()==0) throw MoleExcept(9010302,"IChem::runInts","No ligand found in " + ligand);
        if (cplx.getNumMolecule(MoleType::PROTEIN)==0 || cplx.getMole(MoleType::PROTEIN)->numAtom()==0) throw MoleExcept(9010303,"IChem::runInts","No protein found in "+ protein);

        Grid grid;

        Interactions ins(cplx);

        InterResults results;

        ins.detectInteractions(*cplx.getMole(MoleType::LIGAND),results,true);
        ins.interToMOL2(results,false,false,true,false);


        VolSite volsite(*cplx.getMole(MoleType::PROTEIN),grid);
        volsite.proj(55-20,false);
        // Do Volsite after calcints to have 4.5 grid in calcints and 1.5 in volsite.
        Molecule &intmole = results.Ints;
        vector<double> FGPCount;
        vector<vector<double> > FGPBuried;
        for (size_t i=0;i<8;++i)
        {
            FGPCount.push_back(0);
            if (i==0)continue;
            vector<double> tmp;
            for (size_t j=0;j<10;++j)tmp.push_back(0); FGPBuried.push_back(tmp);
        }


        for (size_t iAtm=0; iAtm < intmole.numAtom();++iAtm)
        {
            const Atom& atm= intmole.getAtom(iAtm);
            Box* box = volsite.getGrid().getBox(atm.fixpos);
            //cout << atm.getName()<<"\t"<<box->getId()<< " " << atm.fixpos.calcDist(box->fixpos)<< " " << volsite.getProjValue(*box)<<endl;
            double pos= floor((volsite.getProjValue(*box)-30)/10);
//            cout << volsite.getProjValue(*box) << "\t col : " << pos+1 << endl;
            if (atm.getName()=="CA")      {FGPBuried[0][pos]+=1;FGPCount[1]++;FGPCount[0]++;}
            else if (atm.getName()=="CZ") {FGPBuried[1][pos]+=1;FGPCount[2]++;FGPCount[0]++;}
            else if (atm.getName()=="O")  {FGPBuried[2][pos]+=1;FGPCount[3]++;FGPCount[0]++;}
            else if (atm.getName()=="OD1"){FGPBuried[3][pos]+=1;FGPCount[4]++;FGPCount[0]++;}
            else if (atm.getName()=="N")  {FGPBuried[4][pos]+=1;FGPCount[5]++;FGPCount[0]++;}
            else if (atm.getName()=="NZ") {FGPBuried[5][pos]+=1;FGPCount[6]++;FGPCount[0]++;}
            else if (atm.getName()=="Zn") {FGPBuried[6][pos]+=1;FGPCount[7]++;FGPCount[0]++;}


        }
        for (size_t i=0; i<7;++i)
        {

            for (size_t j=0;j<10;++j)
            {

                if (FGPBuried[i][j]==0)continue;
                FGPBuried[i][j]=FGPBuried[i][j]/FGPCount[i+1]*100.0;

            }
            FGPCount[i+1]=FGPCount[i+1]/FGPCount[0]*100.0;

        }

        cout <<name<<"\t";
        for (size_t i=0;i<=7;++i) cout << FGPCount[i]<< " ";
        for (size_t i=0; i<7;++i)
        {

            for (size_t j=0;j<9;++j)
            {

                cout << FGPBuried[i][j]<< " ";
            }

        }
        cout <<endl;

        /*
     HBOND_PROT=1,
      HBOND_LIG=2,
     IONIC_PROT=3,
      IONIC_LIG=4,
    HYDROPHOBIC=5,
          METAL=6,
     ARFACEFACE=7,
     AREDGEFACE=8,  */


    }catch (MoleExcept &e)
    {

        throw;
    }


}

