//#define ICHEM_DEBUG
#include "headers/ICCalcs/volsite.h"
#include "headers/ICMole/complex.h"
#include "headers/ICMole/box.h"
#include "headers/ICPars/molewriter.h"
#include "headers/ICPars/pharmwriter.h"
#include "headers/ICMole/atom.h"
#include "headers/ICMole/bond.h"
#include "headers/ICSVM/model.h"
#include "headers/ICMole/molecule.h"
#include "headers/ICMole/pharmprop.h"
#include "headers/ICMole/cycle.h"
#include "headers/ICTools/switch.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


using namespace std;
using namespace ICMole;

double Angl_H = M_PI;        //180
double AngT_H = M_PI/3;      //60
double Angl_D = M_PI/1.65;   //109
double D120         = M_PI/1.5;    //120
double D180         = M_PI;
double angl_aro     = M_PI/2;  // 90

VolSite::VolSite(Complex &cp, const double &boxStep, const int proj_lens,bool def_alternative) throw(MoleExcept):
    grid(cp.genGrid(boxStep)),
    minProj(55),
    proteinOnly(true),
    min_size(70),
    max_size(400),
    drug_value(-1),
    name(""),
    L1(0.1),L2(1),L3(1.5),L4(2),L5(2.5),L6(4),L7(8)
{
    if (boxStep <= 0)
        throw MoleExcept(3010101,
                         "VolSite::VolSite",
                         "Box Step cannot be below or equal to 0"
                         );
    if (cp.getMole(MoleType::PROTEIN)==(Molecule*)NULL)
        throw MoleExcept(3010102,
                         "VolSite::VolSite",
                         "No protein found in complex"
                         );
    try
    {
        grid.setProjLength(proj_lens);
        if (grid.AllBoxes.empty())
        {
            grid.createComplexBoxes(cp,boxStep);
            for (size_t iMole=0; iMole < cp.getNumMolecule();++iMole)
            {
                Molecule &mole = **(cp.firstMole()+iMole);
                grid.colorCube(mole,true);
                defProps(mole);
            }
        }
        else
        {
            Coords center;double NAtm=0;
            for (size_t iMole=0; iMole < cp.getNumMolecule();++iMole)
            {
                Molecule &mole = **(cp.firstMole()+iMole);
                grid.rotateMole(mole,false);
                for (size_t iAtm=0; iAtm < mole.numAtom();++iAtm)
                {
                    const Atom& atm = mole.getAtom(iAtm);
                    if (!atm.isUsed())continue;
                    center+=atm.fixpos;
                    NAtm++;
                }
            }
            center/=NAtm;
            //  Box*box = grid.getBox(center);


            //cout << "NEW CENTER:"<<box->fixpos.toString();
            Coords max(-10000,-10000,-10000),min(10000,10000,10000);
            for (size_t iMole=0; iMole < cp.getNumMolecule();++iMole)
            {
                Molecule &mole = **(cp.firstMole()+iMole);
                grid.rotateMole(mole,false);
                for (size_t iAtm=0; iAtm < mole.numAtom();++iAtm)
                {
                    const Atom& atm = mole.getAtom(iAtm);
                    if (!atm.isUsed())continue;
                    if (atm.rotpos.x > max.x)max.x=atm.rotpos.x;
                    if (atm.rotpos.y > max.y)max.y=atm.rotpos.y;
                    if (atm.rotpos.z > max.z)max.z=atm.rotpos.z;
                    if (atm.rotpos.x < min.x)min.x=atm.rotpos.x;
                    if (atm.rotpos.y < min.y)min.y=atm.rotpos.y;
                    if (atm.rotpos.z < min.z)min.z=atm.rotpos.z;


                }
            }
            const int BxMargin=static_cast<int>(floor(8/boxStep));;
            const int maxX=(int)(ceil(max.x)+BxMargin);//+8 for projection
            const int maxY=(int)(ceil(max.y)+BxMargin);
            const int maxZ=(int)(ceil(max.z)+BxMargin);
            const int minX=(int)(ceil(min.x)-BxMargin);
            const int minY=(int)(ceil(min.y)-BxMargin);
            const int minZ=(int)(ceil(min.z)-BxMargin);
            const int& numPtRangeK= grid.numPtRangeK;
            const int& numPtRangeI= grid.numPtRangeI;
            const int& numPtRangeJ= grid.numPtRangeJ;

            for (int k=0; k < numPtRangeK; k++)
                for (int j=0; j <  numPtRangeJ; j++)
                    for (int i=0; i <  numPtRangeI; i++)
                    {

                        Box &bx= *grid.AllBoxes.at(i+j*numPtRangeI+k*numPtRangeI*numPtRangeJ);
                        bx.Props.clear();
                        //cout <<bx.getId()<<" " << bx.rotpos.toString();
                        if (bx.rotpos.x >= minX-BxMargin && bx.rotpos.x <=maxX+BxMargin
                                &&bx.rotpos.y >= minY-BxMargin && bx.rotpos.y <=maxY+BxMargin
                                &&bx.rotpos.z >= minZ-BxMargin && bx.rotpos.z <=maxZ+BxMargin)
                        {
                            bx.setUse(true);
                            if (bx.rotpos.x >= minX && bx.rotpos.x <=maxX
                                    &&bx.rotpos.y >= minY && bx.rotpos.y <=maxY
                                    &&bx.rotpos.z >= minZ && bx.rotpos.z <=maxZ)
                            {


                                bx.Props.setNoInfo(true);
                                bx.setMargin(false);

                            }
                            else
                            {
                                bx.setMargin(true);

                            }
                        }
                        else
                        {
                            //cout << "OUT"<<endl;
                            bx.setUse(false);
                            // cout << bx.rotpos.toString()<<endl;
                        }


                    }

            //grid.printInFile("_BETS",grid.AllBoxes,false,true);
            for (size_t iMole=0; iMole < cp.getNumMolecule();++iMole)
            {
                Molecule &mole = **(cp.firstMole()+iMole);
                grid.colorCube(mole,false,true);
                defProps(mole);
            }
        }





    }
    catch (MoleExcept &e)
    {
        e.addTrace("VolSite::VolSite");
        throw;
    }


}
VolSite::VolSite(Molecule &protein,
                 Molecule &ligand,
                 Grid &grid,
                 const double &boxStep,
                 const Coords& box_size, const int proj_lens, bool def_alternative) throw(MoleExcept):
    grid(grid),
    minProj(55),
    proteinOnly(false),
    min_size(70),
    max_size(400),
    drug_value(-1),
    name(""),L1(0.1),L2(1),L3(1.5),L4(2),L5(2.5),L6(4),L7(8),clusterPHA_succes(false)
{
    try
    {
        grid.setProjLength(proj_lens);
        grid.createComplexBoxes(protein,ligand,boxStep,box_size);
        grid.colorCube(protein,true,true);
        grid.colorCube(ligand);
        defProps(protein,def_alternative);
    }
    catch (MoleExcept &e)
    {
        e.addTrace("VolSite::VolSite");
        throw;
    }


}

VolSite::VolSite(Complex &cp, Molecule &ligand, const double &boxStep,
                 const Coords& box_size, const int proj_lens, bool def_alternative) throw(MoleExcept):
    grid(cp.genGrid(boxStep)),
    minProj(55),
    proteinOnly(false),
    min_size(70),
    max_size(400),
    drug_value(-1),
    name(""),
    L1(0.1),L2(1),L3(1.5),L4(2),L5(2.5),L6(4),L7(8)
{
    /* grid.setProjLength(8);
    grid.createComplexBoxes(cp,ligand,boxStep,box_size);
   */


    grid.setProjLength(proj_lens);
    grid.rotateMole(ligand);
    // STEP 1 - Getting the ligand center:
    const Coords& fixpos = ligand.getFixpos();
    /*for (size_t i=0;i < ligand.numAtom();++i)
    {

        Box* bx = ligand.getAtom(i).getBox(&grid);
     //   cout << bx->rotpos.toString()<< " " << ligand.getAtom(i).fixpos.toString()<<"\t"<<ligand.getAtom(i).fixpos.calcDist(bx->fixpos)<<endl;
    }*/
    Box* box = grid.getBox(fixpos);
    const Coords& numCubes(box_size/(boxStep*2));


    if (box !=(Box*)NULL)
    {
        const int BxMargin=static_cast<int>(floor(8/boxStep));;
        const int maxX=(int)(box->rotpos.x+ceil(numCubes.x));//+8 for projection
        const int maxY=(int)(box->rotpos.y+ceil(numCubes.y));
        const int maxZ=(int)(box->rotpos.z+ceil(numCubes.z));
        const int minX=(int)(box->rotpos.x-ceil(numCubes.x));
        const int minY=(int)(box->rotpos.y-ceil(numCubes.y));
        const int minZ=(int)(box->rotpos.z-ceil(numCubes.z));
        const int& numPtRangeK= grid.numPtRangeK;
        const int& numPtRangeI= grid.numPtRangeI;
        const int& numPtRangeJ= grid.numPtRangeJ;

        for (int k=0; k < numPtRangeK; k++)
            for (int j=0; j <  numPtRangeJ; j++)
                for (int i=0; i <  numPtRangeI; i++)
                {

                    Box &bx= *grid.AllBoxes.at(i+j*numPtRangeI+k*numPtRangeI*numPtRangeJ);
                    bx.Props.clear();
                    if (bx.rotpos.x >= minX-BxMargin && bx.rotpos.x <=maxX+BxMargin
                            &&bx.rotpos.y >= minY-BxMargin && bx.rotpos.y <=maxY+BxMargin
                            &&bx.rotpos.z >= minZ-BxMargin && bx.rotpos.z <=maxZ+BxMargin)
                    {
                        bx.setUse(true);
                        if (bx.rotpos.x >= minX && bx.rotpos.x <=maxX
                                &&bx.rotpos.y >= minY && bx.rotpos.y <=maxY
                                &&bx.rotpos.z >= minZ && bx.rotpos.z <=maxZ)
                        {

                            bx.Props.clear();
                            bx.Props.setNoInfo(true);
                            bx.setMargin(false);
                        }
                        else
                        {
                            bx.Props.clear();
                            bx.setMargin(true);
                        }
                    }
                    else
                    {
                        bx.setUse(false);
                        // cout << bx.rotpos.toString()<<endl;
                    }


                }

        //   grid.printInFile("testSel",grid.AllBoxes,false,true);

    }
    for (ItCMole it = cp.firstMole();
         it !=cp.lastMole();
         it++)
    {
        Molecule &mole=**it;
        if (&mole==&ligand)
        {
            grid.colorCube(ligand,false);
        }
        else
        {
            grid.colorCube(mole,false);
            defProps(mole, def_alternative);
        }
    }
    // grid.colorCube(molecule,true,true);
    //grid.printInFile(ligand.getName()+"_INI",grid.AllBoxes,false,true);
}


VolSite::VolSite(Molecule &molecule,
                 Grid&grid,
                 const double &boxStep, const int proj_lens, bool def_alternative) throw(MoleExcept):
    grid(grid),
    minProj(55),
    proteinOnly(true),
    min_size(70),
    max_size(400),
    drug_value(-1),
    name(""),
    L1(0.1),L2(1),L3(1.5),L4(2),L5(2.5),L6(4),L7(8)
{
    try
    {
        grid.setProjLength(proj_lens);
        grid.createMolecularBoxes(molecule,boxStep);
        grid.colorCube(molecule,true,true);
        defProps(molecule,def_alternative);
    }
    catch (MoleExcept &e)
    {
        e.addTrace("VolSite::VolSite");
        throw;
    }


}




void VolSite::proj(const int& min, const bool &force_stop)
{
    minProj=min;
#ifdef ICHEM_DEBUG
    cout << "## Calculating projections "<<endl;
    cout << "    |-> threshold : " << minProj<<endl;
#endif
    // FETCHING DATA FROM GRID :
    const int &numPtI = grid.numPtRangeI;
    const int &numPtJ = grid.numPtRangeJ;
    const int &numPtK = grid.numPtRangeK;
    const double &boxStep= grid.boxStep;
    const int &maxCubeNum = grid.maxCubeNum;
    vector<Box*> &AllBoxes = grid.AllBoxes;
    const Coords &X8 = grid.X8;
    const Coords &vect_i=grid.vect_i;
    const Coords &vect_j=grid.vect_j;
    const Coords &vect_k=grid.vect_k;
    const double &projLength = grid.projLength;
    //  vector<Coords> listCubeProjs;


    int  proj=0;
    Box* bx_prj;
    bool calcs=false;

    const Coords vectU(boxStep/2,boxStep/2,boxStep/2);
    Coords p_(0,0,0), S1(0,0,0),tmp(0,0,0),S1xyz(0,0,0);
    float mat[9];
    float A = 0.0;
    float B = 0.0;
    float C = 0.0;
    float D = 0.0;
    float E = 0.0;
    float F = 0.0;
    float AD = 0.0;
    float BD = 0.0;

    int alpha_ = 0, beta_ = 0, gamma_ = 0;  int cube_pos=0;
    vector<vector<Coords> > CubePROJS;

    if ((size_t)maxCubeNum >= projs.max_size())
        throw MoleExcept(3010201,
                         "VolSite::proj",
                         "Unable to reserve projection vector. Too much cubes");
    projs.clear();
    projs.reserve(maxCubeNum);
    for ( int i=0; i < maxCubeNum; i++) projs.push_back(0);

    calcs=true;

    // Getting the grid center :
    cube_pos =(int)floor((double)numPtI/2)+
            (int)floor((double)numPtJ/2)*numPtI+
            (int)floor((double)numPtK/2)*numPtI*numPtJ;
    Box &bxr = *(AllBoxes.at(cube_pos));
#ifdef ICHEM_DEBUG
    cout << "REFERENCE BOX: "<< bxr.fixpos.toString()
         <<"\t"<< bxr.rotpos.toString()<<endl;
#endif
    for( alpha_ = 0; alpha_ < 360; alpha_ += R){
        for(  beta_ = 0;  beta_ < 360;  beta_ += R){
            for( gamma_ = 0; gamma_ < 360; gamma_ += R){
                vector<Coords> OneProj;
#ifdef ICHEM_DEBUG
                cout << "new block "<<endl;
#endif
                // Getting angle in rad:
                A = cos (alpha_*PI/180); B = sin (alpha_*PI/180);
                C = cos (beta_*PI/180);  D = sin (beta_*PI/180);
                E = cos (gamma_*PI/180); F = sin (gamma_*PI/180);
                AD = A*D; BD = B*D;

                // Creating rotation matrix :
                mat[0] = C*E;	      mat[1] = -C*F;	      mat[2] = -D;
                mat[3] = -BD*E + A*F; mat[4] = BD*F + A*E;    mat[5] = -B*C;
                mat[6] = AD*E + B*F;  mat[7] = -AD*F + B*E;   mat[8] = A*C;

                // The p values represents the translation in the i j k reference
                p_.x= ( vectU.x*mat[0] + vectU.y*mat[1] + vectU.z*mat[2]);
                p_.y= ( vectU.x*mat[3] + vectU.y*mat[4] + vectU.z*mat[5]);
                p_.z= ( vectU.x*mat[6] + vectU.y*mat[7] + vectU.z*mat[8]);

                // S1 will then be the value of the transfered point from the looked box, still in i j k reference
                S1.setCoords(bxr.rotpos.x + p_.x, bxr.rotpos.y+p_.y, bxr.rotpos.z+p_.z) ;
                // S1fixpos is the conversion of S1 in the x y z reference.
                //cout << S1.ToString();

                S1xyz.setCoords(X8.x+boxStep*(S1.x*vect_i.x+S1.y*vect_j.x+S1.z*vect_k.x),
                                X8.y+boxStep*(S1.x*vect_i.y+S1.y*vect_j.y+S1.z*vect_k.y),
                                X8.z+boxStep*(S1.x*vect_i.z+S1.y*vect_j.z+S1.z*vect_k.z));
                //cout << S1fixpos-bxr->pos()<<endl;
                bx_prj=(Box*)NULL;
                for ( int i=floor(S1.x); i<=ceil(S1.x);i++){if (bx_prj!=NULL)break;if (i<0)continue;
                    for ( int j=floor(S1.y); j<=ceil(S1.y);j++){if (bx_prj!=NULL)break;if (j<0)continue;
                        for ( int k=floor(S1.z); k<=ceil(S1.z);k++)
                        {
                            if (k<0)continue;
                            cube_pos =i+j*grid.numPtRangeI+k*grid.numPtRangeI*grid.numPtRangeJ;

                            if (cube_pos >=maxCubeNum || cube_pos < 0)continue;
                            Box &bx = *(AllBoxes.at(cube_pos));
                            //  cout << "testA -> "<< bx->rotpos.toString()<<"\t"<< bx->fixpos.calcDist(S1xyz)<<endl;
                            if (S1xyz.x >= bx.fixpos.x-boxStep/2 &&S1xyz.x < bx.fixpos.x+boxStep/2
                                    &&S1xyz.y >= bx.fixpos.y-boxStep/2 &&S1xyz.y < bx.fixpos.y+boxStep/2
                                    &&S1xyz.z >= bx.fixpos.z-boxStep/2 &&S1xyz.z < bx.fixpos.z+boxStep/2){ bx_prj=&bx;break;}
                        }}}

                if (bx_prj != (Box*)NULL){ OneProj.push_back(bx_prj->rotpos-bxr.rotpos);
#ifdef ICHEM_DEBUG
                    cout << bx_prj->rotpos.toString()<<"\t"<<bxr.fixpos.calcDist(S1xyz)<<"\t"<<(bx_prj->rotpos-bxr.rotpos).toString()<<endl;
#endif
                }
                //Second try on the projection vector:

                tmp = S1+p_*0.5*U;

                S1xyz.setCoords(X8.x+boxStep*(tmp.x*vect_i.x+tmp.y*vect_j.x+tmp.z*vect_k.x),
                                X8.y+boxStep*(tmp.x*vect_i.y+tmp.y*vect_j.y+tmp.z*vect_k.y),
                                X8.z+boxStep*(tmp.x*vect_i.z+tmp.y*vect_j.z+tmp.z*vect_k.z));
                while(bxr.fixpos.calcDist(S1xyz)<projLength)
                {

                    bx_prj=(Box*)NULL;
                    for ( int i=floor(tmp.x)-2; i<=ceil(tmp.x)+2;i++){if (bx_prj!=(Box*)NULL)break;
                        for ( int j=floor(tmp.y)-2; j<=ceil(tmp.y)+2;j++){if (bx_prj!=(Box*)NULL)break;
                            for ( int k=floor(tmp.z)-2; k<=ceil(tmp.z)+2;k++)
                            {
                                cube_pos =i+j*numPtI+k*numPtI*numPtJ;

                                if (cube_pos >=maxCubeNum || cube_pos< 0)continue;
                                Box &bx = *(AllBoxes.at(cube_pos));
                                //   cout << "testB -> "<< bx->rotpos.toString()<<"\t"<< bx->fixpos.calcDist(S1xyz)<<endl;
                                if (S1xyz.x >= bx.fixpos.x-boxStep/2 &&S1xyz.x < bx.fixpos.x+boxStep/2
                                        &&S1xyz.y >= bx.fixpos.y-boxStep/2 &&S1xyz.y < bx.fixpos.y+boxStep/2
                                        &&S1xyz.z >= bx.fixpos.z-boxStep/2 &&S1xyz.z < bx.fixpos.z+boxStep/2){ bx_prj=&bx;break;}
                            }}}

                    if( bx_prj != (Box*)NULL)
                    {
#ifdef ICHEM_DEBUG

                        cout << bx_prj->rotpos.toString()<<"\t"<<bxr.fixpos.calcDist(S1xyz)<<"\t"<<(bx_prj->rotpos-bxr.rotpos).toString()<<endl;
#endif
                        OneProj.push_back(bx_prj->rotpos-bxr.rotpos);
                    }
                    if (bxr.fixpos.calcDist(S1xyz)>=projLength)
                    {
                        break;
                    }
                    tmp += p_*0.5*U;
                    S1xyz.setCoords(X8.x+boxStep*(tmp.x*vect_i.x+tmp.y*vect_j.x+tmp.z*vect_k.x),
                                    X8.y+boxStep*(tmp.x*vect_i.y+tmp.y*vect_j.y+tmp.z*vect_k.y),
                                    X8.z+boxStep*(tmp.x*vect_i.z+tmp.y*vect_j.z+tmp.z*vect_k.z));
                }
                const vector<Coords>::iterator itUn = unique(OneProj.begin(), OneProj.end());
                OneProj.resize(itUn-OneProj.begin());
                //listCubeProjs.insert(listCubeProjs.end(),OneProj.begin(),OneProj.end());
                CubePROJS.push_back(OneProj);
            }
            //  break;
        }
        //break;
    }
    /*const vector<Coords>::iterator itUn = unique(listCubeProjs.begin(), listCubeProjs.end());
    listCubeProjs.resize(itUn-listCubeProjs.begin());

    Molecule moleprojs;Residu& res=moleprojs.addResidu("A",1,"ALA");
    for (vector<Coords>::const_iterator itUn=listCubeProjs.begin(); itUn!= listCubeProjs.end();++itUn)
    {
        moleprojs.addAtom("C",*itUn,"C","C.4",&res);
    }
    Atom &ini=moleprojs.addAtom("O",Coords(0,0,0),"O","O.2",&res);
    size_t former_pos;
    for (vector<vector<Coords> >::iterator ITProjs = CubePROJS.begin() ; ITProjs != CubePROJS.end(); ITProjs++)
    {

        for (size_t i=0;i< (*ITProjs).size();++i)
        {
            const vector<Coords>::iterator it=find(listCubeProjs.begin(),listCubeProjs.end(),(*ITProjs).at(i));
            if (it == listCubeProjs.end()) continue;
            const size_t pos =std::distance(listCubeProjs.begin(),it);

            Atom& atom=const_cast<Atom&>(moleprojs.getAtom(pos));
            atom.setBFactor(atom.getBFactor()+1);
            if (i==0){
                if (!atom.hasBondWith(ini)) moleprojs.addBond(ini,atom,BondType::SINGLE);
                former_pos=pos; continue;}
            Atom& atom2=const_cast<Atom&>(moleprojs.getAtom(former_pos));
            if (!atom2.hasBondWith(atom)) moleprojs.addBond(atom,atom2,BondType::SINGLE);
            former_pos=pos;

        }
    }

    AtomList todo,newtodo;
    for (ItCBox bxr = AllBoxes.begin(); bxr != AllBoxes.end(); bxr++)
      {
          proj =0;
          const Box &bx_look = **bxr;
          if (!bx_look.getUse())continue;
          //cout << bx_look.fixpos.toString()<< " " << bx_look.rotpos.toString() <<endl;
          if (bx_look.isMargin())continue;// Don't do projection value on margin
            todo.clear();
          todo.insert(todo.end(),ini.firstAtom(),ini.lastAtom());
moleprojs.setUse(true);
cout <<"#######"<< bx_look.getId()<<endl;
          while (!todo.empty())
          {
                newtodo.clear();
                cout <<"####NEW LEVEL : "<< todo.size()<<endl;
                for (size_t i=0;i< todo.size();++i)
                {
                     Atom& atm=*todo.at(i);
                    cout << "#####"<<atm.fixpos.toString()<<" " << atm.isUsed()<<endl;
                    if (!atm.isUsed())continue;
                    const Coords& at=atm.fixpos;
                    cube_pos = (bx_look.rotpos.x+at.x)+
                               (bx_look.rotpos.y+at.y)*numPtI+
                               (bx_look.rotpos.z+at.z)*numPtI*numPtJ;
                    if (cube_pos >= maxCubeNum || cube_pos < 0) continue;
                    atm.setUse(false,false,false);
                    const Box &bx_prj = *(AllBoxes.at(cube_pos));
                    //cout << "  ||--> "<< bx_prj.rotpos.toString()<<" " << bx_prj.Props.isInProtein() <<endl;
                    if (bx_prj.getUse() && bx_prj.Props.isInProtein())
                    {
                        proj+=atm.getBFactor();
                        cout << "NEW PROJ "<< proj<<endl;
                    }
                    else
                    {
                        for (size_t j=0;j< atm.getNumBond();j++)
                        {
                            Atom& atm2=atm.getBond(j)->getOtherAtom(atm);
                            if (!atm2.isUsed())continue;
                            newtodo.push_back(&atm2);
                        }
                    }
                    if(proj >= NB_PROJ) {break;}
                    if (force_stop && proj > minProj) break;
                }
                if (newtodo.empty())break;
                todo=newtodo;
          }


          Box &bx_update = **bxr;

          projs[bx_update.getId()] = proj;
          if (!(bx_update.Props.isNoInfo() || bx_update.Props.isInLigand())) continue;

          bx_update.Props.clear();
          if (proj >= minProj ) bx_update.Props.setInCavity(true);
          else                  bx_update.Props.setCloseProtein(true);

break;

    }
*/

    /*  for (vector<vector<Coords> >::iterator ITProjs = CubePROJS.begin() ; ITProjs != CubePROJS.end(); ITProjs++)
        {
            cout <<"#####"<<endl;
            for (vector<Coords>::iterator itP = (*ITProjs).begin(); itP != (*ITProjs).end(); itP++)
            {
                cout << (*itP).toString()<<endl;
              }
          }*/


    for (ItCBox bxr = AllBoxes.begin(); bxr != AllBoxes.end(); bxr++)
    {
        proj =0;
        const Box &bx_look = **bxr;
        if (!bx_look.getUse())continue;
        //cout << bx_look.fixpos.toString()<< " " << bx_look.rotpos.toString() <<endl;
        if (bx_look.isMargin())continue;// Don't do projection value on margin
        ;
        for (vector<vector<Coords> >::const_iterator ITProjs = CubePROJS.begin() ; ITProjs != CubePROJS.end(); ITProjs++)
        {
            bool found=false;

            for (vector<Coords>::const_iterator itP = (*ITProjs).begin(); itP != (*ITProjs).end(); itP++)
            {
                //Coords new_Co = bx_look.rotpos+(*itP);

                //cube_pos =new_Co.x+new_Co.y*numPtI+new_Co.z*numPtI*numPtJ;
                cube_pos = (bx_look.rotpos.x+(*itP).x)+
                        (bx_look.rotpos.y+(*itP).y)*numPtI+
                        (bx_look.rotpos.z+(*itP).z)*numPtI*numPtJ;
                if (cube_pos >= maxCubeNum || cube_pos < 0) continue;
                const Box &bx_prj = *(AllBoxes.at(cube_pos));
                if (!bx_prj.getUse())continue;
                //cout << "  ||--> "<< bx_prj.rotpos.toString()<<" " << bx_prj.Props.isInProtein() <<endl;
                if (bx_prj.Props.isInProtein()) {found=true;break;}

            }
            if (found) proj++;
            if(proj >= NB_PROJ) {break;}
            if (force_stop && proj > minProj) break;
        }
        //  cout << "  ||--> "<< proj<<endl;

        Box &bx_update = **bxr;

        projs[bx_update.getId()] = proj;
        if (!(bx_update.Props.isNoInfo() || bx_update.Props.isInLigand())) continue;

        bx_update.Props.clear();
        if (proj >= minProj ) bx_update.Props.setInCavity(true);
        else                  bx_update.Props.setCloseProtein(true);



    }

#ifdef ICHEM_DEBUG
    printProj("proj");
    grid.printInFile("proj_C",AllBoxes);
#endif
}

#ifdef ICHEM_DEBUG
void VolSite::printProj(const string& fname) const
{

    Molecule molen(MoleType::OTHER);

    const int &maxCubeNum = grid.maxCubeNum;
    const vector<Box*> &AllBoxes = grid.AllBoxes;

    string name,mol2;
    ostringstream oss;
    unsigned int pos=0;
    map<unsigned int,Residu*> residid;
    Residu* res=(Residu*)NULL;
    map<unsigned int,Residu*>::iterator residit=residid.end();
    for (unsigned int i = 0; i <= 13; i++)
    {
        oss.str("");oss<<"C"<<i;
        residid.insert(pair<unsigned int,Residu*>(i,&molen.addResidu("A",i,oss.str())));
    }

    try
    {

        ofstream foutput((fname+".csv").c_str(),ios::out);
        foutput << "ID\tX\tY\tZ\tPROJ"<<endl;
        for ( int i=0; i < maxCubeNum; i++)
        {
            const Box & bx = *(AllBoxes.at(i));
            int proj= projs[i];

            if (proj <= 10) continue;
            foutput << bx.getId()<<"\t"<<bx.fixpos.x <<"\t"<<bx.fixpos.y<<"\t"<<bx.fixpos.z<<"\t"<<proj<<endl;
            pos= floor(proj)/NB_PROJ*(NB_MOL2-1)/5;
            mol2 = MOL2_SPEC[pos].mol2type;
            name = AtomData[MOL2_SPEC[pos].atomic_num-1].nom;
            residit= residid.find(pos);
            res = (*residit).second;
            if (res == (Residu*)NULL) continue;
            //        cout << (*it).second <<"\t"<< (*it).first->getId()<<endl;


            molen.addAtom(name,bx.fixpos,name,mol2,res);

        }

        cout << "PROJECTION ATOM SCALE :: "<<endl;
        for (unsigned int Ivar1=0; Ivar1< NB_PROJ; Ivar1+=5)
        {
            pos= floor(Ivar1)/NB_PROJ*(NB_MOL2-1)/5;
            cout << Ivar1<<"-"<<Ivar1+4<< "=>"<< AtomData[MOL2_SPEC[pos].atomic_num-1].nom<<" " << MOL2_SPEC[pos].mol2type << "\t"<<(*residid.find(pos)).second->getName()<<endl;
        }

        MoleWriter writer(fname+".mol2",FileFormat::MOL2);
        writer.writeMOL2(&molen);
        foutput.close();
    }
    catch (MoleExcept &e)
    {
        cout << e.getCode()<<"\t"<<e.getData()<<endl;
    }
}
#endif



const unsigned int& VolSite::getProjValue(const Box& bx)const
{
    return projs.at(bx.getId());
}


void VolSite::treatments(unsigned int minNeighbors, const double dist_prot, bool siteout) throw(MoleExcept)
{
    if (grid.AllBoxes.size() ==0) {throw MoleExcept(1160801,"Grid::traitements","No cube found");}
    vector<Box*> &AllBoxes= grid.AllBoxes;

    // IN_SIT : With enough burying
    // SIT : Are cavities



    unsigned int count=0;
    unsigned int nb_voisins = 0;
    //    unsigned int minNeighbors=3;
    unsigned int maxNeighbors=9;

    // When wanting to keep "ligand" cubes, even without enough burying :
    if (!proteinOnly)
    {
        //        cout << "je rentre dans proteinOnly "<< endl;
        for (ItCBox itBx = AllBoxes.begin(); itBx != AllBoxes.end(); itBx++)
        {
            Box &bx = **itBx;
            if (!bx.getUse()) continue;
            if (!bx.Props.isInLigand()) {continue;}
            bx.Props.clear(); bx.Props.setInCavity(true);
            //            cout << "j'ai mis une bx dans cavity "<<endl;

        }
    }

    count=0;
    // Setting SIT to IN_SIT cubes with enough neighbors
    BoxList closeBox;
    for (ItCBox itBx = AllBoxes.begin(); itBx != AllBoxes.end(); itBx++)
    {
        Box &bx = **itBx;
        if (!bx.getUse())continue;
        nb_voisins=0;
        if (bx.Props.isInCavity()){
            grid.getAdjacency(closeBox,&bx,grid.boxStep-0.1);
            for (ItCBox itBa= closeBox.begin(); itBa != closeBox.end(); ++itBa)
            {
                if (!(*itBa)->getUse())continue;
                if ((*itBa)->Props.isInCavity()) nb_voisins++;
            }
            // Looking at each cube point around this one, closer than max_dist being considered as in_sit
            if(nb_voisins < maxNeighbors )
            {
                bx.Props.clear();
                bx.Props.setCloseProtein(true);
                count++;
            }
        }

    }
#ifdef ICHEM_DEBUG
    cout <<"## Cube Removed (not enough neighbors) : "<<count<<endl;
    grid.printInFile("TREAT_1",AllBoxes);
#endif

    // And recovering thoses with enough neighbors but not enough burying
    BoxList tmpList;
    count=0;
    for (ItCBox itBx = AllBoxes.begin(); itBx != AllBoxes.end(); itBx++)
    {
        Box &bx = **itBx;
        if (!bx.getUse())continue;
        if (bx.Props.isNoInfo() || bx.Props.isCloseProtein()){
            nb_voisins=0;
            grid.getAdjacency(closeBox,&bx,grid.boxStep-0.1);
            for (ItCBox itBa= closeBox.begin(); itBa != closeBox.end(); ++itBa)
            {
                if (!(*itBa)->getUse())continue;
                if ((*itBa)->Props.isInCavity()) nb_voisins++;
            }
            if(nb_voisins >= 9 ){ tmpList.push_back(&bx); count++;}
        }
    }
    for (ItBox bx = tmpList.begin(); bx != tmpList.end(); bx++)
    {
        (*bx)->Props.clear(); (*bx)->Props.setInCavity(true);
    }

#ifdef ICHEM_DEBUG
    cout <<"## Cube recovered by their neighbors : "<<count<<endl;
    grid.printInFile("TREAT_2",AllBoxes);
#endif

    double  bestdist=0;
    Atom *bestatm=(Atom*)NULL;
    AtomList listAtoms;
    cubes_cavite.clear();

    const double thres_dProt=dist_prot; //%toto

    for (ItCBox itBx = AllBoxes.begin(); itBx != AllBoxes.end(); itBx++)
    {
        Box &bx = **itBx;
        if (!bx.getUse())continue;
        if (!bx.Props.isInCavity())continue;
        grid.getAdjacentAtoms(listAtoms,bx,2*thres_dProt);
        // in old version bestdist was not reset to 1000. So if we don't find a bestatm bestdist is the same that for the previous box.
        bestdist=1000;
        bestatm=(Atom*)NULL;

        for (ItCAtom itAtm= listAtoms.begin(); itAtm != listAtoms.end(); ++itAtm)
        {
            Atom &atm = **itAtm;
            if (atm.isHydrogen()  ||atm.getMOL2Type()=="Du") continue;
            Molecule* const parent= &atm.getParent();
            if (parent != (Molecule*)NULL)
            {
                if (proteinOnly && atm.getParent().getMoleType() != MoleType::PROTEIN)continue;
                if (!proteinOnly && (atm.getParent().getMoleType() != MoleType::PROTEIN &&
                                     atm.getParent().getMoleType() == MoleType::LIGAND))continue;

            }

            const double dist=atm.fixpos.calcDist(bx.fixpos);
            if (dist > bestdist) continue;
            bestdist=dist;
            bestatm=&atm;
        }

        // write a file with the site
        if (siteout){
            ofstream siteout("Usefull_Site.cav", ios::app);
            if (bestatm != (Atom*)NULL) siteout << bestatm->getIdentifier()<<endl;
            //                cout << "best dist : " << bestdist << endl;
            siteout.close();
        }

        if (bestdist < thres_dProt || bestdist==1000)
        {
#ifdef ICHEM_DEBUG
            cout << bx.getId() << " too close" <<endl;
#endif
            bx.Props.clear(); bx.Props.setOutCavity(true);

        }
        else if (bestdist >=thres_dProt) {cubes_cavite.push_back(&bx);}
    }
#ifdef ICHEM_DEBUG
    grid.printInFile("TREAT_3",AllBoxes);
#endif

    // Deleting some small cavities that might be existing around the main cavity.
    // Theses empty zones are juste spaces between residues which are not big enough to put two atoms in it.
    for (ItBox itBx = cubes_cavite.begin(); itBx != cubes_cavite.end(); itBx++)
    {
        const Box &bx = **itBx;
        nb_voisins=0;
        for (ItCBox itBv =  cubes_cavite.begin(); itBv != cubes_cavite.end(); itBv++)
            if (	abs((int) (bx.rotpos.x-(*itBv)->rotpos.x)) <= 1
                    && abs((int)(bx.rotpos.y-(*itBv)->rotpos.y)) <= 1
                    && abs((int)(bx.rotpos.z-(*itBv)->rotpos.z))<=1
                    && *itBv != &bx) nb_voisins++;
        if( (nb_voisins < minNeighbors ) || (nb_voisins < 9 && projs[bx.getId()] < 90) ){
#ifdef ICHEM_DEBUG
            cout << "not enough neighbors : " << bx.getId()<<"::"<<bx.rotpos.x<<"-"<<bx.rotpos.y<<"-"<<bx.rotpos.z<<endl;
#endif
            ItBox itBD = itBx;
            itBx--;
            cubes_cavite.erase(itBD);}

    }


#ifdef ICHEM_DEBUG

    grid.printInFile("fcav_prop",AllBoxes);
#endif

}




struct Cluster
{
    set<const Box*> listBox;
    unsigned int idClus;

};


void VolSite::clusterPHA(int i, int j, int k){
    Box *bx=grid.getBox(i,j,k);
    if (bx->isUsed()) {
        return;
    }else {
        if (bx->Props.isInCavity() || bx->Props.isInLigand()){
            bx->setphaUsed(true);
            setPhaSuccess(true);
            clusterPHA(i+1,j,k);
            clusterPHA(i,j+1,k);
            clusterPHA(i,j,k+1);
            clusterPHA(i-1,j,k);
            clusterPHA(i,j-1,k);
            clusterPHA(i,j,k-1);
        }
    }
    //            cout << "no cavity" << endl;
}

double randf(double m)
{
    return m * rand() / (RAND_MAX - 1.);
}

void VolSite::KmeansPHA(int nb_clust,bool prop){
    if (cubes_cavite.size() == 0)throw MoleExcept(1160901,"Grid::cavProp","No cube found");
#ifdef ICHEM_DEBUG
    cout << "###################" << endl <<
            "Kmeans Pha begin" << endl <<
            "###################" << endl;
    cout << "## Cavity cubes = " << cubes_cavite.size() << endl;
    MoleWriter mw;
#endif



    BoxList &AllBoxes = grid.AllBoxes;
    int changed =0;

    for (int i=0; i<nb_clust ; i++){
        int alea = rand()% cubes_cavite.size();
        ClusterPha cluster; cluster.idClus=0;
        do {
            alea = rand()% cubes_cavite.size();
            cluster.Centr.setCoords(cubes_cavite[alea]->fixpos);

        } while (!cubes_cavite[alea]->isUsed());
        clusterList.push_back(cluster);
    }
    if (prop) {
        for (int i=0; i < nb_clust; i++){ //mettre une variable au nombre de kmeans
            clusterList[i].pha = 2;
        }

    }
    int pos =0;
    //initialisation des cluster kmeans
    for (ItCBox itpos = cubes_cavite.begin(); itpos != cubes_cavite.end(); itpos++)
    {

        Box &bx=**itpos;
        if (!bx.getUse() ){ continue;}
        int rand = pos % nb_clust;
        bx.setClust(rand);
        pos++;

        double best_dist = 10000, temp_dist = 10000;
        int best_clus = 1000, temp_clust = 0;
        for (int i=0; i<nb_clust ; i++,temp_clust++){
            clusterList[i].idClus = temp_clust;
            temp_dist = bx.fixpos.calcDist(clusterList[i].Centr);
            if (temp_dist < best_dist) {
                best_dist = temp_dist;
                best_clus = temp_clust;
            }
        }
        if (best_clus != bx.getClust()) {
            changed++;
            bx.setClust(best_clus);
            cout << best_clus << endl;
            //                cout << "transfert vers cluster n "<< best_clus << endl;
        }

    }
    for (int i=0; i<nb_clust ; i++){
        clusterList[i].Centr.setCoords(0,0,0);
        clusterList[i].idClus = 0;
        clusterList[i].pha = i%3;
    }
    do {
        /* group element for centroids are used as counters */
        for (int i=0; i<nb_clust ; i++){
            clusterList[i].Centr.setCoords(0,0,0);
            clusterList[i].idClus = 0;
        }

        for (ItCBox itpos = cubes_cavite.begin(); itpos != cubes_cavite.end(); itpos++)
        {

            const Box &bx=**itpos;
            if (!bx.getUse()) continue;
            Coords transit = clusterList[bx.getClust()].Centr + bx.fixpos;
            clusterList[bx.getClust()].Centr.setCoords(transit);
            clusterList[bx.getClust()].idClus++;

        }//cout << endl<< "Moyenne des position" << endl;
        for (int i=0; i<nb_clust ; i++){
            if (clusterList[i].Centr == Coords(0,0,0)){ continue;}
            //                        cout << clusterList[i].Centr.toString() << "  //  ";
            clusterList[i].Centr.setCoords(clusterList[i].Centr/clusterList[i].idClus);
            //                        cout << clusterList[i].Centr.toString() << endl;

        }
        // moyenne de toutes les coordonées

        changed = 0;
        //        /* find closest centroid of each point */
        for (ItCBox itpos = cubes_cavite.begin(); itpos != cubes_cavite.end(); itpos++)
        {
            Box &bx=**itpos;
            if (!bx.getUse()) continue;
            // calcul le cluster le plus proche
            double best_dist = 10000, temp_dist = 10000;
            int best_clus = 1000, temp_clust = 0;
            for (int i=0; i<nb_clust ; i++,temp_clust++){
                clusterList[i].idClus = temp_clust;
                temp_dist = bx.fixpos.calcDist(clusterList[i].Centr);
                if (temp_dist < best_dist) {
                    best_dist = temp_dist;
                    best_clus = temp_clust;
                }
            }
            if (best_clus != bx.getClust()) {
                changed++;
                bx.setClust(best_clus);
                //                cout << "transfert vers cluster n "<< best_clus << endl;
            }
            //        int toto =
        }
        //        cout << changed << " > " << (cubes_cavite.size() >> 10) << endl;

        //        cout << endl;
    } while (changed > (cubes_cavite.size() >> 10)); /* stop when 99.9% of points are good */



}
void VolSite::clusterCAV(const int &min_proj)
{

    if (cubes_cavite.size() == 0)throw MoleExcept(1160901,"Grid::cavProp","No cube found");
#ifdef ICHEM_DEBUG
    cout << "###################" << endl <<
            "Cluster CAV begin" << endl <<
            "###################" << endl;
    cout << "## Cavity cubes = " << cubes_cavite.size() << endl;
    MoleWriter mw;
#endif

    const size_t max_size=2500;

    bool Checked[max_size][max_size];
    // Sort cubes by decreasing projection values :
    multimap< int,const Box*> boxDecreasProj;
    multimap< int,const Box*>::reverse_iterator ITBDP;
    // Set to know how many adjacent cubes of 1 cube can be within a cluster i.e. not protein or out
    vector<double> NSitCount(grid.AllBoxes.size(),0);
    vector<int> boxAssocClus(grid.AllBoxes.size(),-1);

    double MaxPrj=0;
    const double THRES=40;
    const double distThres=grid.boxStep+1.1;
    unsigned int num_cav=0;
    bool modif;
    unsigned int RId=0,CId=0;
    double RCDist=0,RCCompact=0, NAdj=0,NClose=0, NMerge=0,dist=0,SMDist=0,RAProj=0,CAProj=0,NNProjR=0,NNProjC=0,BestRCC=0,RAprojmax=0,CAprojmax=0;
    const double diagBxstep = sqrt(2)*grid.boxStep+0.1;
    BoxList adjacentBoxes;
    const double& boxStep = grid.boxStep;
    //
    vector<Cluster> clusterList;
    vector<Cluster>::iterator BestCl;
    BoxList &AllBoxes = grid.AllBoxes;




    // Initialization :
    for (size_t I=0; I<max_size; I++)
        for (size_t J=0; J<max_size; J++) Checked[I][J]=false;


    for (ItCBox itpos = cubes_cavite.begin(); itpos != cubes_cavite.end(); itpos++)
    {
        const Box &bx=**itpos;
        boxDecreasProj.insert(pair<unsigned int,const Box* >
                              (projs[bx.getId()],&bx));
        NAdj=0;
        grid.getAdjacency(adjacentBoxes,&bx,grid.boxStep-0.1);
        //cout << "#######"<<(*itpos)->getId()<< "\t "<< adjacentBoxes.size()<<endl;
        for (ItCBox itbb = adjacentBoxes.begin();
             itbb !=adjacentBoxes.end();
             itbb++)
        {
            //  cout << (*itbb)->getId()<<"\t"<<(*itbb)->Props.toCavString();
            // Virer isVid pour créer deséquilibre coté outside cav / inside cav
            if ((*itbb)->Props.isNoInfo() ||
                    (*itbb)->Props.isCloseProtein() ||
                    (*itbb)->Props.isInCavity())
            {
                //  cout << "\tSEL";
                NAdj++;
            }
            //cout <<endl;
        }
        // cout << (*itpos)->getId()<< " " << projs[(*itpos)->getId()]<< " " << NAdj<<endl;
        NSitCount.at((*itpos)->getId())=NAdj;

    }




    // Putting the most buried box as the first cluster
    ITBDP = boxDecreasProj.rbegin();
    Cluster cluster; cluster.idClus=0;
    cluster.listBox.insert((*ITBDP).second);
    clusterList.push_back(cluster);
    boxAssocClus.at((*ITBDP).second->getId())=0;
    // Setting the maximum projection value:
    MaxPrj = (*ITBDP).first;
    MaxPrj = (int)ceil((double)MaxPrj/5.0)*5.0;

    // Deleting this most buried box from the one we treat:
    boxDecreasProj.erase(--ITBDP.base());



    for (; MaxPrj >=THRES; MaxPrj-=5)
    {
#ifdef ICHEM_DEBUG
        cout << "## PROJ MAX : "<< MaxPrj <<"\tPROJ MIN "<< (MaxPrj-5)<<endl;
#endif
        //if (MaxPrj <= 100)break;
        // No other box to treat => we stop
        if (boxDecreasProj.empty())break;

        // For the next step, we need to know if adjacent boxes of
        // our boxes are within a cluster(s) or not
        // For this, we need a table to associate the count of adj box within cluster
        // and of the size of the number of cavities :
        vector<unsigned int> assocBoxCls(clusterList.size(),0);
        //for (size_t i=0; i < clusterList.size();++i) assocBoxCls.push_back(0);

        // We fetch all cubes  that are in the range [MaxPrj;MaxPrj-5[
        for (ITBDP = boxDecreasProj.rbegin();ITBDP != boxDecreasProj.rend(); ITBDP++) //nop
        {

            if ((*ITBDP).first < MaxPrj-5) break;
            const Box &bx = *(*ITBDP).second;

            // We initialize the associated vector :
            for (size_t i=0; i < clusterList.size();++i) assocBoxCls.at(i)=0;

            // We scan each adjacent box :
            grid.getAdjacency(adjacentBoxes,&bx,grid.boxStep-0.1);

            for (vector<Box*>::const_iterator bxL =adjacentBoxes.begin(); bxL != adjacentBoxes.end(); bxL++)
            {
                // And force that theses adjacent box are not on the huge diagonal
                // i.e we keep thoses with at least 1 face or 1 edge common
                // with bx, but not adjacent with only one point
                //
                //cout << (bx.rotpos-(*bxL)->rotpos).toString()<<"\t"<< (fabs((bx.rotpos-(*bxL)->rotpos).x)+fabs((bx.rotpos-(*bxL)->rotpos).y)+fabs((bx.rotpos-(*bxL)->rotpos).z))<<"\t"<<bx.fixpos.calcDist((*bxL)->fixpos,distThres+0.1)<<endl;
                if (bx.fixpos.calcDist((*bxL)->fixpos,distThres+0.1) >distThres)continue;

                // We search for a cluster associated to this adj box :
                int &bxL_clsID=boxAssocClus.at((*bxL)->getId());


                if (bxL_clsID== -1)continue;

                // Now we search if this cluster has already been
                assocBoxCls.at(bxL_clsID)++;

            }
            unsigned int MaxClusADJ=0;
            int bestCluster=-1;
            // Now we search for the closest cluster :
            for (size_t i=0; i < clusterList.size();++i)
            {

                if (assocBoxCls.at(i)==0)continue;
                if (assocBoxCls.at(i) > MaxClusADJ)
                {
                    bestCluster=i;
                    MaxClusADJ  = assocBoxCls.at(i);
                }

            }
            if (( MaxClusADJ >= 11 && min_proj > MaxPrj-1)  ||  (MaxClusADJ > 2 && min_proj < MaxPrj)  ||  (MaxClusADJ > 0 && MaxPrj > 75))
            {
                clusterList.at(bestCluster).listBox.insert(&bx);
                boxAssocClus.at(bx.getId())=bestCluster;
#ifdef ICHEM_DEBUG
                cout << bx.getId()<<"\t"<< bx.fixpos.toString()<<"\t"<< (*ITBDP).first<<"\t"<<MaxClusADJ<<"\tFOUND : "<< bestCluster<<endl;
#endif

            }
            // If we found one, we add it
            else if (min_proj < MaxPrj)
            {
                num_cav++;
                Cluster cluster; cluster.idClus=num_cav;
                cluster.listBox.insert(&bx);
                clusterList.push_back(cluster);
                boxAssocClus.at(bx.getId())=clusterList.size()-1;
                assocBoxCls.push_back(0);
#ifdef ICHEM_DEBUG
                cout << bx.getId()<<"\t"<< bx.fixpos.toString()<<"\t"<< (*ITBDP).first<<"\t"<<MaxClusADJ<<"\tNEW : "<< num_cav<<endl;
#endif

            }
            boxDecreasProj.erase(--ITBDP.base());
            if (boxDecreasProj.empty())break;
            ITBDP = boxDecreasProj.rbegin();



        }//ITBDP = boxDecreasProj.rbegin()

        if (num_cav >= max_size) {cerr << "Too many clusters"<<endl; return;}
        // Cleaning booleans :

        for (size_t I =0; I <= num_cav; I++)
            for (size_t J =0; J <= num_cav; J++)Checked[I][J]=false;


        Coords RACent;bool tomerge;
#ifdef ICHEM_DEBUG
        cout << "Number of clusters : "<< num_cav<< " " << clusterList.size()<<endl;
#endif
        do
        {
            modif=false;
            for (vector<Cluster>::iterator ita= clusterList.begin();
                 ita!=clusterList.end();
                 ++ita)
            {

                RId = (*ita).idClus;
                const set<const Box*> &RCLT = (*ita).listBox;
                BestRCC=0;
                BestCl=clusterList.end();
                vector<Cluster>::iterator itb(ita);++itb;
                for (; itb != clusterList.end();itb++)
                {

                    CId = (*itb).idClus;

                    if (RId >= CId || Checked[RId][CId]==true) continue;

                    const set<const Box*> &CCLT = (*itb).listBox;
                    Checked[RId][CId]=true;
                    Checked[CId][RId]=true;
                    RCDist=0;RCCompact=0; NAdj=0;NClose=0; NMerge=0;dist=0;SMDist=100;RAProj=0;CAProj=0;NNProjR=0;NNProjC=0;RAprojmax=0;CAprojmax=0;
                    RACent.clear();
                    if (RCLT.size() < CCLT.size())
                    {
                        for (set<const Box*>::const_iterator itba = RCLT.begin();itba != RCLT.end();itba++){NAdj=0;
                            for (set<const Box*>::const_iterator itbb = CCLT.begin();itbb != CCLT.end();itbb++)
                            {
                                if (projs[(*itba)->getId()]> RAprojmax) RAprojmax=projs[(*itba)->getId()];
                                if (projs[(*itbb)->getId()]> CAprojmax) CAprojmax=projs[(*itbb)->getId()];
                                dist=(*itba)->fixpos.calcDist((*itbb)->fixpos,SMDist+0.1);
                                if (dist < SMDist){SMDist=dist;}
                                if ((*itba)->isAdjacent((*itbb))) {
                                    NAdj++;
                                    if (dist <=diagBxstep+0.1)
                                    {
                                        RCCompact+=1/NSitCount.at((*itba)->getId())+1/NSitCount.at((*itbb)->getId());
#ifdef ICHEM_DEBUG
                                        cout <<"A::"<< RCCompact
                                            << "\t"<<(*itba)->getId()
                                            <<"\t"<<NSitCount.at((*itba)->getId())
                                           <<"\t"<<(*itbb)->getId()
                                          <<"\t"<<NSitCount.at((*itbb)->getId())<<endl;
#endif
                                    }
                                    RAProj+=projs[(*itba)->getId()];NNProjR++;
                                    CAProj+=projs[(*itbb)->getId()];NNProjC++;
                                    RACent+=(*itba)->fixpos+(*itbb)->fixpos;
                                }

                            }

                            if (NAdj > 0) NMerge++;
                            if (NAdj >= 2)NClose++;
                        }
                    }
                    else
                    {
                        for (set<const Box*>::const_iterator itba = CCLT.begin();itba != CCLT.end();itba++){NAdj=0;
                            for (set<const Box*>::const_iterator itbb = RCLT.begin();itbb !=RCLT.end();itbb++)

                            {
                                if (projs[(*itba)->getId()]> RAprojmax) RAprojmax=projs[(*itba)->getId()];
                                if (projs[(*itbb)->getId()]> CAprojmax) CAprojmax=projs[(*itbb)->getId()];
                                dist=(*itba)->fixpos.calcDist((*itbb)->fixpos,SMDist+0.1);
                                if (dist < SMDist){SMDist=dist;}

                                if ((*itba)->isAdjacent((*itbb))) {
                                    NAdj++;
                                    RAProj+=projs[(*itba)->getId()];
                                    if (dist <= diagBxstep+0.1)
                                    {

                                        RCCompact+=1/NSitCount.at((*itba)->getId())+1/NSitCount.at((*itbb)->getId());
#ifdef ICHEM_DEBUG
                                        cout <<"B::"<< RCCompact
                                            << "\t"<<(*itba)->getId()
                                            <<"\t"<<NSitCount.at((*itba)->getId())
                                           <<"\t"<<(*itbb)->getId()
                                          <<"\t"<<NSitCount.at((*itbb)->getId())<<endl;
#endif
                                    }
                                    NNProjR++;CAProj+=projs[(*itbb)->getId()];NNProjC++;
                                    RACent+=(*itba)->fixpos+(*itbb)->fixpos;
                                }

                            }

                            if (NAdj > 0) NMerge++;
                            if (NAdj >= 2)NClose++;

                        }

                    }
                    if (NNProjR == 0) RAProj=0; else RAProj/=NNProjR;
                    if (NNProjC == 0) CAProj=0; else CAProj/=NNProjC;
                    tomerge=false;
                    //    cout << "ST : "<< L1<< " " << L2 << " " << L3 << " " << L4 << " " << L5 << " " << L6<< L7 <<endl;
                    //L1=0.1;L2=1;L3=1.5;L4=2;L5=2.5;L6=4.5,L7=8;
                    if((RCLT.size() == 1 || CCLT.size()==1) && SMDist == boxStep) tomerge=true;
                    else if(((RCLT.size() ==1  && CCLT.size() <= 5 ) ||(RCLT.size() <=5 && CCLT.size()   == 1)) && SMDist <= diagBxstep  && RCCompact > L1) tomerge=true;
                    else if(((RCLT.size() <=15 && CCLT.size() <= 5 ) ||(RCLT.size() <=5 && CCLT.size()  <= 15)) && SMDist <= boxStep+0.1 && RCCompact > L2) tomerge=true;
                    else if(((RCLT.size() <=15 && CCLT.size() <= 15) ||(RCLT.size() <=15 && CCLT.size() <= 15)) && SMDist <= boxStep+0.1 && RCCompact > L3) tomerge=true;
                    else if(((RCLT.size() <=30 && CCLT.size() <= 15) ||(RCLT.size() <=15 && CCLT.size() <= 30)) && SMDist <= boxStep+0.1 && RCCompact > L4) tomerge=true;
                    else if(((RCLT.size() <=30 && CCLT.size() <= 30) ||(RCLT.size() <=30 && CCLT.size() <= 30)) && SMDist <= boxStep+0.1 && RCCompact > L5) tomerge=true;
                    else if ((RCLT.size() <=70 || CCLT.size() <=70) && SMDist <= boxStep+0.1 && RCCompact > L6) tomerge=true;
                    else if (RCCompact >L7) tomerge=true;

                    if (tomerge && RCCompact > BestRCC)
                    {
                        BestRCC = RCCompact;
                        BestCl=itb;
                    }
#ifdef ICHEM_DEBUG
                    if (SMDist <= 5){
                        std::cout.precision(3);
                        cout <<"CLT\t"<<((tomerge==true) ? "1":"0")<<"\t"
                            <<RId<<"::"<<CId<<"\t"
                           <<RCLT.size()<<"::"<<CCLT.size()<<"\t"
                          <<NClose<<"::"<<NMerge<<"\t"
                         <<RCDist<<"::"<<SMDist<<"\t"
                        <<RAProj<<"::"<<CAProj<<"\t"
                        << RCCompact<<"\t\t" << RAprojmax<<"::"<<CAprojmax
                        <<endl;
                    }
#endif
                }
                if (BestCl != clusterList.end())
                {
                    modif=true;
                    CId = (*BestCl).idClus;
                    const set<const Box*> &CCLT=(*BestCl).listBox;
#ifdef ICHEM_DEBUG
                    cout << "MERGE "<<CId << " in " << RId << endl;
#endif
                    for (unsigned int I =0; I < num_cav; I++){ Checked[I][RId]=false; Checked[RId][I]=false;Checked[I][CId]=false; Checked[CId][I]=false;}

                    for (set<const Box*>::const_iterator itbb = CCLT.begin();itbb != CCLT.end();itbb++)
                    {
                        (*ita).listBox.insert(*itbb);
                        boxAssocClus.at((*itbb)->getId()) = distance(clusterList.begin(),ita);
                    }
                    const int di=std::distance(clusterList.begin(),BestCl);
                    clusterList.erase(BestCl);
                    ita=clusterList.begin();

                    for (size_t i=0; i < AllBoxes.size();i++)
                    {
                        if (boxAssocClus.at(i)>di) boxAssocClus.at(i)--;

                    }

                }
            }
        } while(modif);
#ifdef ICHEM_DEBUG
        for (vector<Cluster>::iterator ita= clusterList.begin();
             ita!=clusterList.end();
             ++ita)
        {
            if (ita==clusterList.end())break;
            set<const Box*> &RCLT = (*ita).listBox;
            //cout << "DEV : "<< (*ita).first<<" " <<devproj<<endl;
            //if (devproj < 10 || RCLT->Size < 20)continue;
            Molecule *ml = new Molecule();
            for (set<const Box*>::iterator itba = RCLT.begin();itba !=RCLT.end();itba++)
            {
                ml->addAtom("C",(*itba)->fixpos,"CA","C.3");
            }

            ostringstream oss; oss << "F_"<<MaxPrj<<"_"<<(*ita).idClus<<".mol2";
            mw.newFile(oss.str());
            mw.writeMOL2(ml);
            delete ml;

        }
#endif


    }// END MaxPrj >=THRES; MaxPrj-=5

    // SELECTING CAVITIES :
    num_cav=0;
#ifdef ICHEM_DEBUG
    cout << "NB CLUSTERS : "<< clusterList.size()<<endl;
#endif
    BoxList closeBox;
    for (vector<Cluster>::const_iterator ita= clusterList.begin();
         ita!=clusterList.end();
         ++ita)
    {
        if (ita==clusterList.end())break;
        const set<const Box*> &RCLT = (*ita).listBox;
        double density=0, N_L=0;
        for (set<const Box*>::iterator itba = RCLT.begin();itba !=RCLT.end();itba++)
        {
            N_L=0;
            grid.getAdjacency(closeBox,*itba,grid.boxStep-0.1);
            for (ItCBox itbb = closeBox.begin(); itbb != closeBox.end(); itbb++)
            {
                if (RCLT.find(const_cast<Box*>(*itbb)) != RCLT.end()) N_L++;
            }
            density += N_L/NSitCount.at((*itba)->getId());
        }
        density/=(double)RCLT.size();
#ifdef ICHEM_DEBUG
        cout << (*ita).idClus <<"\t"<< RCLT.size() <<"\t"<<density<<endl;
#endif
        //        if (density > 0.80){
        //            cout << "size & density :"<< RCLT.size() << " " << density << endl;
        //        }
        if (density < 0.5 || (RCLT.size() <  min_size && density < 0.90) || RCLT.size() < 50 ) continue;
        if (RCLT.size() >  max_size) continue;
#ifdef ICHEM_DEBUG
        cout << "KEPT"<<endl;
#endif
        cavities.insert(pair<unsigned int,set<const Box*> >(num_cav,RCLT));
        //Work here

#ifdef ICHEM_DEBUG
        Molecule *ml = new Molecule();
        for (set<const Box*>::iterator itba = RCLT.begin();itba !=RCLT.end();itba++)
        {
            ml->addAtom("C",(*itba)->fixpos,"CA","C.3");
        }
        ostringstream oss; oss <<"KEPT_"<<(*ita).idClus<<"_"<<num_cav<<".mol2";
        mw.newFile(oss.str());
        mw.writeMOL2(ml);
        delete ml;
#endif
        num_cav++;
    }
    numCav = num_cav;
#ifdef ICHEM_DEBUG
    cout << "NB CAV : " << num_cav << endl;
#endif

}

//void VolSite::clusterPHA(const int &min_proj)
//{

//    if (cubes_cavite.size() == 0)throw MoleExcept(1160901,"Grid::cavProp","No cube found");
//#ifdef ICHEM_DEBUG
//    cout << "###################" << endl <<
//            "Cluster CAV begin" << endl <<
//            "###################" << endl;
//    cout << "## Cavity cubes = " << cubes_cavite.size() << endl;
//    MoleWriter mw;
//#endif

//    const size_t max_size=2500;

//    bool Checked[max_size][max_size];
//    // Sort cubes by decreasing projection values :
//    multimap< int,const Box*> boxDecreasProj;
//    multimap< int,const Box*>::reverse_iterator ITBDP;
//    // Set to know how many adjacent cubes of 1 cube can be within a cluster i.e. not protein or out
//    vector<double> NSitCount(grid.AllBoxes.size(),0);
//    vector<int> boxAssocClus(grid.AllBoxes.size(),-1);

//    double MaxPrj=0;
//    const double THRES=40;
//    const double distThres=grid.boxStep+1.1;
//    unsigned int num_cav=0;
//    bool modif;
//    unsigned int RId=0,CId=0;
//    double RCDist=0,RCCompact=0, NAdj=0,NClose=0, NMerge=0,dist=0,SMDist=0,RAProj=0,CAProj=0,NNProjR=0,NNProjC=0,BestRCC=0,RAprojmax=0,CAprojmax=0;
//    const double diagBxstep = sqrt(2)*grid.boxStep+0.1;
//    BoxList adjacentBoxes;
//    const double& boxStep = grid.boxStep;
//    //
//    vector<Cluster> clusterList;
//    vector<Cluster>::iterator BestCl;
//    BoxList &AllBoxes = grid.AllBoxes;




//    // Initialization :
//    for (size_t I=0; I<max_size; I++)
//        for (size_t J=0; J<max_size; J++) Checked[I][J]=false;


//    for (ItCBox itpos = cubes_cavite.begin(); itpos != cubes_cavite.end(); itpos++)
//    {
//        const Box &bx=**itpos;
//        boxDecreasProj.insert(pair<unsigned int,const Box* >
//                              (projs[bx.getId()],&bx));
//        NAdj=0;
//        grid.getAdjacency(adjacentBoxes,&bx,grid.boxStep-0.1);
//        //cout << "#######"<<(*itpos)->getId()<< "\t "<< adjacentBoxes.size()<<endl;
//        for (ItCBox itbb = adjacentBoxes.begin();
//             itbb !=adjacentBoxes.end();
//             itbb++)
//        {
//            //  cout << (*itbb)->getId()<<"\t"<<(*itbb)->Props.toCavString();
//            // Virer isVid pour créer deséquilibre coté outside cav / inside cav
//            if ((*itbb)->Props.isNoInfo() ||
//                    (*itbb)->Props.isCloseProtein() ||
//                    (*itbb)->Props.isInCavity())
//            {
//                //  cout << "\tSEL";
//                NAdj++;
//            }
//            //cout <<endl;
//        }
//        // cout << (*itpos)->getId()<< " " << projs[(*itpos)->getId()]<< " " << NAdj<<endl;
//        NSitCount.at((*itpos)->getId())=NAdj;

//    }




//    // Putting the most buried box as the first cluster
//    ITBDP = boxDecreasProj.rbegin();
//    Cluster cluster; cluster.idClus=0;
//    cluster.listBox.insert((*ITBDP).second);
//    clusterList.push_back(cluster);
//    boxAssocClus.at((*ITBDP).second->getId())=0;
//    // Setting the maximum projection value:
//    MaxPrj = (*ITBDP).first;
//    MaxPrj = (int)ceil((double)MaxPrj/5.0)*5.0;

//    // Deleting this most buried box from the one we treat:
//    boxDecreasProj.erase(--ITBDP.base());



//    for (; MaxPrj >=THRES; MaxPrj-=5)
//    {
//#ifdef ICHEM_DEBUG
//        cout << "## PROJ MAX : "<< MaxPrj <<"\tPROJ MIN "<< (MaxPrj-5)<<endl;
//#endif
//        //if (MaxPrj <= 100)break;
//        // No other box to treat => we stop
//        if (boxDecreasProj.empty())break;

//        // For the next step, we need to know if adjacent boxes of
//        // our boxes are within a cluster(s) or not
//        // For this, we need a table to associate the count of adj box within cluster
//        // and of the size of the number of cavities :
//        vector<unsigned int> assocBoxCls(clusterList.size(),0);
//        //for (size_t i=0; i < clusterList.size();++i) assocBoxCls.push_back(0);

//        // We fetch all cubes  that are in the range [MaxPrj;MaxPrj-5[
//        for (ITBDP = boxDecreasProj.rbegin();ITBDP != boxDecreasProj.rend(); ITBDP++) //nop
//        {

//            if ((*ITBDP).first < MaxPrj-5) break;
//            const Box &bx = *(*ITBDP).second;

//            // We initialize the associated vector :
//            for (size_t i=0; i < clusterList.size();++i) assocBoxCls.at(i)=0;

//            // We scan each adjacent box :
//            grid.getAdjacency(adjacentBoxes,&bx,grid.boxStep-0.1);

//            for (vector<Box*>::const_iterator bxL =adjacentBoxes.begin(); bxL != adjacentBoxes.end(); bxL++)
//            {
//                // And force that theses adjacent box are not on the huge diagonal
//                // i.e we keep thoses with at least 1 face or 1 edge common
//                // with bx, but not adjacent with only one point
//                //
//                //cout << (bx.rotpos-(*bxL)->rotpos).toString()<<"\t"<< (fabs((bx.rotpos-(*bxL)->rotpos).x)+fabs((bx.rotpos-(*bxL)->rotpos).y)+fabs((bx.rotpos-(*bxL)->rotpos).z))<<"\t"<<bx.fixpos.calcDist((*bxL)->fixpos,distThres+0.1)<<endl;
//                if (bx.fixpos.calcDist((*bxL)->fixpos,distThres+0.1) >distThres)continue;

//                // We search for a cluster associated to this adj box :
//                int &bxL_clsID=boxAssocClus.at((*bxL)->getId());


//                if (bxL_clsID== -1)continue;

//                // Now we search if this cluster has already been
//                assocBoxCls.at(bxL_clsID)++;

//            }
//            unsigned int MaxClusADJ=0;
//            int bestCluster=-1;
//            // Now we search for the closest cluster :
//            for (size_t i=0; i < clusterList.size();++i)
//            {

//                if (assocBoxCls.at(i)==0)continue;
//                if (assocBoxCls.at(i) > MaxClusADJ)
//                {
//                    bestCluster=i;
//                    MaxClusADJ  = assocBoxCls.at(i);
//                }

//            }
//            if (( MaxClusADJ >= 11 && min_proj > MaxPrj-1)  ||  (MaxClusADJ > 2 && min_proj < MaxPrj)  ||  (MaxClusADJ > 0 && MaxPrj > 75))
//            {
//                clusterList.at(bestCluster).listBox.insert(&bx);
//                boxAssocClus.at(bx.getId())=bestCluster;
//#ifdef ICHEM_DEBUG
//                cout << bx.getId()<<"\t"<< bx.fixpos.toString()<<"\t"<< (*ITBDP).first<<"\t"<<MaxClusADJ<<"\tFOUND : "<< bestCluster<<endl;
//#endif

//            }
//            // If we found one, we add it
//            else if (min_proj < MaxPrj)
//            {
//                num_cav++;
//                Cluster cluster; cluster.idClus=num_cav;
//                cluster.listBox.insert(&bx);
//                clusterList.push_back(cluster);
//                boxAssocClus.at(bx.getId())=clusterList.size()-1;
//                assocBoxCls.push_back(0);
//#ifdef ICHEM_DEBUG
//                cout << bx.getId()<<"\t"<< bx.fixpos.toString()<<"\t"<< (*ITBDP).first<<"\t"<<MaxClusADJ<<"\tNEW : "<< num_cav<<endl;
//#endif

//            }
//            boxDecreasProj.erase(--ITBDP.base());
//            if (boxDecreasProj.empty())break;
//            ITBDP = boxDecreasProj.rbegin();



//        }//ITBDP = boxDecreasProj.rbegin()

//        if (num_cav >= max_size) {cerr << "Too many clusters"<<endl; return;}
//        // Cleaning booleans :

//        for (size_t I =0; I <= num_cav; I++)
//            for (size_t J =0; J <= num_cav; J++)Checked[I][J]=false;


//        Coords RACent;bool tomerge;
//#ifdef ICHEM_DEBUG
//        cout << "Number of clusters : "<< num_cav<< " " << clusterList.size()<<endl;
//#endif
//        do
//        {
//            modif=false;
//            for (vector<Cluster>::iterator ita= clusterList.begin();
//                 ita!=clusterList.end();
//                 ++ita)
//            {

//                RId = (*ita).idClus;
//                const set<const Box*> &RCLT = (*ita).listBox;
//                BestRCC=0;
//                BestCl=clusterList.end();
//                vector<Cluster>::iterator itb(ita);++itb;
//                for (; itb != clusterList.end();itb++)
//                {

//                    CId = (*itb).idClus;

//                    if (RId >= CId || Checked[RId][CId]==true) continue;

//                    const set<const Box*> &CCLT = (*itb).listBox;
//                    Checked[RId][CId]=true;
//                    Checked[CId][RId]=true;
//                    RCDist=0;RCCompact=0; NAdj=0;NClose=0; NMerge=0;dist=0;SMDist=100;RAProj=0;CAProj=0;NNProjR=0;NNProjC=0;RAprojmax=0;CAprojmax=0;
//                    RACent.clear();
//                    if (RCLT.size() < CCLT.size())
//                    {
//                        for (set<const Box*>::const_iterator itba = RCLT.begin();itba != RCLT.end();itba++){NAdj=0;
//                            for (set<const Box*>::const_iterator itbb = CCLT.begin();itbb != CCLT.end();itbb++)
//                            {
//                                if (projs[(*itba)->getId()]> RAprojmax) RAprojmax=projs[(*itba)->getId()];
//                                if (projs[(*itbb)->getId()]> CAprojmax) CAprojmax=projs[(*itbb)->getId()];
//                                dist=(*itba)->fixpos.calcDist((*itbb)->fixpos,SMDist+0.1);
//                                if (dist < SMDist){SMDist=dist;}
//                                if ((*itba)->isAdjacent((*itbb))) {
//                                    NAdj++;
//                                    if (dist <=diagBxstep+0.1)
//                                    {
//                                        RCCompact+=1/NSitCount.at((*itba)->getId())+1/NSitCount.at((*itbb)->getId());
//#ifdef ICHEM_DEBUG
//                                        cout <<"A::"<< RCCompact
//                                            << "\t"<<(*itba)->getId()
//                                            <<"\t"<<NSitCount.at((*itba)->getId())
//                                           <<"\t"<<(*itbb)->getId()
//                                          <<"\t"<<NSitCount.at((*itbb)->getId())<<endl;
//#endif
//                                    }
//                                    RAProj+=projs[(*itba)->getId()];NNProjR++;
//                                    CAProj+=projs[(*itbb)->getId()];NNProjC++;
//                                    RACent+=(*itba)->fixpos+(*itbb)->fixpos;
//                                }

//                            }

//                            if (NAdj > 0) NMerge++;
//                            if (NAdj >= 2)NClose++;
//                        }
//                    }
//                    else
//                    {
//                        for (set<const Box*>::const_iterator itba = CCLT.begin();itba != CCLT.end();itba++){NAdj=0;
//                            for (set<const Box*>::const_iterator itbb = RCLT.begin();itbb !=RCLT.end();itbb++)

//                            {
//                                if (projs[(*itba)->getId()]> RAprojmax) RAprojmax=projs[(*itba)->getId()];
//                                if (projs[(*itbb)->getId()]> CAprojmax) CAprojmax=projs[(*itbb)->getId()];
//                                dist=(*itba)->fixpos.calcDist((*itbb)->fixpos,SMDist+0.1);
//                                if (dist < SMDist){SMDist=dist;}

//                                if ((*itba)->isAdjacent((*itbb))) {
//                                    NAdj++;
//                                    RAProj+=projs[(*itba)->getId()];
//                                    if (dist <= diagBxstep+0.1)
//                                    {

//                                        RCCompact+=1/NSitCount.at((*itba)->getId())+1/NSitCount.at((*itbb)->getId());
//#ifdef ICHEM_DEBUG
//                                        cout <<"B::"<< RCCompact
//                                            << "\t"<<(*itba)->getId()
//                                            <<"\t"<<NSitCount.at((*itba)->getId())
//                                           <<"\t"<<(*itbb)->getId()
//                                          <<"\t"<<NSitCount.at((*itbb)->getId())<<endl;
//#endif
//                                    }
//                                    NNProjR++;CAProj+=projs[(*itbb)->getId()];NNProjC++;
//                                    RACent+=(*itba)->fixpos+(*itbb)->fixpos;
//                                }

//                            }

//                            if (NAdj > 0) NMerge++;
//                            if (NAdj >= 2)NClose++;

//                        }

//                    }
//                    if (NNProjR == 0) RAProj=0; else RAProj/=NNProjR;
//                    if (NNProjC == 0) CAProj=0; else CAProj/=NNProjC;
//                    tomerge=false;
//                    //    cout << "ST : "<< L1<< " " << L2 << " " << L3 << " " << L4 << " " << L5 << " " << L6<< L7 <<endl;
//                    //L1=0.1;L2=1;L3=1.5;L4=2;L5=2.5;L6=4.5,L7=8;
//                    if((RCLT.size() == 1 || CCLT.size()==1) && SMDist == boxStep) tomerge=true;
//                    else if(((RCLT.size() ==1  && CCLT.size() <= 5 ) ||(RCLT.size() <=5 && CCLT.size()   == 1)) && SMDist <= diagBxstep  && RCCompact > L1) tomerge=true;
//                    else if(((RCLT.size() <=15 && CCLT.size() <= 5 ) ||(RCLT.size() <=5 && CCLT.size()  <= 15)) && SMDist <= boxStep+0.1 && RCCompact > L2) tomerge=true;
//                    else if(((RCLT.size() <=15 && CCLT.size() <= 15) ||(RCLT.size() <=15 && CCLT.size() <= 15)) && SMDist <= boxStep+0.1 && RCCompact > L3) tomerge=true;
//                    else if(((RCLT.size() <=30 && CCLT.size() <= 15) ||(RCLT.size() <=15 && CCLT.size() <= 30)) && SMDist <= boxStep+0.1 && RCCompact > L4) tomerge=true;
//                    else if(((RCLT.size() <=30 && CCLT.size() <= 30) ||(RCLT.size() <=30 && CCLT.size() <= 30)) && SMDist <= boxStep+0.1 && RCCompact > L5) tomerge=true;
//                    else if ((RCLT.size() <=70 || CCLT.size() <=70) && SMDist <= boxStep+0.1 && RCCompact > L6) tomerge=true;
//                    else if (RCCompact >L7) tomerge=true;

//                    if (tomerge && RCCompact > BestRCC)
//                    {
//                        BestRCC = RCCompact;
//                        BestCl=itb;
//                    }
//#ifdef ICHEM_DEBUG
//                    if (SMDist <= 5){
//                        std::cout.precision(3);
//                        cout <<"CLT\t"<<((tomerge==true) ? "1":"0")<<"\t"
//                            <<RId<<"::"<<CId<<"\t"
//                           <<RCLT.size()<<"::"<<CCLT.size()<<"\t"
//                          <<NClose<<"::"<<NMerge<<"\t"
//                         <<RCDist<<"::"<<SMDist<<"\t"
//                        <<RAProj<<"::"<<CAProj<<"\t"
//                        << RCCompact<<"\t\t" << RAprojmax<<"::"<<CAprojmax
//                        <<endl;
//                    }
//#endif
//                }
//                if (BestCl != clusterList.end())
//                {
//                    modif=true;
//                    CId = (*BestCl).idClus;
//                    const set<const Box*> &CCLT=(*BestCl).listBox;
//#ifdef ICHEM_DEBUG
//                    cout << "MERGE "<<CId << " in " << RId << endl;
//#endif
//                    for (unsigned int I =0; I < num_cav; I++){ Checked[I][RId]=false; Checked[RId][I]=false;Checked[I][CId]=false; Checked[CId][I]=false;}

//                    for (set<const Box*>::const_iterator itbb = CCLT.begin();itbb != CCLT.end();itbb++)
//                    {
//                        (*ita).listBox.insert(*itbb);
//                        boxAssocClus.at((*itbb)->getId()) = distance(clusterList.begin(),ita);
//                    }
//                    const int di=std::distance(clusterList.begin(),BestCl);
//                    clusterList.erase(BestCl);
//                    ita=clusterList.begin();

//                    for (size_t i=0; i < AllBoxes.size();i++)
//                    {
//                        if (boxAssocClus.at(i)>di) boxAssocClus.at(i)--;

//                    }

//                }
//            }
//        } while(modif);
//#ifdef ICHEM_DEBUG
//        for (vector<Cluster>::iterator ita= clusterList.begin();
//             ita!=clusterList.end();
//             ++ita)
//        {
//            if (ita==clusterList.end())break;
//            set<const Box*> &RCLT = (*ita).listBox;
//            //cout << "DEV : "<< (*ita).first<<" " <<devproj<<endl;
//            //if (devproj < 10 || RCLT->Size < 20)continue;
//            Molecule *ml = new Molecule();
//            for (set<const Box*>::iterator itba = RCLT.begin();itba !=RCLT.end();itba++)
//            {
//                ml->addAtom("C",(*itba)->fixpos,"CA","C.3");
//            }

//            ostringstream oss; oss << "F_"<<MaxPrj<<"_"<<(*ita).idClus<<".mol2";
//            mw.newFile(oss.str());
//            mw.writeMOL2(ml);
//            delete ml;

//        }
//#endif


//    }// END MaxPrj >=THRES; MaxPrj-=5

//    // SELECTING CAVITIES :
//    num_cav=0;
//#ifdef ICHEM_DEBUG
//    cout << "NB CLUSTERS : "<< clusterList.size()<<endl;
//#endif
//    BoxList closeBox;
//    for (vector<Cluster>::const_iterator ita= clusterList.begin();
//         ita!=clusterList.end();
//         ++ita)
//    {
//        if (ita==clusterList.end())break;
//        const set<const Box*> &RCLT = (*ita).listBox;
//        double density=0, N_L=0;
//        for (set<const Box*>::iterator itba = RCLT.begin();itba !=RCLT.end();itba++)
//        {
//            N_L=0;
//            grid.getAdjacency(closeBox,*itba,grid.boxStep-0.1);
//            for (ItCBox itbb = closeBox.begin(); itbb != closeBox.end(); itbb++)
//            {
//                if (RCLT.find(const_cast<Box*>(*itbb)) != RCLT.end()) N_L++;
//            }
//            density += N_L/NSitCount.at((*itba)->getId());
//        }
//        density/=(double)RCLT.size();
//#ifdef ICHEM_DEBUG
//        cout << (*ita).idClus <<"\t"<< RCLT.size() <<"\t"<<density<<endl;
//#endif
//        //        if (density > 0.80){
//        //            cout << "size & density :"<< RCLT.size() << " " << density << endl;
//        //        }
//        if (density < 0.5 || (RCLT.size() <  min_size && density < 0.80) || RCLT.size() < 30 ) continue;
//#ifdef ICHEM_DEBUG
//        cout << "KEPT"<<endl;
//#endif
//        cavities.insert(pair<unsigned int,set<const Box*> >(num_cav,RCLT));
//        //Work here

//#ifdef ICHEM_DEBUG
//        Molecule *ml = new Molecule();
//        for (set<const Box*>::iterator itba = RCLT.begin();itba !=RCLT.end();itba++)
//        {
//            ml->addAtom("C",(*itba)->fixpos,"CA","C.3");
//        }
//        ostringstream oss; oss <<"KEPT_"<<(*ita).idClus<<"_"<<num_cav<<".mol2";
//        mw.newFile(oss.str());
//        mw.writeMOL2(ml);
//        delete ml;
//#endif
//        num_cav++;
//    }
//    numCav = num_cav;
//#ifdef ICHEM_DEBUG
//    cout << "NB CAV : " << num_cav << endl;
//#endif

//}

void VolSite::cavToMOL(Molecule& molen,  vector<const Box*> &boxlist, const string& molename, const string& size) const
{

    int MoleType=0;
    if (size == "4") MoleType=MoleType::CAV_4;
    else if (size=="6")MoleType=MoleType::CAV_6;
    else if (size=="8")MoleType=MoleType::CAV_8;
    else if (size=="12")MoleType=MoleType::CAV_12;
    else if (size=="ALL")MoleType=MoleType::CAV_ALL;
    //    Molecule molen= new Molecule(MoleType);

    molen.setName(molename);
    molen.clear();
    Residu &resSER = molen.addResidu("X",1,"SER");
    Residu &resALA = molen.addResidu("X",2,"ALA");
    Residu &resASP = molen.addResidu("X",3,"ASP");
    Residu &resLYS = molen.addResidu("X",4,"LYS");
    Residu &resPHE = molen.addResidu("X",5,"PHE");
    Residu &resGLY = molen.addResidu("X",6,"GLY");
    Residu &resCUB = molen.addResidu("X",7,"CUB");

    //    Box* bx=(Box*)NULL;
    int cube_pos;
    Atom * bestatm = (Atom*)NULL;
    double bestdist;

    AtomList listAtoms;
    for (vector<const Box*>::iterator itBx = boxlist.begin(); itBx != boxlist.end(); itBx++)
    {
        const Box &bx = **itBx;
        cube_pos=bx.getId();
        bestatm=(Atom*)NULL;
#ifdef ICHEM_DEBUG
        cout << bx.getId()<<" ";
#endif

        bestdist = 1000;
        bestatm= (Atom*)NULL;

        if (!bx.Props.isInCavity())continue;

        grid.getAdjacentAtoms(listAtoms,bx,5);

        for (ItCAtom itAtm= listAtoms.begin(); itAtm != listAtoms.end(); ++itAtm)
        {
            Atom &atm = **itAtm;
            //            atm = (*itAtm).second;
            if (atm.getParent().getMoleType() != MoleType::PROTEIN) continue;
            if (atm.isHydrogen()
                    || atm.getMOL2Type()=="Du"
                    || (!atm.props.isAromatic()
                        && !atm.props.isAcceptor()
                        && !atm.props.isDonor()
                        && !atm.props.isCation()
                        && !atm.props.isAnion()
                        && !atm.props.isHydrophobic()
                        && !atm.props.isMetal()        ) ) continue;
            // Calcul de distance
            const double dist=atm.fixpos.calcDist(bx.fixpos);

            if (dist > bestdist) continue;

            bestdist=dist;
            bestatm=&atm;
        }
#ifdef ICHEM_DEBUG
        cout << bestatm<<" "<< bestdist<<endl;
#endif


        if (bestatm != (Atom*)NULL){
            //            bx.Props=bestatm->props;
#ifdef ICHEM_DEBUG
            cout << "PROPS::"<<cube_pos<<"::"<<bx.rotpos.x<<":"<<bx.rotpos.y<<":"<<bx.rotpos.z<<"=>""=>("<<bx.fixpos.x<<":"<<bx.fixpos.y<<":"<<bx.fixpos.z<<")=>"<<bestatm->fixpos.calcDist(bx.fixpos)<<"=>"<<bestatm->getIdentifier()<<"\t"<< bestatm->props.toString()<<endl;
#endif
        }
        else continue;
        if (bestatm == (Atom*)NULL || bestdist >= 4)                        {molen.addAtom("H",bx.fixpos,"DU","H"   , &resCUB).setFormalCharge((double)projs[bx.getId()]); }
        else if (bestatm->props.isCation())                                 {molen.addAtom("O",bx.fixpos,"OD1","O.co2", &resASP).setFormalCharge((double)projs[bx.getId()]);}
        else if (bestatm->props.isAnion())                                  {molen.addAtom("N",bx.fixpos,"NZ","N.4" , &resLYS).setFormalCharge((double)projs[bx.getId()]);}
        else if (bestatm->props.isDonor()&& bestatm->props.isAcceptor())    {molen.addAtom("O",bx.fixpos,"OG","O.3" , &resSER).setFormalCharge((double)projs[bx.getId()]);}
        else if (bestatm->props.isAcceptor())                               {molen.addAtom("N",bx.fixpos,"N","N.am"  , &resALA).setFormalCharge((double)projs[bx.getId()]);}
        else if (bestatm->props.isDonor())                                  {molen.addAtom("O",bx.fixpos,"O","O.2"  , &resALA).setFormalCharge((double)projs[bx.getId()]);}
        else if (bestatm->props.isAromatic())                               {molen.addAtom("C",bx.fixpos,"CZ","C.ar", &resPHE).setFormalCharge((double)projs[bx.getId()]);}
        else if (bestatm->props.isHydrophobic())                            {molen.addAtom("C",bx.fixpos,"CA","C.3" , &resGLY).setFormalCharge((double)projs[bx.getId()]);}
        else {cout << "ATOM NOT RECOGNIZED "<< bestatm->toString()<<endl;continue;}


    }

}void VolSite::cavToMOLB(Molecule& molen,  vector<const Box*> &boxlist, const string& molename, const string& size,const int bestb) const
{

    int MoleType=0;
    if (size == "4") MoleType=MoleType::CAV_4;
    else if (size=="6")MoleType=MoleType::CAV_6;
    else if (size=="8")MoleType=MoleType::CAV_8;
    else if (size=="12")MoleType=MoleType::CAV_12;
    else if (size=="ALL")MoleType=MoleType::CAV_ALL;
    //    Molecule molen= new Molecule(MoleType);

    molen.setName(molename);
    molen.clear();
    Residu &resSER = molen.addResidu("X",1,"SER");
    Residu &resALA = molen.addResidu("X",2,"ALA");
    Residu &resASP = molen.addResidu("X",3,"ASP");
    Residu &resLYS = molen.addResidu("X",4,"LYS");
    Residu &resPHE = molen.addResidu("X",5,"PHE");
    Residu &resGLY = molen.addResidu("X",6,"GLY");
    Residu &resCUB = molen.addResidu("X",7,"CUB");

    //    Box* bx=(Box*)NULL;
    int cube_pos;
    Atom * bestatm = (Atom*)NULL;
    double bestdist;

    AtomList listAtoms;
    for (vector<const Box*>::iterator itBx = boxlist.begin(); itBx != boxlist.end(); itBx++)
    {
        const Box &bx = **itBx;
        cube_pos=bx.getId();
        bestatm=(Atom*)NULL;
#ifdef ICHEM_DEBUG
        cout << bx.getId()<<" ";
#endif

        bestdist = 1000;
        bestatm= (Atom*)NULL;

        if (!bx.Props.isInCavity())continue;

        grid.getAdjacentAtoms(listAtoms,bx,5);

        for (ItCAtom itAtm= listAtoms.begin(); itAtm != listAtoms.end(); ++itAtm)
        {
            Atom &atm = **itAtm;
            //            atm = (*itAtm).second;
            if (atm.getParent().getMoleType() != MoleType::PROTEIN) continue;
            if (atm.isHydrogen()
                    || atm.getMOL2Type()=="Du"
                    || (!atm.props.isAromatic()
                        && !atm.props.isAcceptor()
                        && !atm.props.isDonor()
                        && !atm.props.isCation()
                        && !atm.props.isAnion()
                        && !atm.props.isHydrophobic()
                        && !atm.props.isMetal()        ) ) continue;
            // Calcul de distance
            const double dist=atm.fixpos.calcDist(bx.fixpos);

            if (dist > bestdist) continue;

            bestdist=dist;
            bestatm=&atm;
        }
#ifdef ICHEM_DEBUG
        cout << bestatm<<" "<< bestdist<<endl;
#endif


        if (bestatm != (Atom*)NULL){
            //            bx.Props=bestatm->props;
#ifdef ICHEM_DEBUG
            cout << "PROPS::"<<cube_pos<<"::"<<bx.rotpos.x<<":"<<bx.rotpos.y<<":"<<bx.rotpos.z<<"=>""=>("<<bx.fixpos.x<<":"<<bx.fixpos.y<<":"<<bx.fixpos.z<<")=>"<<bestatm->fixpos.calcDist(bx.fixpos)<<"=>"<<bestatm->getIdentifier()<<"\t"<< bestatm->props.toString()<<endl;
#endif
        }
        else continue;
        if (projs[bx.getId()] >= bestb){
            if (bestatm == (Atom*)NULL || bestdist >= 4)                        {molen.addAtom("H",bx.fixpos,"DU","H"   , &resCUB).setFormalCharge((double)projs[bx.getId()]); }
            else if (bestatm->props.isCation())                                 {molen.addAtom("O",bx.fixpos,"OD1","O.co2", &resASP).setFormalCharge((double)projs[bx.getId()]);}
            else if (bestatm->props.isAnion())                                  {molen.addAtom("N",bx.fixpos,"NZ","N.4" , &resLYS).setFormalCharge((double)projs[bx.getId()]);}
            else if (bestatm->props.isDonor()&& bestatm->props.isAcceptor())    {molen.addAtom("O",bx.fixpos,"OG","O.3" , &resSER).setFormalCharge((double)projs[bx.getId()]);}
            else if (bestatm->props.isAcceptor())                               {molen.addAtom("N",bx.fixpos,"N","N.am"  , &resALA).setFormalCharge((double)projs[bx.getId()]);}
            else if (bestatm->props.isDonor())                                  {molen.addAtom("O",bx.fixpos,"O","O.2"  , &resALA).setFormalCharge((double)projs[bx.getId()]);}
            else if (bestatm->props.isAromatic())                               {molen.addAtom("C",bx.fixpos,"CZ","C.ar", &resPHE).setFormalCharge((double)projs[bx.getId()]);}
            else if (bestatm->props.isHydrophobic())                            {molen.addAtom("C",bx.fixpos,"CA","C.3" , &resGLY).setFormalCharge((double)projs[bx.getId()]);}
            else {cout << "ATOM NOT RECOGNIZED "<< bestatm->toString()<<endl;continue;}

        }
    }

}


void VolSite::cubeToPHA(Molecule& molen, PharmList& MonPharma , Molecule& ligand, const string& molename, const string& size, const int mod, const bool check,const bool exclu) const
{
    const int sizetot =     cubes_cavite.size();
    const Box *boxliste[40000]  = {NULL} ;
    Atom *atomlist[cubes_cavite.size()];
    float distlist[cubes_cavite.size()]; //0
    bool Usedlist[cubes_cavite.size()]; // Ne faire q'un tableau
    bool PHAlist[cubes_cavite.size()];
    bool cluslist[cubes_cavite.size()];
    bool Arolist[cubes_cavite.size()];
    Coords coorlist[cubes_cavite.size()];
    int tol[cubes_cavite.size()];
    int group[cubes_cavite.size()];
    int nb_enfoui[120] = {0};
    int cbnum = 0, bestenfoui=0, hydro_present=0;

    for (int i = 0;i < cubes_cavite.size(); i++){
        atomlist[i] = NULL;
        distlist[i] = 0;
        Usedlist[i] = false;
        PHAlist[i] = false;
        cluslist[i] = false;
        Arolist[i] = false;
        coorlist[i] = Coords(0,0,0);
        tol[i]= 0;
        group[i]= 0;
    }

    int MoleType=0;
    if (size == "4") MoleType=MoleType::CAV_4;
    else if (size=="6")MoleType=MoleType::CAV_6;
    else if (size=="8")MoleType=MoleType::CAV_8;
    else if (size=="12")MoleType=MoleType::CAV_12;
    else if (size=="ALL")MoleType=MoleType::CAV_ALL;

    //    Molecule molen= new Molecule(MoleType);
    molen.setName(molename);
    molen.clear();
    Residu &resSER = molen.addResidu("X",1,"SER");
    Residu &resALA = molen.addResidu("X",2,"ALA");
    Residu &resASP = molen.addResidu("X",3,"ASP");
    Residu &resLYS = molen.addResidu("X",4,"LYS");
    Residu &resPHE = molen.addResidu("X",5,"PHE");
    Residu &resGLY = molen.addResidu("X",6,"GLY");
    Residu &resTRP = molen.addResidu("X",7,"TRP");
    Residu &resCUB = molen.addResidu("X",8,"CUB");

    //    Box* bx=(Box*)NULL;
    int cube_pos = 0, envi = 0 , envi_H=0;
    Atom * bestatm = (Atom*)NULL;
    Atom * bestatm_h = (Atom*)NULL;
    Atom * bestatm_h1 = (Atom*)NULL;
    Atom * bestatm_h2 = (Atom*)NULL;

    double bestdist;
    double dist_hyd1 = 10000;
    double dist_hyd2 = 10000;

    AtomList listAtoms;
    AtomList listAtomsH;
    AtomList listatomsAro;
    ///% reviens here


    //    free(v); free(c);
    for (ItCBox itpos = cubes_cavite.begin(); itpos != cubes_cavite.end(); itpos++)
    {
        const Box &bx=**itpos;
        double angle;
        if (check){
            if (!bx.isUsed()){
                continue;
            }
        }
        double dist=1000;
        double tmpdist = 0;
        for (ItCAtom atom = ligand.firstAtom();atom != ligand.lastAtom();++atom)
        {
            tmpdist = bx.fixpos.calcDist((*atom)->fixpos);
            if (tmpdist < dist) dist = tmpdist;
        }
        if (dist > 2) {
            continue;
        }


        cube_pos=bx.getId();
        bestatm = (Atom*)NULL;
        bestatm_h = (Atom*)NULL;
        bestatm_h1 = (Atom*)NULL;
        bestatm_h2 = (Atom*)NULL;

#ifdef ICHEM_DEBUG
        cout << bx.getId()<<" ";
#endif

        bestdist = 10000;
        dist_hyd1 = 10000;
        dist_hyd2 = 10000;
        double best_hyd = 10000;
        bool hyd1 = false, hyd2 = false;
        int best_res;

        if (!bx.Props.isInCavity())continue;
        grid.getAdjacentAtoms(listAtoms,bx,4.5);
        envi = 0;
        envi_H=0;
        int in_ring=0;
        for (ItCAtom itAtm = listAtoms.begin(); itAtm != listAtoms.end(); ++itAtm)
        {

            Atom &atm = **itAtm;

            int res_temp = atm.getResidu()->getNum();
            envi++;
            //            atm = (*itAtm).second;
            if (atm.getParent().getMoleType() != MoleType::PROTEIN) continue;
            if (atm.isHydrogen()) { /*cout << "hydro passe a 1"  << endl;*/ hydro_present = 1; continue;}
            if (atm.getMOL2Type()=="Du"
                    || (!atm.props.isAromatic()
                        && !atm.props.isAcceptor()
                        && !atm.props.isDonor()
                        && !atm.props.isCation()
                        && !atm.props.isAnion()
                        && !atm.props.isHydrophobic()
                        && !atm.props.isMetal()        ) ) continue;
            // Calcul de distance
            double dist = 10000;
            double dist_temp = 10000;

            dist_temp = atm.fixpos.calcDist(bx.fixpos);

            //            if (Arolist[cbnum] == true ) continue;
            if (atm.props.isMetal()){
                if ( dist_temp < 2.8) {
                    dist = atm.fixpos.calcDist(bx.fixpos);
                }
            }
            else if (atm.props.isDonor() && atm.props.isAcceptor() && dist_temp < 3.5 ){
                if (hydro_present && atm.getAtomicName() =="N") {
                    atm.checkMOL2type();
                    atm.props.setAromatic(true);
                    atm.props.setRing(true);
                }
                //                cout << atm.toString() << endl;
                dist = atm.fixpos.calcDist(bx.fixpos);

            }
            else if (atm.props.isDonor() && dist_temp < 3.5){
                if (hydro_present) {

                    for (size_t i=0; i< atm.getNumBond();++i)
                    {
                        const Atom &atm2 = atm.getAtomLinked(i);
                        if (!atm2.isHydrogen())  hydro_present =1;   continue;
                        angle = atm2.fixpos.calcAngle(atm.fixpos,  bx.fixpos);
                        if (angle >= Angl_H-AngT_H && angle <= Angl_H+AngT_H )  {
                            dist = dist_temp;
                        }
                    }
                }else {
                    dist = dist_temp;
                }

            }
            else if ( atm.props.isAcceptor() && dist_temp < 3.5){
                if (hydro_present) {
                    for (size_t i=0; i< atm.getNumBond();++i)
                    {
                        const Atom &atm2 = atm.getAtomLinked(i);
                        if (atm2.isHydrogen())     continue;
                        angle = atm.fixpos.calcAngle(atm2.fixpos,  bx.fixpos);
                        if (angle >= Angl_D-AngT_H && angle <= Angl_D+AngT_H )  {
                            dist = dist_temp;
                            //                        cout << " acceptor correct" << endl;
                        }

                    }
                }else {
                    dist = dist_temp;

                }
            }
            else if (atm.props.isAnion() && dist_temp < 4 ){
                dist = dist_temp;

            }

            else if (atm.props.isCation() && dist_temp < 4){
                dist = dist_temp;
            }


            else if ( atm.props.isAcceptor() && dist_temp < 3.5){
                for (size_t i=0; i< atm.getNumBond();++i)
                {
                    const Atom &atm2 = atm.getAtomLinked(i);
                    if (atm2.isHydrogen())     continue;
                    angle = atm.fixpos.calcAngle(atm2.fixpos,  bx.fixpos);
                    if (angle >= Angl_D-AngT_H && angle <= Angl_D+AngT_H )  {
                        dist = dist_temp;
                        //                        cout << " acceptor correct" << endl;
                    }
                }
            }
            if (atm.props.isAromatic() && dist_temp < 4.5 && dist_temp < bestdist){
                Atom * keep = (Atom*)NULL;
                grid.getAdjacentAtoms(listatomsAro,*(atm.getBox(&grid)),6);
                int max=0;
                Coords centroid;
                centroid.setCoords(atm.fixpos);
                int weight = 1;
                in_ring=1;
                for (ItCAtom itAtm= listatomsAro.begin(); itAtm != listatomsAro.end(); ++itAtm)
                {
                    Atom &atmaro = **itAtm;
                    if (atmaro.props.isAromatic() && atmaro.getResidu() == atm.getResidu()){
                        if (&atmaro != &atm){
                            keep = &atmaro;
                        }
                        if (atmaro.fixpos.calcDist(atm.fixpos)<8){
                            centroid.setCoords( (centroid.x+atmaro.fixpos.x),(centroid.y+atmaro.fixpos.y),(centroid.z+atmaro.fixpos.z));
                            weight++;
                            in_ring++;
                        }
                    }
                }
                //                cout << " polop " << in_ring << " " << weight<< endl;
                centroid.setCoords(centroid.x/weight,centroid.y/weight,centroid.z/weight);

                Coords vecto = centroid.getNormal(atm.fixpos,keep->fixpos)-centroid;
                //                cout << centroid.toString() << " " << atm.fixpos.toString() << " " << keep->fixpos.toString() << vecto.toString() << endl;

                Coords distant = centroid+vecto+vecto+vecto+vecto;
                Coords distantneg = centroid-vecto-vecto-vecto-vecto;
                //                cout << distant.toString() << endl;
                //                cout << "distance " << bx.fixpos.calcDist(distant)<< endl;
                if (bx.fixpos.calcDist(distant) < 1.5 || bx.fixpos.calcDist(distantneg) < 1.5){
                    dist = dist_temp;
                    Arolist[cbnum] = true;
                }


            }
            if (atm.props.isHydrophobic()){
                envi_H++;
            }

            if (atm.props.isHydrophobic() && dist_temp < 4.5 && hyd1 && res_temp != best_res){
                if (dist_temp < dist_hyd2){
                    dist_hyd2 = dist_temp;
                    bestatm_h2 = &atm;
                }
                hyd2=true;

            }

            else if (atm.props.isHydrophobic() && dist_temp < 4.5 && !hyd1){
                dist_hyd1 = dist_temp;
                hyd1 = true;
                best_res = res_temp;
                bestatm_h1 = &atm;
            }

            /////////////////////////HERE
            if (dist >= bestdist ) {

                continue;
            }
            //            cout << " dista assigment \t" << bestdist << "\t " << cbnum << endl;
            //            cout << dist << " : " << bestdist << endl;
            bestdist = dist;
            bestatm = &atm;


        }

        if (hyd2){

            if (dist_hyd1 < dist_hyd2){
                best_hyd = dist_hyd1;
                bestatm_h =  bestatm_h1;
            }else {
                best_hyd = dist_hyd2;
                bestatm_h =  bestatm_h2;
            }
            if (100*envi_H*envi >=50){
                //                if (bestdist <= best_hyd ){
                //                    continue;
                //                }
                //                else
                if (bestdist > dist_hyd1 && bestdist > dist_hyd2){
                    bestdist = best_hyd;
                    bestatm = bestatm_h;
                    Arolist[cbnum] = false;
                }
                //                else {
                //                    //                    cout << " between hydrophobe (>50%)  " <<  bestatm->props.toString() << "   \t" << dist_hyd1 << "\t||\t" <<  bestdist  << "\t||\t" << dist_hyd2 << "  " << bestatm->getIdentifier() << endl;
                //                    bestdist = best_hyd;
                //                    bestatm = bestatm_h;
                //                    Arolist[cbnum] = false;
                //                }

            }
        }

        //    for (ItCCycle itC = cyclelist.begin()    ; itC != cyclelist.end();itC++)
        //    {
        //    }
#ifdef ICHEM_DEBUG
        cout << "cubetopha "<< bestdist << " " <<  bestatm->toString() <<endl;
#endif
        //        cout << "end" << endl;
        boxliste[cbnum] = &bx;
        atomlist[cbnum] = bestatm;
        distlist[cbnum] = bestdist;

        cbnum ++;
        for (int enf=0; enf <projs[bx.getId()];enf ++)
        {
            nb_enfoui[enf]++;
        }
        //        if (bestatm->props.isAromatic()){
        //            cout << "Aromatique proche de " << in_ring << endl;
        //        }
    }
    int pourcent=0,ccube=0;
    for (int enf=0; enf<120;enf++)
    {
        if (( nb_enfoui[enf] <cavitykeepC || 100*nb_enfoui[enf]/ cbnum < cavitykeepP) && bestenfoui==0){
            bestenfoui = enf;
            pourcent= 100*nb_enfoui[enf]/ cbnum;
            ccube=nb_enfoui[enf];
        }
    }
    //    cout << endl<< endl<< endl<< endl<< endl<< endl ;
//    if (mod > 2){
//        if (bestenfoui > 55){
//            cout << "best b : -b" << bestenfoui<<" || "<< ccube <<" Cavity cubes || Initial CAV : "<< pourcent <<"%"<< endl;
//        }else
//            cout <<"best b : default" << endl;
//    }

    // mod finaux qui prennent en compte la distance au centre du cube general
    if (mod == 22 || mod == 23 || mod == 24 || mod == 25 || mod == 26 || mod == 27 ) {
        double dist_noH=0;
        double dist_h=0;
        double dist_center = 8;
        if ( mod == 22 ) {
            dist_noH=7;
            dist_h = 3.1;
        }else if (mod == 23) {
            dist_noH=5.5;
            dist_h = 3.1;
        }else if (mod == 24 || mod == 25 || mod == 26 || mod == 27 ) {
            dist_noH=3.1;
            dist_h = 3.1;
        }
        for ( int itPHA =0 ;itPHA<cbnum;itPHA++){
            if (atomlist[itPHA] == (Atom*)NULL || Usedlist[itPHA] ) { continue;}

            if (boxliste[itPHA]->fixpos.calcDist(grid.getcenter()) > dist_center){continue;}
            if (coorlist[itPHA].calcDist(Coords(0,0,0)) == 0) { coorlist[itPHA]=boxliste[itPHA]->fixpos;} // barycentre des cluster coordonées des cubes de la cavité
            if (mod == 27 ){ if ((double)projs[boxliste[itPHA]->getId()] < 65) { continue;}}
            else if ((double)projs[boxliste[itPHA]->getId()] < 80) { continue;}
            //            PHAlist[itPHA]=false;


            PHAlist[itPHA]=true;
            //            if (atomlist[itPHA]->props.isMetal()){
            //                int weight =1;
            //                if (atomlist[itPHA2] == (Atom*)NULL  ){ continue;}
            //                if (atomlist[itPHA] == atomlist[itPHA2] && !Usedlist[itPHA2]){
            //                    coorlist[itPHA].setCoords( coorlist[itPHA].x +boxliste[itPHA2]->fixpos.x , coorlist[itPHA].y +boxliste[itPHA2]->fixpos.y ,coorlist[itPHA].z +boxliste[itPHA2]->fixpos.z );
            //                    Usedlist[itPHA2] = true;
            //                    PHAlist[itPHA2] = false;
            //                    weight++;
            //                    PHAlist[itPHA] = true;
            //                }
            //            }
            //            Usedlist[itPHA] = true;
            //            tol[itPHA]=weight;
            //            coorlist[itPHA].setCoords(coorlist[itPHA].x/weight,coorlist[itPHA].y/weight,coorlist[itPHA].z/weight);

            //            }
            if ( Arolist[itPHA] == true ){
                //                cout << " un aro " << endl;
                int weight =1;
                for ( int itPHA2 =itPHA+1 ;itPHA2<cbnum;itPHA2++){
                    if (atomlist[itPHA2] == (Atom*)NULL  ){ continue;}
                    if (atomlist[itPHA]->getResidu()->getNum() == atomlist[itPHA2]->getResidu()->getNum() && Arolist[itPHA2] == true && !Usedlist[itPHA2]){
                        coorlist[itPHA].setCoords( coorlist[itPHA].x +boxliste[itPHA2]->fixpos.x , coorlist[itPHA].y +boxliste[itPHA2]->fixpos.y ,coorlist[itPHA].z +boxliste[itPHA2]->fixpos.z );
                        Usedlist[itPHA2] = true;
                        PHAlist[itPHA2] = false;
                        weight++;
                        PHAlist[itPHA] = true;
                        cluslist[itPHA] = true;
                    }
                }
                Usedlist[itPHA] = true;
                tol[itPHA]=weight;
                coorlist[itPHA].setCoords(coorlist[itPHA].x/weight,coorlist[itPHA].y/weight,coorlist[itPHA].z/weight);



            }else if ( mod == 25 && atomlist[itPHA]->props.isAcceptor() && atomlist[itPHA]->props.isDonor()){
                int weight =1;
                if (atomlist[itPHA]->getResiduName() == "HOH"){
                    bool d_ok = false;
                    bool a_ok = false;
                    for (size_t i=0; i< atomlist[itPHA]->getNumBond();++i)
                    {
                        const Atom &atm2 = atomlist[itPHA]->getAtomLinked(i);
                        if (!atm2.isHydrogen())     continue;
                        double angle_d = atm2.fixpos.calcAngle(atomlist[itPHA]->fixpos,  coorlist[itPHA]);
                        double angle_a = atomlist[itPHA]->fixpos.calcAngle(atm2.fixpos,  coorlist[itPHA]);

                        if (      angle_d >= Angl_H-AngT_H && angle_d <= Angl_H+AngT_H )  {
                            d_ok = true;
                        }else if (angle_a >= Angl_H-AngT_H && angle_a <= Angl_H+AngT_H )  {
                            a_ok = true;
                        }
                    }
                    if (d_ok){
                        for ( int itPHA2 =itPHA+1 ;itPHA2<cbnum;itPHA2++){
                            if (atomlist[itPHA2] == (Atom*)NULL  ){ continue;}
                            if ( atomlist[itPHA]!=atomlist[itPHA2] || Usedlist[itPHA2] ){continue;}
                            if (boxliste[itPHA2]->fixpos.calcDist(grid.getcenter()) > dist_center || coorlist[itPHA].calcDist(boxliste[itPHA2]->fixpos) > dist_noH){continue;}

                            for (size_t i=0; i< atomlist[itPHA]->getNumBond();++i)
                            {
                                const Atom &atm2 = atomlist[itPHA]->getAtomLinked(i);
                                if (!atm2.isHydrogen())     continue;
                                double angle_d = atm2.fixpos.calcAngle(atomlist[itPHA]->fixpos,  coorlist[itPHA]);
                                if (      angle_d >= Angl_H-AngT_H && angle_d <= Angl_H+AngT_H )  {
                                    PHAlist[itPHA] = true;
                                    PHAlist[itPHA2] = false;
                                    coorlist[itPHA].setCoords((weight*coorlist[itPHA].x +boxliste[itPHA2]->fixpos.x)/(weight+1), (weight*coorlist[itPHA].y +boxliste[itPHA2]->fixpos.y)/(weight+1), (weight*coorlist[itPHA].z +boxliste[itPHA2]->fixpos.z)/(weight+1));
                                    Usedlist[itPHA2] = true;
                                    weight++;
                                }
                            }
                        }
                    }
                    else if (a_ok){
                        for ( int itPHA2 =itPHA+1 ;itPHA2<cbnum;itPHA2++){
                            if (atomlist[itPHA2] == (Atom*)NULL  ){ continue;}
                            if ( atomlist[itPHA]!=atomlist[itPHA2] || Usedlist[itPHA2] ){continue;}
                            if (boxliste[itPHA2]->fixpos.calcDist(grid.getcenter()) > dist_center || coorlist[itPHA].calcDist(boxliste[itPHA2]->fixpos) > dist_noH){continue;}

                            for (size_t i=0; i< atomlist[itPHA]->getNumBond();++i)
                            {
                                const Atom &atm2 = atomlist[itPHA]->getAtomLinked(i);
                                if (!atm2.isHydrogen())     continue;
                                double angle_a = atomlist[itPHA]->fixpos.calcAngle(atm2.fixpos,  coorlist[itPHA]);
                                if (      angle_a >= Angl_H-AngT_H && angle_a <= Angl_H+AngT_H )  {
                                    PHAlist[itPHA] = true;
                                    PHAlist[itPHA2] = false;
                                    coorlist[itPHA].setCoords((weight*coorlist[itPHA].x +boxliste[itPHA2]->fixpos.x)/(weight+1), (weight*coorlist[itPHA].y +boxliste[itPHA2]->fixpos.y)/(weight+1), (weight*coorlist[itPHA].z +boxliste[itPHA2]->fixpos.z)/(weight+1));
                                    Usedlist[itPHA2] = true;
                                    weight++;
                                }
                            }
                        }
                    }
                    tol[itPHA] = weight;
                    Usedlist[itPHA] = true;
                } else if (atomlist[itPHA]->getAtomicName() == "O"){
                    bool d_ok = false;
                    bool a_ok = false;
                    for (size_t i=0; i< atomlist[itPHA]->getNumBond();++i)
                    {
                        const Atom &atm2 = atomlist[itPHA]->getAtomLinked(i);
                        if (!atm2.isHydrogen())     continue;
                        double angle_d = atm2.fixpos.calcAngle(atomlist[itPHA]->fixpos,  boxliste[itPHA]->fixpos);
                        double angle_a = atomlist[itPHA]->fixpos.calcAngle(atm2.fixpos,  boxliste[itPHA]->fixpos);
                        if (angle_d >= D180-AngT_H && angle_d <= D180+AngT_H )  {
                            d_ok = true;
                        }else if (angle_a >= D120-AngT_H && angle_a <= D120+AngT_H )  {
                            a_ok = true;
                        }
                    }
                    if (d_ok){
                        for ( int itPHA2 =itPHA+1 ;itPHA2<cbnum;itPHA2++){
                            if (atomlist[itPHA2] == (Atom*)NULL  ){ continue;}
                            if ( atomlist[itPHA]!=atomlist[itPHA2] || Usedlist[itPHA2] ){continue;}
                            if (boxliste[itPHA2]->fixpos.calcDist(grid.getcenter()) > dist_center || coorlist[itPHA].calcDist(boxliste[itPHA2]->fixpos) > dist_noH){continue;}

                            for (size_t i=0; i< atomlist[itPHA]->getNumBond();++i)
                            {
                                const Atom &atm2 = atomlist[itPHA]->getAtomLinked(i);
                                if (!atm2.isHydrogen())     continue;
                                double angle_d = atm2.fixpos.calcAngle(atomlist[itPHA]->fixpos,  boxliste[itPHA]->fixpos);
                                if (angle_d >= D180-AngT_H && angle_d <= D180+AngT_H )  {
                                    PHAlist[itPHA] = true;
                                    PHAlist[itPHA2] = false;
                                    coorlist[itPHA].setCoords((weight*coorlist[itPHA].x +boxliste[itPHA2]->fixpos.x)/(weight+1), (weight*coorlist[itPHA].y +boxliste[itPHA2]->fixpos.y)/(weight+1), (weight*coorlist[itPHA].z +boxliste[itPHA2]->fixpos.z)/(weight+1));
                                    Usedlist[itPHA2] = true;
                                    weight++;
                                }
                            }
                        }
                    }
                    else if (a_ok){
                        for ( int itPHA2 =itPHA+1 ;itPHA2<cbnum;itPHA2++){
                            if (atomlist[itPHA2] == (Atom*)NULL  ){ continue;}
                            if ( atomlist[itPHA]!=atomlist[itPHA2] || Usedlist[itPHA2] ){continue;}
                            if (boxliste[itPHA2]->fixpos.calcDist(grid.getcenter()) > dist_center || coorlist[itPHA].calcDist(boxliste[itPHA2]->fixpos) > dist_noH){continue;}

                            for (size_t i=0; i< atomlist[itPHA]->getNumBond();++i)
                            {
                                const Atom &atm2 = atomlist[itPHA]->getAtomLinked(i);
                                if (!atm2.isHydrogen())     continue;
                                double angle_a = atomlist[itPHA]->fixpos.calcAngle(atm2.fixpos,  boxliste[itPHA]->fixpos);
                                if (      angle_a >= D120-AngT_H && angle_a <= D120+AngT_H)  {
                                    PHAlist[itPHA] = true;
                                    PHAlist[itPHA2] = false;
                                    coorlist[itPHA].setCoords((weight*coorlist[itPHA].x +boxliste[itPHA2]->fixpos.x)/(weight+1), (weight*coorlist[itPHA].y +boxliste[itPHA2]->fixpos.y)/(weight+1), (weight*coorlist[itPHA].z +boxliste[itPHA2]->fixpos.z)/(weight+1));
                                    Usedlist[itPHA2] = true;
                                    weight++;
                                }
                            }
                        }
                    } else {

                        for ( int itPHA2 =itPHA+1 ;itPHA2<cbnum;itPHA2++){
                            if (atomlist[itPHA2] == (Atom*)NULL  ){ continue;}
                            if ( atomlist[itPHA]!=atomlist[itPHA2] || Usedlist[itPHA2] ){continue;}
                            if (boxliste[itPHA2]->fixpos.calcDist(grid.getcenter()) > dist_center || coorlist[itPHA].calcDist(boxliste[itPHA2]->fixpos) > dist_noH){continue;}

                            for (size_t i=0; i< atomlist[itPHA]->getNumBond();++i)
                            {
                                const Atom &atm2 = atomlist[itPHA]->getAtomLinked(i);
                                if (!atm2.isHydrogen())     continue;
                                double angle_a = atomlist[itPHA]->fixpos.calcAngle(atm2.fixpos,  boxliste[itPHA]->fixpos);
                                double angle_d = atm2.fixpos.calcAngle(atomlist[itPHA]->fixpos,  boxliste[itPHA]->fixpos);
                                if ( !(angle_a >= D120-AngT_H && angle_a <= D120+AngT_H)  && !(angle_d >= D180-AngT_H && angle_d <= D180+AngT_H))  {
                                    PHAlist[itPHA] = true;
                                    PHAlist[itPHA2] = false;
                                    coorlist[itPHA].setCoords((weight*coorlist[itPHA].x +boxliste[itPHA2]->fixpos.x)/(weight+1), (weight*coorlist[itPHA].y +boxliste[itPHA2]->fixpos.y)/(weight+1), (weight*coorlist[itPHA].z +boxliste[itPHA2]->fixpos.z)/(weight+1));
                                    Usedlist[itPHA2] = true;
                                    weight++;
                                }
                            }
                        }
                    }


                }else {
                    int weight =1;
                    for ( int itPHA2 =itPHA+1 ;itPHA2<cbnum;itPHA2++){
                        if (atomlist[itPHA2] == (Atom*)NULL  ){ continue;}
                        if (atomlist[itPHA]!=atomlist[itPHA2] || Usedlist[itPHA2] ){continue;}
                        if (boxliste[itPHA2]->fixpos.calcDist(grid.getcenter()) > dist_center){continue;}

                        if ( atomlist[itPHA]==atomlist[itPHA2] && !Usedlist[itPHA2] )
                        {
                            //                            if (atomlist[itPHA]->props.isAcceptor() && atomlist[itPHA]->props.isDonor()){
                            //                                cout << "A/D\t" << coorlist[itPHA].calcDist(boxliste[itPHA2]->fixpos) << "\t" << atomlist[itPHA]->props.toString() << "\tvs  " << atomlist[itPHA2]->props.toString()<< endl;
                            //                            }
                            if ( coorlist[itPHA].calcDist(boxliste[itPHA2]->fixpos) < dist_noH ) // cluster de XA pour les non hydrophobes
                            {
                                PHAlist[itPHA] = true;
                                PHAlist[itPHA2] = false;
                                coorlist[itPHA].setCoords((weight*coorlist[itPHA].x +boxliste[itPHA2]->fixpos.x)/(weight+1), (weight*coorlist[itPHA].y +boxliste[itPHA2]->fixpos.y)/(weight+1), (weight*coorlist[itPHA].z +boxliste[itPHA2]->fixpos.z)/(weight+1));
                                Usedlist[itPHA2] = true;
                                weight++;
                                //                                if (atomlist[itPHA]->props.isAcceptor() && atomlist[itPHA]->props.isDonor()){
                                //                                    cout << "FuuuGL_VENDOR:   NVIDIA Corporationusion "<< endl;
                                //                                }
                            }
                        }
                    }
                    tol[itPHA] = weight;
                    Usedlist[itPHA] = true;
                }

            }
            else if (!atomlist[itPHA]->props.isHydrophobic() )
            {
                int weight =1;
                for ( int itPHA2 =itPHA+1 ;itPHA2<cbnum;itPHA2++){
                    if (atomlist[itPHA2] == (Atom*)NULL  ){ continue;}
                    if (boxliste[itPHA2]->fixpos.calcDist(grid.getcenter()) > dist_center){continue;}

                    if ( atomlist[itPHA]==atomlist[itPHA2] && !Usedlist[itPHA2] )
                    {
                        //                            if (atomlist[itPHA]->props.isAcceptor() && atomlist[itPHA]->props.isDonor()){
                        //                                cout << "A/D\t" << coorlist[itPHA].calcDist(boxliste[itPHA2]->fixpos) << "\t" << atomlist[itPHA]->props.toString() << "\tvs  " << atomlist[itPHA2]->props.toString()<< endl;
                        //                            }
                        if ( coorlist[itPHA].calcDist(boxliste[itPHA2]->fixpos) < dist_noH ) // cluster de XA pour les non hydrophobes
                        {
                            PHAlist[itPHA] = true;
                            cluslist[itPHA] =true;
                            PHAlist[itPHA2] = false;
                            coorlist[itPHA].setCoords((weight*coorlist[itPHA].x +boxliste[itPHA2]->fixpos.x)/(weight+1), (weight*coorlist[itPHA].y +boxliste[itPHA2]->fixpos.y)/(weight+1), (weight*coorlist[itPHA].z +boxliste[itPHA2]->fixpos.z)/(weight+1));
                            Usedlist[itPHA2] = true;
                            weight++;
                            //                                if (atomlist[itPHA]->props.isAcceptor() && atomlist[itPHA]->props.isDonor()){
                            //                                    cout << "FuuuGL_VENDOR:   NVIDIA Corporationusion "<< endl;
                            //                                }
                        }
                    }
                }
                tol[itPHA] = weight;
                Usedlist[itPHA] = true;

            }
            else if( atomlist[itPHA]->props.isHydrophobic())
            {

                int weight =1;
                for ( int itPHA2 =itPHA+1 ;itPHA2<cbnum;itPHA2++){
                    if (atomlist[itPHA2] == (Atom*)NULL || Usedlist[itPHA2]  ){ continue;}

                    if (atomlist[itPHA]->getResidu()->getNum() == atomlist[itPHA2]->getResidu()->getNum() && coorlist[itPHA].calcDist(boxliste[itPHA2]->fixpos) < dist_h) {
                        PHAlist[itPHA] = true;
                        cluslist[itPHA] = true;
                        PHAlist[itPHA2]=false;
                        coorlist[itPHA].setCoords((weight*coorlist[itPHA].x +boxliste[itPHA2]->fixpos.x)/(weight+1), (weight*coorlist[itPHA].y +boxliste[itPHA2]->fixpos.y)/(weight+1), (weight*coorlist[itPHA].z +boxliste[itPHA2]->fixpos.z)/(weight+1));

                        //                        coorlist[itPHA].setCoords( coorlist[itPHA].x +boxliste[itPHA2]->fixpos.x ,coorlist[itPHA].y +boxliste[itPHA2]->fixpos.y ,coorlist[itPHA].z +boxliste[itPHA2]->fixpos.z );
                        Usedlist[itPHA2] = true;
                        weight++;
                    }

                }
                tol[itPHA] = weight;
                //                coorlist[itPHA].setCoords(coorlist[itPHA].x/weight,coorlist[itPHA].y/weight,coorlist[itPHA].z/weight);
                Usedlist[itPHA]= true;
            }
            //                cout << itPHA << "\t" << atomlist[itPHA]->props.toString() << "\t" << PHAlist[itPHA]<< endl;
            else { PHAlist[itPHA]= false;}


        }

        // Verifie les points hydrophobe restant pour etre sur qu'il y est toujours au moins deux résidus hydrophobes autour de chaque points hydrophobes
        //        for ( int itPHA =0 ;itPHA<cbnum;itPHA++){
        //            int keepRes =100000;
        //            int Res=0;
        //            if (projs[boxliste[itPHA]->getId()]<bestenfoui) {  continue;}
        //            if (atomlist[itPHA] != (Atom*)NULL && PHAlist[itPHA] && atomlist[itPHA]->props.isHydrophobic()){
        //                //                cout << endl<<"hydrophobe" <<endl;
        //                AtomList listAtoms;
        //                Box *mbx= grid.getBox(coorlist[itPHA]);
        //                grid.getAdjacentAtoms(listAtoms,*mbx,4);

        //                for (ItCAtom itAtm= listAtoms.begin(); itAtm != listAtoms.end(); ++itAtm)
        //                {
        //                    Atom &atm = **itAtm;
        //                    //            atm = (*itAtm).second;
        //                    if (atm.getParent().getMoleType() != MoleType::PROTEIN) continue;
        //                    if (atm.isHydrogen() || atm.getMOL2Type()=="Du" || !atm.props.isHydrophobic() ) continue;
        //                    Res= atm.getResidu()->getNum();
        //                    if (keepRes==100000){
        //                        keepRes=Res;
        //                    }
        //                    if (Res != keepRes){
        //                        PHAlist[itPHA]=true;
        //                        //                        cout << "on conserve" << endl;
        //                        break;
        //                    }else {
        //                        PHAlist[itPHA]=false;
        //                    }
        //                }
        //            }
        //        }

        // Fin mod
    }else if (mod == 0  || mod == 8 || mod == 9  || mod == 10 || mod == 19 || mod == 20) {
        for ( int itPHA =0 ;itPHA<cbnum;itPHA++)
        {
            if (atomlist[itPHA] == (Atom*)NULL ) { continue;}

            Coords vectorD = atomlist[itPHA]->fixpos - boxliste[itPHA]->fixpos;
            double norme = sqrt(pow(vectorD.x,2)+pow(vectorD.y,2)+pow(vectorD.z,2));
            Coords Vector1 = Coords(vectorD.x/norme,vectorD.y/norme,vectorD.z/norme);
            double norme2 = sqrt(pow(Vector1.x,2)+pow(Vector1.y,2)+pow(Vector1.z,2));
            Coords pos1 = boxliste[itPHA]->fixpos+Vector1;
            Coords pos2 = pos1+Vector1;
            Coords pos3 = pos2+Vector1;


            //            if (mod == 8){
            //                if ((double)projs[boxliste[itPHA]->getId()] < 70) { continue;}
            //            }else
            if (mod == 9 || mod == 20){
                if ((double)projs[boxliste[itPHA]->getId()] < 80) { continue;}
            }else if (mod == 10){
                if ((double)projs[boxliste[itPHA]->getId()] < 90) { continue;}
            }
            if (mod == 19 || mod == 20){
                if ( boxliste[itPHA]->fixpos.calcDist(grid.getcenter()) > 8) { continue;}
            }



            if ( Arolist[itPHA] == true)    {
                //                cout << "arommm" << endl;
                molen.addAtom("C",boxliste[itPHA]->fixpos,"CZ","C.ar", &resPHE);
                grid.getAdjacentAtoms(listAtoms,*(boxliste[itPHA]),10);
                int max=0;
                Coords centroid;
                int weight = 0;
                for (ItCAtom itAtm= listAtoms.begin(); itAtm != listAtoms.end(); ++itAtm)
                {
                    Atom &atm = **itAtm;
                    if (atm.props.isAromatic() && atm.getResidu() == atomlist[itPHA]->getResidu()){
                        if (boxliste[itPHA]->fixpos.calcDist(atm.fixpos)<8){
                            centroid.setCoords( (centroid.x+atm.fixpos.x),(centroid.y+atm.fixpos.y),(centroid.z+atm.fixpos.z));
                            weight++;
                            //                                cout << "Add \t"<< atm.fixpos.toString() << endl;
                        }
                    }
                }
                centroid.setCoords(centroid.x/weight,centroid.y/weight,centroid.z/weight); // On se sert du centroïde pour les aromatiques
                vectorD = centroid - boxliste[itPHA]->fixpos;
                norme = sqrt(pow(vectorD.x,2)+pow(vectorD.y,2)+pow(vectorD.z,2));
                Vector1 = Coords(vectorD.x/norme,vectorD.y/norme,vectorD.z/norme);
                pos1 = boxliste[itPHA]->fixpos+Vector1;
                pos3 = boxliste[itPHA]->fixpos+Vector1+Vector1+Vector1;
                //                molen.addAtom("P",centroid,"P","P",&resPHE).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                PharmProp *newPharm = new PharmProp(PharType::AROMATIC,boxliste[itPHA]->fixpos,pos1,pos3,0);
                MonPharma.push_back(newPharm);


            }
            if(atomlist[itPHA]->props.isMetal()){
                molen.addAtom("O",boxliste[itPHA]->fixpos,"Zn","O.co2", &resTRP).setPartialCharge(-1);
                PharmProp* newPharm = new PharmProp(PharType::METAL,boxliste[itPHA]->fixpos,pos1,pos3,0);
                MonPharma.push_back(newPharm);
            }
            else if (atomlist[itPHA]->props.isCation())
            {
                //                molen.addAtom("O",boxliste[itPHA]->fixpos,"OD1","O.3", &resASP).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                molen.addAtom("O",boxliste[itPHA]->fixpos,"OD1","O.co2", &resASP).setPartialCharge(-1);
                PharmProp* newPharm = new PharmProp(PharType::HBACCEPTOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                MonPharma.push_back(newPharm);

            }
            else if (atomlist[itPHA]->props.isAnion())                                  {
                //                molen.addAtom("N",boxliste[itPHA]->fixpos,"NZ","N.2" , &resLYS).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                molen.addAtom("N",boxliste[itPHA]->fixpos,"NZ","N.4" , &resLYS).setPartialCharge(1);
                PharmProp* newPharm = new PharmProp(PharType::HBDONOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                MonPharma.push_back(newPharm);
            }
            else if (atomlist[itPHA]->props.isDonor() && atomlist[itPHA]->props.isAcceptor())
            {
                if (atomlist[itPHA]->getResiduName() == "HOH"){
                    bool d_ok = false;
                    bool a_ok = false;
                    for (size_t i=0; i< atomlist[itPHA]->getNumBond();++i)
                    {
                        const Atom &atm2 = atomlist[itPHA]->getAtomLinked(i);
                        if (!atm2.isHydrogen())     continue;
                        double angle_d = atm2.fixpos.calcAngle(atomlist[itPHA]->fixpos,  boxliste[itPHA]->fixpos);
                        double angle_a = atomlist[itPHA]->fixpos.calcAngle(atm2.fixpos,  boxliste[itPHA]->fixpos);

                        if (angle_d >= Angl_H-AngT_H && angle_d <= Angl_H+AngT_H )  {
                            d_ok = true;
                        }else if (angle_a >= D120-AngT_H && angle_a <= D120+AngT_H )  {
                            a_ok = true;
                            cout << "angle accepteur ok entre : "<< atomlist[itPHA]->getIdentifier() << " || " << atm2.getIdentifier() << " || " << boxliste[itPHA]->getId() << endl;
                        }
                    }
                    if (d_ok){
                        molen.addAtom("O",boxliste[itPHA]->fixpos,"O","O.2"  , &resALA);
                        PharmProp* newPharm = new PharmProp(PharType::HBDONOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                        MonPharma.push_back(newPharm);

                    }
                    else if (a_ok){
                        molen.addAtom("N",boxliste[itPHA]->fixpos,"N","N.am"  , &resALA);
                        PharmProp* newPharm2 = new PharmProp(PharType::HBACCEPTOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                        MonPharma.push_back(newPharm2);
                    }
                }
                else if (atomlist[itPHA]->getAtomicName() == "O")
                {
                    bool d_ok = false;
                    bool a_ok = false;
                    for (size_t i=0; i< atomlist[itPHA]->getNumBond();++i)
                    {
                        const Atom &atm2 = atomlist[itPHA]->getAtomLinked(i);
                        if (!atm2.isHydrogen())     continue;
                        double angle_d = atm2.fixpos.calcAngle(atomlist[itPHA]->fixpos,  boxliste[itPHA]->fixpos);
                        double angle_a = atomlist[itPHA]->fixpos.calcAngle(atm2.fixpos,  boxliste[itPHA]->fixpos);
                        if (angle_d >= D180-AngT_H && angle_d <= D180+AngT_H )  {
                            d_ok = true;
                        }else if (angle_a >= D120-AngT_H && angle_a <= D120+AngT_H )  {
                            cout << "angle accepteur ok entre : "<< atomlist[itPHA]->getIdentifier() << " || " << atm2.getIdentifier() << " || " << boxliste[itPHA]->getId() << endl;
                            a_ok = true;
                        }
                    }
                    if (d_ok){
                        molen.addAtom("O",boxliste[itPHA]->fixpos,"O","O.2"  , &resALA);
                        PharmProp* newPharm = new PharmProp(PharType::HBDONOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                        MonPharma.push_back(newPharm);

                    }else if (a_ok){
                        molen.addAtom("N",boxliste[itPHA]->fixpos,"N","N.am"  , &resALA);
                        PharmProp* newPharm2 = new PharmProp(PharType::HBACCEPTOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                        MonPharma.push_back(newPharm2);
                    }
                    else
                    {
                        molen.addAtom("O",boxliste[itPHA]->fixpos,"OG","O.3" , &resSER);
                        //                molen.addAtom("P",atomlist[itPHA]->fixpos,"P","P",&resALA).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);  // Creer le vecteur en faisant la différence des coordonées
                        PharmProp* newPharm = new PharmProp(PharType::HBDONOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                        MonPharma.push_back(newPharm);
                        PharmProp* newPharm2 = new PharmProp(PharType::HBACCEPTOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                        MonPharma.push_back(newPharm2);
                    }


                }
                else
                {
                    bool d_ok = false;
                    bool a_ok = false;
                    for (size_t i=0; i< atomlist[itPHA]->getNumBond();++i)
                    {
                        const Atom &atm2 = atomlist[itPHA]->getAtomLinked(i);
                        if (!atm2.isHydrogen())     continue;
                        double angle_d = atm2.fixpos.calcAngle(atomlist[itPHA]->fixpos,  boxliste[itPHA]->fixpos);
                        double angle_a = atomlist[itPHA]->fixpos.calcAngle(atm2.fixpos,  boxliste[itPHA]->fixpos);
                        if (angle_d >= D180-AngT_H && angle_d <= D180+AngT_H )  {
                            d_ok = true;

                        }else if (angle_a >= D120-AngT_H && angle_a <= D120+AngT_H )  {
                            cout << "angle accepteur ok entre : "<< atomlist[itPHA]->getIdentifier() << " || " << atm2.getIdentifier() << " || " << boxliste[itPHA]->getId() << endl;
                            a_ok = true;
                        }
                    }
                    if (d_ok){
                        molen.addAtom("O",boxliste[itPHA]->fixpos,"O","O.2"  , &resALA);
                        PharmProp* newPharm = new PharmProp(PharType::HBDONOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                        MonPharma.push_back(newPharm);

                    }else if (a_ok){
                        molen.addAtom("N",boxliste[itPHA]->fixpos,"N","N.am"  , &resALA);
                        PharmProp* newPharm2 = new PharmProp(PharType::HBACCEPTOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                        MonPharma.push_back(newPharm2);
                    }
                    else
                    {
                        molen.addAtom("O",boxliste[itPHA]->fixpos,"OG","O.3" , &resSER);
                        //                molen.addAtom("P",atomlist[itPHA]->fixpos,"P","P",&resALA).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);  // Creer le vecteur en faisant la différence des coordonées
                        PharmProp* newPharm = new PharmProp(PharType::HBDONOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                        MonPharma.push_back(newPharm);
                        PharmProp* newPharm2 = new PharmProp(PharType::HBACCEPTOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                        MonPharma.push_back(newPharm2);
                    }

                }
            }
            else if (atomlist[itPHA]->props.isAcceptor())                               {
                molen.addAtom("N",boxliste[itPHA]->fixpos,"N","N.am"  , &resALA);
                //                molen.addAtom("P",atomlist[itPHA]->fixpos,"P","P",&resALA).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);  // Creer le vecteur en faisant la différence des coordonées
                PharmProp* newPharm = new PharmProp(PharType::HBDONOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                MonPharma.push_back(newPharm);
            }
            else if (atomlist[itPHA]->props.isDonor())                                  {
                molen.addAtom("O",boxliste[itPHA]->fixpos,"O","O.2"  , &resALA);
                //                molen.addAtom("P",atomlist[itPHA]->fixpos,"P","P",&resALA).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);  // Creer le vecteur en faisant la différence des coordonées
                PharmProp* newPharm = new PharmProp(PharType::HBACCEPTOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                MonPharma.push_back(newPharm);

            }

            else if (atomlist[itPHA]->props.isHydrophobic())                            {
                molen.addAtom("C",boxliste[itPHA]->fixpos,"CA","C.3" , &resGLY);
                PharmProp *newPharm = new PharmProp(PharType::HYDROPHOBIC,boxliste[itPHA]->fixpos,pos1,pos3,0);
                MonPharma.push_back(newPharm);
            }
            //            else {cout << "ATOM NOT RECOGNIZED "<< atomlist[itPHA]->toString()<<endl;continue;}

        }

    }


    if (mod == 8 || mod == 22 || mod == 23 || mod == 24 || mod == 25  || mod == 26  || mod == 27) {
        vector<int> Resused;
        vector<int> in;
        vector<int> close;
        double plp_atm =0;
        for ( int itPHA =0 ;itPHA<cbnum;itPHA++)
        {
            if (PHAlist[itPHA]){
                //                Coords vectorD = Coords(atomlist[itPHA]->fixpos.x-coorlist[itPHA].x,atomlist[itPHA]->fixpos.y-coorlist[itPHA].y,atomlist[itPHA]->fixpos.z-coorlist[itPHA].z) ;


                float plpX,plpY,plpZ;
                double best_plp = 10000;
                Coords bestcoord;
//                cout << "Pharmacophor : " << itPHA << ", prop inverse: " << atomlist[itPHA]->props.toString() << endl;
                for (plpX=-0.5;plpX < 0.6; plpX+=0.5){
                    for (plpY=-0.5;plpY < 0.6; plpY+=0.5){
                        for (plpZ=-0.5;plpZ < 0.6; plpZ+=0.5){
                            Coords localbox = coorlist[itPHA]+Coords(plpX,plpY,plpZ);
                            Atom &atomL= *atomlist[itPHA];
                            double plp,plpA,plpB,plpC,plpD,plpE;
                            double dist;
                            plp_atm =0;
                            if (atomL.isHydrogen() || (atomL.getName()=="DuCy") || !atomL.isUsed()) continue;
                            //                if (atomL.getBox(&grid) == (Box*)NULL) continue;
                            AtomList atmlist;
                            grid.getAdjacentAtoms(atmlist,*atomlist[itPHA],6);
                            for (ItCAtom itA = atmlist.begin(); itA != atmlist.end(); itA++){
                                Atom &atomP = **itA;
                                dist = 0;
                                if (&atomP.getParent() == &ligand || !atomP.isUsed() || atomP.isHydrogen())continue;
                                dist = localbox.calcDist(atomP.fixpos);
                                plp=0;plpA=0;plpB=0;plpC=0;plpD=0;plpE=0;
                                if( atomL.props.isMetal() && atomP.props.isMetal()){
                                    plpA = 1.4; plpB = 2.2; plpC = 3.1; plpD = 3.4;plpE = -5;
                                }
                                else if ( ((atomL.props.isAcceptor() || atomL.props.isAnion()) && (atomP.props.isAcceptor() || atomP.props.isAnion()))
                                          || ((atomL.props.isDonor() || atomL.props.isCation()) && (atomP.props.isDonor() || atomP.props.isCation()))) {
                                    plpA = 2.3; plpB = 2.6; plpC = 3.1; plpD = 3.4; plpE = -2;

                                }
                                else if (atomL.props.isApolar() && atomP.props.isApolar() ){
                                    plpA = 3.4; plpB = 3.6; plpC = 4.5; plpD = 5.5; plpE = -0.4;
                                }
                                if (dist < plpA && plpA != 0){
                                    plp = 20*(plpA-dist)/plpA;
                                }
                                else if (dist >= plpA && dist < plpB){
                                    plp = plpE*(dist-plpA)/(plpB-plpA);
                                }
                                else if (dist >= plpB && dist < plpC){
                                    plp = plpE;
                                }
                                else if (dist >= plpC && dist < plpD){
                                    plp = plpE*(plpD-dist)/(plpD-plpC);
                                }
                                else if (dist >= plpD){
                                    plp = 0;
                                }
                                plp_atm += plp;
                                if (plp != 0) {
                                cout << " "<<plpX<<" " <<plpY<<" " <<plpZ<<"\t : " << atomP.getIdentifier() << " | "<< atomP.props.toString() << " ||\t  dist : " << dist << "||\t plp : " << plp << endl;
                                }
                            }
//                            cout << "\t|| Somme de la pose : " << plp_atm << endl;
                            if (best_plp > plp_atm || (best_plp == plp_atm && (plpX ==0 && plpY ==0 && plpZ ==0))){
                                best_plp = plp_atm;
                                bestcoord = localbox;
                            }
                        }
                    }



                }

                cout << "prop inverse " << atomlist[itPHA]->props.toString() << " || plp : " << best_plp << endl;
                coorlist[itPHA] = bestcoord;

                Coords vectorD = atomlist[itPHA]->fixpos - coorlist[itPHA];
                double norme = sqrt(pow(vectorD.x,2)+pow(vectorD.y,2)+pow(vectorD.z,2));
                Coords Vector1 = Coords(vectorD.x/norme,vectorD.y/norme,vectorD.z/norme);
                Coords pos1 = coorlist[itPHA]+Vector1;
                Coords pos2 = pos1+Vector1;
                Coords pos3 = pos2+Vector1;
                double tolS = 0;
                int tolEx = 100;
                if (tolerance){
                    if (tol[itPHA] == 1){
                        tolS -= 30;
                    }else if (tol[itPHA]> 10){
                        tolS = 30;
                    }
                }

                AtomList atmlist;
                grid.getAdjacentAtoms(atmlist,*atomlist[itPHA],5);
                for (ItCAtom itA = atmlist.begin(); itA != atmlist.end(); itA++){
                    Atom &atmR = **itA;

                    if (atmR.getParent().getMoleType() != MoleType::PROTEIN) { continue;}
                    if (atmR.isHydrogen()) {
                        hydro_present = 1;
                    }
                    if (find(Resused.begin(), Resused.end(), atmR.getResidu()->getNum()) == Resused.end()) {
                        Coords cent;
                        int nbatm = 0;
                        cent.setCoords(Coords(0,0,0));
                        for (ItCAtom itB = atmR.getResidu()->firstAtom();itB != atmR.getResidu()->lastAtom(); itB++ )
                        {
                            Atom &atm = **itB;
                            bool exclu[2] = {false};
                            for ( int itPHAE =0 ;itPHAE<cbnum;itPHAE++){
                                if ( atm.fixpos.calcDist(coorlist[itPHAE]) > 4 && atm.fixpos.calcDist(coorlist[itPHAE]) < 5 && atm.fixpos != atomlist[itPHA]->fixpos){
                                    exclu[0]=true;
                                }else if ( atm.fixpos.calcDist(coorlist[itPHAE]) < 4 && atm.fixpos != atomlist[itPHA]->fixpos){
                                    exclu[1]=true;
                                }

                            }
                            if (exclu[0] && !exclu[1]){
                                //                                cout << atm.getIdentifier() << endl; ;
                                Resused.push_back(atm.getResidu()->getNum());
                                cent.setCoords(cent+atm.fixpos);
                                nbatm++;
                            }
                        }
                        if (nbatm!=0){
                            cent.setCoords(cent/nbatm);
                            switch (nbatm)
                            {
                            case 0: tolEx = 100; break;
                            case 1: tolEx = 115; break;
                            case 2: tolEx = 125; break;
                            case 3: tolEx = 135; break;
                            case 4: tolEx = 145; break;
                            case 5: tolEx = 150; break;
                            case 6: tolEx = 155; break;
                            case 7: tolEx = 160; break;
                            case 8: tolEx = 165; break;
                            default: tolEx = 170;break;
                            }
                            if (exclu) {
                                PharmProp* newPharm = new PharmProp(PharType::EXCLUSION,cent,Coords(0,0,0),Coords(0,0,0),tolEx);
                                molen.addAtom("H",cent,"Du","H",&resCUB);
                                MonPharma.push_back(newPharm);
                            }
                        }

                    }

                }

                if (Arolist[itPHA] == true)    {
//                    if (best_plp > -1.6) { continue;}
                    //                cout << "arommm" << endl;
                    //                    cout << "here" <<endl;
                    //                    cout << "res : " << atomlist[itPHA]->props.toString() << endl;
                    //                    cout << "init : "<< boxliste[itPHA]->fixpos.toString() << endl;
                    AtomList listAtoms;
                    molen.addAtom("C",coorlist[itPHA],"CZ","C.ar", &resPHE);
                    grid.getAdjacentAtoms(listAtoms,*(boxliste[itPHA]),10);
                    int max=0;
                    Coords centroid;
                    int weight = 0;
                    for (ItCAtom itAtm= listAtoms.begin(); itAtm != listAtoms.end(); ++itAtm)
                    {
                        Atom &atm = **itAtm;
                        if (atm.props.isAromatic() && atm.getResidu() == atomlist[itPHA]->getResidu()){
                            if (coorlist[itPHA].calcDist(atm.fixpos)<8){
                                centroid.setCoords( (centroid.x+atm.fixpos.x),(centroid.y+atm.fixpos.y),(centroid.z+atm.fixpos.z));
                                weight++;
                                //                                cout << weight << endl;
                                //                                                                cout << "Add \t"<< atm.fixpos.toString() << endl;
                            }
                        }
                    }
                    centroid.setCoords(centroid.x/weight,centroid.y/weight,centroid.z/weight); // On se sert du centroïde pour les aromatiques
                    vectorD = centroid - coorlist[itPHA];
                    norme = sqrt(pow(vectorD.x,2)+pow(vectorD.y,2)+pow(vectorD.z,2));
                    Vector1 = Coords(vectorD.x/norme,vectorD.y/norme,vectorD.z/norme);
                    pos1 = coorlist[itPHA]+Vector1;
                    pos3 = coorlist[itPHA]+Vector1+Vector1+Vector1;
                    //                                    molen.addAtom("P",centroid,"P","P",&resPHE).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);

                    PharmProp *newPharm = new PharmProp(PharType::AROMATIC,coorlist[itPHA],pos1,pos3,0);
                    MonPharma.push_back(newPharm);
                    //                    cout << "all : "<<centroid.toString() << " " << Vector1.toString() << " " << pos1.toString() <<" " <<pos3.toString() << endl;


                }

                if(atomlist[itPHA]->props.isMetal()){
//                    if (best_plp > 0) { continue;}
                    molen.addAtom("O",coorlist[itPHA],"Zn","O.co2", &resTRP).setPartialCharge(-1);
                    PharmProp* newPharm = new PharmProp(PharType::METAL,coorlist[itPHA],pos1,pos3,0);
                    MonPharma.push_back(newPharm);
                }
                else if (atomlist[itPHA]->props.isCation())                                 {
//                    if (best_plp > 1.9) { continue;}
                    molen.addAtom("O",coorlist[itPHA],"OD1","O.co2", &resASP).setPartialCharge(-1);

                    PharmProp* newPharm = new PharmProp(PharType::HBACCEPTOR,coorlist[itPHA],pos1,pos3,tolS);
                    MonPharma.push_back(newPharm);

                }
                else if (atomlist[itPHA]->props.isAnion())                                  {
//                    if (best_plp > 0.7) { continue;}
                    molen.addAtom("N",coorlist[itPHA],"NZ","N.4" , &resLYS).setPartialCharge(1);
                    PharmProp* newPharm = new PharmProp(PharType::HBDONOR,coorlist[itPHA],pos1,pos3,tolS);
                    MonPharma.push_back(newPharm);

                }
                else if (atomlist[itPHA]->props.isDonor() && atomlist[itPHA]->props.isAcceptor())    {
//                    if (best_plp > 1.9) { continue;}
                    if (atomlist[itPHA]->getResiduName() == "HOH"){
                        bool d_ok = false;
                        bool a_ok = false;
                        for (size_t i=0; i< atomlist[itPHA]->getNumBond();++i)
                        {
                            const Atom &atm2 = atomlist[itPHA]->getAtomLinked(i);
                            if (!atm2.isHydrogen())     continue;
                            double angle_d = atm2.fixpos.calcAngle(atomlist[itPHA]->fixpos,  coorlist[itPHA]);
                            double angle_a = atomlist[itPHA]->fixpos.calcAngle(atm2.fixpos,  coorlist[itPHA]);

                            if (angle_d >= Angl_H-AngT_H && angle_d <= Angl_H+AngT_H )  {
                                d_ok = true;
//                                cout << "angle donneur EAU ok "<< angle_d << " entre : "<< atomlist[itPHA]->getIdentifier() << " || " << atm2.getIdentifier() << " || " << boxliste[itPHA]->getId() << endl;
                            }else if (angle_a >= Angl_H-AngT_H && angle_a <= Angl_H+AngT_H )  {
                                a_ok = true;

                            }
                        }
                        if (d_ok){
                            molen.addAtom("O",coorlist[itPHA],"O","O.2"  , &resALA);
                            PharmProp* newPharm = new PharmProp(PharType::HBDONOR,coorlist[itPHA],pos1,pos3,0);
                            MonPharma.push_back(newPharm);

                        }else if (a_ok){
                            molen.addAtom("N",coorlist[itPHA],"N","N.am"  , &resALA);
                            PharmProp* newPharm2 = new PharmProp(PharType::HBACCEPTOR,coorlist[itPHA],pos1,pos3,0);
                            MonPharma.push_back(newPharm2);
                        }
                    }else if (atomlist[itPHA]->getAtomicName() == "O"){
                        bool d_ok = false;
                        bool a_ok = false;
                        for (size_t i=0; i< atomlist[itPHA]->getNumBond();++i)
                        {
                            const Atom &atm2 = atomlist[itPHA]->getAtomLinked(i);
                            if (!atm2.isHydrogen())     continue;
                            double angle_d = atm2.fixpos.calcAngle(atomlist[itPHA]->fixpos,  coorlist[itPHA]);
                            double angle_a = atomlist[itPHA]->fixpos.calcAngle(atm2.fixpos,  coorlist[itPHA]);
                            if (angle_d >= Angl_H-AngT_H && angle_d <= Angl_H+AngT_H )  {
                                d_ok = true;
//                                cout << Angl_H-AngT_H << " < " << angle_d << " < " << Angl_H+AngT_H << endl;
//                                cout << "angle donneur O   ok "<< angle_d << " entre : "<< atomlist[itPHA]->getIdentifier() << " || " << atm2.getIdentifier() << " || " << boxliste[itPHA]->getId() << endl;
                            }else if (angle_a >= D120-AngT_H && angle_a <= D120+AngT_H )  {
                                a_ok = true;

                            }
                        }
                        if (d_ok){
                            molen.addAtom("O",coorlist[itPHA],"O","O.2"  , &resALA);
                            PharmProp* newPharm = new PharmProp(PharType::HBDONOR,coorlist[itPHA],pos1,pos3,0);
                            MonPharma.push_back(newPharm);

                        }else if (a_ok){
                            molen.addAtom("N",coorlist[itPHA],"N","N.am"  , &resALA);
                            PharmProp* newPharm2 = new PharmProp(PharType::HBACCEPTOR,coorlist[itPHA],pos1,pos3,0);
                            MonPharma.push_back(newPharm2);
                        }
                        else {

                            molen.addAtom("O",coorlist[itPHA],"OG","O.3" , &resSER);
                            //                molen.addAtom("P",atomlist[itPHA]->fixpos,"P","P",&resALA).setPartialCharge((double)projs[boxliste[itPHA]->getId()]);  // Creer le vecteur en faisant la différence des coordonées
                            PharmProp* newPharm = new PharmProp(PharType::HBDONOR,coorlist[itPHA],pos1,pos3,0);
                            MonPharma.push_back(newPharm);
                            PharmProp* newPharm2 = new PharmProp(PharType::HBACCEPTOR,coorlist[itPHA],pos1,pos3,0);
                            MonPharma.push_back(newPharm2);
                        }



                    } else {
                        bool d_ok = false;
                        bool a_ok = false;
                        for (size_t i=0; i< atomlist[itPHA]->getNumBond();++i)
                        {
                            if (hydro_present == 1){
                                const Atom &atm2 = atomlist[itPHA]->getAtomLinked(i);
                                if (!atm2.isHydrogen())     continue;
                                double angle_d = atm2.fixpos.calcAngle(atomlist[itPHA]->fixpos,  boxliste[itPHA]->fixpos);
                                double angle_a = atomlist[itPHA]->fixpos.calcAngle(atm2.fixpos,  boxliste[itPHA]->fixpos);
                                if (angle_d >= Angl_H-AngT_H && angle_d <= Angl_H+AngT_H )  {
                                    d_ok = true;
//                                    cout << "angle donneur N   ok "<< angle_d << " entre : "<< atomlist[itPHA]->getIdentifier() << " || " << atm2.getIdentifier() << " || " << boxliste[itPHA]->getId() << endl;

                                }
                                else if (angle_a >= D120-AngT_H && angle_a <= D120+AngT_H )  {
                                    a_ok = true;

                                }
                                if (d_ok){
                                    molen.addAtom("O",boxliste[itPHA]->fixpos,"O","O.2"  , &resALA);
                                    PharmProp* newPharm = new PharmProp(PharType::HBDONOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                                    MonPharma.push_back(newPharm);

                                }else {
                                    molen.addAtom("N",coorlist[itPHA],"N","N.am"  , &resALA);
                                    PharmProp* newPharm2 = new PharmProp(PharType::HBACCEPTOR,coorlist[itPHA],pos1,pos3,0);
                                    MonPharma.push_back(newPharm2);
                                }
                            }
                            else
                            {
                                molen.addAtom("O",boxliste[itPHA]->fixpos,"OG","O.3" , &resSER);
                                //                molen.addAtom("P",atomlist[itPHA]->fixpos,"P","P",&resALA).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);  // Creer le vecteur en faisant la différence des coordonées
                                PharmProp* newPharm = new PharmProp(PharType::HBDONOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                                MonPharma.push_back(newPharm);
                                PharmProp* newPharm2 = new PharmProp(PharType::HBACCEPTOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                                MonPharma.push_back(newPharm2);
                            }
                        }


                    }
                }
                //                else if (atomlist[itPHA]->props.isDonor() && atomlist[itPHA]->props.isAcceptor())    {
                //                    molen.addAtom("O",coorlist[itPHA],"OG","O.2" , &resSER).setPartialCharge((double)projs[boxliste[itPHA]->getId()]);
                //                    molen.addAtom("P",atomlist[itPHA]->fixpos,"P","P",&resALA).setPartialCharge((double)projs[boxliste[itPHA]->getId()]);
                //                    PharmProp* newPharm = new PharmProp(PharType::HBDONOR,coorlist[itPHA],pos1,pos3,tolS);
                //                    MonPharma.push_back(newPharm);
                //                    PharmProp* newPharm2 = new PharmProp(PharType::HBACCEPTOR,coorlist[itPHA],pos1,pos3,tolS);
                //                    MonPharma.push_back(newPharm2);

                //                }
                else if (atomlist[itPHA]->props.isAcceptor())                               {
//                    if (best_plp > 0.7) { continue;}
                    molen.addAtom("N",coorlist[itPHA],"N","N.am"  , &resALA);
                    //                    molen.addAtom("P",atomlist[itPHA]->fixpos,"P","P",&resALA).setPartialCharge((double)projs[boxliste[itPHA]->getId()]);
                    PharmProp* newPharm = new PharmProp(PharType::HBDONOR,coorlist[itPHA],pos1,pos3,tolS);
                    MonPharma.push_back(newPharm);

                }
                else if (atomlist[itPHA]->props.isDonor())                                  {

//                    if (best_plp > 1.9) { continue;}
                    molen.addAtom("O",coorlist[itPHA],"O","O.2"  , &resALA);
                    //                    molen.addAtom("P",atomlist[itPHA]->fixpos,"P","P",&resALA).setPartialCharge((double)projs[boxliste[itPHA]->getId()]);

                    PharmProp* newPharm = new PharmProp(PharType::HBACCEPTOR,coorlist[itPHA],pos1,pos3,tolS);
                    MonPharma.push_back(newPharm);

                }
                //                else if (atomlist[itPHA]->props.isAromatic())                               {
                //                    molen.addAtom("C",coorlist[itPHA],"CZ","C.ar", &resPHE).setPartialCharge((double)projs[boxliste[itPHA]->getId()]);
                //                    grid.getAdjacentAtoms(listAtoms,*(boxliste[itPHA]),10);
                //                    int max=0;
                //                    Coords centroid;
                //                    int weight = 0;
                //                    for (ItCAtom itAtm= listAtoms.begin(); itAtm != listAtoms.end(); ++itAtm)
                //                    {
                //                        Atom &atm = **itAtm;
                //                        if (atm.props.isAromatic() && atm.getResidu() == atomlist[itPHA]->getResidu()){
                //                            if (boxliste[itPHA]->fixpos.calcDist(atm.fixpos)<8){
                //                                centroid.setCoords( (centroid.x+atm.fixpos.x),(centroid.y+atm.fixpos.y),(centroid.z+atm.fixpos.z));
                //                                weight++;
                //                            }
                //                        }
                //                    }
                //                    centroid.setCoords(centroid.x/weight,centroid.y/weight,centroid.z/weight); // On se sert du centroïde pour les aromatiques
                //                    vectorD = centroid - coorlist[itPHA];
                //                    norme = sqrt(pow(vectorD.x,2)+pow(vectorD.y,2)+pow(vectorD.z,2));
                //                    Vector1 = Coords(vectorD.x/norme,vectorD.y/norme,vectorD.z/norme);
                //                    pos1 = coorlist[itPHA]+Vector1;
                //                    pos3 = coorlist[itPHA]+Vector1+Vector1+Vector1;
                //                    //                    molen.addAtom("P",centroid,"P","P",&resPHE).setPartialCharge((double)projs[boxliste[itPHA]->getId()]);
                //                    PharmProp *newPharm = new PharmProp(PharType::AROMATIC,coorlist[itPHA],pos1,pos3,tolS);
                //                    MonPharma.push_back(newPharm);

                //                }
                else if (atomlist[itPHA]->props.isHydrophobic())                            {
//                    if (best_plp > -1.6) { continue;}
                    molen.addAtom("C",coorlist[itPHA],"CA","C.3" , &resGLY);
                    PharmProp *newPharm = new PharmProp(PharType::HYDROPHOBIC,coorlist[itPHA],pos1,pos3,tolS);
                    MonPharma.push_back(newPharm);
                    //                                        PharmProp* newPharm2 = new PharmProp(PharType::EXCLUSION,atomlist[itPHA]->fixpos,Coords(0,0,0),Coords(0,0,0),tolS);
                    //                                        MonPharma.push_back(newPharm2);

                }
                //                else {cout << "ATOM NOT RECOGNIZED Cube"<< atomlist[itPHA]->toString()<<endl;continue;}

            }
        }
    }

}


void VolSite::cavToPHA(Molecule& molen, PharmList& MonPharma , vector<const Box*> &boxlist, const string& molename, const string& size, const int mod, int& bestb) const
{
    const Box *boxliste[10000]  = {NULL} ;
    Atom *atomlist[10000] = {NULL};
    float distlist[10000]= {0};
    bool Usedlist[10000] = {false};
    bool PHAlist[10000] ={false};
    Coords coorlist[10000];
    int tol[10000]={0};
    int nb_enfoui[120] = {0};
    int cbnum = 0, bestenfoui=0;

    int MoleType=0;
    if (size == "4") MoleType=MoleType::CAV_4;
    else if (size=="6")MoleType=MoleType::CAV_6;
    else if (size=="8")MoleType=MoleType::CAV_8;
    else if (size=="12")MoleType=MoleType::CAV_12;
    else if (size=="ALL")MoleType=MoleType::CAV_ALL;

    //    Molecule molen= new Molecule(MoleType);
    molen.setName(molename);
    molen.clear();
    Residu &resSER = molen.addResidu("X",1,"SER");
    Residu &resALA = molen.addResidu("X",2,"ALA");
    Residu &resASP = molen.addResidu("X",3,"ASP");
    Residu &resLYS = molen.addResidu("X",4,"LYS");
    Residu &resPHE = molen.addResidu("X",5,"PHE");
    Residu &resGLY = molen.addResidu("X",6,"GLY");
    Residu &resCUB = molen.addResidu("X",8,"CUB");

    //    Box* bx=(Box*)NULL;
    int cube_pos;
    Atom * bestatm = (Atom*)NULL;
    double bestdist;

    ///% reviens here
    AtomList listAtoms;

    for (vector<const Box*>::iterator itBx = boxlist.begin(); itBx != boxlist.end(); itBx++)
    {
        const Box &bx = **itBx;
        cube_pos=bx.getId();
        bestatm=(Atom*)NULL;
#ifdef ICHEM_DEBUG
        cout << bx.getId()<<" ";
#endif

        bestdist = 1000;
        bestatm= (Atom*)NULL;

        if (!bx.Props.isInCavity())continue;

        grid.getAdjacentAtoms(listAtoms,bx,5);

        for (ItCAtom itAtm= listAtoms.begin(); itAtm != listAtoms.end(); ++itAtm)
        {
            Atom &atm = **itAtm;
            //            atm = (*itAtm).second;
            if (atm.getParent().getMoleType() != MoleType::PROTEIN) continue;
            if (atm.isHydrogen()
                    || atm.getMOL2Type()=="Du"
                    || (!atm.props.isAromatic()
                        && !atm.props.isAcceptor()
                        && !atm.props.isDonor()
                        && !atm.props.isCation()
                        && !atm.props.isAnion()
                        && !atm.props.isHydrophobic()
                        && !atm.props.isMetal()        ) ) continue;
            // Calcul de distance
            const double dist=atm.fixpos.calcDist(bx.fixpos);

            if (dist > bestdist) continue;

            bestdist=dist;
            bestatm=&atm;
        }

#ifdef ICHEM_DEBUG
        cout << bestatm<<" "<< bestdist<<endl;
#endif
        boxliste[cbnum]  = &bx;
        atomlist[cbnum] = bestatm;
        distlist[cbnum] = bestdist;
        cbnum ++;
        for (int enf=0; enf <projs[bx.getId()];enf ++)
        {
            nb_enfoui[enf]++;
        }
    }
    int pourcent=0,ccube=0;
    for (int enf=0; enf<120;enf++)
    {
        if (( nb_enfoui[enf] <cavitykeepC || 100*nb_enfoui[enf]/ cbnum < cavitykeepP) && bestenfoui==0){
            bestenfoui = enf;
            pourcent= 100*nb_enfoui[enf]/ cbnum;
            ccube=nb_enfoui[enf];
        }
    }
//    if (mod >2){
//        if (bestenfoui > 55){
//            cout << "best b : -b" << bestenfoui<<" || "<< ccube <<" Cavity cubes || Initial CAV : "<< pourcent <<"%"<< endl;
//        }else
//            cout <<"best b : default" << endl;
//    }

    for ( int itPHA =0 ;itPHA<cbnum;itPHA++){
        if (atomlist[itPHA] != (Atom*)NULL && !Usedlist[itPHA] ){
            //            bx.Props=bestatm->props;
            //#ifdef ICHEM_DEBUG
            //            cout << "PROPS::"<<cube_pos<<"::"<<bx.rotpos.x<<":"<<bx.rotpos.y<<":"<<bx.rotpos.z<<"=>""=>("<<bx.fixpos.x<<":"<<bx.fixpos.y<<":"<<bx.fixpos.z<<")=>"<<bestatm->fixpos.calcDist(bx.fixpos)<<"=>"<<bestatm->getIdentifier()<<"\t"<< bestatm->props.toString()<<endl;
            //#endif
        }
        else continue;
        if (coorlist[itPHA].calcDist(Coords(0,0,0)) ==0) { coorlist[itPHA]=boxliste[itPHA]->fixpos;}

        // Calcul for the cube point

        ///////////////////////////////////// Comment here
        if (mod ==1)
        {
            PHAlist[itPHA]=true;
            if (!atomlist[itPHA]->props.isHydrophobic() && !atomlist[itPHA]->props.isAromatic())
            {
                for ( int itPHA2 =itPHA ;itPHA2<cbnum;itPHA2++){
                    if (atomlist[itPHA]->getNum() == atomlist[itPHA2]->getNum()) {
                        //                    cout << "doneur/accepteur Atom : " << itPHA <<" " <<itPHA2 <<endl;
                        //                    cout << "\t\t\t\tdistance vs : "<<  distlist[itPHA] << " vs "<<  distlist[itPHA2]<< endl;
                        if (distlist[itPHA] > distlist[itPHA2]){
                            PHAlist[itPHA2]=false;
                            PHAlist[itPHA]=false;
                            break;
                        }else { PHAlist[itPHA]= true;
                            Usedlist[itPHA2]  = true;
                        }
                    }
                }
            }else if( atomlist[itPHA]->props.isAromatic())
            {
                for ( int itPHA2 =itPHA ;itPHA2<cbnum;itPHA2++){
                    if (atomlist[itPHA]->getNum() == atomlist[itPHA2]->getNum() && atomlist[itPHA2]->props.isAromatic()) {
                        //                    cout << "doneur/accepteur Atom : " << itPHA <<" " <<itPHA2 <<endl;
                        //                    cout << "\t\t\t\tdistance vs : "<<  distlist[itPHA] << " vs "<<  distlist[itPHA2]<< endl;
                        if (distlist[itPHA] > distlist[itPHA2]){
                            PHAlist[itPHA2]=false;
                            PHAlist[itPHA]=false;
                            break;
                        }else { PHAlist[itPHA]= true;
                            Usedlist[itPHA2]  = true;
                        }
                    }
                }
            }else if( atomlist[itPHA]->props.isHydrophobic())
            {
                for ( int itPHA2 =itPHA ;itPHA2<cbnum;itPHA2++){
                    if (atomlist[itPHA2]->props.isHydrophobic() && !Usedlist[itPHA2] )
                    {
                        if (atomlist[itPHA]->getResidu()->getNum() == atomlist[itPHA2]->getResidu()->getNum()) {
                            //                        cout << "Hydro Residu : " << atomlist[itPHA]->getResidu()->getNum() <<" " <<atomlist[itPHA2]->getResidu()->getNum() <<endl;
                            //                        cout << "\t\t\t\tdistance vs : "<<  distlist[itPHA] << " vs "<<  distlist[itPHA2]<< endl;
                            if (distlist[itPHA]> distlist[itPHA2]){
                                PHAlist[itPHA2]=true;
                                PHAlist[itPHA]=false;
                                break;
                            }else { PHAlist[itPHA]= true;
                                Usedlist[itPHA2]  = true;

                            }
                        }
                    }
                }
                Usedlist[itPHA]= true;
            }
        }else if (mod ==4){
            if (projs[boxliste[itPHA]->getId()]<bestenfoui) {  continue;}
            PHAlist[itPHA]=true;
            if (!atomlist[itPHA]->props.isHydrophobic() && !atomlist[itPHA]->props.isAromatic())
            {
                for ( int itPHA2 =itPHA ;itPHA2<cbnum;itPHA2++){
                    if (atomlist[itPHA]->getNum() == atomlist[itPHA2]->getNum()) {
                        //                    cout << "doneur/accepteur Atom : " << itPHA <<" " <<itPHA2 <<endl;
                        //                    cout << "\t\t\t\tdistance vs : "<<  distlist[itPHA] << " vs "<<  distlist[itPHA2]<< endl;
                        if (distlist[itPHA] > distlist[itPHA2]){
                            PHAlist[itPHA2]=false;
                            PHAlist[itPHA]=false;
                            break;
                        }else { PHAlist[itPHA]= true;
                            Usedlist[itPHA2]  = true;
                        }
                    }
                }
            }else if( atomlist[itPHA]->props.isAromatic())
            {
                for ( int itPHA2 =itPHA ;itPHA2<cbnum;itPHA2++){
                    if (atomlist[itPHA]->getResidu()->getNum() == atomlist[itPHA2]->getResidu()->getNum() && atomlist[itPHA2]->props.isAromatic()) {
                        //                    cout << "doneur/accepteur Atom : " << itPHA <<" " <<itPHA2 <<endl;
                        //                    cout << "\t\t\t\tdistance vs : "<<  distlist[itPHA] << " vs "<<  distlist[itPHA2]<< endl;
                        if (distlist[itPHA] > distlist[itPHA2]){
                            PHAlist[itPHA2]=false;
                            PHAlist[itPHA]=false;
                            break;
                        }else { PHAlist[itPHA]= true;
                            Usedlist[itPHA2]  = true;
                        }
                    }
                }
            }else if( atomlist[itPHA]->props.isHydrophobic())
            {
                for ( int itPHA2 =itPHA ;itPHA2<cbnum;itPHA2++){
                    if (atomlist[itPHA2]->props.isHydrophobic() && !Usedlist[itPHA2] )
                    {
                        if (atomlist[itPHA]->getResidu()->getNum() == atomlist[itPHA2]->getResidu()->getNum()) {
                            //                        cout << "Hydro Residu : " << atomlist[itPHA]->getResidu()->getNum() <<" " <<atomlist[itPHA2]->getResidu()->getNum() <<endl;
                            //                        cout << "\t\t\t\tdistance vs : "<<  distlist[itPHA] << " vs "<<  distlist[itPHA2]<< endl;
                            if (distlist[itPHA]> distlist[itPHA2]){
                                PHAlist[itPHA2]=true;
                                PHAlist[itPHA]=false;
                                break;
                            }else { PHAlist[itPHA]= true;
                                Usedlist[itPHA2]  = true;

                            }
                        }
                    }
                }
                Usedlist[itPHA]= true;
            }


        }else if (mod == 2 || mod == 8 ){
            ////////////////////////////// to here
            // Calcul for the median
            PHAlist[itPHA]=true;
            if (!atomlist[itPHA]->props.isHydrophobic() && !atomlist[itPHA]->props.isAromatic())
            {
                int weight =1;
                //                cout << "coord\t" << coorlist[itPHA].toString() << endl;
                for ( int itPHA2 =itPHA+1 ;itPHA2<cbnum;itPHA2++){
                    if (atomlist[itPHA]->getNum() == atomlist[itPHA2]->getNum()) {
                        coorlist[itPHA].setCoords( coorlist[itPHA].x +boxliste[itPHA2]->fixpos.x , coorlist[itPHA].y +boxliste[itPHA2]->fixpos.y ,coorlist[itPHA].z +boxliste[itPHA2]->fixpos.z );
                        Usedlist[itPHA2]  = true;
                        PHAlist[itPHA]    = true;
                        PHAlist[itPHA2]   =false;
                        weight++;
                    }
                }
                tol[itPHA]=weight;
                coorlist[itPHA].setCoords(coorlist[itPHA].x/weight,coorlist[itPHA].y/weight,coorlist[itPHA].z/weight);

            }else if( atomlist[itPHA]->props.isAromatic())
            {
                int weight =1;
                for ( int itPHA2 =itPHA+1 ;itPHA2<cbnum;itPHA2++){
                    if (atomlist[itPHA]->getResidu()->getNum() == atomlist[itPHA2]->getResidu()->getNum() && atomlist[itPHA2]->props.isAromatic()) {
                        coorlist[itPHA].setCoords( coorlist[itPHA].x +boxliste[itPHA2]->fixpos.x , coorlist[itPHA].y +boxliste[itPHA2]->fixpos.y ,coorlist[itPHA].z +boxliste[itPHA2]->fixpos.z );
                        //                        int transit = projs[boxliste[itPHA]->getId()];
                        //                        projs[boxliste[itPHA]->getId()]= (transit*weight +projs[boxliste[itPHA2]->getId()] )/3 ;
                        Usedlist[itPHA2]  = true;
                        PHAlist[itPHA]= true;
                        PHAlist[itPHA2]=false;
                        weight++;
                    }
                }
                tol[itPHA]=weight;
                coorlist[itPHA].setCoords(coorlist[itPHA].x/weight,coorlist[itPHA].y/weight,coorlist[itPHA].z/weight);


            }else if( atomlist[itPHA]->props.isHydrophobic())
            {
                int weight =1;
                for ( int itPHA2 =itPHA+1 ;itPHA2<cbnum;itPHA2++){
                    if (mod == 2) {
                        if (atomlist[itPHA]->getResidu()->getNum() == atomlist[itPHA2]->getResidu()->getNum()) {
                            PHAlist[itPHA]= true;
                            PHAlist[itPHA2]=false;
                            coorlist[itPHA].setCoords( coorlist[itPHA].x +boxliste[itPHA2]->fixpos.x ,coorlist[itPHA].y +boxliste[itPHA2]->fixpos.y ,coorlist[itPHA].z +boxliste[itPHA2]->fixpos.z );
                            Usedlist[itPHA2]  = true;
                            weight++;
                        }
                    } else if (mod == 8){

                        if ( boxliste[itPHA]->fixpos.calcDist(boxliste[itPHA2]->fixpos) <5 || coorlist[itPHA].calcDist(boxliste[itPHA2]->fixpos) <5)
                        {
                            PHAlist[itPHA]= true;
                            PHAlist[itPHA2]=false;
                            coorlist[itPHA].setCoords( coorlist[itPHA].x +boxliste[itPHA2]->fixpos.x ,coorlist[itPHA].y +boxliste[itPHA2]->fixpos.y ,coorlist[itPHA].z +boxliste[itPHA2]->fixpos.z );
                            Usedlist[itPHA2]  = true;
                            weight++;
                        }
                    }
                }
                tol[itPHA]=weight;
                coorlist[itPHA].setCoords(coorlist[itPHA].x/weight,coorlist[itPHA].y/weight,coorlist[itPHA].z/weight);
                Usedlist[itPHA]= true;
            }

        }else if (mod ==3 || mod > 4 ){
            ////////////////////////////// to here
            // Calcul for the median
            if (projs[boxliste[itPHA]->getId()]<bestenfoui) {  continue;}
            PHAlist[itPHA]=true;
            if (!atomlist[itPHA]->props.isHydrophobic() && !atomlist[itPHA]->props.isAromatic())
            {
                int weight =1;
                //                cout << "coord\t" << coorlist[itPHA].toString() << endl;
                for ( int itPHA2 =itPHA+1 ;itPHA2<cbnum;itPHA2++){
                    if (atomlist[itPHA]->getNum() == atomlist[itPHA2]->getNum()) {
                        coorlist[itPHA].setCoords( coorlist[itPHA].x + boxliste[itPHA2]->fixpos.x , coorlist[itPHA].y +boxliste[itPHA2]->fixpos.y ,coorlist[itPHA].z +boxliste[itPHA2]->fixpos.z );
                        Usedlist[itPHA2]  = true;
                        PHAlist[itPHA]    = true;
                        PHAlist[itPHA2]   = false;
                        weight++;
                    }
                }
                tol[itPHA]=weight;
                coorlist[itPHA].setCoords(coorlist[itPHA].x/weight,coorlist[itPHA].y/weight,coorlist[itPHA].z/weight);


            }else if( atomlist[itPHA]->props.isAromatic())
            {
                int weight =1;
                for ( int itPHA2 =itPHA+1 ;itPHA2<cbnum;itPHA2++){
                    if (atomlist[itPHA]->getResidu()->getNum() == atomlist[itPHA2]->getResidu()->getNum() && atomlist[itPHA2]->props.isAromatic()) {
                        coorlist[itPHA].setCoords( coorlist[itPHA].x +boxliste[itPHA2]->fixpos.x , coorlist[itPHA].y +boxliste[itPHA2]->fixpos.y ,coorlist[itPHA].z +boxliste[itPHA2]->fixpos.z );
                        //                        int transit = projs[boxliste[itPHA]->getId()];
                        //                        projs[boxliste[itPHA]->getId()]= (transit*weight +projs[boxliste[itPHA2]->getId()] )/3 ;
                        Usedlist[itPHA2]  = true;
                        PHAlist[itPHA]= true;
                        PHAlist[itPHA2]=false;
                        weight++;
                    }
                }
                tol[itPHA]=weight;
                coorlist[itPHA].setCoords(coorlist[itPHA].x/weight,coorlist[itPHA].y/weight,coorlist[itPHA].z/weight);

            }else if( atomlist[itPHA]->props.isHydrophobic())
            {

                //                PHAlist[itPHA]= false; // Décommenter pour enlever les hydrophobe alone
                int weight =1;
                for ( int itPHA2 =itPHA+1 ;itPHA2<cbnum;itPHA2++){
                    if (atomlist[itPHA2]->props.isHydrophobic() && !Usedlist[itPHA2] )
                    {

                        if (mod == 3 || mod == 5 || mod == 7 ) {
                            if (atomlist[itPHA]->getResidu()->getNum() == atomlist[itPHA2]->getResidu()->getNum()) {
                                PHAlist[itPHA]= true;
                                PHAlist[itPHA2]=false;
                                coorlist[itPHA].setCoords( coorlist[itPHA].x +boxliste[itPHA2]->fixpos.x ,coorlist[itPHA].y +boxliste[itPHA2]->fixpos.y ,coorlist[itPHA].z +boxliste[itPHA2]->fixpos.z );
                                Usedlist[itPHA2]  = true;
                                weight++;
                            }
                        } else if (mod == 6){

                            if ( boxliste[itPHA]->fixpos.calcDist(boxliste[itPHA2]->fixpos) <5 || coorlist[itPHA].calcDist(boxliste[itPHA2]->fixpos) <5)
                            {
                                PHAlist[itPHA]= true;
                                PHAlist[itPHA2]=false;
                                coorlist[itPHA].setCoords( coorlist[itPHA].x +boxliste[itPHA2]->fixpos.x ,coorlist[itPHA].y +boxliste[itPHA2]->fixpos.y ,coorlist[itPHA].z +boxliste[itPHA2]->fixpos.z );
                                Usedlist[itPHA2]  = true;
                                weight++;
                            }
                        }

                    }
                }
                tol[itPHA]=weight;
                coorlist[itPHA].setCoords(coorlist[itPHA].x/weight,coorlist[itPHA].y/weight,coorlist[itPHA].z/weight);
                Usedlist[itPHA]= true;
            }
            Usedlist[itPHA]= true;
        }
    }
    if (mod ==5) {
        for ( int itPHA =0 ;itPHA<cbnum;itPHA++){
            if (projs[boxliste[itPHA]->getId()]<bestenfoui) {  continue;}
            if (atomlist[itPHA] != (Atom*)NULL && PHAlist[itPHA] ){
                double distref=100000,dist =100000;
                int keep=itPHA;
                for ( int itPHA2 =0 ;itPHA2<cbnum;itPHA2++){
                    if ( (atomlist[itPHA]->getNum() == atomlist[itPHA2]->getNum()) ||  ( ( atomlist[itPHA]->props.isHydrophobic() || atomlist[itPHA]->props.isAromatic()) &&  atomlist[itPHA]->getResidu()->getNum() == atomlist[itPHA2]->getResidu()->getNum())) {
                        if (projs[boxliste[itPHA2]->getId()]<bestenfoui) {  continue;}
                        dist = coorlist[itPHA].calcDist(atomlist[itPHA2]->fixpos);
                        if (dist <  distref){
                            distref =dist;
                            keep=itPHA2;
                        }

                    }
                }
                PHAlist[itPHA] =false;
                PHAlist[keep] = true;
            }
        }
    }else if (mod == 7) {
        for ( int itPHA =0 ;itPHA<cbnum;itPHA++){
            int keepRes =100000;
            int Res=0;
            if (projs[boxliste[itPHA]->getId()]<bestenfoui) {  continue;}
            if (atomlist[itPHA] != (Atom*)NULL && PHAlist[itPHA] && atomlist[itPHA]->props.isHydrophobic()){
                //                cout << endl<<"hydrophobe" <<endl;
                AtomList listAtoms;
                Box *mbx= grid.getBox(coorlist[itPHA]);
                grid.getAdjacentAtoms(listAtoms,*mbx,4);

                for (ItCAtom itAtm= listAtoms.begin(); itAtm != listAtoms.end(); ++itAtm)
                {
                    Atom &atm = **itAtm;
                    //            atm = (*itAtm).second;
                    if (atm.getParent().getMoleType() != MoleType::PROTEIN) continue;
                    if (atm.isHydrogen() || atm.getMOL2Type()=="Du" || !atm.props.isHydrophobic() ) continue;
                    Res= atm.getResidu()->getNum();
                    if (keepRes==100000){
                        keepRes=Res;
                    }
                    if (Res != keepRes){
                        PHAlist[itPHA]=true;
                        //                        cout << "on conserve" << endl;
                        break;
                    }else {
                        PHAlist[itPHA]=false;
                    }
                }
            }
        }
    }

    if (mod ==1 || mod == 4 || mod == 5 ){
        for ( int itPHA =0 ;itPHA<cbnum;itPHA++)
        {
            if (PHAlist[itPHA]){
                Coords vectorD = atomlist[itPHA]->fixpos - boxliste[itPHA]->fixpos;
                double norme = sqrt(pow(vectorD.x,2)+pow(vectorD.y,2)+pow(vectorD.z,2));
                Coords Vector1 = Coords(vectorD.x/norme,vectorD.y/norme,vectorD.z/norme);
                double norme2 = sqrt(pow(Vector1.x,2)+pow(Vector1.y,2)+pow(Vector1.z,2));
                Coords pos1 = boxliste[itPHA]->fixpos+Vector1;
                Coords pos2 = pos1+Vector1;
                Coords pos3 = pos2+Vector1;



                if (atomlist[itPHA]->props.isCation())                                 {
                    molen.addAtom("O",boxliste[itPHA]->fixpos,"OD1","O.co2", &resASP).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                    PharmProp* newPharm = new PharmProp(PharType::HBACCEPTOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                    MonPharma.push_back(newPharm);

                }
                else if (atomlist[itPHA]->props.isAnion())                                  {
                    molen.addAtom("N",boxliste[itPHA]->fixpos,"NZ","N.4" , &resLYS).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                    PharmProp* newPharm = new PharmProp(PharType::HBDONOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                    MonPharma.push_back(newPharm);
                }
                else if (atomlist[itPHA]->props.isDonor() && atomlist[itPHA]->props.isAcceptor())    {
                    molen.addAtom("O",boxliste[itPHA]->fixpos,"OG","O.3" , &resSER).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                    molen.addAtom("P",atomlist[itPHA]->fixpos,"P","P",&resALA).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);  // Creer le vecteur en faisant la différence des coordonées
                    PharmProp* newPharm = new PharmProp(PharType::HBDONOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                    MonPharma.push_back(newPharm);
                    PharmProp* newPharm2 = new PharmProp(PharType::HBACCEPTOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                    MonPharma.push_back(newPharm2);
                }
                else if (atomlist[itPHA]->props.isAcceptor())                               {
                    molen.addAtom("N",boxliste[itPHA]->fixpos,"N","N.am"  , &resALA).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                    molen.addAtom("P",atomlist[itPHA]->fixpos,"P","P",&resALA).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);  // Creer le vecteur en faisant la différence des coordonées
                    PharmProp* newPharm = new PharmProp(PharType::HBDONOR,boxliste[itPHA]->fixpos,pos1,pos3,0);
                    MonPharma.push_back(newPharm);
                }
                else if (atomlist[itPHA]->props.isDonor())                                  {
                    molen.addAtom("O",boxliste[itPHA]->fixpos,"O","O.2"  , &resALA).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                    molen.addAtom("P",atomlist[itPHA]->fixpos,"P","P",&resALA).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);  // Creer le vecteur en faisant la différence des coordonées
                    PharmProp* newPharm = new PharmProp(PharType::HBACCEPTOR,coorlist[itPHA],pos1,pos3,0);
                    MonPharma.push_back(newPharm);

                }
                else if (atomlist[itPHA]->props.isAromatic())                               {
                    molen.addAtom("C",boxliste[itPHA]->fixpos,"CZ","C.ar", &resPHE).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                    grid.getAdjacentAtoms(listAtoms,*(boxliste[itPHA]),10);
                    int max=0;
                    Coords centroid;
                    int weight = 0;
                    for (ItCAtom itAtm= listAtoms.begin(); itAtm != listAtoms.end(); ++itAtm)
                    {
                        Atom &atm = **itAtm;
                        if (atm.props.isAromatic() && atm.getResidu() == atomlist[itPHA]->getResidu()){
                            if (boxliste[itPHA]->fixpos.calcDist(atm.fixpos)<8){
                                centroid.setCoords( (centroid.x+atm.fixpos.x),(centroid.y+atm.fixpos.y),(centroid.z+atm.fixpos.z));
                                weight++;
                                //                                cout << "Add \t"<< atm.fixpos.toString() << endl;
                            }
                        }
                    }
                    centroid.setCoords(centroid.x/weight,centroid.y/weight,centroid.z/weight); // On se sert du centroïde pour les aromatiques
                    vectorD = centroid - coorlist[itPHA];
                    norme = sqrt(pow(vectorD.x,2)+pow(vectorD.y,2)+pow(vectorD.z,2));
                    Vector1 = Coords(vectorD.x/norme,vectorD.y/norme,vectorD.z/norme);
                    pos1 = coorlist[itPHA]+Vector1;
                    pos3 = coorlist[itPHA]+Vector1+Vector1+Vector1;
                    molen.addAtom("P",centroid,"P","P",&resPHE).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                    PharmProp *newPharm = new PharmProp(PharType::AROMATIC,boxliste[itPHA]->fixpos,pos1,pos3,0);
                    MonPharma.push_back(newPharm);


                }
                else if (atomlist[itPHA]->props.isHydrophobic())                            {
                    molen.addAtom("C",boxliste[itPHA]->fixpos,"CA","C.3" , &resGLY).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                    PharmProp *newPharm = new PharmProp(PharType::HYDROPHOBIC,boxliste[itPHA]->fixpos,pos1,pos3,0);
                    MonPharma.push_back(newPharm);
                }
                else {cout << "ATOM NOT RECOGNIZED "<< bestatm->toString()<<endl;continue;}

            }
        }
    }else if (mod == 2 || mod == 3 || mod == 6 || mod == 7 || mod == 8) {

        for ( int itPHA =0 ;itPHA<cbnum;itPHA++)
        {
            if (PHAlist[itPHA]){
                //                Coords vectorD = Coords(atomlist[itPHA]->fixpos.x-coorlist[itPHA].x,atomlist[itPHA]->fixpos.y-coorlist[itPHA].y,atomlist[itPHA]->fixpos.z-coorlist[itPHA].z) ;

                Coords vectorD = atomlist[itPHA]->fixpos - coorlist[itPHA];
                double norme = sqrt(pow(vectorD.x,2)+pow(vectorD.y,2)+pow(vectorD.z,2));
                Coords Vector1 = Coords(vectorD.x/norme,vectorD.y/norme,vectorD.z/norme);
                Coords pos1 = coorlist[itPHA]+Vector1;
                Coords pos2 = pos1+Vector1;
                Coords pos3 = pos2+Vector1;
                double tolS=0;
                if (tolerance){
                    if (tol[itPHA]==1){
                        tolS-=30;
                    }else if (tol[itPHA]> 10){
                        tolS=30;
                    }
                }
                if (atomlist[itPHA]->props.isCation())                                 {
                    molen.addAtom("O",coorlist[itPHA],"OD1","O.co2", &resASP).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);

                    PharmProp* newPharm = new PharmProp(PharType::HBACCEPTOR,coorlist[itPHA],pos1,pos3,tolS);
                    MonPharma.push_back(newPharm);

                }
                else if (atomlist[itPHA]->props.isAnion())                                  {
                    molen.addAtom("N",coorlist[itPHA],"NZ","N.4" , &resLYS).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                    PharmProp* newPharm = new PharmProp(PharType::HBDONOR,coorlist[itPHA],pos1,pos3,tolS);
                    MonPharma.push_back(newPharm);

                }
                else if (atomlist[itPHA]->props.isDonor() && atomlist[itPHA]->props.isAcceptor())    {
                    molen.addAtom("O",coorlist[itPHA],"OG","O.3" , &resSER).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                    molen.addAtom("P",atomlist[itPHA]->fixpos,"P","P",&resALA).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                    PharmProp* newPharm = new PharmProp(PharType::HBDONOR,coorlist[itPHA],pos1,pos3,tolS);
                    MonPharma.push_back(newPharm);
                    PharmProp* newPharm2 = new PharmProp(PharType::HBACCEPTOR,coorlist[itPHA],pos1,pos3,tolS);
                    MonPharma.push_back(newPharm2);

                }
                else if (atomlist[itPHA]->props.isAcceptor())                               {
                    molen.addAtom("N",coorlist[itPHA],"N","N.am"  , &resALA).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                    molen.addAtom("P",atomlist[itPHA]->fixpos,"P","P",&resALA).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                    PharmProp* newPharm = new PharmProp(PharType::HBDONOR,coorlist[itPHA],pos1,pos3,tolS);
                    MonPharma.push_back(newPharm);

                }
                else if (atomlist[itPHA]->props.isDonor())                                  {
                    molen.addAtom("O",coorlist[itPHA],"O","O.2"  , &resALA).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                    molen.addAtom("P",atomlist[itPHA]->fixpos,"P","P",&resALA).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                    PharmProp* newPharm = new PharmProp(PharType::HBACCEPTOR,coorlist[itPHA],pos1,pos3,tolS);
                    MonPharma.push_back(newPharm);

                }
                else if (atomlist[itPHA]->props.isAromatic())                               {
                    molen.addAtom("C",coorlist[itPHA],"CZ","C.ar", &resPHE).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                    grid.getAdjacentAtoms(listAtoms,*(boxliste[itPHA]),10);
                    int max=0;
                    Coords centroid;
                    int weight = 0;
                    for (ItCAtom itAtm= listAtoms.begin(); itAtm != listAtoms.end(); ++itAtm)
                    {
                        Atom &atm = **itAtm;
                        if (atm.props.isAromatic() && atm.getResidu() == atomlist[itPHA]->getResidu()){
                            if (boxliste[itPHA]->fixpos.calcDist(atm.fixpos)<8){
                                centroid.setCoords( (centroid.x+atm.fixpos.x),(centroid.y+atm.fixpos.y),(centroid.z+atm.fixpos.z));
                                weight++;
                            }
                        }
                    }
                    centroid.setCoords(centroid.x/weight,centroid.y/weight,centroid.z/weight); // On se sert du centroïde pour les aromatiques
                    vectorD = centroid - coorlist[itPHA];
                    norme = sqrt(pow(vectorD.x,2)+pow(vectorD.y,2)+pow(vectorD.z,2));
                    Vector1 = Coords(vectorD.x/norme,vectorD.y/norme,vectorD.z/norme);
                    pos1 = coorlist[itPHA]+Vector1;
                    pos3 = coorlist[itPHA]+Vector1+Vector1+Vector1;
                    molen.addAtom("P",centroid,"P","P",&resPHE).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                    PharmProp *newPharm = new PharmProp(PharType::AROMATIC,coorlist[itPHA],pos1,pos3,tolS);
                    MonPharma.push_back(newPharm);

                }
                else if (atomlist[itPHA]->props.isHydrophobic())                            {
                    molen.addAtom("C",coorlist[itPHA],"CA","C.3" , &resGLY).setFormalCharge((double)projs[boxliste[itPHA]->getId()]);
                    PharmProp *newPharm = new PharmProp(PharType::HYDROPHOBIC,coorlist[itPHA],pos1,pos3,tolS);
                    MonPharma.push_back(newPharm);
                    //                    PharmProp* newPharm2 = new PharmProp(PharType::EXCLUSION,atomlist[itPHA]->fixpos,Coords(0,0,0),Coords(0,0,0),tolS);
                    //                    MonPharma.push_back(newPharm2);

                }
                else {cout << "ATOM NOT RECOGNIZED "<< bestatm->toString()<<endl;continue;}

            }
        }
    }

    //return of best b
    bestb=bestenfoui;
}




void VolSite::structCavLig(VolSiteResult  result[10],const Molecule &ligand, const string& name,const  bool& svm)
{

    ostringstream ossname;
    unsigned int num_cav=1,ns=0;
    double dist,tmpdist;
    const double volcube = grid.boxStep*grid.boxStep*grid.boxStep;
    set<const Box*> *bxlst=NULL;
    size_cavs.clear();
    // Choix de la cavité du ligand ar distance
    for (map<unsigned int, set<const Box*>  >::const_iterator itRef = cavities.begin(); itRef != cavities.end();itRef++)
    {
        const set<const Box*> &listbx=(*itRef).second;
        bool cavlig=false;
        ns=0;
        for (set<const Box*>::iterator  itComp = listbx.begin(); itComp != listbx.end(); itComp++){
            dist=1000;
            for (ItCAtom atom = ligand.firstAtom();atom != ligand.lastAtom();++atom)
            {
                tmpdist=(*itComp)->fixpos.calcDist((*atom)->fixpos);
                if (tmpdist < dist) dist =tmpdist;
            }
            if (dist <= 6){  ns++;cavlig=true;}
        }
        if (cavlig) size_cavs.insert(pair<unsigned int,unsigned int>(ns,(*itRef).first));
    }

    num_cav=1;
    vector<const Box*> cav[5];

    for (map<unsigned int,unsigned int>::const_reverse_iterator
         itsize = size_cavs.rbegin();
         itsize!=size_cavs.rend();
         itsize++)
    {
        for (size_t i=0;i<5;++i) cav[i].clear();

        bxlst = &((cavities.at((*itsize).second)));
        int i=0;
        for (set<const Box*>::const_iterator  itb = bxlst->begin(); itb != bxlst->end(); itb++){
            dist=1000;
            cav[4].push_back(*itb);

            for (ItCAtom atm=ligand.firstAtom();atm!=ligand.lastAtom();++atm){
                if (!grid.wHydrogen && (*atm)->isHydrogen())continue;
                tmpdist = (*itb)->fixpos.calcDist((*atm)->fixpos);
                if (tmpdist < dist) dist =tmpdist;
            }

            if (dist <= 12)  { cav[3].push_back(*itb);
                if (dist <= 8) { cav[2].push_back(*itb);
                    if (dist <= 6) { cav[1].push_back(*itb);
                        if (dist <= 4) { cav[0].push_back(*itb);

                            i++;
                        }   }  }// END dist <12
            }
        }

        if (cav[4].size() < 10) continue;

        bool check_name=false;
        if (name.length() >= 4) check_name=true;


        for (size_t iSize=0; iSize< 5;++iSize)
        {
            ossname.str("");
            if (check_name) ossname << name << "_CAVITY_N"<<num_cav;
            else ossname << "CAVITY_N"<<num_cav;
            string size="";
            switch (iSize)
            {
            case 0:result[0].type=MoleType::CAV_4;size="4";break;
            case 1:result[1].type=MoleType::CAV_6;size="6";break;
            case 2:result[2].type=MoleType::CAV_8;size="8";break;
            case 3:result[3].type=MoleType::CAV_12;size="12";break;
            case 4:result[4].type=MoleType::CAV_ALL;size="ALL";break;
            }
            result[iSize].nbcube=cav[iSize].size();
            result[iSize].volume=cav[iSize].size()*volcube;

            ossname << "_"<<size;
            cavToMOL(result[iSize].moleRes, cav[iSize] , ossname.str(),size);
            result[iSize].moleRes.setMoleType(result[iSize].type);
            result[iSize].ligRecov=0;
            result[iSize].cavRecov=0;
            for (ItCAtom itAtm= result[iSize].moleRes.firstAtom();
                 itAtm!= result[iSize].moleRes.lastAtom();
                 ++itAtm)
            {
                const Atom& atomole = **itAtm;
                if      (atomole.getName()=="CA") result[iSize].CA++;
                else if (atomole.getName()=="CZ") result[iSize].CZ++;
                else if (atomole.getName()=="O")  result[iSize].O++;
                else if (atomole.getName()=="OD1")result[iSize].OD1++;
                else if (atomole.getName()=="OG") result[iSize].OG++;
                else if (atomole.getName()=="N")  result[iSize].N++;
                else if (atomole.getName()=="NZ") result[iSize].NZ++;
                else if (atomole.getName()=="DU") result[iSize].DU++;
                for (ItCAtom atmL = ligand.firstAtom();atmL != ligand.lastAtom();++atmL)
                {
                    Atom &atmLigand = **atmL;
                    if (atmLigand.isHydrogen() || atmLigand.getMOL2Type()=="Du") continue;
                    if (atmLigand.fixpos.calcDist(atomole.fixpos,atmLigand.getVdWRadius()+0.1) < atmLigand.getVdWRadius())
                    {
                        result[iSize].cavRecov++;break;
                    }
                }
            }
            for (ItCAtom atmL = ligand.firstAtom();atmL != ligand.lastAtom();++atmL)
            {
                Atom &atmLigand = **atmL;
                if (atmLigand.isHydrogen() || atmLigand.getMOL2Type()=="Du") continue;
                for (ItCAtom itAtm= result[iSize].moleRes.firstAtom();
                     itAtm!= result[iSize].moleRes.lastAtom();
                     ++itAtm)
                {
                    const Atom& atomole = **itAtm;
                    if (atmLigand.fixpos.calcDist(atomole.fixpos,atmLigand.getVdWRadius()+0.1) < atmLigand.getVdWRadius()){result[iSize].ligRecov++;break;}
                }
            }


            if (iSize != 4)continue;
            if (num_cav == 1 && svm)
            {
                getDruggability(result[4].drugg,result[4].moleRes);
            }
        }
        num_cav++;
        break;
    }

}



void VolSite::structCavAll(VolSiteResult result[],const std::string& name, const bool &svm)
{

#ifdef ICHEM_DEBUG
    cout << "###################" << endl <<
            "Structure Cav All" << endl <<
            "###################" << endl;
#endif
    ostringstream ossname;
    // Testing all cavities :
    unsigned int num_cav=1;
    set<const Box*> *bxlst=NULL;
    size_cavs.clear();
    // SORTING CAVITIES  :
    for (map<unsigned int, set<const Box*> >::const_iterator
         ita = cavities.begin();
         ita != cavities.end();
         ita++)
        size_cavs.insert(pair<unsigned int,unsigned int>((*ita).second.size(),(*ita).first));

#ifdef ICHEM_DEBUG
    cout << name.substr(0,4) << "|CAVITY|NBCAV \t" << size_cavs.size() << endl;
    cout << name.substr(0,4) << "|CAVITY|DESCAV\tName\tNCav\tSize\tNPts\tVolume\tCA\tCZ\tO\tOG\tOD1\tN\tNZ\tDU\tDrugg"<<endl;
#endif

    num_cav=1;
    const double volcube = grid.boxStep*grid.boxStep*grid.boxStep;
    for (multimap<unsigned int,unsigned int>::const_reverse_iterator
         itsize = size_cavs.rbegin();
         itsize!=size_cavs.rend();
         itsize++)
    {
        vector<const Box*> cavAll;
        cavAll.clear();
        bxlst = &((cavities.at((*itsize).second)));
        for (set<const Box*>::const_iterator  itb = bxlst->begin(); itb != bxlst->end(); itb++){
            cavAll.push_back(*itb);

        }

        bool check_name=false;
        if (name.length() >= 4) check_name=true;


        ossname.str("");
        if (check_name) ossname << name << "_CAVITY";
        else ossname << "CAVITY";
        ossname<<  "_N"<<num_cav<< "_ALL.mol2";
        cavToMOL(result[num_cav-1].moleRes, cavAll, ossname.str(),"ALL");
        result[num_cav-1].NumCav=num_cav;
        result[num_cav-1].type=MoleType::CAV_ALL;
        result[num_cav-1].nbcube=cavAll.size();
        result[num_cav-1].volume=cavAll.size()*volcube;

        double drugg;
        if (svm) {

            getDruggability(drugg,result[num_cav-1].moleRes,drug_value) ;
        }else
        {
            getDruggability(drugg,result[num_cav-1].moleRes);
        }
        result[num_cav-1].drugg=drugg;
        for (ItCAtom atm= result[num_cav-1].moleRes.firstAtom();
             atm!= result[num_cav-1].moleRes.lastAtom();
             ++atm)
        {
            Atom& atomole = **atm;
            if      (atomole.getName()=="CA") result[num_cav-1].CA++;
            else if (atomole.getName()=="CZ") result[num_cav-1].CZ++;
            else if (atomole.getName()=="O")  result[num_cav-1].O++;
            else if (atomole.getName()=="OD1")result[num_cav-1].OD1++;
            else if (atomole.getName()=="OG") result[num_cav-1].OG++;
            else if (atomole.getName()=="N")  result[num_cav-1].N++;
            else if (atomole.getName()=="NZ") result[num_cav-1].NZ++;
            else if (atomole.getName()=="DU") result[num_cav-1].DU++;

        }

#ifdef ICHEM_DEBUG
        //cout<< name.substr(0,4) << "|CAVITY|DESCAV\t" << ((check_name)? name:"/")<<"\t"<< num_cav <<"\tALL\t"<<bxlst->size()*volcube<<"\t"<<bxlst->size()<<"\t"<<CA<<"\t"<<CZ<<"\t"<<O<<"\t"<<OG<<"\t"<<OD1<<"\t"<<N<<"\t"<<NZ<<"\t"<<DU<<"\t"<<drugg<<endl;
#endif
        num_cav++;


    }
#ifdef ICHEM_DEBUG
    cout << name.substr(0,4) << "|CAVITY|END" << endl;
#endif
}



void  VolSite::output_cavs_lig(const Molecule &ligand, const string& name,const  bool& svm, double &drug, bool Pharma, int PHA_mod)
{

#ifdef ICHEM_DEBUG
    cout << "###################" << endl <<
            "Ouptut Cav Lig" << endl <<
            "###################" << endl;
#endif

    ostringstream oss,ossname;
    // Testing all cavities :
    unsigned int num_cav=1,ns=0;
    double dist,tmpdist;
    double drugg=-100000;
    double druggNew=-100000;
    double recovery=0;
    set<const Box*> *bxlst=NULL;
    size_cavs.clear();

    //mole_cavs.clear();
    for (map<unsigned int, set<const Box*>  >::const_iterator
         itRef = cavities.begin();
         itRef != cavities.end();
         itRef++)
    {

        bool cavlig=false;
        ns=0;
        const set<const Box*> &listbx=(*itRef).second;
        for (set<const Box*>::iterator  itComp = listbx.begin(); itComp != listbx.end(); itComp++){
            dist=1000;
            for (ItCAtom atom = ligand.firstAtom();atom != ligand.lastAtom();++atom)
            {
                tmpdist=(*itComp)->fixpos.calcDist((*atom)->fixpos);
                if (tmpdist < dist) dist = tmpdist;
            }
            if (dist <= 6){  ns++;cavlig=true;}
        }
        if (cavlig) size_cavs.insert(pair<unsigned int,unsigned int>(ns,(*itRef).first));
    }

    ifstream ifs("VolSite_Stat.csv",ios::in);
    ofstream ofs;
    if (!ifs.is_open())
    {
        ofs.open("VolSite_Stat.csv",ios::out);
        ofs << "Name\tNCav\tSize\tNPts\tVolume\tCA\tCZ\tO\tOG\tOD1\tN\tNZ\tDU\tRecovery\tDrugg"<<endl;
    }else
    {
        ifs.close();
        ofs.open("VolSite_Stat.csv",ios::out|ios::app);
    }


    num_cav=1;
    const double volcube = grid.boxStep*grid.boxStep*grid.boxStep;
    vector<const Box*> cav[5];
    unsigned int CA=0,CZ=0,O=0,OD1=0,N=0,NZ=0,OG=0,DU=0;
    MoleWriter mw;
    for (map<unsigned int,unsigned int>::const_reverse_iterator itsize = size_cavs.rbegin(); itsize!=size_cavs.rend(); itsize++)
    {

        for (size_t i=0;i<5;++i) cav[i].clear();
        bxlst = &((cavities.at((*itsize).second)));
        for (set<const Box*>::const_iterator  itb = bxlst->begin(); itb != bxlst->end(); itb++){
            dist=1000;
            cav[4].push_back(*itb);

            for (ItCAtom atm=ligand.firstAtom();atm!=ligand.lastAtom();++atm){
                if (!grid.wHydrogen && (*atm)->isHydrogen())continue;
                tmpdist = (*itb)->fixpos.calcDist((*atm)->fixpos);
                if (tmpdist < dist) dist =tmpdist;
            }

            if (dist <= 12)  { cav[3].push_back(*itb);
                if (dist <= 8) { cav[2].push_back(*itb);
                    if (dist <= 6) { cav[1].push_back(*itb);
                        if (dist <= 4) {cav[0].push_back(*itb);}
                    }   }  }
        }
        if (cav[4].size() < 10) continue;

        bool check_name=false;
        if (name.length() >= 4) check_name=true;

        for (size_t iSize=0;iSize < 5;++iSize)
        {

            Molecule mole;
            CA=0,CZ=0,O=0,OD1=0,N=0,NZ=0,OG=0,DU=0;recovery=0;
            ossname.str("");
            if (check_name) ossname << name << "_CAVITY_N"<<num_cav;
            else ossname << "CAVITY_N"<<num_cav;
            string size="";
            switch (iSize)
            {
            case 0:size="4";break;
            case 1:size="6";break;
            case 2:size="8";break;
            case 3:size="12";break;
            case 4:size="ALL";break;
            }
            ossname << "_"<<size<<".mol2";
            cavToMOL(mole, cav[iSize] , ossname.str(),size);

            for (ItCAtom itAtm= mole.firstAtom();itAtm!= mole.lastAtom();++itAtm)
            {
                const Atom& atomole = **itAtm;
                if (atomole.getName()=="CA")CA++;
                else if (atomole.getName()=="CZ")CZ++;
                else if (atomole.getName()=="O")O++;
                else if (atomole.getName()=="OD1")OD1++;
                else if (atomole.getName()=="OG")OG++;
                else if (atomole.getName()=="N")N++;
                else if (atomole.getName()=="NZ")NZ++;
                else if (atomole.getName()=="DU")DU++;
                for (ItCAtom atmL = ligand.firstAtom();atmL != ligand.lastAtom();++atmL)
                {
                    Atom &atmLigand = **atmL;
                    if (atmLigand.isHydrogen() || atmLigand.getMOL2Type()=="Du") continue;
                    if (atmLigand.fixpos.calcDist(atomole.fixpos,atmLigand.getVdWRadius()+0.1) < atmLigand.getVdWRadius()){recovery++;break;}
                }
            }
            mw.newFile(ossname.str(),FileFormat::MOL2); mw.writeMOL2(&mole);


            if (iSize!=4)ofs<<((check_name)? name:ligand.getName())<<"\t"<< num_cav <<"\t"<<size<<"\t"<<cav[iSize].size()<<"\t"<<cav[iSize].size()*volcube<<"\t"<<CA<<"\t"<<CZ<<"\t"<<O<<"\t"<<OG<<"\t"<<OD1<<"\t"<<N<<"\t"<<NZ<<"\t"<<DU<<"\t"<<recovery/(double)cav[iSize].size()*100
                           <<"\t"<< "\t/" <<endl;

            if (iSize==4)
            {
                if (num_cav == 1 )
                {
                    drugg = 0;
                    getDruggability(drugg,mole,svm,0,num_cav);
                    //                    getDruggability(druggNew,mole,svm,1);
                    //% rajouter ici pour duomode
                    ofs<<((check_name)? name:ligand.getName())<<"\t"<< num_cav <<"\t"<<size<<"\t"<<cav[iSize].size()<<"\t"<<cav[iSize].size()*volcube<<"\t"<<CA<<"\t"<<CZ<<"\t"<<O<<"\t"<<OG<<"\t"<<OD1<<"\t"<<N<<"\t"<<NZ<<"\t"<<DU<<"\t"<<recovery/(double)cav[iSize].size()*100<<"\t\t"<< drugg<<endl;
                }else{
                    double druggtemp=0;
                    getDruggability(druggtemp,mole,svm,0,num_cav);
                    //                    getDruggability(druggNew,mole,svm,1);
                    ofs<<((check_name)? name:ligand.getName())<<"\t"<< num_cav <<"\t"<<size<<"\t"<<cav[iSize].size()<<"\t"<<cav[iSize].size()*volcube<<"\t"<<CA<<"\t"<<CZ<<"\t"<<O<<"\t"<<OG<<"\t"<<OD1<<"\t"<<N<<"\t"<<NZ<<"\t"<<DU<<"\t"<<recovery/(double)cav[iSize].size()*100<<"\t\t"<< druggtemp <<endl;
                }
                //            else {ofs<<((check_name)? name:ligand.getName())<<"\t"<< num_cav <<"\t"<<size<<"\t"<<cav[iSize].size()<<"\t"<<cav[iSize].size()*volcube<<"\t"<<CA<<"\t"<<CZ<<"\t"<<O<<"\t"<<OG<<"\t"<<OD1<<"\t"<<N<<"\t"<<NZ<<"\t"<<DU<<"\t"<<recovery/(double)cav[iSize].size()*100<<"\t"<< "\t/" <<endl;
                //}
            }//else if (num_cav == 1 && svm) getDruggability(drugg,moleAll,1);
            //        else ofs<<((check_name)? name:ligand.getName())<<"\t"<< num_cav <<"\t"<<size<<"\t"<<cav[4].size()*volcube<<"\t"<<cav[4].size()<<"\t"<<CA<<"\t"<<CZ<<"\t"<<O<<"\t"<<OG<<"\t"<<OD1<<"\t"<<N<<"\t"<<NZ<<"\t"<<DU<<"\t";


        }
        if (Pharma)
        {


            //////////////////////////////////////////////////
            // Pharma .mol2
            int bestb = 0;

            Molecule mole_pha;
            PharmList pharmaco;
            ossname.str("");
            string size="ALL";
            if (check_name) ossname << name << "_PHARMA_N"<<num_cav;
            else ossname << "PHARMA_N"<<num_cav;

            ossname << ".mol2";
            cavToPHA(mole_pha,pharmaco, cav[4] , ossname.str(),size,PHA_mod,bestb);
            mw.newFile(ossname.str(),FileFormat::MOL2); mw.writeMOL2(&mole_pha);

            ///////////// Cavity used to buil pharmacophore
            Molecule mole_cav;
            ossname.str("");
            if (check_name) ossname << name << "_CAVITY_B"<<num_cav;
            else ossname << "CAVITY_B"<<num_cav;

            ossname << "_" <<size <<".mol2";
            cavToMOLB(mole_cav, cav[4], ossname.str(),size,bestb);
            mw.newFile(ossname.str(),FileFormat::MOL2); mw.writeMOL2(&mole_cav);


            ///////////////////////////////////////////////////////
            //pharma .chm

            ossname.str("");
            if (check_name) ossname << name << "_PHARMA_N"<<num_cav;
            else ossname << "PHARMA_N"<<num_cav;
            ossname<<".chm";
            PharmWriter pw(ossname.str(),FileFormat::CHM);
            string toto =ossname.str();
            pw.writeCHM(pharmaco,toto.c_str(),Bweight);

        }
        num_cav++;


    }

    ofs.close();
    drug = drugg;
}
void VolSite::output_cavs_all(const std::string& name, const bool& svm, bool Pharma, int PHA_mod, vector<string> Input)
{

#ifdef ICHEM_DEBUG
    cout << "###################" << endl <<
            "Ouptut Cav All" << endl <<
            "###################" << endl;
#endif
    ostringstream oss,ossname;
    // Testing all cavities :
    unsigned int num_cav=1;
    unsigned int CA=0,CZ=0,O=0,OD1=0,N=0,NZ=0,OG=0,DU=0;
    set<const Box*> *bxlst=NULL;
    size_cavs.clear();
    double avgBuried=0;

    //mole_cavs.clear();
    //    cout << "CAVITIES : " <<cavities.size()<<endl;
    // SORTING CAVITIES  :
    for (map<unsigned int, set<const Box*> >::iterator ita = cavities.begin(); ita != cavities.end();ita++) size_cavs.insert(pair<unsigned int,unsigned int>((*ita).second.size(),(*ita).first));

    //    cout << "SIZE CAV : "<<size_cavs.size()<<endl;

    ifstream ifs("VolSite_Stat.csv",ios::in);
    ofstream ofs;

#ifdef ICHEM_DEBUG
    cout << name.substr(0,4) << "|CAVITY|NBCAV \t" << size_cavs.size() << endl;
#endif

    if (!ifs.is_open())
    {
        ofs.open("VolSite_Stat.csv",ios::out);
        ofs << "Name\tNCav\tSize\tBuriedness\tNPts\tVolume\tCA\tCZ\tO\tOG\tOD1\tN\tNZ\tDU\tDrugg"<<endl;
    }else
    {
        ifs.close();
        ofs.open("VolSite_Stat.csv",ios::out|ios::app);
    }
#ifdef ICHEM_DEBUG
    cout << name.substr(0,4) << "|CAVITY|DESCAV\tName\tNCav\tSize\tNPts\tVolume\tCA\tCZ\tO\tOG\tOD1\tN\tNZ\tDU\tDrugg"<<endl;
#endif

    /*  oss << "********************************     VOLSITE    *******************************"<<endl<<endl
                    << "-------------------------------------------------------------------------------"<<endl
                    << "-------------------------------------------------------------------------------"<<endl
                    << "CAVITY \t ID \t NCube \t Volume \t Druggability \n"<<endl;
               */
    //unsigned int
    num_cav=1;
    const double volcube = grid.boxStep*grid.boxStep*grid.boxStep;
    // changer le tri sur la distance a l'interface.//%
    for (multimap<unsigned int,unsigned int>::reverse_iterator itsize = size_cavs.rbegin(); itsize!=size_cavs.rend(); itsize++)
    {
        vector<const Box*> cavAll;
        cavAll.clear();
        bxlst = &((cavities.at((*itsize).second)));
        for (set<const Box*>::iterator  itb = bxlst->begin(); itb != bxlst->end(); itb++){
            cavAll.push_back(*itb);
            avgBuried += projs[(*itb)->getId()];
        }
        avgBuried /= cavAll.size();
        bool check_name=false;
        if (name.length() >= 4) check_name=true;

        Molecule mole;

        ossname.str("");
        if (check_name) ossname << "CAVITY";
        else ossname << "CAVITY";
        ossname<<  "_N"<<num_cav<< "_ALL.mol2";
        cavToMOL(mole, cavAll, ossname.str(),"ALL");
        //        cavToPHA(mole, cavAll, ossname.str(),"ALL",1);
        double drugg;
        double druggNew;
        getDruggability(drugg,mole,svm,0,num_cav) ;
        //        getDruggability(druggNew,mole,svm,1) ;

        //     oss << drugg<<endl;
        MoleWriter  mw(ossname.str(),FileFormat::MOL2);
        mw.writeMOL2(&mole);

        if (Pharma)
        {
            vector<string> input;
            input.push_back("--pharm");
            input.push_back("volsite");

            input.push_back(Input[0]);
            input.push_back(ossname.str());
//            for (int ii = 0; ii<input.size();ii++ ){
//                cout << input[ii] << endl;
//            }
//            cout << Input. << endl;
            IChemSwitch(input.size(),input);
        }

        CA=0,CZ=0,O=0,OD1=0,N=0,NZ=0,OG=0,DU=0;

        for (ItCAtom atm= mole.firstAtom();atm!= mole.lastAtom();++atm)
        {
            //            cout << atm.getPointer()<<endl;
            //            cout << atm->toString()<<endl;
            Atom& atomole = **atm;
            if (atomole.getName()=="CA")CA++;
            else if (atomole.getName()=="CZ")CZ++;
            else if (atomole.getName()=="O")O++;
            else if (atomole.getName()=="OD1")OD1++;
            else if (atomole.getName()=="OG")OG++;
            else if (atomole.getName()=="N")N++;
            else if (atomole.getName()=="NZ")NZ++;
            else if (atomole.getName()=="DU")DU++;

        }
        //       MoleWriter  mw(ossname.str(),FileFormat::MOL2);
        //mw.writeMOL2(&mole);

        // for Duomode
        //ofs<<((check_name)? name:"/")<<"\t"<< num_cav <<"\tALL\t"<<bxlst->size()<<"\t"<<bxlst->size()*volcube<<"\t"<<CA<<"\t"<<CZ<<"\t"<<O<<"\t"<<OG<<"\t"<<OD1<<"\t"<<N<<"\t"<<NZ<<"\t"<<DU<<"\t"<<drugg<<"/"<<druggNew<<endl;
        ofs<<((check_name)? name:"/")<<"\t"<< num_cav <<"\tALL\t"<<avgBuried<<"\t"<<bxlst->size()<<"\t"<<bxlst->size()*volcube<<"\t"<<CA<<"\t"<<CZ<<"\t"<<O<<"\t"<<OG<<"\t"<<OD1<<"\t"<<N<<"\t"<<NZ<<"\t"<<DU<<"\t"<<drugg<<endl;
#ifdef ICHEM_DEBUG
        cout<< name.substr(0,4) << "|CAVITY|DESCAV\t" << ((check_name)? name:"/")<<"\t"<< num_cav <<"\tALL\t"<<bxlst->size()*volcube<<"\t"<<bxlst->size()<<"\t"<<CA<<"\t"<<CZ<<"\t"<<O<<"\t"<<OG<<"\t"<<OD1<<"\t"<<N<<"\t"<<NZ<<"\t"<<DU<<"\t"<<drugg<<endl;
#endif
        //        mole.clear();
        //else mole_cavs.push_back(mole);
        num_cav++;


    }
#ifdef ICHEM_DEBUG
    cout << name.substr(0,4) << "|CAVITY|END" << endl;
    cout << oss.str() << endl;
#endif

    ofs << "END" << endl;
    ofs.close();
    cout << num_cav-1 << " Cavities found" << endl;
}


void VolSite::output_cavs_ppi(const std::string& name, const bool& svm)
{

#ifdef ICHEM_DEBUG
    cout << "###################" << endl <<
            "Ouptut Cav All" << endl <<
            "###################" << endl;
#endif
    ostringstream oss,ossname;
    // Testing all cavities :
    unsigned int num_cav=1;
    unsigned int CA=0,CZ=0,O=0,OD1=0,N=0,NZ=0,OG=0,DU=0;
    set<const Box*> *bxlst=NULL;
    size_cavs.clear();
    double avgBuried=0;
    //mole_cavs.clear();
    //    cout << "CAVITIES : " <<cavities.size()<<endl;
    // SORTING CAVITIES  :
    for (map<unsigned int, set<const Box*> >::iterator ita = cavities.begin(); ita != cavities.end();ita++) size_cavs.insert(pair<unsigned int,unsigned int>((*ita).second.size(),(*ita).first));

    //    cout << "SIZE CAV : "<<size_cavs.size()<<endl;

    ifstream ifs("VolSite_Stat.csv",ios::in);
    ofstream ofs;

#ifdef ICHEM_DEBUG
    cout << name.substr(0,4) << "|CAVITY|NBCAV \t" << size_cavs.size() << endl;
#endif

    if (!ifs.is_open())
    {
        ofs.open("VolSite_Stat.csv",ios::out);
        ofs << "Name\tNCav\tSize\tBuriedness\tNPts\tVolume\tCA\tCZ\tO\tOG\tOD1\tN\tNZ\tDU\tDrugg"<<endl;
    }else
    {
        ifs.close();
        ofs.open("VolSite_Stat.csv",ios::out|ios::app);
    }
#ifdef ICHEM_DEBUG
    cout << name.substr(0,4) << "|CAVITY|DESCAV\tName\tNCav\tSize\tNPts\tVolume\tCA\tCZ\tO\tOG\tOD1\tN\tNZ\tDU\tDrugg"<<endl;
#endif

    /*  oss << "********************************     VOLSITE    *******************************"<<endl<<endl
                    << "-------------------------------------------------------------------------------"<<endl
                    << "-------------------------------------------------------------------------------"<<endl
                    << "CAVITY \t ID \t NCube \t Volume \t Druggability \n"<<endl;
               */
    //unsigned int
    num_cav=1;
    const double volcube = grid.boxStep*grid.boxStep*grid.boxStep;
    ofs << "Cavity for "<< name << endl;
    for (multimap<unsigned int,unsigned int>::reverse_iterator itsize = size_cavs.rbegin(); itsize!=size_cavs.rend(); itsize++)
    {
        vector<const Box*> cavAll;
        cavAll.clear();
        bxlst = &((cavities.at((*itsize).second)));
        for (set<const Box*>::iterator  itb = bxlst->begin(); itb != bxlst->end(); itb++){
            cavAll.push_back(*itb);
            avgBuried += projs[(*itb)->getId()];
        }
        avgBuried /= cavAll.size();
        //if (cavAll.size()*volcube <min_size ) continue;

        /*    oss.precision(5);
                    oss << "----------------------------\t-----\t-----\t-----\t-----\t-----"<<endl
                        << "CAVITY\t"<< num_cav <<"\t"<<bxlst->size()<<"\t"<<bxlst->size()*volcube<<"\t";
                */
        bool check_name=false;
        if (name.length() >= 4) check_name=true;


        Molecule mole;

        ossname.str("");
        if (check_name) ossname << "CAVITY";
        else ossname << "CAVITY";
        ossname <<"_" << name;
        ossname<<  "_N"<<num_cav<< "_ALL.mol2";
        cavToMOL(mole, cavAll, ossname.str(),"ALL");
        //        cout << "cavité distante de :" << intscoor.calcDist(mole.getFixpos()) << " Ang" << endl;
        double drugg;
        double druggNew;
        getDruggability(drugg,mole,svm,0,num_cav) ;
        //        getDruggability(druggNew,mole,svm,1) ;

        //     oss << drugg<<endl;
        MoleWriter  mw(ossname.str(),FileFormat::MOL2);
        mw.writeMOL2(&mole);
        CA=0,CZ=0,O=0,OD1=0,N=0,NZ=0,OG=0,DU=0;

        for (ItCAtom atm= mole.firstAtom();atm!= mole.lastAtom();++atm)
        {
            //            cout << atm.getPointer()<<endl;
            //            cout << atm->toString()<<endl;
            Atom& atomole = **atm;
            if (atomole.getName()=="CA")CA++;
            else if (atomole.getName()=="CZ")CZ++;
            else if (atomole.getName()=="O")O++;
            else if (atomole.getName()=="OD1")OD1++;
            else if (atomole.getName()=="OG")OG++;
            else if (atomole.getName()=="N")N++;
            else if (atomole.getName()=="NZ")NZ++;
            else if (atomole.getName()=="DU")DU++;

        }
        //       MoleWriter  mw(ossname.str(),FileFormat::MOL2);
        //mw.writeMOL2(&mole);

        // for Duomode
        //ofs<<((check_name)? name:"/")<<"\t"<< num_cav <<"\tALL\t"<<bxlst->size()<<"\t"<<bxlst->size()*volcube<<"\t"<<CA<<"\t"<<CZ<<"\t"<<O<<"\t"<<OG<<"\t"<<OD1<<"\t"<<N<<"\t"<<NZ<<"\t"<<DU<<"\t"<<drugg<<"/"<<druggNew<<endl;
        ofs<<((check_name)? name:"/")<<"\t"<< num_cav <<"\tALL\t"<<avgBuried<<"\t"<<bxlst->size()<<"\t"<<bxlst->size()*volcube<<"\t"<<CA<<"\t"<<CZ<<"\t"<<O<<"\t"<<OG<<"\t"<<OD1<<"\t"<<N<<"\t"<<NZ<<"\t"<<DU<<"\t"<<drugg<<endl;
#ifdef ICHEM_DEBUG
        cout<< name.substr(0,4) << "|CAVITY|DESCAV\t" << ((check_name)? name:"/")<<"\t"<< num_cav <<"\tALL\t"<<bxlst->size()*volcube<<"\t"<<bxlst->size()<<"\t"<<CA<<"\t"<<CZ<<"\t"<<O<<"\t"<<OG<<"\t"<<OD1<<"\t"<<N<<"\t"<<NZ<<"\t"<<DU<<"\t"<<drugg<<endl;
#endif
        //        mole.clear();plp
        //else mole_cavs.push_back(mole);
        num_cav++;


    }
#ifdef ICHEM_DEBUG
    cout << name.substr(0,4) << "|CAVITY|END" << endl;
    cout << oss.str() << endl;
#endif

    ofs << "END" << endl;
    ofs.close();
    cout << num_cav-1 << " Cavities found" << endl;
}


void VolSite::getDruggability(double& dist, const Molecule &mole, const bool& svm, int combo,int numcav)
{
    dist = 0;
    double projsval[8][10];
    double propsval[9];
    unsigned int tproj;
    size_t tpos=0;
    const double volcube = grid.boxStep*grid.boxStep*grid.boxStep;
    std::ofstream trainsvm;
    //% emplacement du fichier training que je n'aime pas.
    if (svm && combo==0){
        string descname;
        descname = "./"+ name +"_descriptor.txt";
        trainsvm.open(descname.c_str(), ios::out|ios::app);
    }
    // Initialize values :
    for(unsigned int i=0; i < 8; i++){propsval[i]=0.0;
        for (unsigned int j=0;j<10; j++){projsval[i][j]=0.0;}}
    propsval[8]=0.0;
    for (ItCAtom itAtom = mole.firstAtom();itAtom!= mole.lastAtom() ;++itAtom)
    {
        const Atom &atm = **itAtom;

        int atmproj=(int)atm.getFormalCharge();
        if (atmproj <minProj)continue;
        tproj =static_cast<unsigned int>(floor(((int)atm.getFormalCharge()-minProj)/10));
        propsval[0]+=volcube;

        if (atm.getName()=="DU"   ){tpos=7;}
        else if (atm.getName()=="OD1" ){tpos=3;}
        else if (atm.getName()=="NZ" ){tpos=6;}
        else if (atm.getName()=="OG" ){tpos=4;}
        else if (atm.getName()=="N" ){tpos=5;}
        else if (atm.getName()=="O" ){tpos=2;}
        else if (atm.getName()=="CA"){tpos=1;}
        else if (atm.getName()=="CZ" ){tpos=0;}
        else cerr << "not found"<<atm.toString()<<endl;

        propsval[tpos+1]+= volcube;
        if (tpos != 7 || (tpos == 7 && tproj < 10))   projsval[tpos][tproj]+=volcube;
    }
    ostringstream oss;

    for(unsigned int i=0; i< 8;i++)
        for (unsigned int j=0;j<10; j++) projsval[i][j] = projsval[i][j]/propsval[0]*100;
    for (unsigned int i=1; i<9; i++)     propsval[i]=propsval[i]/propsval[0]*100;



    if (svm&& combo==0){
        //trainsvm << drug_value <<" "; //only for svm file
        oss << drug_value << " ";
    } else
    {oss << "0 ";
    }

    unsigned int pos=1;

    for (unsigned int i=0; i<9; i++)
    {

        if (svm && combo ==0){
            trainsvm << propsval[i]<<" ";

            /* switch (i)
                    {
                    case 0: cout << "VOLUME:";break;
                    case 1: cout << "CZ:";break;
                    case 2: cout << "CA:";break;
                    case 3: cout << "O:";break;
                    case 4: cout << "OD1:";break;
                    case 5: cout << "OG:";break;
                    case 6: cout << "N:";break;
                    case 7: cout << "NZ:";break;
                    case 8: cout << "DU:";break;
                    }*/
        }

        //cout << propsval[i]<<endl;
        oss << pos<<":"<<propsval[i]<<" ";pos++;
    }
    for(unsigned int i=0; i< 8;i++) for (unsigned int j=0;j<10; j++)
    {
        if (svm && combo==0){
            trainsvm << projsval[i][j] << " ";
        }
        oss << pos<<":"<< projsval[i][j]<<" ";pos++;

    }

    // cout << oss.str()<<endl;

    //% Need SVM

    WORD *words;DOC *doc;
    char *comment;
    long max_words_doc=100000,queryid,slackid, wnum;
    double doc_label,costfactor;

    words = (WORD *)my_malloc(sizeof(WORD)*(170));


    MODEL model;
    if (combo==0){
        read_VolSitemodel(model,const_cast<char*>(modelsvm.c_str()));
    }else{
        read_VolSitemodel(model,const_cast<char*>(modelsvm2.c_str()));
    }

    parse_document((char*)oss.str().c_str(),words,&doc_label,&queryid,&slackid,&costfactor,&wnum,
                   max_words_doc,&comment);

    doc = create_example(-1,0,0,0.0,create_svector(words,comment,1.0));
    dist=classify_example(&model,doc);
    free_example(doc,1);

    if (svm&& combo==0){
        trainsvm << "#"<< name <<"_"<< numcav  <<endl;
        trainsvm.close();
    }
    free(words);

    //    free_model(&model,1);

    //    return dist;
}


void VolSite::read_VolSitemodel(MODEL &modelin, const std::string msvm )
{
    vector<std::string> lignes;
    //string msvm=tu;
    tokenStr(msvm,lignes,"|_|");
    vector<std::string>::iterator ligne=lignes.begin();
    string tmp = *ligne;
    //   modelin = (MODEL *)my_malloc(sizeof(MODEL));
    char version_buffer[100];
    sscanf(tmp.c_str(),"SVM-light Version %s\n",version_buffer);
    if(strcmp(version_buffer,VERSION)) {
        perror ("Version of model-file does not match version of svm_classify!");
        exit (1);
    }
    ligne++; tmp = *ligne;sscanf(tmp.c_str(),"%ld%*[^\n]\n", &modelin.kernel_parm.kernel_type);
    ligne++; tmp = *ligne;sscanf(tmp.c_str(),"%ld%*[^\n]\n", &modelin.kernel_parm.poly_degree);
    ligne++; tmp = *ligne;sscanf(tmp.c_str(),"%lf%*[^\n]\n", &modelin.kernel_parm.rbf_gamma);
    ligne++; tmp = *ligne;sscanf(tmp.c_str(),"%lf%*[^\n]\n", &modelin.kernel_parm.coef_lin);
    ligne++; tmp = *ligne;sscanf(tmp.c_str(),"%lf%*[^\n]\n", &modelin.kernel_parm.coef_const);
    ligne++; tmp = *ligne;sscanf(tmp.c_str(),"%[^#]%*[^\n]\n", modelin.kernel_parm.custom);
    ligne++; tmp = *ligne;sscanf(tmp.c_str(),"%ld%*[^\n]\n", &modelin.totwords);
    ligne++; tmp = *ligne;sscanf(tmp.c_str(),"%ld%*[^\n]\n", &modelin.totdoc);
    ligne++; tmp = *ligne;sscanf(tmp.c_str(),"%ld%*[^\n]\n", &modelin.sv_num);
    ligne++; tmp = *ligne;sscanf(tmp.c_str(),"%lf%*[^\n]\n", &modelin.b);

    modelin.supvec = (DOC **)my_malloc(sizeof(DOC *)*modelin.sv_num);
    modelin.alpha = (double *)my_malloc(sizeof(double)*modelin.sv_num);
    modelin.index=NULL;
    modelin.lin_weights=NULL;
    long max_words = 160;
    WORD *words;long i,queryid,slackid,wpos; double costfactor;char *comment;
    words = (WORD *)my_malloc(sizeof(WORD)*(max_words+10));
    char t[100];
    comment =t;

    for(i=1;i<modelin.sv_num;i++) {
        ligne++; tmp = *ligne;
        if (tmp.length()==0) continue;

        if(!parse_document((char*)tmp.c_str(),words,&(modelin.alpha[i]),&queryid,&slackid,
                           &costfactor,&wpos,max_words,&comment)) {
            printf("\nParsing error while reading model file in SV %ld!\n%s",
                   i,tmp.c_str());
            exit(1);
        }
        modelin.supvec[i] = create_example(-1,
                                           0,0,
                                           0.0,
                                           create_svector(words,comment,1.0));
    }
    free(words);
    //    return(model);
}


bool appartient( std::string name, const std::string& liste)
{
    if (name.find(" ,") == string::npos && liste.find(name) != string::npos) return true;
    std::string::size_type lastPos = liste.find_first_not_of(" ,", 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = liste.find_first_of(" ,", lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        if(name == liste.substr(lastPos, pos - lastPos))
            return true;
        // Skip delimiters.  Note the "not_of"
        lastPos = liste.find_first_not_of(" ,", pos);
        // Find next "non-delimiter"
        pos = liste.find_first_of(" ,", lastPos);
    }

    return false;
}

void VolSite::defProps(Molecule& protein, bool pharma)
{
    string  resname;
    //    cout <<" je suis dans defprops " << protein.getName() << endl;

    bool modif = false;
    for (ItCAtom itA = protein.firstAtom(); itA != protein.lastAtom(); ++itA)
    {
        Atom &atm = **itA;
        const std::string& name = atm.getName();
        modif=false;
        PhysProp p = atm.props;
        //
        if (atm.getResidu() != (Residu*)NULL) resname = atm.getResidu()->getName();else resname="";

        if (!pharma){
            atm.props.clear();
            if (atm.isHydrogen()) continue;
            if (name.compare("C")==0)                               {atm.props.setHydrophobic(true);continue;}//modifier
            if (appartient(name,"CA, CB") && resname != "CA")                       {atm.props.setHydrophobic(true);continue;}
            if (name.compare("CD") == 0 && appartient(resname,"ARG, LYS, PRO"))     {atm.props.setHydrophobic(true);continue;}
            if (name.compare("CD1")==0){
                if (appartient(resname,"ILE, LEU"))                     {atm.props.setHydrophobic(true);continue;}
                else if (appartient(resname, "PHE, TRP, TYR"))                  {atm.props.setAromatic(true);continue;}}
            if (name.compare("CD2")==0){
                if (resname.compare("LEU")==0)                          {atm.props.setHydrophobic(true);continue;}
                else if (appartient(resname, "PHE, TRP, TYR"))                  {atm.props.setAromatic(true);continue;}
                else if (appartient(resname, "HIS,HID,HIE"))                    {atm.props.setAromatic(true);
                    atm.props.setHydrophobic(false);continue;}}
            if (name.compare("CE")==0 && appartient(resname,"LYS, MET"))            {atm.props.setHydrophobic(true);continue;}
            if (name.compare("CE1")==0){
                if (appartient(resname,"HIS,HID,HIE"))                      {atm.props.setAromatic(true);
                    atm.props.setHydrophobic(false);continue;}
                else if (appartient(resname,"PHE, TYR"))                    {atm.props.setAromatic(true);continue;}}
            if (name.compare("CE2")==0 && appartient(resname, "PHE, TRP, TYR"))     {atm.props.setAromatic(true);continue;}
            if (name.compare("CE3")==0 && resname.compare("TRP")==0)                {atm.props.setAromatic(true);continue;}
            if (name.compare("CG")==0)
            {    if (appartient(resname,"ARG, GLN, GLU, LEU, LYS, MET, PRO")){ atm.props.setHydrophobic(true);continue;}
                else if (appartient(resname,"HIS, HID, HIE, PHE, TRP, TYR")){ atm.props.setAromatic(true);continue;}}
            if ((name.compare("CG1")==0 && appartient(resname,"ILE, VAL"))
                    || (name.compare("CG2")==0 && appartient(resname,"ILE, THR, VAL"))){ atm.props.setHydrophobic(true);continue;}
            if (name.compare("CH2")==0 && resname.compare("TRP")==0) {atm.props.setAromatic(true);continue;}
            if (name.compare("CZ")==0){
                if (resname.compare("ARG")==0)                          { atm.props.setHydrophobic(false);continue;}
                else if (appartient(resname,"PHE, TYR"))                    { atm.props.setAromatic(true);continue;}}
            if (appartient(name, "CZ2, CZ3") && resname.compare("TRP")==0)          { atm.props.setAromatic(true);continue;}
            if (name.compare("N")==0)                               { atm.props.setDonor(true);continue;}
            if (appartient(name,"ND1, ND2") && appartient(resname,"HIS, HID, HIE"))          { atm.props.setDonor(true);atm.props.setAcceptor(true);continue;}
            if (name.compare("ND2")==0 && resname.compare("ASN")==0)                { atm.props.setDonor(true);continue;}
            if (name.compare("NE")==0 && appartient(resname,"ARG"))            { atm.props.setCation(true);continue;}
            if (name.compare("NE1")==0 && resname.compare("TRP")==0)                { atm.props.setDonor(true);continue;}
            if (name.compare("NE2")==0  ){
                if (resname.compare("GLN")==0)                          {atm.props.setDonor(true);continue;}
                else if (appartient(resname,"HIS, HID, HIE"))                   { atm.props.setDonor(true);
                    atm.props.setAcceptor(true);continue;}}
            if ((appartient(name,"NH1, NH2") && resname.compare("ARG")==0)
                    || (name.compare("NZ")==0 && resname.compare("LYS")==0))            {atm.props.setCation(true);continue;}
            if (name.compare("O")==0 && resname.compare("HOH")==0)                               {
                atm.props.setAcceptor(true);atm.props.setDonor(true);continue;}
            else if (name.compare("O")==0){
                atm.props.setAcceptor(true);continue;}
            if (appartient(name,"OD1, OD2") && resname.compare("ASP")==0)           {atm.props.setAnion(true);continue;}
            if (name.compare("OD1")==0 && resname.compare("ASN")==0)                {atm.props.setAcceptor(true);continue;}
            if (name.compare("OE1")==0 && resname.compare("GLN")==0)                {atm.props.setAcceptor(true);continue;}
            if (appartient(name,"OE1, OE2")  && resname.compare("GLU")==0)          {atm.props.setAnion(true);continue;}
            if ((name.compare("OG") ==0 && resname.compare("SER")==0)||
                    (name.compare("OG1")==0 && resname.compare("THR")==0)||
                    (name.compare("OH") ==0 && resname.compare("TYR")==0))              {atm.props.setDonor(true);
                atm.props.setAcceptor(true);continue;}
            if (name.compare("OXT")==0)                             {atm.props.setAnion(true);continue;}
            if (appartient(name, "SD, SG")){
                if (appartient(resname, "CYS, CYX")) {atm.props.setAcceptor(true);continue;}
                else if (resname.compare("MET")==0) {atm.props.setHydrophobic(true);continue;}}
            if (appartient(resname, "FE, MG, MG, MN, ZN, CO, CA, NA")) {atm.props.setCation(true),atm.props.setMetal(true);continue;}
            atm.props = p;
        }
        else {
            atm.props.clear();
            if (atm.isHydrogen()) continue;
            if (appartient(name,"CB") && resname != "CA")                       {
                if (appartient(resname,"ALA, VAL, ILE, LEU, PRO, CYS, MET, PHE, TYR, TRP")) {atm.props.setHydrophobic(true);continue;}}
            if (name.compare("CD") == 0 && resname.compare("PRO"))     {atm.props.setHydrophobic(true);continue;}
            if (name.compare("CD1")==0){
                if (appartient(resname,"ILE, LEU"))                     {atm.props.setHydrophobic(true);continue;}
                else if (appartient(resname, "PHE, TRP, TYR"))                  {atm.props.setAromatic(true);continue;}}
            if (name.compare("CD2")==0){
                if (resname.compare("LEU")==0)                          {atm.props.setHydrophobic(true);continue;}
                else if (appartient(resname, "PHE, TRP, TYR"))                  {atm.props.setAromatic(true);continue;}
                else if (appartient(resname, "HIS,HID,HIE"))                    {atm.props.setAromatic(true);
                    atm.props.setHydrophobic(false);continue;}}
            if (name.compare("CE")==0 && resname.compare("MET"))            {atm.props.setHydrophobic(true);continue;}
            if (name.compare("CE1")==0){
                if (appartient(resname,"HIS,HID,HIE"))                      {atm.props.setAromatic(true);
                    atm.props.setHydrophobic(false);continue;}
                else if (appartient(resname,"PHE, TYR"))                    {atm.props.setAromatic(true);continue;}}
            if (name.compare("CE2")==0 && appartient(resname, "PHE, TRP, TYR"))     {atm.props.setAromatic(true);continue;}
            if (name.compare("CE3")==0 && resname.compare("TRP")==0)                {atm.props.setAromatic(true);continue;}
            if (name.compare("CG")==0)
            {    if (appartient(resname,"LEU, MET, PRO")){ atm.props.setHydrophobic(true);continue;}
                else if (appartient(resname,"HIS, HID, HIE, PHE, TRP, TYR")){ atm.props.setAromatic(true);continue;}}
            if ((name.compare("CG1")==0 && appartient(resname,"ILE, VAL"))
                    || (name.compare("CG2")==0 && appartient(resname,"ILE, VAL"))){ atm.props.setHydrophobic(true);continue;}
            if (name.compare("CH2")==0 && resname.compare("TRP")==0) {atm.props.setAromatic(true);continue;}
            if (name.compare("CZ")==0){
                if (appartient(resname,"PHE, TYR"))                    { atm.props.setAromatic(true);continue;}}
            if (appartient(name, "CZ2, CZ3") && resname.compare("TRP")==0)          { atm.props.setAromatic(true);continue;}
            if (name.compare("N")==0)                               { atm.props.setDonor(true);continue;}
            if (appartient(name,"ND1, ND2") && appartient(resname,"HIS, HID, HIE"))          { atm.props.setDonor(true);atm.props.setAcceptor(true);continue;}
            if (name.compare("ND2")==0 && resname.compare("ASN")==0)                { atm.props.setDonor(true);continue;}
            if (name.compare("NE")==0 && appartient(resname,"ARG"))            { atm.props.setCation(true);continue;}
            if (name.compare("NE1")==0 && resname.compare("TRP")==0)                { atm.props.setDonor(true);continue;}
            if (name.compare("NE2")==0  ){
                if (resname.compare("GLN")==0)                          {atm.props.setDonor(true);continue;}
                else if (appartient(resname,"HIS, HID, HIE"))                   { atm.props.setDonor(true);
                    atm.props.setAcceptor(true);continue;}}
            if ((appartient(name,"NH1, NH2") && resname.compare("ARG")==0)
                    || (name.compare("NZ")==0 && resname.compare("LYS")==0))            {atm.props.setCation(true);continue;}
            if (name.compare("O")==0 && resname.compare("HOH")==0)                               {
                atm.props.setAcceptor(true);atm.props.setDonor(true);continue;}
            else if (name.compare("O")==0){
                atm.props.setAcceptor(true);continue;}

            if (appartient(name,"OD1, OD2") && resname.compare("ASP")==0)           {atm.props.setAnion(true);continue;}
            if (name.compare("OD1")==0 && resname.compare("ASN")==0)                {atm.props.setAcceptor(true);continue;}
            if (name.compare("OE1")==0 && resname.compare("GLN")==0)                {atm.props.setAcceptor(true);continue;}
            if (appartient(name,"OE1, OE2")  && resname.compare("GLU")==0)          {atm.props.setAnion(true);continue;}
            if ((name.compare("OG") ==0 && resname.compare("SER")==0)||
                    (name.compare("OG1")==0 && resname.compare("THR")==0)||
                    (name.compare("OH") ==0 && resname.compare("TYR")==0))              {atm.props.setDonor(true);
                atm.props.setAcceptor(true);continue;}
            if (name.compare("OXT")==0)                             {atm.props.setAnion(true);continue;}
            if (appartient(name, "SD, SG")){
                if (appartient(resname, "CYS, CYX")) {atm.props.setAcceptor(true);continue;}
                else if (resname.compare("MET")==0) {atm.props.setHydrophobic(true);continue;}}
            if (appartient(resname, "FE, MG, MG, MN, ZN, CO, CA, NA")) {atm.props.setCation(true),atm.props.setMetal(true);continue;}
            atm.props = p;
        }
    }

}
