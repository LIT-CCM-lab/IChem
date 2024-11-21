#include <iostream>
#include "headers/ICMole/box.h"
#include "headers/ICMole/complex.h"
#include "headers/ICMole/molecule.h"
#include "headers/ICMole/atom.h"
#include "headers/ICPars/molewriter.h"
//#undef ICHEM_DEBUG
//#define ICHEM_DEBUG
using namespace std;
using namespace ICMole;


Grid::Grid()
{
    setIniValues();


}

Grid::~Grid()
{

    if (AllBoxes.empty())return;
    for(size_t i = 0; i < AllBoxes.size(); ++i) delete AllBoxes[i];
}

void Grid::clear()
{
    for(size_t i = 0; i < AllBoxes.size(); ++i) delete AllBoxes[i];
    AllBoxes.clear();
}

void Grid::setIniValues()
{
    Residu::loadRules();
    Molecule::loadRules();

    for (unsigned int I=0; I< 9; I++)
    {
        rotMatrix[I]=0;
    }

    projLength=8;
    numPtRangeI=0;
    numPtRangeJ=0;
    numPtRangeK=0;
    maxCubeNum=0;
    adjLevel=5;
    boxStep=1.5;
    box_size=Coords(0,0,0);
    vect_i=Coords(0,0,0);
    vect_j=Coords(0,0,0);
    vect_k=Coords(0,0,0);
    X8=Coords(0,0,0);
    center=Coords(0,0,0);
    wHydrogen=false;
    wVolSite=false;
    status=0;


}


/**
 * @brief Grid::rotating
 * @param ini : coordinates that will be transformed
 *
 * Rotation Coords ini according to the generated rotation matrix:
 */
void Grid::rotating(Coords &ini) const
{
    ini.setCoords(round(100*(rotMatrix[0]*ini.x+rotMatrix[1]*ini.y+rotMatrix[2]*ini.z))/100,
                  round(100*(rotMatrix[3]*ini.x+rotMatrix[4]*ini.y+rotMatrix[5]*ini.z))/100,
                  round(100*(rotMatrix[6]*ini.x+rotMatrix[7]*ini.y+rotMatrix[8]*ini.z))/100);
//    ini.setCoords(rotMatrix[0]*ini.x+rotMatrix[1]*ini.y+rotMatrix[2]*ini.z,
//            rotMatrix[3]*ini.x+rotMatrix[4]*ini.y+rotMatrix[5]*ini.z,
//            rotMatrix[6]*ini.x+rotMatrix[7]*ini.y+rotMatrix[8]*ini.z);
}



/**
 * @brief Grid::prepareMatrix : Generate axis and rotation matrix
 * @param center : Coordinate of the grid center
 * @param pt1: Coordinates of the first point used to generate new axis
 * @param pt2: Coordinates of the second point used to generate new axis
 *
 *
 * The first step is to generate the new i,j,k axis centered on O' (center param).
 * To do so, it generate the first axis by calculating the cross-product
 * between pt1,pt2 and center. The second axis is made bt calculating the cross
 * product between the first axis, pt1 and the center. At last, the third axis
 * is the cross product between the center, the first and the second axis.
 *
 * To get the final i,j,k axis, we simply calculates their vectors by retrieving
 * the first, second and this axis from the center, leading the our new coordinates
 * system, vect_i,vect_j,vect_k
 *
 * To generate the rotation matrix from (O,x,y,z) to (O',i,j,k), we calculate
 * the inverse matrix of the matrix M=[vect_i,vect_j,vect_k].
 * http://fr.wikipedia.org/wiki/Matrice_inversible#Inversion_des_matrices_3_x_3
 *
 *
 */
void Grid::prepareMatrix(const Coords& center,
                         const Coords& pt1,
                         const Coords& pt2)
{


    //--------------------------------------------------------------------------
    //------------------------ STEP 1 - Getting the normal vector --------------
    //--------------------------------------------------------------------------

    // The normal vector between the (protein-Mole) barycenters vector
    // and the (inertial-barycenter) Mole vector will give us the first axis.
    //center.getNormal(&pt1,&pt2,&x_axis); // y axes
    const Coords x_axis= center.getNormal(pt1,pt2);
    // The second axis is given by the normal vector
    // between the x axis and the barycenter of the protein.
    // center.getNormal(&pt1,&x_axis,&y_axis); // z axes
    const Coords y_axis=center.getNormal(pt1,x_axis);
    // The last axis is the normal vector with x and y axis.
    //center.getNormal(&x_axis,&y_axis,&z_axis); // x axes
    const Coords z_axis =center.getNormal(x_axis,y_axis);
    //--------------------------------------------------------------------------
    //----------------------- STEP 2 - Getting the new O',i,j,k ----------------
    //--------------------------------------------------------------------------
//    cout << "X_AXIS " << x_axis.toString() << endl;
//    cout << "Y_AXIS " << y_axis.toString() << endl;
//    cout << "Z_AXIS " << z_axis.toString() << endl;

    vect_i.setCoords(x_axis.x-center.x,x_axis.y-center.y,x_axis.z-center.z);
    vect_j.setCoords(y_axis.x-center.x,y_axis.y-center.y,y_axis.z-center.z);
    vect_k.setCoords(z_axis.x-center.x,z_axis.y-center.y,z_axis.z-center.z);

    const double
            A = vect_i.x,
            D = vect_i.y,
            G = vect_i.z,
            B = vect_j.x,
            E = vect_j.y,
            H = vect_j.z,
            C = vect_k.x,
            F = vect_k.y,
            I = vect_k.z;
    const double detA = A*E*I+B*F*G+C*D*H-C*E*G-F*H*A-I*B*D;
    rotMatrix[0] = (E*I-F*H)/detA;
    rotMatrix[1] = (C*H-B*I)/detA;
    rotMatrix[2] = (B*F-C*E)/detA;
    rotMatrix[3] = (F*G-D*I)/detA;
    rotMatrix[4] = (A*I-C*G)/detA;
    rotMatrix[5] = (C*D-A*F)/detA;
    rotMatrix[6] = (D*H-E*G)/detA;
    rotMatrix[7] = (B*G-A*H)/detA;
    rotMatrix[8] = (A*E-B*D)/detA;
}





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// GRID CREATION ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


/**
 * @brief Grid::createComplexBoxes: scan the complex and generates the cube dimension
 * @param complex : Complex to create the grid over
 * @param box_Step : Length of each box
 */
void  Grid::createComplexBoxes (Complex& complex, const double& box_Step) throw(MoleExcept)
{

    if (box_Step <= 0)
        throw MoleExcept(1120101,
                         "Grid::createComplexBoxes",
                         "Box step is null or negative. MUST BE POSITIVE VALUE");

    // INITIALIZATION :

    boxStep = box_Step;
    // Maximal distance on negative x axis from grid center
    double dist1 = 0;
    // Maximal distance on positive x axis from grid center
    double dist2 = 0;
    // Maximal distance on negative y axis from grid center
    double dist3 = 0;
    // Maximal distance on positive y axis from grid center
    double dist4 = 0;
    // Maximal distance on negative z axis from grid center
    double dist5 = 0;
    // Maximal distance on positive z axis from grid center
    double dist6 = 0;

    double disttmp=0;

    // Cleaning the mass center of the complex :
    center.clear();

    // Coordinates of the inertial moment of the first molecule in the complex
    // Used to generate i,j,k axis
    Coords pt1;

    // Coordinates of the first atom of the first molecule in the complex
    // Used to generate i,j,k axis
    Coords pt2;

    // Number of molecules within the complex:
    double nMole=0;



    //--------------------------------------------------------------------------
    //--------- STEP 1 - SCANNING MOLECULES WITHIN COMPLEX - FILTERING ---------
    //--------------------------------------------------------------------------

    for (ItCMole it=complex.firstMole();it!= complex.lastMole();++it)
    {
        Molecule &molecule=**it;
        if (Molecule::Rules[molecule.getMoleType()]==MoleType::UNDEFINED)continue;
        if (molecule.numAtom() == 0)
            throw MoleExcept(1120102,
                             "Grid::createComplexBoxes",
                             "Molecule in complex has no atoms");
        try{
        center+= molecule.getFixpos();

        nMole+=1;
        if (nMole==1)
        {
            pt1=  molecule.calcInertialMoment();
            pt2= molecule.getAtom(0).fixpos;
        }
        }catch (MoleExcept &e)
        {
            std::cerr << e.getData() << " | " <<molecule.getName() << molecule.getResidu(0).getIdentifier()<<std::endl;
        }
    }

    if (nMole==0)
        throw MoleExcept(1120103,
                         "Grid::createComplexBoxes",
                         "No molecule in complex");
    center /= nMole;


#ifdef ICHEM_DEBUG
    cout << "CENTER : " << center.toString()<<endl
         << "POINT1 : " << pt1.toString()<<endl
         << "POINT2 : " << pt2.toString()<<endl;
#endif


    //----------------------------------------------------------------------
    //--------------- STEP 2 - GETTING MAIN GRID COORDINATES ---------------
    //----------------------------------------------------------------------

    prepareMatrix(center,pt1,pt2);

    Coords temp, nul(0,0,0);

    // Scanning each molecular type:
    for (size_t mtype=0; mtype < NB_MOLETYPE;mtype++)
    {

        // Avoiding unselected one :
        if (Molecule::Rules[mtype]==MoleType::UNDEFINED)continue;

        // Scanning each molecule for this molecular type :
        for (ItCMole itM = complex.firstMole(mtype);
             itM!= complex.lastMole(mtype);
             itM++)
        {
            Molecule &molecule=**itM;

            // Scanning each residu for this molecule :
            for (ItCRes itR = molecule.firstResidu();
                 itR!= molecule.lastResidu();
                 ++itR)
            {
                Residu &residu=**itR;
                // Avoiding unselected one :

                if (Residu::Rules[mtype][residu.getResType()]==MoleType::UNDEFINED)continue;
                // Scanning each atom of this residu :
                for (ItCAtom itA = residu.firstAtom();
                     itA!= residu.lastAtom();
                     itA++)
                {
                    Atom &atom =**itA;
                    // Avoiding ununsed one :
                    if (!atom.isUsed())continue;

                    // Calculating coordinates on new axis :
                    atom.rotpos=atom.fixpos;
                    atom.rotpos-=center;
                    rotating(atom.rotpos);

                    // To get the projected distance onto an axis between the atom
                    // and the grid center, we need to calculate the projection,
                    // For x axis : x value for the atom * vect_i
                    // And calculate the distance of this position against the center
                    // Since we need the extrema min and max for each axis, we
                    // check that for every coordinates :

#ifdef ICHEM_DEBUG
                    cout << atom.getIdentifier()<<" :: "<<atom.rotpos.toString()<<endl;
#endif
                    temp=(vect_i*atom.rotpos.x);
                    disttmp = temp.calcDist(nul);
                    if (atom.rotpos.x > 0 && disttmp > dist1)
                    { dist1 = disttmp;
#ifdef ICHEM_DEBUG
                        cout <<"  ||===> dist1: "<< disttmp<<endl;
#endif
                    }
                    if (atom.rotpos.x < 0 && disttmp > dist2){
                        dist2 = disttmp;
#ifdef ICHEM_DEBUG
                        cout <<"  ||===> dist2: "<< disttmp<<endl;
#endif
                    }
                    temp=vect_j*atom.rotpos.y;
                    disttmp = temp.calcDist(nul);
                    if (atom.rotpos.y> 0 && disttmp > dist3)
                    { dist3 = disttmp;
#ifdef ICHEM_DEBUG
                        cout <<"  ||===> dist3: "<< disttmp<<endl;
#endif
                    }
                    if (atom.rotpos.y< 0 && disttmp > dist4)
                    { dist4 = disttmp;
#ifdef ICHEM_DEBUG
                        cout <<"  ||===> dist4: "<< disttmp<<endl;
#endif
                    }
                    temp=vect_k*atom.rotpos.z;
                    disttmp = temp.calcDist(nul);
                    if (atom.rotpos.z > 0 && disttmp > dist5)
                    { dist5 = disttmp;
#ifdef ICHEM_DEBUG
                        cout <<"  ||===> dist5: "<< disttmp<<endl;
#endif
                    }
                    if (atom.rotpos.z <0 && disttmp > dist6)
                    { dist6 = disttmp;
#ifdef ICHEM_DEBUG
                        cout <<"  ||===> dist6: "<< disttmp<<endl;
#endif
                    }
                }


            }// END RESIDU
        }// END MOLECULE
    }// END MTYPE


    //Adding margin :
    dist1 += projLength;
    dist2 += projLength;
    dist3 += projLength;
    dist4 += projLength;
    dist5 += projLength;
    dist6 += projLength;

    if ((dist1+dist2)/box_Step < 3
            ||(dist3+dist4)/box_Step < 3
            ||(dist5+dist6)/box_Step < 3)
        throw MoleExcept(1120104,
                         "Grid::createComplexBoxes",
                         "Number of small boxes less than 3.\
                         Must increase box_size or reduce box_step");




                         try{
                             genBoxes( dist1,dist2, dist3,dist4,dist5,dist6);
                         }catch (MoleExcept &e)
    {
                             e.addTrace("Grid::createComplexBoxes");
                             throw;
                         }
}

void  Grid::createComplexBoxes  (Complex&  protein,
                                 Molecule & ligand,
                                 const double & box_Step,
                                 const Coords& box_size)
throw(MoleExcept)
{
    boxStep=box_Step;

    if (box_size.x/box_Step < 3
            ||box_size.y/box_Step < 3
            ||box_size.z/box_Step < 3)
        throw MoleExcept(1120205,
                         "Grid::createComplexBoxes",
                         "Number of small boxes less than 3. Must increase box_size or reduce box_step");


    Coords pt1;double nMole=0;
    ligand.setUse(false);
    for (ItCMole it = protein.firstMole();
         it !=protein.lastMole();
         it++)
    {
        Molecule &mole=**it;

        try
        {
            pt1+=mole.getFixpos();
            nMole++;
        }catch (MoleExcept &e)
        {

        }

    }
    pt1/=nMole;
    ligand.setUse(true);


    //---------------------------------------------------------------------------
    //----------- STEP 1 - SELECTION OF 3 POINTS FOR NEW AXIS CALCULATION -------
    //---------------------------------------------------------------------------

    center           = ligand.getFixpos();
    const Coords pt2 = ligand.calcInertialMoment(false);
    if (pt1 == center) pt1+=1;

    // Maximal distance on negative x axis from grid center
    const double dist1 = box_size.x/2+projLength;
    // Maximal distance on positive x axis from grid center
    const double dist2 = box_size.x/2+projLength;
    // Maximal distance on negative y axis from grid center
    const double dist3 = box_size.y/2+projLength;
    // Maximal distance on positive y axis from grid center
    const double dist4 = box_size.y/2+projLength;
    // Maximal distance on negative z axis from grid center
    const double dist5 = box_size.z/2+projLength;
    // Maximal distance on positive z axis from grid center
    const double dist6 = box_size.z/2+projLength;
#ifdef ICHEM_DEBUG
    cout << "CENTER : " << center.toString()<<endl
         << "POINT1 : " << pt1.toString()<<endl
         << "POINT2 : " << pt2.toString()<<endl;
    cout << "CURRENT DISTANCES : "<< dist1<< " "<<dist2<< " "<<dist3<< " "
         << dist4<< " "<<dist5<< " "<<dist6<<endl
         << "PROJ LENGTH : "<< projLength<<endl;
#endif
    prepareMatrix(center,pt1,pt2);

    try{
        genBoxes(dist1,dist2, dist3,dist4,dist5,dist6);
    }catch (MoleExcept &e)
    {
        e.addTrace("Grid::createComplexBoxes");
        throw;
    }
}

/**
 * @brief Grid::createComplexBoxes
 * @param protein = Protein molecule used for creating the grid
 * @param ligand : Ligand molecule
 * @param box_Step: width of each small boxes (Default 1.5Angstroem)
 * @param box_size : width of the main boxes
 * @throw 1120201 - When no protein given or no atoms in it
 * @throw 1120202 - When no ligand given or no atoms in it
 * @throw 1120203 - Box step is below 0 or equals 0
 * @throw 1120204 - Not enough boxes on one axis. Minimal value is 5.
 * The function uses both protein and ligand coordinates to create a grid,
 * centered at the ligand barycenter. Since grid must be coordinate orientation
 * independant, the grid class creates a new coordinate reference (O',i,j,k)
 * implemented with the use of rotpos coords.The first step is to select 3
 * points for the new axis. The first is the ligand barycenter, the second is
 * the protein barycenter, and the last one is the inertial moment of the ligand.
 * From this, the function calculate normal vector to get i, j and k axis.\n
 * The function automatically generates a cubic grid, where the edge length on
 * each i,j and k axis is defined by box_size Coords.\n At last, the function
 * calls the private function genBoxes() to generate small boxes.
 */
void  Grid::createComplexBoxes  (Molecule&  protein,
                                 Molecule & ligand,
                                 const double & box_Step,
                                 const Coords& box_size) throw(MoleExcept)
{
    if (protein.numAtom() == 0)
        throw MoleExcept(1120201,
                         "Grid::createComplexBoxes",
                         "No protein given");
    if (ligand.numAtom() == 0)
        throw MoleExcept(1120202,
                         "Grid::createComplexBoxes",
                         "No ligand given");
    if (box_Step <= 0)
        throw MoleExcept(1120204,
                         "Grid::createComplexBoxes"
                         ,"Box step is null or negative. MUST BE POSITIVE VALUE");
    boxStep=box_Step;

    if (box_size.x/box_Step < 3
            ||box_size.y/box_Step < 3
            ||box_size.z/box_Step < 3)
        throw MoleExcept(1120205,
                         "Grid::createComplexBoxes",
                         "Number of small boxes less than 3. Must increase box_size or reduce box_step");
    //---------------------------------------------------------------------------
    //----------- STEP 1 - SELECTION OF 3 POINTS FOR NEW AXIS CALCULATION -------
    //---------------------------------------------------------------------------

    center           = ligand.getFixpos();
    Coords pt1       = protein.getFixpos();
    const Coords pt2 = ligand.calcInertialMoment(false);
    if (pt1 == center) pt1+=1;

    // Maximal distance on negative x axis from grid center
    const double dist1 = box_size.x/2+projLength;
    // Maximal distance on positive x axis from grid center
    const double dist2 = box_size.x/2+projLength;
    // Maximal distance on negative y axis from grid center
    const double dist3 = box_size.y/2+projLength;
    // Maximal distance on positive y axis from grid center
    const double dist4 = box_size.y/2+projLength;
    // Maximal distance on negative z axis from grid center
    const double dist5 = box_size.z/2+projLength;
    // Maximal distance on positive z axis from grid center
    const double dist6 = box_size.z/2+projLength;
#ifdef ICHEM_DEBUG
    cout << "CENTER : " << center.toString()<<endl
         << "POINT1 : " << pt1.toString()<<endl
         << "POINT2 : " << pt2.toString()<<endl;
    cout << "CURRENT DISTANCES : "<< dist1<< " "<<dist2<< " "<<dist3<< " "
         << dist4<< " "<<dist5<< " "<<dist6<<endl
         << "PROJ LENGTH : "<< projLength<<endl;
#endif
    prepareMatrix(center,pt1,pt2);

    try{
        genBoxes(dist1,dist2, dist3,dist4,dist5,dist6);
    }catch (MoleExcept &e)
    {
        e.addTrace("Grid::createComplexBoxes");
        throw;
    }
}


/**
 * @brief Grid::createMolecularBoxes
 * @param molecule : Molecule to enclose inside a grid
 * @param box_Step :width of each small box (Default 1.5Angstroem)
 * The function creates a grid around the given molecule Mole. Since grid must
 * be coordinate orientation independant, the grid class creates a new
 * coordinate reference (O',i,j,k) implemented with the use of rotpos coords.
 * The first step is to select 3 points for the new axis. The first is the
 * molecule barycenter, the second is the molecule inertial moment, and the last
 * one is the coordinates of the first atom. From this, the second step is to
 * calculate normal vector to get i, j and k axis.\n
 * The next step is to look over coordinates to search for minima and maxima
 * on each i, j and k axis.\n
 * A margin is added to theses values (Default 3 boxes of size box_Step).\n
 * At last, the function calls the private function genBoxes() to generate
 * small boxes.\n Therefore, the grid might not be cubic.
 *  @throw 1120301 - When no molecule given or no atoms in it
 *  @throw 1120302 - Box step is below 0 or equals 0
 *  @throw 1120303 - Selection rules for this molecule is set to 'ignore'
 *  @throw 1120303 - Not enough boxes on one axis. Minimal value is 5.
 */
void  Grid::createMolecularBoxes(Molecule& molecule,const double& box_Step) throw(MoleExcept)
{
    if (molecule.numAtom() == 0)
        throw MoleExcept(1120301,
                         "Grid::createMolecularBoxes",
                         "No Mole given");
    if (box_Step <= 0)
        throw MoleExcept(1120302,
                         "Grid::createMolecularBoxes",
                         "Box step is null or negative. MUST BE POSITIVE VALUE");

    if (Molecule::Rules[molecule.getMoleType()]==MoleType::UNDEFINED)
        throw MoleExcept(1120303,
                         "Grid::createMolecularBoxes",
                         "Selection rules for this molecule is set to 'ignore'");

    boxStep=box_Step;
    // Maximal distance on negative x axis from grid center
    double dist1 = 0;
    // Maximal distance on positive x axis from grid center
    double dist2 = 0;
    // Maximal distance on negative y axis from grid center
    double dist3 = 0;
    // Maximal distance on positive y axis from grid center
    double dist4 = 0;
    // Maximal distance on negative z axis from grid center
    double dist5 = 0;
    // Maximal distance on positive z axis from grid center
    double dist6 = 0;
    double disttmp=0;

    //--------------------------------------------------------------------------
    //------ STEP 1 - SELECTION OF 3 POINTS FOR NEW AXIS CALCULATION -----------
    //--------------------------------------------------------------------------

    center = molecule.getFixpos();
//    Coords pt1    = molecule.calcInertialMoment();
//    center = molecule.getAtom(50).fixpos;
//    Coords pt2 = molecule.getalpha(10);
//    Coords pt1 = molecule.getalpha(0);
    Coords pt1    = molecule.calcInertialMoment();
//    center = molecule.getAtom(50).fixpos;
    Coords pt2 = molecule.getAtom(0).fixpos;
    if (molecule.SizeRes()>9){
        pt2 = molecule.getalpha(10);
        pt1 = molecule.getalpha(0);
    }
//    molecule.getextrem();
//    Coords pt1    = molecule.extremH;
//    Coords pt2    = molecule.extremL;
//    round all coordinate trying to minimise differnece between grid
//    Coords pt1 = molecule.getAtom(0).fixpos;
//    Coords pt2 = molecule.getAtom(100).fixpos;
//    double CX = round(1000*center.x)/1000;  double CY = round(1000*center.y)/1000;  double CZ = round(1000*center.z)/1000;
//    double p1X = round(1000*pt1.x)/1000;    double p1Y = round(1000*pt1.y)/1000;    double p1Z = round(1000*pt1.z)/1000;
//    double p2X = round(1000*pt2.x)/1000;    double p2Y = round(1000*pt2.y)/1000;    double p2Z = round(1000*pt2.z)/1000;

//    center.setCoords(CX,CY,CZ);
//    pt1.setCoords(p1X,p1Y,p1Z);
//    pt2.setCoords(p2X,p2Y,p2Z);
#ifdef ICHEM_DEBUG
    cout << "C-PT1 "<< center.calcDist(pt1)<< endl
         << "C-pt2 " << center.calcDist(pt2) <<endl
         <<"pt1-pt2" << pt1.calcDist(pt2)<< endl;
    cout << "CENTER : " << center.toString()<<endl
         << "POINT1 : " << pt1.toString()<<endl
         << "POINT2 : " << pt2.toString()<<endl;
#endif
    prepareMatrix(center,pt1,pt2);

    Coords temp;
    const Coords nul(0,0,0);
    for (ItCAtom itA = molecule.firstAtom(); itA != molecule.lastAtom();++itA)
    {
        Atom &atom =**itA;
        if (!atom.isUsed())continue;
        if (atom.isHydrogen())continue;
        // Avoiding unselected one :
        if (Residu::Rules[molecule.getMoleType()]
                [atom.getResidu()->getResType()]==MoleType::UNDEFINED)continue;

        // Calculating coordinates on new axis :
        atom.rotpos=atom.fixpos;
        atom.rotpos-=center;
        rotating(atom.rotpos);

        // To get the projected distance onto an axis between the atom
        // and the grid center, we need to calculate the projection,
        // For x axis : x value for the atom * vect_i
        // And calculate the distance of this position against the center
        // Since we need the extrema min and max for each axis, we
        // check that for every coordinates :
        temp=(vect_i*atom.rotpos.x);
        disttmp = temp.calcDist(nul);
        if (atom.rotpos.x > 0 && disttmp > dist1)
        { dist1 = disttmp;
#ifdef ICHEM_DEBUG
            cout <<atom.getIdentifier()<<"::"
                <<atom.rotpos.toString()<<"\t  ||===> dist1: "<< disttmp<<endl;
#endif
        }
        if (atom.rotpos.x < 0 && disttmp > dist2){
            dist2 = disttmp;
#ifdef ICHEM_DEBUG
            cout <<atom.getIdentifier()<<"::"
                <<atom.rotpos.toString()<<"\t  ||===> dist2: "<< disttmp<<endl;
#endif
        }
        temp=vect_j*atom.rotpos.y;
        disttmp = temp.calcDist(nul);
        if (atom.rotpos.y> 0 && disttmp > dist3)
        { dist3 = disttmp;
#ifdef ICHEM_DEBUG
            cout <<atom.getIdentifier()<<"::"
                <<atom.rotpos.toString()<<"\t  ||===> dist3: "<< disttmp<<endl;
#endif
        }
        if (atom.rotpos.y< 0 && disttmp > dist4)
        { dist4 = disttmp;
#ifdef ICHEM_DEBUG
            cout <<atom.getIdentifier()<<"::"
                <<atom.rotpos.toString()<<"\t  ||===> dist4: "<< disttmp<<endl;
#endif
        }
        temp=vect_k*atom.rotpos.z;
        disttmp = temp.calcDist(nul);
        if (atom.rotpos.z > 0 && disttmp > dist5)
        { dist5 = disttmp;
#ifdef ICHEM_DEBUG
            cout <<atom.getIdentifier()<<"::"
                <<atom.rotpos.toString()<<"\t  ||===> dist5: "<< disttmp<<endl;
#endif
        }
        if (atom.rotpos.z <0 && disttmp > dist6)
        { dist6 = disttmp;
#ifdef ICHEM_DEBUG
            cout <<atom.getIdentifier()<<"::"
                <<atom.rotpos.toString()<<"\t  ||===> dist6: "<< disttmp<<endl;
#endif
        }
    }

    dist1 = round(dist1*1000)/1000;
    dist2 = round(dist2*1000)/1000;
    dist3 = round(dist3*1000)/1000;
    dist4 = round(dist4*1000)/1000;
    dist5 = round(dist5*1000)/1000;
    dist6 = round(dist6*1000)/1000;
#ifdef ICHEM_DEBUG
    cout << "CURRENT DISTANCES : "<< dist1<< " "<<dist2<< " "<<dist3<< " "
         << dist4<< " "<<dist5<< " "<<dist6<<endl
         << "PROJ LENGTH : "<< projLength<<endl;
#endif

    //Adding margin :

    dist1 += projLength+3;
    dist2 += projLength+3;
    dist3 += projLength+3;
    dist4 += projLength+3;
    dist5 += projLength+3;
    dist6 += projLength+3;

    if ((dist1+dist2)/box_Step < 3
            ||(dist3+dist4)/box_Step < 3
            ||(dist5+dist6)/box_Step < 3)
        throw MoleExcept(1120303,
                         "Grid::createMolecularBoxes",
                         "Number of small boxes less than 3. \
                         Must increase box_size or reduce box_step");




                         try{
                             genBoxes(dist1,dist2, dist3,dist4,dist5,dist6);
                         }catch (MoleExcept &e)
    {
                             e.addTrace("Grid::createMolecularBoxes");
                             throw;
                         }
}

/**
 * @brief Grid::genBoxes: Generates the grid boxes:
 * @param dist1: Maximal distance on negative x axis from grid center
 * @param dist2: Maximal distance on positive x axis from grid center
 * @param dist3: Maximal distance on negative y axis from grid center
 * @param dist4: Maximal distance on positive y axis from grid center
 * @param dist5: Maximal distance on negative z axis from grid center
 * @param dist6: Maximal distance on positive z axis from grid center
 * @throw 1120401 : Bad allocation exception - exiting
 */
void Grid::genBoxes(const double& dist1,
                    const double& dist2,
                    const double& dist3,
                    const double& dist4,
                    const double& dist5,
                    const double& dist6 )
{
    // Defines grid length over all axis :
    box_size.setCoords(dist1+dist2,dist3+dist4,dist5+dist6);

    // X8 is the lower coordinates of the grid that we will use as box starter :
    X8 = center -vect_i*dist2-vect_j*dist4-vect_k*dist6;

    const int BxMargin = static_cast<int>(floor(projLength/boxStep));

    // Number of cubes over each axis :
    numPtRangeI = static_cast<int>(ceil(box_size.x/boxStep));
    numPtRangeJ = static_cast<int>(ceil(box_size.y/boxStep));
    numPtRangeK = static_cast<int>(ceil(box_size.z/boxStep));

    // Number of boxes in the grid :
    maxCubeNum = numPtRangeI*numPtRangeJ*numPtRangeK;
    AllBoxes.reserve(maxCubeNum);

#ifdef ICHEM_DEBUG
    const Coords X7 = center -vect_i*dist2-vect_j*dist4+vect_k*dist5;
    const Coords X6 = center -vect_i*dist2+vect_j*dist3-vect_k*dist6;
    const Coords X5 = center -vect_i*dist2+vect_j*dist3+vect_k*dist5;
    const Coords X4 = center +vect_i*dist1-vect_j*dist4-vect_k*dist6;
    const Coords X3 = center +vect_i*dist1-vect_j*dist4+vect_k*dist5;
    const Coords X2 = center +vect_i*dist1+vect_j*dist3-vect_k*dist6;
    const Coords X1 = center +vect_i*dist1+vect_j*dist3+vect_k*dist5;
    cout << "## BOX STEP : "<<boxStep<<endl;
    cout << "## BOX SIZE : "<<box_size.toString()<<endl;
    cout << "## Num cubes per axes : (" << numPtRangeI << " " << numPtRangeJ << " " << numPtRangeK << ")"<<endl;
    cout << "## Max num cubes : " << maxCubeNum<<endl;


    Molecule boxmol(8,12,MoleType::OTHER);
    boxmol.setName("BOX");
    Atom &atm  = boxmol.addAtom("C",X1,             "CA","C.3");
    Atom &atm1 = boxmol.addAtom("C",X2,             "CB","C.3");
    Atom &atm2 = boxmol.addAtom("C",X3,             "CC","C.3");
    Atom &atm3 = boxmol.addAtom("C",X4,             "CD","C.3");
    Atom &atm4 = boxmol.addAtom("C",X5,             "CE","C.3");
    Atom &atm5 = boxmol.addAtom("C",X6,             "CF","C.3");
    Atom &atm6 = boxmol.addAtom("C",X7,             "CG","C.3");
    Atom &atm7 = boxmol.addAtom("C",X8,             "CH","C.3");
    Atom &atm8 = boxmol.addAtom("C",center,         "CL","C.3");
    Atom &atm11= boxmol.addAtom("C",(center+vect_i),"CI","C.3");
    Atom &atm12= boxmol.addAtom("C",(center+vect_j),"CJ","C.3");
    Atom &atm13= boxmol.addAtom("C",(center+vect_k),"CK","C.3");
    boxmol.addBond(atm ,atm1,BondType::SINGLE); boxmol.addBond(atm ,atm2,BondType::SINGLE); boxmol.addBond(atm ,atm4,BondType::SINGLE);
    boxmol.addBond(atm1,atm3,BondType::SINGLE); boxmol.addBond(atm1,atm5,BondType::SINGLE);
    boxmol.addBond(atm2,atm6,BondType::SINGLE); boxmol.addBond(atm2,atm3,BondType::SINGLE);
    boxmol.addBond(atm5,atm7,BondType::SINGLE); boxmol.addBond(atm3,atm7,BondType::SINGLE);boxmol.addBond(atm6,atm7,BondType::SINGLE);
    boxmol.addBond(atm4,atm5,BondType::SINGLE); boxmol.addBond(atm4,atm6,BondType::SINGLE);
    boxmol.addBond(atm8,atm11,BondType::TRIPLE);
    boxmol.addBond(atm8,atm12,BondType::TRIPLE);boxmol.addBond(atm8,atm13,BondType::TRIPLE);
    MoleWriter mw("box.mol2",FileFormat::MOL2);
    mw.writeMOL2(&boxmol);


    cout << "DEBUG - Generating boxes";
#endif

    status=1;


    // Now we have the limits of the main box. We are now able to discretize the main box space in smaller ones.
    // X8 is the "smallest" point in the cube, so we will use it as reference for initial position.
    // The triple loop is here to generate the cubes in all the 3D directions.
    // i represent the number of the cube regarding the x' axis.
    // j represent the number of the cube regarding the y' axis.
    // k represent the number of the cube regarding the z' axis.

    Box *bx=NULL;
    try{
        for (int k=0; k < numPtRangeK; k++)
            for (int j=0; j < numPtRangeJ; j++)
                for (int i=0; i < numPtRangeI; i++)
                {

                    bx = new Box(floor(0.5+100*(X8.x+i*boxStep*vect_i.x+j*boxStep*vect_j.x+k*boxStep*vect_k.x))/100,
                                 floor(0.5+100*(X8.y+i*boxStep*vect_i.y+j*boxStep*vect_j.y+k*boxStep*vect_k.y))/100,
                                 floor(0.5+100*(X8.z+i*boxStep*vect_i.z+j*boxStep*vect_j.z+k*boxStep*vect_k.z))/100);

                    bx->rotpos.setCoords(i,j,k);
                    bx->Props.clear();
                    if ((k < BxMargin || k > numPtRangeK-BxMargin)||
                            (i < BxMargin || i > numPtRangeI-BxMargin)||
                            (j < BxMargin || j > numPtRangeJ-BxMargin)
                            ) bx->setMargin(true);
                    else bx->Props.setNoInfo(true);
                    bx->setId(i+j*numPtRangeI+k*numPtRangeI*numPtRangeJ);
                    AllBoxes.push_back(bx)  ;
                }
    }catch (std::bad_alloc &e)
    {
        throw MoleExcept(1120401,
                         "Grid::genBoxes",
                         "Bad allocation exception - exiting ");
    }
#ifdef ICHEM_DEBUG

    printInFile("GENERATION",AllBoxes);
#endif

}







void Grid::createParentGrid(const Grid& grid, const double&box_Step) throw(MoleExcept)
{

    if (fmod(box_Step,grid.boxStep) != 0)
        throw MoleExcept(1120501,
                         "Grid::createParentGrid",
                         "Given boxStep is not a factor of given grid boxstep");
    if (grid.AllBoxes.empty() || grid.status==0)
        throw MoleExcept(1120502,
                         "Grid::createParentGrid",
                         "Parent grid cannot be generated without child grid box generation");
    if (grid.status==1)
        throw MoleExcept(1120503,
                         "Grid::createParentGrid",
                         "Parent grid color cannot be generated without child grid color");
    for (unsigned int I=0; I< 9; I++)
    {
        rotMatrix[I]=grid.rotMatrix[I];
    }


    projLength=grid.projLength;
    box_size = grid.box_size;
    boxStep = box_Step;
    numPtRangeI = static_cast<int>(ceil(box_size.x/boxStep))+1;
    numPtRangeJ = static_cast<int>(ceil(box_size.y/boxStep))+1;
    numPtRangeK = static_cast<int>(ceil(box_size.z/boxStep))+1;
    maxCubeNum=numPtRangeI*numPtRangeJ*numPtRangeK;
    const int BxMargin = static_cast<int>(floor(projLength/boxStep));
    AllBoxes.reserve(maxCubeNum);
    adjLevel=grid.adjLevel;

    vect_i=grid.vect_i;
    vect_j=grid.vect_j;
    vect_k=grid.vect_k;
    X8=grid.X8;
    center=grid.center;
    wHydrogen=grid.wHydrogen;
    wVolSite=grid.wVolSite;
    status=1;


    // Now we have the limits of the main box. We are now able to discretize the main box space in smaller ones.
    // X8 is the "smallest" point in the cube, so we will use it as reference for initial position.
    // The triple loop is here to generate the cubes in all the 3D directions.
    // i represent the number of the cube regarding the x' axis.
    // j represent the number of the cube regarding the y' axis.
    // k represent the number of the cube regarding the z' axis.

    Box *bx=NULL;
    try{
        for (int k=0; k < numPtRangeK; k++)
            for (int j=0; j < numPtRangeJ; j++)
                for (int i=0; i < numPtRangeI; i++)
                {

                    bx = new Box(X8.x+i*boxStep*vect_i.x+j*boxStep*vect_j.x+k*boxStep*vect_k.x,
                                 X8.y+i*boxStep*vect_i.y+j*boxStep*vect_j.y+k*boxStep*vect_k.y,
                                 X8.z+i*boxStep*vect_i.z+j*boxStep*vect_j.z+k*boxStep*vect_k.z);

                    bx->rotpos.setCoords(i,j,k);
                    bx->Props.clear();
                    if ((k < BxMargin || k > numPtRangeK-BxMargin)||
                            (i < BxMargin || i > numPtRangeI-BxMargin)||
                            (j < BxMargin || j > numPtRangeJ-BxMargin)
                            ) bx->setMargin(true);
                    else bx->Props.setNoInfo(true);
                    bx->setId(i+j*numPtRangeI+k*numPtRangeI*numPtRangeJ);
                    AllBoxes.push_back(bx)  ;
                }
    }catch (std::bad_alloc &e)
    {
        throw MoleExcept(1120401,
                         "Grid::genBoxes",
                         "Bad allocation exception - exiting ");
    }
#ifdef ICHEM_DEBUG

    printInFile("PAR_GEN",AllBoxes);
#endif


    status=1;


    const int diff_box=(int)floor(boxStep/grid.boxStep);

    for (ItCBox itb =  AllBoxes.begin();
         itb != AllBoxes.end();
         ++itb)
    {
        Box &bx = **itb;
#ifdef ICHEM_DEBUG
        cout << bx.getId()<<" " <<bx.fixpos.toString()<< " " << bx.rotpos.toString()<<endl;
#endif
        const int iR=static_cast<int>(bx.rotpos.x*diff_box);
        const int jR=static_cast<int>(bx.rotpos.y*diff_box);
        const int kR=static_cast<int>(bx.rotpos.z*diff_box);

        if (kR > grid.numPtRangeK)continue;
        if (jR > grid.numPtRangeJ)continue;
        if (iR > grid.numPtRangeI)continue;
        for (int k=kR; k <kR+diff_box;++k){if (k <0 || k >= grid.numPtRangeK) continue;
            for (int j=jR; j <jR+diff_box;++j){if (j <0 || j >= grid.numPtRangeJ) continue;
                for (int i=iR; i <iR+diff_box;++i){if (i <0 || i >= grid.numPtRangeI) continue;


                    const int cubePos =  i+j*grid.numPtRangeI+k*grid.numPtRangeI*grid.numPtRangeJ;
                    if (cubePos<0 || cubePos >=grid.maxCubeNum)continue;
                    Box &bxCGrid=*grid.AllBoxes.at(cubePos);
#ifdef ICHEM_DEBUG
                    cout << "  ||--> "<<i<< " " <<j<<" "<<k<<" "<< cubePos <<" " << bxCGrid.getId()<< " " << bxCGrid.rotpos.toString()<< "\t"<< bxCGrid.fixpos.calcDist(bx.fixpos)<<endl;
#endif
                    for (ItCAtom itCA = bxCGrid.firstAtom();
                         itCA!= bxCGrid.lastAtom();
                         itCA++)
                    {

                        (*itCA)->addBox(this,&bx);

#ifdef ICHEM_DEBUG
                        cout <<   "     ||--> "<<(*itCA)->getIdentifier()<< (*itCA)->fixpos.calcDist(bxCGrid.fixpos)<< " " << (*itCA)->fixpos.calcDist(bx.fixpos)<<endl;
#endif
                    }
                    if (bxCGrid.Props.isInCavity()) bx.Props.setInCavity(true);
                    if (bxCGrid.Props.isInLigand()) bx.Props.setInLigand(true);
                    if (bxCGrid.Props.isInProtein())    bx.Props.setInProtein(true);
                    if (bxCGrid.Props.isNoInfo())       bx.Props.setNoInfo(true);
                    if (bxCGrid.Props.isCloseProtein()) bx.Props.setCloseProtein(true);
                    if (bxCGrid.Props.isOutCavity())    bx.Props.setOutCavity(true);



                }}}

    }
    status=2;
#ifdef ICHEM_DEBUG

    printInFile("PAR_COL",AllBoxes);
#endif



}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////// PRINTERS //////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/**  \fn void Grid::printInFile(const string f) const
    \param f : Name of files
    \param liste : liste of cube to be outputed
    \param wBorder : output cube borders
    \brief Convert grid cubes in MOL2 files

    Group cubes according to their properties (i.e. SIT, OUT, PROT, VID, IN_SIT) and output a MOL2 file for each group.


    */
void Grid::printInFile(const string& f, const BoxList& liste, const bool& wBorder, const bool& onlyUse) const
{
    //  Molecule molen7(MoleType::OTHER);
    Molecule molen (MoleType::OTHER);Residu &R_PROT=molen.addResidu("X",1,"PRO");//Residu &X_PROT=molen7.addResidu("X",1,"PRO");
    Molecule molen2(MoleType::OTHER);Residu &R_SIT=molen2.addResidu("X",1,"SIT");//Residu &X_SIT=molen7.addResidu("X",1,"SIT");
    Molecule molen3(MoleType::OTHER);Residu &R_OUT=molen3.addResidu("X",1,"OUT");//Residu &X_OUT=molen7.addResidu("X",1,"OUT");
    Molecule molen4(MoleType::OTHER);Residu &R_VID=molen4.addResidu("X",1,"VID");//Residu &X_VID=molen7.addResidu("X",1,"VID");
    Molecule molen5(MoleType::OTHER);Residu &R_ISIT=molen5.addResidu("X",1,"INS");//Residu &X_ISIT=molen7.addResidu("X",1,"INS");
    Molecule molen6(MoleType::OTHER);Residu &R_DUM=molen6.addResidu("X",1,"DUM"); //Residu &X_DUM=molen7.addResidu("X",1,"DUM");
    Molecule molen8(MoleType::OTHER);Residu &R_IN=molen8.addResidu("X",1,"ISI");  //Residu &X_IN=molen7.addResidu("X",1,"ISI");

    int cube_pos;



    string s;

    for (ItCBox it=liste.begin(); it != liste.end(); it++)
    {
        Box *bx = *it;
        if (onlyUse && !bx->getUse()) continue;
        //if (cube_pos >= Max_Cube_Num)continue;
        cube_pos = bx->getId();

        if (bx->isMargin())               {   molen6.addAtom("H",bx->fixpos,"H","Du", &R_DUM); } //molen7.addAtom("H",bx->fixpos,"H","Du" ,&X_DUM);}
        if (bx->Props.isInProtein())      {    molen.addAtom("O",bx->fixpos,"O","O.2",&R_PROT);}// molen7.addAtom("O",bx->fixpos,"O","O.2",&X_PROT);}
        if (bx->Props.isCloseProtein())   {   molen3.addAtom("C",bx->fixpos,"C","C.3",&R_OUT); } //molen7.addAtom("C",bx->fixpos,"C","C.3",&X_OUT);}
        if (bx->Props.isOutCavity())      {   molen2.addAtom("N",bx->fixpos,"N","N.3",&R_SIT); } //molen7.addAtom("N",bx->fixpos,"N","N.3",&X_SIT);}
        if (bx->Props.isInCavity())       {   molen8.addAtom("S",bx->fixpos,"S","S.2",&R_IN);  } //molen7.addAtom("S",bx->fixpos,"S","S.2",&X_IN);}
        if (bx->Props.isNoInfo())         {   molen4.addAtom("P",bx->fixpos,"P","P.2",&R_VID); } //molen7.addAtom("P",bx->fixpos,"P","P.2",&X_VID);}
        if (bx->Props.isInLigand())       {   molen5.addAtom("I",bx->fixpos,"I","I", &R_ISIT); } //molen7.addAtom("I",bx->fixpos,"I","I",&X_ISIT);}
        if (wBorder)
        {
            /*   Atom *atmX[8];unsigned int iix=0;
        for (double x=-boxStep/2; x<=boxStep/2;x+=boxStep)
        for (double y=-boxStep/2; y<=boxStep/2;y+=boxStep)
        for (double z=-boxStep/2; z<=boxStep/2;z+=boxStep)
        {
        Coords tmp=bx->fixpos+vect_i*x+vect_j*y+vect_k*z;
        atmX[iix]=&molen7.addAtom("Cl",tmp,"Cl","Cl" ,&X_DUM);
        iix++;
        }
        molen7.addBond(atmX[0],atmX[1],BondType::SINGLE);
        molen7.addBond(atmX[0],atmX[2],BondType::SINGLE);
        molen7.addBond(atmX[0],atmX[4],BondType::SINGLE);
        molen7.addBond(atmX[1],atmX[3],BondType::SINGLE);
        molen7.addBond(atmX[1],atmX[5],BondType::SINGLE);

        molen7.addBond(atmX[2],atmX[3],BondType::SINGLE);
        molen7.addBond(atmX[2],atmX[6],BondType::SINGLE);
        molen7.addBond(atmX[3],atmX[7],BondType::SINGLE);
        molen7.addBond(atmX[4],atmX[5],BondType::SINGLE);
        molen7.addBond(atmX[4],atmX[6],BondType::SINGLE);
        molen7.addBond(atmX[5],atmX[7],BondType::SINGLE);
        molen7.addBond(atmX[6],atmX[7],BondType::SINGLE);*/
        }
    }
    /*
    for (ItCBox it=AllBoxes.begin(); it != AllBoxes.end(); it++)
        {
            Box *bx = *it;
            //if (cube_pos >= Max_Cube_Num)continue;
            cube_pos = bx->getId();

            if (bx->isMargin())               {   molen6.addAtom("H",bx->rotpos,"H","Du", &R_DUM);  molen7.addAtom("H",bx->rotpos,"H","Du" ,&X_DUM);}
            if (bx->Props.isInProtein())      {    molen.addAtom("O",bx->rotpos,"O","O.2",&R_PROT); molen7.addAtom("O",bx->rotpos,"O","O.2",&X_PROT);}
            if (bx->Props.isCloseProtein())   {   molen3.addAtom("C",bx->rotpos,"C","C.3",&R_OUT);  molen7.addAtom("C",bx->rotpos,"C","C.3",&X_OUT);}
            if (bx->Props.isOutCavity())      {   molen2.addAtom("N",bx->rotpos,"N","N.3",&R_SIT);  molen7.addAtom("N",bx->rotpos,"N","N.3",&X_SIT);}
            if (bx->Props.isInCavity())       {   molen8.addAtom("F",bx->rotpos,"F","F",  &R_IN);   molen7.addAtom("F",bx->rotpos,"F","F"  ,&X_IN);}
            if (bx->Props.isNoInfo())         {   molen4.addAtom("P",bx->rotpos,"P","P.2",&R_VID);  molen7.addAtom("P",bx->rotpos,"P","P.2",&X_VID);}

            Atom *atmX[8];unsigned int iix=0;
            for (double x=-0.5; x<=0.5;x+=1)
            for (double y=-0.5; y<=0.5;y+=1)
            for (double z=-0.5; z<=0.5;z+=1)
            {
            Coords tmp=bx->rotpos+Coords(x,y,z);
            atmX[iix]=&molen7.addAtom("Cl",tmp,"Cl","Cl" ,&X_DUM);
            iix++;
            }
            molen7.addBond(atmX[0],atmX[1],BondType::SINGLE);
            molen7.addBond(atmX[0],atmX[2],BondType::SINGLE);
            molen7.addBond(atmX[0],atmX[4],BondType::SINGLE);
            molen7.addBond(atmX[1],atmX[3],BondType::SINGLE);
            molen7.addBond(atmX[1],atmX[5],BondType::SINGLE);

            molen7.addBond(atmX[2],atmX[3],BondType::SINGLE);
            molen7.addBond(atmX[2],atmX[6],BondType::SINGLE);
            molen7.addBond(atmX[3],atmX[7],BondType::SINGLE);
            molen7.addBond(atmX[4],atmX[5],BondType::SINGLE);
            molen7.addBond(atmX[4],atmX[6],BondType::SINGLE);
            molen7.addBond(atmX[5],atmX[7],BondType::SINGLE);
            molen7.addBond(atmX[6],atmX[7],BondType::SINGLE);
        }
*/
    // s = f+"_FULL.mol2";     MoleWriter writer(s,FileFormat::MOL2); writer.writeMOL2(&molen7);
    if (molen.numAtom()) {  s = f+"_PROTEIN.mol2"; MoleWriter writer(s,FileFormat::MOL2); writer.writeMOL2(&molen); }
    if (molen2.numAtom()){  s = f+"_OUTCAV.mol2";     MoleWriter writer2(s,FileFormat::MOL2); writer2.writeMOL2(&molen2); }
    if (molen3.numAtom()){  s = f+"_CLOSEPROTEIN.mol2";     MoleWriter writer3(s,FileFormat::MOL2); writer3.writeMOL2(&molen3); }
    if (molen4.numAtom()){  s = f+"_NOINFO.mol2";     MoleWriter writer4(s,FileFormat::MOL2); writer4.writeMOL2(&molen4); }
    if (molen5.numAtom()){  s = f+"_LIGAND.mol2";  MoleWriter writer5(s,FileFormat::MOL2); writer5.writeMOL2(&molen5); }
    if (molen6.numAtom()){  s = f+"_MARGIN.mol2";  MoleWriter writer6(s,FileFormat::MOL2); writer6.writeMOL2(&molen6); }
    if (molen8.numAtom()){  s = f+"_INCAVITY.mol2";  MoleWriter writer7(s,FileFormat::MOL2); writer7.writeMOL2(&molen8); }


}


/**
 * @brief Grid::getAdjacency : return a list of box adjacent to the given box
 * @param liste : list that will be filled with adjacent boxes
 * @param box : given box to seek adjacent box from
 * @param min_level : minimal number of box separating the given box from the adjacent box
 * @param max_level : maximal number of box separating the given box from the adjacent box
 *
 * Return a list of box adjacent to the given box from a number or layer between
 * the min_level and the max_level. A layer is the list of box. The layer 0 is
 * the box itself. The layer 1 is boxes that have at least one edge or vertex in
 * common with the given box. The layer 2 does not touch the box itself, but have
 * edge and vertex in common with the layer 1.
 * Layer 0 is not saved in the list.
 *
 */
void Grid::getAdjacency(multimap<unsigned int,Box*>& liste, const Box* const box,
                        const  int& min_level,
                        const  int& max_level)
{
    if (box == (Box*)NULL) return;
    liste.clear();

    int cube_pos;
    // Explicit conversion from unsinged int to int:
    const int NPPRi = static_cast<int>(numPtRangeI);
    const int NPPRj = static_cast<int>(numPtRangeJ);
    const int NPPRk = static_cast<int>(numPtRangeK);
    const int MCN = static_cast<int>(maxCubeNum);

    Box * adjBx;
    int max;
    // Getting box coordinates over i,j,k axis :
    const int i=(int)box->rotpos.x;
    const int j=(int)box->rotpos.y;
    const int k=(int)box->rotpos.z;

    // Scanning over the range [-max_level;max_level] for each axis.
    // Avoid outside boundaries cube :
    for (int k1 = k-max_level; k1 <= k+max_level; k1++){if (k1 <0 || k1 > NPPRk) continue;
        for (int j1 = j-max_level; j1 <= j+max_level; j1++){if (j1 <0 || j1 > NPPRj) continue;
            for (int i1 = i-max_level; i1 <= i+max_level; i1++){if (i1 <0 || i1 > NPPRi) continue;
                cube_pos = i1+j1*NPPRi+k1*NPPRi*NPPRj;

                //Avoid outside boundaries cube :
                if (cube_pos < 0 || cube_pos >= MCN) continue;

                // Getting the maximal distance betwwen the 3 axis :
                max=abs(k1-k);
                if (abs(j1-j)>max)max=abs(j1-j);
                if (abs(i1-i)>max)max=abs(i1-i);

                if (max > max_level || max< min_level || max==0)continue;
                adjBx=AllBoxes.at(cube_pos);
                liste.insert(pair<unsigned int,Box*>(max,adjBx));


            }}}

}

/**
 * @brief Grid::getAdjacency  : return a list of box adjacent to the given box
 * @param liste : list that will be filled with adjacent boxes
 * @param box : given box to seek adjacent box from
 * @param dist : distance from the box center to consider an adjacent box
 * Return a list of box adjacent to the given box from a number or layer between
 * the min_level and the max_level.
 * min_level is defined by floor(distance/boxStep)
 * max_level is defined by the distance/boxStep+1.
 * A layer is the list of box. The layer 0 is
 * the box itself. The layer 1 is boxes that have at least one edge or vertex in
 * common with the given box. The layer 2 does not touch the box itself, but have
 * edge and vertex in common with the layer 1.
 * Layer 0 is not saved in the list.

 */
void Grid::getAdjacency(vector<Box*>& liste,
                        const Box* const box,
                        const  double &dist,
                        const bool& spherical,
                        const bool& borderonly)
{
    if (box == (Box*)NULL) return;
    liste.clear();

    int cube_pos;
    const int NPPRi = static_cast<int>(numPtRangeI);
    const int NPPRj = static_cast<int>(numPtRangeJ);
    const int NPPRk = static_cast<int>(numPtRangeK);
    const int MCN = static_cast<int>(maxCubeNum);
    const int min_level= (int)floor(dist/boxStep);
    const int max_level= ((int)floor(dist/boxStep))+1;
    const int steps=min_level*min_level+1;
    Box * adjBx;
    int max,min;

    const int    i=(int)box->rotpos.x;
    const int    j=(int)box->rotpos.y;
    const int    k=(int)box->rotpos.z;
    for (int k1 = k-max_level; k1 <= k+max_level; k1++){if (k1 <0 || k1 > NPPRk) continue;
        for (int j1 = j-max_level; j1 <= j+max_level; j1++){if (j1 <0 || j1 > NPPRj) continue;
            for (int i1 = i-max_level; i1 <= i+max_level; i1++){if (i1 <0 || i1 > NPPRi) continue;
                cube_pos = i1+j1*NPPRi+k1*NPPRi*NPPRj;

                if (cube_pos < 0 || cube_pos >= MCN) continue;
                //if (borderonly && max != max_level)continue;// BORDER ONLY
                if (spherical)
                {
                    const int bo=(k1-k)*(k1-k)+(j1-j)*(j1-j)+(i1-i)*(i1-i);
                    if ((borderonly && bo < steps-(min_level-1)*(min_level-1)) || bo >steps)continue;
                }
                else
                {
                    max=abs(k1-k);min=abs(k1-k);
                    if (abs(j1-j)>max)max=abs(j1-j);if (abs(j1-j)<min)min=abs(j1-j);
                    if (abs(i1-i)>max)max=abs(i1-i);if (abs(j1-j)<min)min=abs(j1-j);
                    if (max > max_level|| max==0)continue;
                }



                adjBx=AllBoxes.at(cube_pos);
                if (!adjBx->getUse())continue;
                liste.push_back(adjBx);


            }}}

}


/**
 * @brief Grid::getFullAdjacency
 * @param liste : list that will be filled with adjacent boxes
 * @param box : given box to seek adjacent box from
 * @param dist : distance from the box center to consider an adjacent box
 * Return a list of box adjacent to the given box from a number or layer between
 * 0 and the max_level. max_level is defined by the distance/boxStep+1.
 * A layer is the list of box. The layer 0 is
 * the box itself. The layer 1 is boxes that have at least one edge or vertex in
 * common with the given box. The layer 2 does not touch the box itself, but have
 * edge and vertex in common with the layer 1.
 * Layer 0 is not saved in the list.

 */
void Grid::getFullAdjacency(vector<Box*>& liste, const Box* const box,
                            const  double &dist) const
{
    if (box == (Box*)NULL) return;
    liste.clear();

    int cube_pos;
    const int NPPRi = static_cast<int>(numPtRangeI);
    const int NPPRj = static_cast<int>(numPtRangeJ);
    const int NPPRk = static_cast<int>(numPtRangeK);
    const int MCN = static_cast<int>(maxCubeNum);
    const int min_level=1;
    const int max_level= ((int)floor(dist/boxStep))+1;

    Box * adjBx;
    int max;
    const int    i=(int)box->rotpos.x;
    const int    j=(int)box->rotpos.y;
    const int    k=(int)box->rotpos.z;
    for (int k1 = k-max_level; k1 <= k+max_level; k1++){if (k1 <0 || k1 > NPPRk) continue;
        for (int j1 = j-max_level; j1 <= j+max_level; j1++){if (j1 <0 || j1 > NPPRj) continue;
            for (int i1 = i-max_level; i1 <= i+max_level; i1++){if (i1 <0 || i1 > NPPRi) continue;
                cube_pos = i1+j1*NPPRi+k1*NPPRi*NPPRj;

                if (cube_pos < 0 || cube_pos >= MCN) continue;
                max=abs(k1-k);
                if (abs(j1-j)>max)max=abs(j1-j);
                if (abs(i1-i)>max)max=abs(i1-i);

                if (max > max_level || max< min_level || max==0)continue;
                adjBx=AllBoxes.at(cube_pos);
                if (!adjBx->getUse())continue;
                liste.push_back(adjBx);


            }}}

}

void Grid::getAdjacentAtoms(vector<Atom*>& liste,const Atom&  atom, const double &dist,const bool &wHydrogen) const
{
    if (atom.getBox(this) == (Box*)NULL)return;
    const Box*const box= atom.getBox(this);
    vector<Box*> adjacents;
    liste.clear();
    getFullAdjacency(adjacents,box,dist+boxStep);

    for (ItCBox itC = adjacents.begin(); itC != adjacents.end(); ++itC)
    {
        Box &bx = **itC;
        if (!bx.getUse())continue;
        for (ItCAtom itA = bx.firstAtom();      itA != bx.lastAtom();
             ++itA)
        {
            Atom& atomC = **itA;
            if (!wHydrogen &&atomC.isHydrogen())continue;
            if (atom.fixpos.calcDist(atomC.fixpos) >= dist)continue;
            liste.push_back(&atomC);
        }
    }
    if (!box->getUse())return;
    for (ItCAtom itA = box->firstAtom();      itA != box->lastAtom();
         ++itA)
    {
        Atom& atomC = **itA;
        if (!wHydrogen &&atomC.isHydrogen())continue;
        if (atom.fixpos.calcDist(atomC.fixpos) >= dist)continue;
        liste.push_back(&atomC);
    }
}

void Grid::getAdjacentAtoms(vector<Atom*>& liste,const Box&  box, const double &dist,const bool &wHydrogen) const
{
    vector<Box*> adjacents;
    liste.clear();
    getFullAdjacency(adjacents,&box,dist);
    for (ItCBox itC = adjacents.begin(); itC != adjacents.end(); ++itC)
    {
        Box &bx = **itC;
//        cout << bx.
        if (!bx.getUse())continue;
        for (ItCAtom itA = bx.firstAtom();      itA != bx.lastAtom();
             ++itA)
        {
            Atom& atomC = **itA;
//            if (atomC.getName() =="DuCy") { cout << "Bingo j'ai un centroide" << endl << endl<< endl;}
            if (!wHydrogen &&atomC.isHydrogen())continue;
            if (box.fixpos.calcDist(atomC.fixpos) >= dist)continue;
            liste.push_back(&atomC);
        }
    }
}

Box* Grid::getBox(const Coords& pos)const
{
    Coords rotpos = (pos-X8)/boxStep;
    const int NPPRi = numPtRangeI;
    const int NPPRj = numPtRangeJ;
    const int NPPRk = numPtRangeK;

    // Maximal distance an atom can be from the box center
    // to be considered as within it.
    // Represent a sphere that goes through all vertex of the cube :
    const double maxThres=sqrt(3)*boxStep/2;

    // Threshold for distance calculation :
    const double maxThresEs=maxThres+0.1;

    // Maximal distance where we are sure the atom is within the box:
    const double minThres=boxStep/2;
    // Id of the supposed atom cube :
    int cube_pos;
    // Iterators over i,j,k axis :
    int i,j,k;
    Box *bx;
    double dist;
    double bestdist=1000;
    Box *bestbx;

    rotpos.setCoords(rotpos.x*vect_i.x+rotpos.y*vect_i.y+rotpos.z*vect_i.z,
                     rotpos.x*vect_j.x+rotpos.y*vect_j.y+rotpos.z*vect_j.z,
                     rotpos.x*vect_k.x+rotpos.y*vect_k.y+rotpos.z*vect_k.z);

    i=rotpos.x;
    j=rotpos.y;
    k=rotpos.z;
    cube_pos = i+j*NPPRi+k*NPPRi*NPPRj;
    // Impossible box :
    if (cube_pos >= maxCubeNum || cube_pos < 0) return (Box*)NULL;
    // STEP 2 - Scanning around the position in O',i,j,k to search for the closest box
    for (int k1 = rotpos.z-1; k1 <= rotpos.z+1; k1++){if (k1 <0 || k1 > NPPRk) continue;
        for (int j1 = rotpos.y-1; j1 <= rotpos.y+1; j1++){if (j1 <0 || j1 > NPPRj) continue;
            for (int i1 = rotpos.x-1; i1 <= rotpos.x+1; i1++){if (i1 <0 || i1 > NPPRi) continue;

                cube_pos = i1+j1*NPPRi+k1*NPPRi*NPPRj;

                if (cube_pos >= maxCubeNum || cube_pos < 0) continue;
                bx   = AllBoxes.at(cube_pos);
                dist = bx->fixpos.calcDist(pos,maxThresEs);
                // Outside the sphere passing by all cube vertex, i.e. the sphere
                // that contains the cube, we are sure that the atom is outside the box
                if (dist > maxThres) continue;
                // When within the small sphere, i.e the sphere that is contained
                // within the cube, we are sure that the atom is in the box.
                if (dist < minThres)
                {
                    bestbx=bx;bestdist=0;break;
                }
                else if (dist < bestdist)
                {
                    bestdist=dist; bestbx=bx;
                }
            }
            if (bestdist==0)break;
        }//j1
        if (bestdist==0)break;
    }//k1
    if (bestbx != (Box*)NULL)return bestbx;
    return (Box*)NULL;

}

void Grid::rotateAtoms(AtomList &atomlist,const bool& find_box)const
{
    // Used as rotpos iterator for atoms:
    Coords xpos;

    // Id of the supposed atom cube :
    int cube_pos;

    // Distance between cube and atom :
    double dist;

    // Box where we supposed the atom is in :
    Box *bx=(Box*)NULL;

    // Box where the atom is the closest to
    // This box will change during iteration :
    Box *bestbx=(Box*)NULL;
    // by trying to minimize the distance:
    double bestdist;

    // Box borders :
    const int NPPRi = numPtRangeI;
    const int NPPRj = numPtRangeJ;
    const int NPPRk = numPtRangeK;

    // Maximal distance an atom can be from the box center
    // to be considered as within it.
    // Represent a sphere that goes through all vertex of the cube :
    const double maxThres=sqrt(3)*boxStep/2;

    // Threshold for distance calculation :
    const double maxThresEs=maxThres+0.1;

    // Maximal distance where we are sure the atom is within the box:
    const double minThres=boxStep/2;

    // Iterators over i,j,k axis :
    int i,j,k;



    // Scanning all atoms of the molecule :
    for (ItCAtom itA = atomlist.begin();itA != atomlist.end();
         ++itA)
    {

        Atom *atm= *itA;

        // Calculating coordinates based on O',i,j,k axis :
        atm->rotpos = (atm->fixpos-X8)/boxStep;


        atm->rotpos.setCoords(atm->rotpos.x*vect_i.x+atm->rotpos.y*vect_i.y+atm->rotpos.z*vect_i.z,
                              atm->rotpos.x*vect_j.x+atm->rotpos.y*vect_j.y+atm->rotpos.z*vect_j.z,
                              atm->rotpos.x*vect_k.x+atm->rotpos.y*vect_k.y+atm->rotpos.z*vect_k.z);
#ifdef ICHEM_DEBUG
        cout << atm->getIdentifier()<<"\t";
#endif
#ifdef ICHEM_FULL_DEBUG
        cout  << "  |--> FIX:"<< atm->fixpos.toString()<< endl
              << "  |--> ROT:"<<atm->rotpos.toString()<<endl;
#endif

        // If we just want to rotate and not find the box, we stop here
        if (!find_box || !atm->isUsed())
        {
#ifdef ICHEM_DEBUG
            if (!atm->isUsed())cout << "|ATOM NOT USED";
            else cout << "|NO FINDING BOX";
            cout <<endl;
#endif
            continue;

        }

        bestdist=1000;
        xpos = atm->rotpos;
        i=xpos.x;
        j=xpos.y;
        k=xpos.z;
        cube_pos = i+j*NPPRi+k*NPPRi*NPPRj;
        // Impossible box :
        if (cube_pos >= maxCubeNum || cube_pos < 0) continue;


        bestbx = (Box*)NULL;
        bx = AllBoxes.at(cube_pos);

#ifdef ICHEM_FULL_DEBUG
        cout << "  |--> ROTi" << " "<<i << " " <<j<< " " <<k<<endl
             << "  |--> DIST:"<<bx->rotpos.calcDist(atm->rotpos)<<" "<<bx->fixpos.calcDist(atm->fixpos)<<endl;
#endif

        // STEP 2 - Scanning around the position in O',i,j,k to search for the closest box
        for (int k1 = xpos.z-1; k1 <= xpos.z+1; k1++){if (k1 <0 || k1 > NPPRk) continue;
            for (int j1 = xpos.y-1; j1 <= xpos.y+1; j1++){if (j1 <0 || j1 > NPPRj) continue;
                for (int i1 = xpos.x-1; i1 <= xpos.x+1; i1++){if (i1 <0 || i1 > NPPRi) continue;

                    cube_pos = i1+j1*NPPRi+k1*NPPRi*NPPRj;

                    if (cube_pos >= maxCubeNum || cube_pos < 0) continue;

                    bx   = AllBoxes.at(cube_pos);
                    dist = bx->fixpos.calcDist(atm->fixpos,maxThresEs);
                    // Outside the sphere passing by all cube vertex, i.e. the sphere
                    // that contains the cube, we are sure that the atom is outside the box
                    if (dist > maxThres) continue;
                    // When within the small sphere, i.e the sphere that is contained
                    // within the cube, we are sure that the atom is in the box.
                    if (dist < minThres)
                    {
                        bestbx=bx;bestdist=0;break;
                    }
                    else if (dist < bestdist)
                    {
                        bestdist=dist; bestbx=bx;
                    }
                }
                if (bestdist==0)break;
            }//j1
            if (bestdist==0)break;
        }//k1


        if (bestbx == (Box*)NULL){
#ifdef ICHEM_DEBUG
            cout << "\tCUBE TO ATOM NOT FOUND"<<endl;
#endif
            continue;}
        if (bestbx != (Box*)NULL)
        {
#ifdef ICHEM_DEBUG
            cout << "\tCUBE FOUND :"<<bestbx->getId()<<"\t" << bestbx->fixpos.toString()<<endl;
#endif
            atm->addBox(this,bestbx);
        }

    }
}

/**
 * @brief Grid::rotateMole : generation rotpos for the given molecule
 * @param mole : molecule to rotate
 * @param find_box : find for each atom the box it is enclosed in
 * @param wDummy : also find box for dummy atom such as center of cycle
 *
 */
void Grid::rotateMole(Molecule & mole,
                      const bool& find_box,const bool& withResSelRules) const
{

    // Used as rotpos iterator for atoms:
    Coords xpos;

    // Id of the supposed atom cube :
    int cube_pos;

    // Distance between cube and atom :
    double dist;

    // Box where we supposed the atom is in :
    Box *bx=(Box*)NULL;

    // Box where the atom is the closest to
    // This box will change during iteration :
    Box *bestbx=(Box*)NULL;
    // by trying to minimize the distance:
    double bestdist;

    // Box borders :
    const int NPPRi = numPtRangeI;
    const int NPPRj = numPtRangeJ;
    const int NPPRk = numPtRangeK;

    // Maximal distance an atom can be from the box center
    // to be considered as within it.
    // Represent a sphere that goes through all vertex of the cube :
    const double maxThres=sqrt(3)*boxStep;

    // Threshold for distance calculation :
    const double maxThresEs=maxThres+0.1;

    // Maximal distance where we are sure the atom is within the box:
    const double minThres=boxStep/2;

    // Iterators over i,j,k axis :
    int i,j,k;

    const unsigned int &MType= mole.getMoleType();

    // Scanning all atoms of the molecule :
    for (ItCAtom itA = mole.firstAtom();itA != mole.lastAtom();
         ++itA)
    {

        Atom *atm= *itA;
        if (withResSelRules &&
                Residu::Rules[MType][atm->getResidu()->getResType()] != MoleType::PROTEIN
                && Residu::Rules[MType][atm->getResidu()->getResType()] != MoleType::LIGAND
                ) continue;
        // Calculating coordinates based on O',i,j,k axis :
        atm->rotpos = (atm->fixpos-X8)/boxStep;


        atm->rotpos.setCoords(atm->rotpos.x*vect_i.x+atm->rotpos.y*vect_i.y+atm->rotpos.z*vect_i.z,
                              atm->rotpos.x*vect_j.x+atm->rotpos.y*vect_j.y+atm->rotpos.z*vect_j.z,
                              atm->rotpos.x*vect_k.x+atm->rotpos.y*vect_k.y+atm->rotpos.z*vect_k.z);
#ifdef ICHEM_DEBUG
        cout << atm->getIdentifier()<<"\t";
#endif
#ifdef ICHEM_FULL_DEBUG
        cout  << "  |--> FIX:"<< atm->fixpos.toString()<< endl
              << "  |--> ROT:"<<atm->rotpos.toString()<<endl;
#endif

        // If we just want to rotate and not find the box, we stop here
        if (!find_box || !atm->isUsed())
        {
#ifdef ICHEM_DEBUG
            if (!atm->isUsed())cout << "|ATOM NOT USED";
            else cout << "|NO FINDING BOX";
            cout <<endl;
#endif
            continue;

        }
        if (atm->getBox(this) != (Box*)NULL)
        {
            const Box& tmpbox=*atm->getBox(this);
            if (tmpbox.fixpos.calcDist(atm->fixpos)<minThres) continue;
        }

        bestdist=1000;
        xpos = atm->rotpos;
        i=xpos.x;
        j=xpos.y;
        k=xpos.z;
        cube_pos = i+j*NPPRi+k*NPPRi*NPPRj;



        bestbx = (Box*)NULL;
        bx = (Box*)NULL;//

#ifdef ICHEM_FULL_DEBUG
        cout << "  |--> ROTi" << " "<<i << " " <<j<< " " <<k<<endl
             << "  |--> DIST:"<<bx->rotpos.calcDist(atm->rotpos)<<" "<<bx->fixpos.calcDist(atm->fixpos)<<endl;
#endif

        // STEP 2 - Scanning around the position in O',i,j,k to search for the closest box
        for (int k1 = xpos.z-1; k1 <= xpos.z+1; k1++){if (k1 <0 || k1 > NPPRk) continue;
            for (int j1 = xpos.y-1; j1 <= xpos.y+1; j1++){if (j1 <0 || j1 > NPPRj) continue;
                for (int i1 = xpos.x-1; i1 <= xpos.x+1; i1++){if (i1 <0 || i1 > NPPRi) continue;

                    cube_pos = i1+j1*NPPRi+k1*NPPRi*NPPRj;

                    if (cube_pos >= maxCubeNum || cube_pos < 0) continue;

                    bx   = AllBoxes.at(cube_pos);
                    dist = bx->fixpos.calcDist(atm->fixpos,maxThresEs);
                    // Outside the sphere passing by all cube vertex, i.e. the sphere
                    // that contains the cube, we are sure that the atom is outside the box
                    if (dist > maxThres) continue;
                    // When within the small sphere, i.e the sphere that is contained
                    // within the cube, we are sure that the atom is in the box.
                    if (dist < minThres)
                    {
                        bestbx=bx;bestdist=0;break;
                    }
                    else if (dist < bestdist)
                    {
                        bestdist=dist; bestbx=bx;
                    }
                }
                if (bestdist==0)break;
            }//j1
            if (bestdist==0)break;
        }//k1


        if (bestbx == (Box*)NULL){
#ifdef ICHEM_DEBUG
            cout << "\tCUBE TO ATOM NOT FOUND"<<endl;
#endif
            continue;}
        if (bestbx != (Box*)NULL)
        {
#ifdef ICHEM_DEBUG
            cout << "\tCUBE FOUND :"<<bestbx->getId()<<"\t" << bestbx->fixpos.toString()<<endl;
#endif
            atm->addBox(this,bestbx);
        }

    }

}




void Grid::colorCube(Molecule & mole,
                     const bool& withRotation,
                     const bool& withResSelRules,
                     const bool& withBorder)throw(MoleExcept)
{
    if (AllBoxes.size() ==0)
    {
        throw MoleExcept(1160701,
                         "Grid::colorCube",
                         "No cube found");
    }

#ifdef ICHEM_DEBUG
    cout << "######################"<<endl
         << "### COLORING CUBES ###"<<endl
         << "######################"<<endl;
#endif
    // Generates rotpos coordinates according to this grid :
    if (withRotation)rotateMole(mole);


    Box *bxclose=(Box*)NULL,
            *bxAtm=(Box*)NULL;
    BoxList boxlist;
    double dist;

    const unsigned int MType = mole.getMoleType();
    const unsigned int color = Molecule::Rules[MType];


#ifdef ICHEM_DEBUG
    cout << " MTYPE:"<< mole.getLongMoletype()<<endl;
    cout << "COLORS:"<< color<<endl;
#endif


    // Scanning all atoms of the molecule :
    for (ItCAtom itA=mole.firstAtom();itA != mole.lastAtom();++itA)
    {
        Atom &atom=**itA;
        // HYDROGEN FILTERING
        if (!wHydrogen && atom.isHydrogen()) continue;
        if (!atom.isUsed())continue;
#ifdef ICHEM_DEBUG
        cout << "####################"<<endl
             <<atom.getIdentifier()<<endl
            <<atom.rotpos.toString()<<endl
           <<atom.getBox(this)<<endl;
#endif

#ifdef ICHEM_DEBUG
        cout << "  ||--> RULES : "<<Residu::Rules[MType][atom.getResidu()->getResType()]<<endl;
#endif
        const unsigned int& ruler=Residu::Rules[MType][atom.getResidu()->getResType()];
        // RULE TYPE FILTERING :
        if (withResSelRules &&
                ruler != MoleType::PROTEIN
                &&ruler != MoleType::LIGAND
                ) continue;

        // FETCHING BOX:
        bxAtm=atom.getBox(this);

        if (bxAtm== (Box*)NULL || !bxAtm->getUse())continue;

        // CALCULATING DISTANCES :
        const double &vdRadius = atom.getVdWRadius();
        dist = atom.fixpos.calcDist(bxAtm->fixpos);


#ifdef ICHEM_DEBUG
        cout << " ||--> CURR BOX : "<< bxAtm->getId()<< " " << dist<<endl;
        cout << " ||--> ATM VdW radius : "<< vdRadius<<endl;
#endif

        // CHECKING THAT ATOM IS WITHING THE BOX :
        if (atom.fixpos.x >= bxAtm->fixpos.x-boxStep/2
                &&atom.fixpos.x  < bxAtm->fixpos.x+boxStep/2
                &&atom.fixpos.y >= bxAtm->fixpos.y-boxStep/2
                &&atom.fixpos.y  < bxAtm->fixpos.y+boxStep/2
                &&atom.fixpos.z >= bxAtm->fixpos.z-boxStep/2
                &&atom.fixpos.z  < bxAtm->fixpos.z+boxStep/2){
            bxAtm->Props.clear();
#ifdef ICHEM_DEBUG
            cout << " ||--> CCL : IN BOX"<<endl;
#endif

            if (color == MoleType::PROTEIN)
            {
                bxAtm->Props.setNoInfo(false);
                bxAtm->Props.setInProtein(true);
            }
            else if (color == MoleType::LIGAND)
            {
                bxAtm->Props.setNoInfo(false);
                bxAtm->Props.setInLigand(true);
            }
        }
        // OTHERWISE CHECKING THAT ATOM SPHERE INTERSECT THE BOX :
        else if (dist <= vdRadius )
        {
#ifdef ICHEM_DEBUG
            cout << " ||--> CCL : IN VdW Radius"<<endl;
#endif
            bxAtm->Props.clear();
            if (color == MoleType::PROTEIN)
            {
                bxAtm->Props.setInProtein(true);
                bxAtm->Props.setNoInfo(false);
            }
            else if (color == MoleType::LIGAND) {
                bxAtm->Props.setInLigand(true);
                bxAtm->Props.setNoInfo(false);}

        }
        // ATOM SPHERE OUT BUT CLOSE ENOUGH
        /*   else if (withBorder && dist < 2.5
               &&!bxAtm->Props.isInLigand()
               && !bxAtm->Props.isInProtein())
        {
#ifdef ICHEM_DEBUG
            cout << " ||--> CCL : Close enough"<<endl;
#endif
            bxAtm->Props.clear();
            bxAtm->Props.setCloseProtein(true);

        }*/



        getFullAdjacency(boxlist,bxAtm,boxStep*2);

        for (ItCBox itB =boxlist.begin(); itB != boxlist.end(); itB++)
        {
            bxclose = *itB;
            if (!bxclose->getUse())continue;
            dist = bxclose->fixpos.calcDist(atom.fixpos);


            if (dist <= vdRadius)
            {
#ifdef ICHEM_DEBUG
                cout << " ||--> OTHER BOX : "
                     << bxclose->getId()
                     << " " << bxclose->fixpos.toString()
                     <<"\t"<< dist<<"\t"
                    << "WITHIN VdW Radius"<<endl;
#endif
                bxclose->Props.clear();
                if (color == MoleType::PROTEIN){
                    bxclose->Props.clear();
                    bxclose->Props.setInProtein(true);
                }
                else if (color == MoleType::LIGAND)
                {
                    bxclose->Props.clear();
                    bxclose->Props.setInLigand(true);
                }
            }
            /*  else if (withBorder && dist < 2.5 &&
                     !bxclose->Props.isInLigand() &&
                     !bxclose->Props.isInProtein())
            {
#ifdef ICHEM_DEBUG
                cout << " ||--> OTHER BOX : "
                     << bxclose->getId()
                     << " " << bxclose->fixpos.toString()
                     <<"\t"<< dist<<"\t"
                     << "WITHIN 2.5"<<endl;
#endif
                bxclose->Props.clear();
                bxclose->Props.setCloseProtein(true);
            }*/

        }

#ifdef ICHEM_DEBUG
        cout <<endl;
#endif

    }
#ifdef ICHEM_DEBUG
    /* if (color == MoleType::PROTEIN)printInFile("COLOR_PROT",AllBoxes);
        else if (color == MoleType::LIGAND) printInFile("COLOR_LIG",AllBoxes);
        else printInFile("COLOR_OTH",AllBoxes);*/
#endif

    status=2;
}





void Grid::calcBSA(Molecule &protein, Molecule& ligand, double&bsa, double &volume_ligand)
{
    rotateMole(ligand);
    rotateMole(protein);
    AtomList listProtAtm;
    double dist,bestdist;
    Atom * bestatm;
    bool check;

    for (ItCAtom itPA = protein.firstAtom();
         itPA!= protein.lastAtom();
         ++itPA)
    {
        Atom &atomP=**itPA;
        check=false;
        if (atomP.getResiduName().substr(0,3)=="HOH"
                ||atomP.isHydrogen()
                ||!atomP.isUsed())continue;
        for (ItCAtom itLA = ligand.firstAtom();
             itLA!= ligand.lastAtom();
             ++itLA)
        {
            Atom &atomL=**itLA;
            if (!atomL.isHydrogen()
                    && atomL.calcFixpos(atomP,4.5)<4.5) {check=true;break;}
        }
        if (!check)continue;
        listProtAtm.push_back(&atomP);

    }

    double N_SIT=0, N_SURF=0;
    bool in;BoxList closeBox;


    for (ItCBox itB = AllBoxes.begin();
         itB!= AllBoxes.end();
         ++itB)
    {
        Box &bx=**itB;
        bestatm=(Atom*)NULL;
        bestdist=1000;
        bx.Props.clear();
        for (ItCAtom itPA = listProtAtm.begin();
             itPA!= listProtAtm.end();
             ++itPA)
        {
            Atom& atomP=**itPA;
            dist = atomP.fixpos.calcDist(bx.fixpos);
            if (dist < bestdist)
            {
                bestdist=dist;
                bestatm=&atomP;
            }
        }
        if (bestdist <= bestatm->getVdWRadius()) {bx.Props.setInProtein(true);}
    }

    for (ItCBox itB = AllBoxes.begin();
         itB!= AllBoxes.end();
         ++itB)
    {
        Box& bx =**itB;
        bestatm=(Atom*)NULL;
        bestdist=1000;
        for (ItCAtom itLA = ligand.firstAtom();
             itLA!= ligand.lastAtom();
             ++itLA)
        {
            Atom& atomL=**itLA;
            if (atomL.isHydrogen() || !atomL.isUsed())continue;
            dist = atomL.fixpos.calcDist(bx.fixpos);
            if (dist >= bestdist)continue;
            bestdist=dist;
            bestatm = &atomL;

        }
        if (bestatm == (Atom*)NULL)continue;


        if (bestdist < bestatm->getVdWRadius()-0.2)
        {
            bx.Props.clear();
            bx.Props.setInLigand(true);

        }
        else if (bestdist >= bestatm->getVdWRadius()-0.35
                 && bestdist <= bestatm->getVdWRadius()+0.35
                 &&!bx.Props.isInLigand())
        {
            bx.Props.clear();
            bx.Props.setCloseProtein(true);

        }
        else if (bestdist >= bestatm->getVdWRadius()+1.06
                 && bestdist<= bestatm->getVdWRadius()+1.74
                 && !bx.Props.isInLigand() && !bx.Props.isCloseProtein())
        {
            bx.Props.clear();
            bx.Props.setInCavity(true);

        }


    }

    volume_ligand=0;

    for (ItCBox itB = AllBoxes.begin();
         itB!= AllBoxes.end();
         ++itB)
    {
        const Box &bx=**itB;
        if (bx.Props.isCloseProtein()||bx.Props.isInLigand()) volume_ligand++;
        if (!bx.Props.isInCavity())continue;
        N_SIT++;
        if (bx.Props.isInProtein()){ N_SURF++;continue;}
        in=false;
        closeBox.clear();
        getAdjacency(closeBox,&bx,0.2);

        for (ItCBox itCB = closeBox.begin(); itCB != closeBox.end(); itCB++)
        {
            Box& boxcloser = **itCB;
            if (!boxcloser.Props.isInProtein()) continue;
            in=true;break;
        }
        if (in) N_SURF++;

    }
    bsa= ((N_SURF/N_SIT)*100);

    volume_ligand*=boxStep*boxStep*boxStep;




}


void Grid::calcBSA(Complex& complex, Molecule& ligand, double&bsa, double &volume_ligand)
{
    for (ItCMole it = complex.firstMole();
         it !=complex.lastMole();
         it++)
    {
        Molecule &mole=**it;
        if (&mole==&ligand)continue;
        rotateMole(mole);
    }

    rotateMole(ligand);



    AtomList listProtAtm;
    double dist,bestdist;
    Atom * bestatm;
    bool check;

    for (ItCMole it = complex.firstMole();
         it !=complex.lastMole();
         it++)
    {
        Molecule &mole=**it;
        if (&mole==&ligand)continue;
        for (ItCAtom itPA = mole.firstAtom();
             itPA!= mole.lastAtom();
             ++itPA)
        {
            Atom &atomP=**itPA;
            check=false;
            if (atomP.getResiduName().substr(0,3)=="HOH"
                    ||atomP.isHydrogen()
                    ||!atomP.isUsed())continue;
            for (ItCAtom itLA = ligand.firstAtom();
                 itLA!= ligand.lastAtom();
                 ++itLA)
            {
                Atom &atomL=**itLA;
                if (!atomL.isHydrogen()
                        && atomL.calcFixpos(atomP,4.5)<4.5) {check=true;break;}
            }
            if (!check)continue;
            listProtAtm.push_back(&atomP);

        }
    }

    double N_SIT=0, N_SURF=0;
    bool in;BoxList closeBox;


    for (ItCBox itB = AllBoxes.begin();
         itB!= AllBoxes.end();
         ++itB)
    {
        Box &bx=**itB;
        bestatm=(Atom*)NULL;
        bestdist=1000;
        bx.Props.clear();
        for (ItCAtom itPA = listProtAtm.begin();
             itPA!= listProtAtm.end();
             ++itPA)
        {
            Atom& atomP=**itPA;
            dist = atomP.fixpos.calcDist(bx.fixpos);
            if (dist < bestdist)
            {
                bestdist=dist;
                bestatm=&atomP;
            }
        }
        if (bestatm != (Atom*)NULL&& bestdist <= bestatm->getVdWRadius()) {bx.Props.setInProtein(true);}
    }

    for (ItCBox itB = AllBoxes.begin();
         itB!= AllBoxes.end();
         ++itB)
    {
        Box& bx =**itB;
        bestatm=(Atom*)NULL;
        bestdist=1000;
        for (ItCAtom itLA = ligand.firstAtom();
             itLA!= ligand.lastAtom();
             ++itLA)
        {
            Atom& atomL=**itLA;
            if (atomL.isHydrogen() || !atomL.isUsed())continue;
            dist = atomL.fixpos.calcDist(bx.fixpos);
            if (dist >= bestdist)continue;
            bestdist=dist;
            bestatm = &atomL;

        }
        if (bestatm == (Atom*)NULL)continue;


        if (bestdist < bestatm->getVdWRadius()-0.2)
        {
            bx.Props.clear();
            bx.Props.setInLigand(true);

        }
        else if (bestdist >= bestatm->getVdWRadius()-0.35
                 && bestdist <= bestatm->getVdWRadius()+0.35
                 &&!bx.Props.isInLigand())
        {
            bx.Props.clear();
            bx.Props.setCloseProtein(true);

        }
        else if (bestdist >= bestatm->getVdWRadius()+1.06
                 && bestdist<= bestatm->getVdWRadius()+1.74
                 && !bx.Props.isInLigand() && !bx.Props.isCloseProtein())
        {
            bx.Props.clear();
            bx.Props.setInCavity(true);

        }


    }

    volume_ligand=0;

    for (ItCBox itB = AllBoxes.begin();
         itB!= AllBoxes.end();
         ++itB)
    {
        const Box &bx=**itB;
        if (bx.Props.isCloseProtein()||bx.Props.isInLigand()) volume_ligand++;
        if (!bx.Props.isInCavity())continue;
        N_SIT++;
        if (bx.Props.isInProtein()){ N_SURF++;continue;}
        in=false;
        closeBox.clear();
        getAdjacency(closeBox,&bx,0.2);

        for (ItCBox itCB = closeBox.begin(); itCB != closeBox.end(); itCB++)
        {
            Box& boxcloser = **itCB;
            if (!boxcloser.Props.isInProtein()) continue;
            in=true;break;
        }
        if (in) N_SURF++;

    }
    volume_ligand*=boxStep*boxStep*boxStep;
    bsa= ((N_SURF/N_SIT)*100);




}







