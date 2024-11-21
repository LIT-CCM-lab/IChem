#include <iostream>
#include "headers/ICTools/switch.h"
#include "headers/ICMole/box.h"
#include "headers/ICCalcs/volsite.h"
#include "headers/ICCalcs/interaction.h"
#include "headers/ICPars/molewriter.h"
#include "headers/ICCalcs/pdbconvert.h"
#include "headers/ICMole/complex.h"
#include "headers/ICMole/graphmatch.h"
#include "headers/ICMole/graph.h"
#include "headers/ICMole/similarity.h"
#include "headers/ICCalcs/fragments.h"
using namespace std;
using namespace ICMole;




/*
    clock_t t;
      t = clock();
  t = clock() - t;
  cout << ((float)t)/CLOCKS_PER_SEC<<endl;
*/
int main(int argc, char *argv[])
{
try{
    if (argc == 3)
      {
        const std::string i1 = argv[1];
        if (i1 == "-F")
          {
            ifstream ifs(argv[2]);
            string params;
            vector<string> tokens;
            while (!ifs.eof())
              {
                std::getline(ifs,params);
                if (params.length()==0)continue;
                tokens.clear();
                tokenStr(params,tokens," ");
                try
                {
                  IChemSwitch(tokens.size(),tokens);
                }
                catch(MoleExcept &e)
                {
                  cerr << e.getCode()<<"\t"<<e.getSource()<<"\t"<< e.getData()<<endl;

                }

              }
          }
        else
        {

            IChemSwitch(argc,argv);
        }
      }
    else
    {

    IChemSwitch(argc,argv);
    }
    }
    catch(MoleExcept &e)
    {
    if (e.getCode() == 9010101)
    {

        return 1;}
        cerr << e.getCode()<<"\t"<<e.getSource()<<"\t"<< e.getData()<<endl;
        cerr << "TRACE : "<< endl << e.getTrace()<<endl;

    }
    catch(std::exception &e)
    {
        cerr << e.what()<<endl;
    }
}

void runTest()
{
    cout << "###################"<<endl
         << "###################"<<endl
         << "TESTING COORDINATES"<<endl
         << "###################"<<endl
         << "###################"<<endl;



    unsigned int error=0;

// STEP 1
  // PRINT
    cout << "Testing initialization ...   \t  \t";
  // CODE
    Coords coord1;
    Coords coord2(1,2,3);
    Coords coord3(coord2);

  // TEST
    if (coord1.x != 0){cerr <<endl<< "X coordinate not good for coord1"<<endl;error++;}
    if (coord1.y != 0){cerr <<endl<< "Y coordinate not good for coord1"<<endl;error++;}
    if (coord1.z != 0){cerr <<endl<< "Z coordinate not good for coord1"<<endl;error++;}
    if (coord2.x != 1){cerr <<endl<< "X coordinate not good for coord2"<<endl;error++;}
    if (coord2.y != 2){cerr <<endl<< "Y coordinate not good for coord2"<<endl;error++;}
    if (coord2.z != 3){cerr <<endl<< "Z coordinate not good for coord2"<<endl;error++;}
    if (coord3.x != 1){cerr <<endl<< "X coordinate not good for coord3"<<endl;error++;}
    if (coord3.y != 2){cerr <<endl<< "Y coordinate not good for coord3"<<endl;error++;}
    if (coord3.z != 3){cerr <<endl<< "Z coordinate not good for coord3"<<endl;error++;}
    if (error == 0) {cout << "PASSED"<<endl;}
    error=0;

// STEP 2
  // PRINT
    cout << "Testing assignation ....\t  \t";
  // CODE
    Coords coord4(10,11,12);
    Coords coord5;                      // Values set  to 0
    coord5.setCoords(coord4);           // Should be 10, 11, 12
    coord1=coord2;                      // operator=, should be 1,2,3
    coord3.x=7; coord3.y=8; coord3.z=9; // assigning class parameters
    coord2.setCoords(3,4,5);            // setCoords function
    coord4.clear();                     // Values set to 0
// TEST
    if (coord1.x != 1){cerr <<endl<< "X coordinate not good for coord1"<<endl;error++;}
    if (coord1.y != 2){cerr <<endl<< "Y coordinate not good for coord1"<<endl;error++;}
    if (coord1.z != 3){cerr <<endl<< "Z coordinate not good for coord1"<<endl;error++;}
    if (coord2.x != 3){cerr <<endl<< "X coordinate not good for coord2"<<endl;error++;}
    if (coord2.y != 4){cerr <<endl<< "Y coordinate not good for coord2"<<endl;error++;}
    if (coord2.z != 5){cerr <<endl<< "Z coordinate not good for coord2"<<endl;error++;}
    if (coord3.x != 7){cerr <<endl<< "X coordinate not good for coord3"<<endl;error++;}
    if (coord3.y != 8){cerr <<endl<< "Y coordinate not good for coord3"<<endl;error++;}
    if (coord3.z != 9){cerr <<endl<< "Z coordinate not good for coord3"<<endl;error++;}
    if (coord4.x != 0){cerr <<endl<< "X coordinate not good for coord4"<<endl;error++;}
    if (coord4.y != 0){cerr <<endl<< "Y coordinate not good for coord4"<<endl;error++;}
    if (coord4.z != 0){cerr <<endl<< "Z coordinate not good for coord4"<<endl;error++;}
    if (coord5.x !=10){cerr <<endl<< "X coordinate not good for coord5"<<endl;error++;}
    if (coord5.y !=11){cerr <<endl<< "Y coordinate not good for coord5"<<endl;error++;}
    if (coord5.z !=12){cerr <<endl<< "Z coordinate not good for coord5"<<endl;error++;}
    if (error == 0) {cout << "PASSED"<<endl;}
    error=0;


// STEP 3
  // PRINT
    cout << "Swapping values ....   \t  \t";
  // CODE
    coord1.swap(coord5);
  // TEST
    if (coord1.x !=10){cerr <<endl<< "X coordinate not good for coord1"<<endl;error++;}
    if (coord1.y !=11){cerr <<endl<< "Y coordinate not good for coord1"<<endl;error++;}
    if (coord1.z !=12){cerr <<endl<< "Z coordinate not good for coord1"<<endl;error++;}
    if (coord5.x != 1){cerr <<endl<< "X coordinate not good for coord5"<<endl;error++;}
    if (coord5.y != 2){cerr <<endl<< "Y coordinate not good for coord5"<<endl;error++;}
    if (coord5.z != 3){cerr <<endl<< "Z coordinate not good for coord5"<<endl;error++;}
    if (error == 0) {cout << "PASSED"<<endl;}
    error=0;

// STEP 4
  // PRINT
    cout << "Translate coordinates ....   \t  \t";
  // CODE
    Coords coord6(1,2,3);
    coord1.translate(3,2,1);        // Translate x axis to 3, y axis to 2, z axis to 1
    coord2.translateX(2);           // Translate X by 2
    coord3.translateY(2);           // Translate X by 2
    coord4.translateZ(2);           // Translate X by 2
    coord5+=coord1;                 // Translate coord5 by adding coord1 (warning: coord1 has been translated just before)
    coord6-=coord2;                 // Translate coord6 by removing coord2 (warning: coord2 has been translated just before)
  // TEST
    if (coord1.x !=13){cerr <<endl<< "X coordinate not good for coord1"<<endl;error++;}
    if (coord1.y !=13){cerr <<endl<< "Y coordinate not good for coord1"<<endl;error++;}
    if (coord1.z !=13){cerr <<endl<< "Z coordinate not good for coord1"<<endl;error++;}
    if (coord2.x != 5){cerr <<endl<< "X coordinate not good for coord2"<<endl;error++;}
    if (coord2.y != 4){cerr <<endl<< "Y coordinate not good for coord2"<<endl;error++;}
    if (coord2.z != 5){cerr <<endl<< "Z coordinate not good for coord2"<<endl;error++;}
    if (coord3.x != 7){cerr <<endl<< "X coordinate not good for coord3"<<endl;error++;}
    if (coord3.y !=10){cerr <<endl<< "Y coordinate not good for coord3"<<endl;error++;}
    if (coord3.z != 9){cerr <<endl<< "Z coordinate not good for coord3"<<endl;error++;}
    if (coord4.x != 0){cerr <<endl<< "X coordinate not good for coord4"<<endl;error++;}
    if (coord4.y != 0){cerr <<endl<< "Y coordinate not good for coord4"<<endl;error++;}
    if (coord4.z != 2){cerr <<endl<< "Z coordinate not good for coord4"<<endl;error++;}
    if (coord5.x !=14){cerr <<endl<< "X coordinate not good for coord5"<<endl;error++;}
    if (coord5.y !=15){cerr <<endl<< "Y coordinate not good for coord5"<<endl;error++;}
    if (coord5.z !=16){cerr <<endl<< "Z coordinate not good for coord5"<<endl;error++;}
    if (coord6.x !=-4){cerr <<endl<< "X coordinate not good for coord6"<<endl;error++;}
    if (coord6.y !=-2){cerr <<endl<< "Y coordinate not good for coord6"<<endl;error++;}
    if (coord6.z !=-2){cerr <<endl<< "Z coordinate not good for coord6"<<endl;error++;}
    if (error == 0) {cout << "PASSED"<<endl;}
    error=0;


// STEP 5:
  // PRINT
   cout << "Arithmetic operations ....   \t  \t";

  // CODE
   coord4.setCoords(1,2,3); // We reset coord4 to 1,2,3 because making multiplication and divison on null value is pointless
   coord1*=coord2;      // X1=13*5 ; Y1=13*4; Z1=13*5
   coord3/=coord5;      // X3=7/14 ; Y3=10/15; Z3=9/16
   coord4*=2;           // X4=1*2  ; Y4=2*2 ;  Z4=3*2
   coord6/=3;           // X6=-4/3 ; Y5=-2/3 ; Z5=-2/3

 // TEST
   // Note : For division, since there are floating numbers, we have to compare the difference in the substraction of the results
   if (coord1.x !=13*5){cerr <<endl<< "X coordinate not good for coord1"<<endl;error++;}
   if (coord1.y !=13*4){cerr <<endl<< "Y coordinate not good for coord1"<<endl;error++;}
   if (coord1.z !=13*5){cerr <<endl<< "Z coordinate not good for coord1"<<endl;error++;}
   if (fabs(coord3.x - 7/14)< 0.01){cerr <<endl<< "X coordinate not good for coord3"<<endl;error++;}
   if (fabs(coord3.y -10/15)< 0.01){cerr <<endl<< "Y coordinate not good for coord3"<<endl;error++;}
   if (fabs(coord3.z - 9/16)< 0.01){cerr <<endl<< "Z coordinate not good for coord3"<<endl;error++;}
   if (coord4.x != 1*2){cerr <<endl<< "X coordinate not good for coord4"<<endl;error++;}
   if (coord4.y != 2*2){cerr <<endl<< "Y coordinate not good for coord4"<<endl;error++;}
   if (coord4.z != 3*2){cerr <<endl<< "Z coordinate not good for coord4"<<endl;error++;}
   if (fabs(coord6.x +4/3)< 0.01){cerr <<endl<< "X coordinate not good for coord6"<<endl;error++;}
   if (fabs(coord6.y +2/3)< 0.01){cerr <<endl<< "Y coordinate not good for coord6"<<endl;error++;}
   if (fabs(coord6.z +2/3)< 0.01){cerr <<endl<< "Z coordinate not good for coord6"<<endl;error++;}
   if (error == 0) {cout << "PASSED"<<endl;}
   error=0;

   // What if we divide by 0 :
   cout << "Division by 0 ....   \t  \t"<<endl;
   coord4.clear();
   try
   {
   coord1/=coord4;
   error++;
   }catch (MoleExcept &e)
   {
       cerr << "This EXPECTED exception is thrown because of division by 0 => "<<e.getCode()<< " " << e.getSource()<< " " << e.getData()<<endl;
   }
   try
   {
   coord1/=0;
   error++;
   }catch (MoleExcept &e)
   {
       cerr << "This EXPECTED exception is thrown because of division by 0 => "<<e.getCode()<< " " << e.getSource()<< " " << e.getData()<<endl;
   }
   if (error == 0) {cout << "PASSED"<<endl;}
   error=0;



}


