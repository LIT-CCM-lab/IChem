//
// Created by cjacquemard on 7/13/22.
//

#include "headers/ICTools/switch.h"
#include "headers/ICMole/atom.h"
#include "headers/ICCalcs/interaction.h"
#include "headers/ICPars/molereader.h"
#include "headers/ICMole/box.h"

using namespace std;
using namespace ICMole;

void IChemSwitch::helpTest()      const
{
    cout
            << "TEST:" <<endl
            << "\tRun tests to check program integrity" <<endl
            << std::string(75, '#') <<endl
            << endl;

}

void IChemSwitch::test()      const throw(ICMole::MoleExcept)
{
    cout << "Running IChem test procedures..." << endl;
    cout << "Atom" << endl;
}