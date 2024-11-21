#include <iostream>
#include "headers/ICMole/global.h"
using namespace std;
using namespace ICMole;

std::string ICMole::get_IChem_LIB_Path(const string alternate)
{
      char * pPath;
      const string var_env = (alternate.length() ==0) ? "ICHEM_LIB": alternate;
      pPath=getenv (var_env.c_str());
      if (pPath==NULL){
          throw MoleExcept(1110201,
                           "GLOBAL::getIchemLibPath",
                           "No ICHEM_LIB DEFINED");
      }
      return static_cast<string>(pPath);
}




void ICMole::tokenStr(const std::string& str,
                      std::vector<std::string>& tokens,
                      const std::string& delimiters)
    {
        // Skip delimiters at beginning.
        string::size_type lastPos = str.find_first_not_of(delimiters, 0);
        // Find first "non-delimiter".
        string::size_type pos     = str.find_first_of(delimiters, lastPos);

        while (string::npos != pos || string::npos != lastPos)
        {
            // Found a token, add it to the vector.
            tokens.push_back(str.substr(lastPos, pos - lastPos));
            // Skip delimiters.  Note the "not_of"
            lastPos = str.find_first_not_of(delimiters, pos);
            // Find next "non-delimiter"
            pos = str.find_first_of(delimiters, lastPos);
        }
    }


void ICMole::tokenize(const std::string& str,
                      std::vector<unsigned int>& tokens,
                      const std::string& delimiters)
    {
        // Skip delimiters at beginning.
        string::size_type lastPos = str.find_first_not_of(delimiters, 0);
        // Find first "non-delimiter".
        string::size_type pos     = str.find_first_of(delimiters, lastPos);

        while (string::npos != pos || string::npos != lastPos)
        {
            // Found a token, add it to the vector.
            tokens.push_back(atoi(str.substr(lastPos, pos - lastPos).c_str()));
            // Skip delimiters.  Note the "not_of"
            lastPos = str.find_first_not_of(delimiters, pos);
            // Find next "non-delimiter"
            pos = str.find_first_of(delimiters, lastPos);
        }
        sort(tokens.begin(),tokens.end());
        vector<unsigned int>::iterator it= unique(tokens.begin(),tokens.end());
        tokens.resize(it-tokens.begin());
    }



double getplp(double plpA, double plpB, double plpC, double plpD, double plpE, double distance_plp) {
    double plp = 0;
    int plpF = 20;
    if (distance_plp < plpA && plpA != 0){
        plp = plpF*(plpA-distance_plp)/plpA;
    }
    else if (distance_plp >= plpA && distance_plp < plpB){
        plp = plpE*(distance_plp-plpA)/(plpB-plpA);
    }
    else if (distance_plp >= plpB && distance_plp < plpC){
        plp = plpE;
    }
    else if (distance_plp >= plpC && distance_plp < plpD){
        plp = plpE*(plpD-distance_plp)/(plpD-plpC);
    }
    else if (distance_plp >= plpD){
        plp = 0;
    }
    return plp;
}



