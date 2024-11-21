#ifndef SVM_MAIN_H
#define SVM_MAIN_H

#include "headers/ICMole/global.h"
#include "headers/ICSVM/svm_common.h"
#include "headers/ICSVM/svm_learn.h"
//#include "headers/ICSvm/statsdrug.h"
//#include "headers/ICSvm/statsregr.h"



MODEL* genModel (double g ,int c , int t, char *file ,std::string modelfileIN, int writemod,std::string&,int mod);

MODEL* genModelRegressionloop(double g, int c, int t, std::string modelfileIN,int writemod,std::string& modelchar,DOC **docs,double *target,long totwords,long totdoc);

MODEL* genModelloop (double g ,int c , int t, std::string modelfileIN, int writemod,std::string&,DOC **docs, double *label,  long int totwords, long int totdoc);



std::string model_to_char(MODEL *model);

//statsdrug classification(std::string filename, model *models, bool verbose, const int nbword);

//statsregr classificationRegression(std::string filename, model *models,bool verbose,const int nbword);

#endif // SVM_MAIN_H
