#include "headers/ICSVM/svm_main.h"
#include <fstream>
#include <iostream>
#include <math.h>

using namespace std;

MODEL* genModel(double g, int c, int t, char *docfile, string modelfileIN,int writemod,string& modelchar,int mod){
    DOC **docs;  /* training examples */
    long totwords,totdoc,i;
    double *target;
    double *alpha_in=NULL;
    KERNEL_CACHE *kernel_cache;
    LEARN_PARM learn_parm;
    KERNEL_PARM kernel_parm;
    MODEL *model=(MODEL *)my_malloc(sizeof(MODEL));
    char type[100];

    //  char * modelfile =(char *) malloc(modelfileIN.size());
    //  strcpy (modelfile, modelfileIN.c_str());

    /* set default */
    strcpy (learn_parm.predfile, "trans_predictions");
    strcpy (learn_parm.alphafile, "");
    learn_parm.biased_hyperplane=1;
    learn_parm.sharedslack=0;
    learn_parm.remove_inconsistent=0;
    learn_parm.skip_final_opt_check=0;
    learn_parm.svm_maxqpsize=10;
    learn_parm.svm_newvarsinqp=0;
    learn_parm.svm_iter_to_shrink=-9999;
    learn_parm.maxiter=100000;
    learn_parm.kernel_cache_size=40;
    learn_parm.eps=0.1;
    learn_parm.transduction_posratio=-1.0;
    learn_parm.svm_costratio=1.0;
    learn_parm.svm_costratio_unlab=1.0;
    learn_parm.svm_unlabbound=1E-5;
    learn_parm.epsilon_crit=0.001;
    learn_parm.epsilon_a=1E-15;
    learn_parm.compute_loo=0;
    learn_parm.rho=1.0;
    learn_parm.xa_depth=0;
    kernel_parm.poly_degree=3;
    kernel_parm.coef_lin=1;
    kernel_parm.coef_const=1;
    strcpy(kernel_parm.custom,"empty");
    strcpy(type,"c");


    //variable parameters
    kernel_parm.kernel_type=t;
    learn_parm.svm_c=c;
    kernel_parm.rbf_gamma=g;
    learn_parm.type=CLASSIFICATION;
    kernel_parm.kernel_type = RBF;


    read_documents(docfile,&docs,&target,&totwords,&totdoc);

    /* Always get a new kernel cache. It is not possible to use the
       same cache for two different training runs */
    kernel_cache=kernel_cache_init(totdoc,learn_parm.kernel_cache_size);

    if (mod==1){
    svm_learn_classification(docs,target,totdoc,totwords,&learn_parm,
                             &kernel_parm,kernel_cache,model,alpha_in);
    }else{
    svm_learn_regression(docs,target,totdoc,totwords,&learn_parm,
                             &kernel_parm,&kernel_cache,model);
    }
    if(kernel_cache) {
        /* Free the memory used for the cache. */
        kernel_cache_cleanup(kernel_cache);
    }
    /* Warning: The model contains references to the original data 'docs'.
       If you want to free the original data, and only keep the model, you
       have to make a deep copy of 'model'. */
    /* deep_copy_of_model=copy_model(model); */
    if (writemod == 1 ){
        write_model(const_cast<char*>(modelfileIN.c_str()),model);
    }
    else{
        modelchar=model_to_char( model);
        //        printf("%s",modelchar);
    }



    free(alpha_in);
    for(i=0;i<totdoc;i++) free_example(docs[i],1);
    free(docs);
    free(target);
    return model;
}


MODEL* genModelRegressionloop(double g, int c, int t, string modelfileIN,int writemod,string& modelchar,DOC **docs,double *target,long totwords,long totdoc){
//    DOC **docs;  /* training examples */
//    long totwords,totdoc,i;

//    double *target;
    KERNEL_CACHE *kernel_cache;
    LEARN_PARM learn_parm;
    KERNEL_PARM kernel_parm;
    MODEL *model=(MODEL *)my_malloc(sizeof(MODEL));
    char type[100];

    //  char * modelfile =(char *) malloc(modelfileIN.size());
    //  strcpy (modelfile, modelfileIN.c_str());

    /* set default */
    strcpy (learn_parm.predfile, "trans_predictions");
    strcpy (learn_parm.alphafile, "");
    learn_parm.biased_hyperplane=1;
    learn_parm.sharedslack=0;
    learn_parm.remove_inconsistent=0;
    learn_parm.skip_final_opt_check=0;
    learn_parm.svm_maxqpsize=10;
    learn_parm.svm_newvarsinqp=0;
    learn_parm.svm_iter_to_shrink=-9999;
    learn_parm.maxiter=100000;
    learn_parm.kernel_cache_size=40;
    learn_parm.eps=0.1;
    learn_parm.transduction_posratio=-1.0;
    learn_parm.svm_costratio=1.0;
    learn_parm.svm_costratio_unlab=1.0;
    learn_parm.svm_unlabbound=1E-5;
    learn_parm.epsilon_crit=0.001;
    learn_parm.epsilon_a=1E-15;
    learn_parm.compute_loo=0;
    learn_parm.rho=1.0;
    learn_parm.xa_depth=0;
    kernel_parm.poly_degree=3;
    kernel_parm.coef_lin=1;
    kernel_parm.coef_const=1;
    strcpy(kernel_parm.custom,"empty");
    strcpy(type,"c");


    //variable parameters
    kernel_parm.kernel_type=t;
    learn_parm.svm_c=c;
    kernel_parm.rbf_gamma=g;
    learn_parm.type=CLASSIFICATION;
    kernel_parm.kernel_type = RBF;


//    read_documents(docfile,&docs,&target,&totwords,&totdoc);

    /* Always get a new kernel cache. It is not possible to use the
       same cache for two different training runs */
    kernel_cache=kernel_cache_init(totdoc,learn_parm.kernel_cache_size);

    svm_learn_regression(docs,target,totdoc,totwords,&learn_parm,
                             &kernel_parm,&kernel_cache,model);

    if(kernel_cache) {
        /* Free the memory used for the cache. */
        kernel_cache_cleanup(kernel_cache);
    }
    /* Warning: The model contains references to the original data 'docs'.
       If you want to free the original data, and only keep the model, you
       have to make a deep copy of 'model'. */
    /* deep_copy_of_model=copy_model(model); */
    if (writemod == 1 ){
        write_model(const_cast<char*>(modelfileIN.c_str()),model);
    }
    else{
        modelchar=model_to_char( model);
        //        printf("%s",modelchar);
    }

    return model;
}


MODEL* genModelloop(double g, int c, int t, string modelfileIN,int writemod,string& modelchar,DOC **docs,double *target,long totwords,long totdoc){
//    DOC **docs;  /* training examples */
//    long totwords,totdoc,i;

//    double *target;
    double *alpha_in=NULL;
    KERNEL_CACHE *kernel_cache;
    LEARN_PARM learn_parm;
    KERNEL_PARM kernel_parm;
    MODEL *model=(MODEL *)my_malloc(sizeof(MODEL));
    char type[100];

    //  char * modelfile =(char *) malloc(modelfileIN.size());
    //  strcpy (modelfile, modelfileIN.c_str());

    /* set default */
    strcpy (learn_parm.predfile, "trans_predictions");
    strcpy (learn_parm.alphafile, "");
    learn_parm.biased_hyperplane=1;
    learn_parm.sharedslack=0;
    learn_parm.remove_inconsistent=0;
    learn_parm.skip_final_opt_check=0;
    learn_parm.svm_maxqpsize=10;
    learn_parm.svm_newvarsinqp=0;
    learn_parm.svm_iter_to_shrink=-9999;
    learn_parm.maxiter=100000;
    learn_parm.kernel_cache_size=40;
    learn_parm.eps=0.1;
    learn_parm.transduction_posratio=-1.0;
    learn_parm.svm_costratio=1.0;
    learn_parm.svm_costratio_unlab=1.0;
    learn_parm.svm_unlabbound=1E-5;
    learn_parm.epsilon_crit=0.001;
    learn_parm.epsilon_a=1E-15;
    learn_parm.compute_loo=0;
    learn_parm.rho=1.0;
    learn_parm.xa_depth=0;
    kernel_parm.poly_degree=3;
    kernel_parm.coef_lin=1;
    kernel_parm.coef_const=1;
    strcpy(kernel_parm.custom,"empty");
    strcpy(type,"c");


    //variable parameters
    kernel_parm.kernel_type=t;
    learn_parm.svm_c=c;
    kernel_parm.rbf_gamma=g;
    learn_parm.type=CLASSIFICATION;
    kernel_parm.kernel_type = RBF;


//    read_documents(docfile,&docs,&target,&totwords,&totdoc);

    /* Always get a new kernel cache. It is not possible to use the
       same cache for two different training runs */
    kernel_cache=kernel_cache_init(totdoc,learn_parm.kernel_cache_size);

    svm_learn_classification(docs,target,totdoc,totwords,&learn_parm,
                             &kernel_parm,kernel_cache,model,alpha_in);

    if(kernel_cache) {
        /* Free the memory used for the cache. */
        kernel_cache_cleanup(kernel_cache);
    }
    /* Warning: The model contains references to the original data 'docs'.
       If you want to free the original data, and only keep the model, you
       have to make a deep copy of 'model'. */
    /* deep_copy_of_model=copy_model(model); */
    if (writemod == 1 ){
        write_model(const_cast<char*>(modelfileIN.c_str()),model);
    }
    else{
        modelchar=model_to_char( model);
        //        printf("%s",modelchar);
    }



    free(alpha_in);
//    for(i=0;i<totdoc;i++) free_example(docs[i],1);
//    free(docs);
//    free(target);
    return model;
}

//% manque stat regr

//statsregr classificationRegression(string filename, model *models,bool verbose,const int nbword){

//    WORD *words;DOC *doc =(DOC*)NULL;
//    statsregr stat;
//    stat.devZero=0;
//    char *comment;
//    long max_words_doc=1000000,queryid=0,slackid=0, wnum=0;
//    double dist=0,doc_label,costfactor,TP=0, TN=0,FP=0,FN=0;
//    double myints[2000]; //nomber of maximal read line
//    string line;
//    ifstream svm1((char*)filename.c_str() , ios::in );
//    for (int i=0; i < 2000; i++) myints[i]=0.;
//    unsigned int pos=0;
//    double RpUp=0, Rpxy=0, Rpyx=0,distdroit=0;
//    while(getline(svm1, line))
//    {
//        words = (WORD *)my_malloc(sizeof(WORD)*(nbword));
//        parse_document((char*)line.c_str(),words,&doc_label,&queryid,&slackid,&costfactor,&wnum,
//                       max_words_doc,&comment);

//        doc = create_example(doc_label,queryid,slackid,costfactor,create_svector(words,comment,1.0));
//        dist=classify_example(models,doc);
//        myints[pos]=dist;
//        pos++;
//        if (verbose) cout << "set: " << comment << "  ||  real: " << doc_label << "  ||  pred:"<< dist << endl;
//        free_example (doc,1);
//        free(words);

//        stat.rmse+= (doc_label-dist)*(doc_label-dist);
//        distdroit= fabs(doc_label-dist);
//        RpUp+= ((dist-distdroit)*(doc_label-distdroit));
//        Rpxy+=pow((dist-distdroit),2);
//        Rpyx+=pow((doc_label-distdroit),2);

//    }
//    stat.Rpears= RpUp/(sqrt(Rpxy)*sqrt(Rpyx));
//    stat.devZero=stat.rmse;
//    stat.rmse=sqrt(stat.rmse/pos);
//    std::vector<double> myvector (myints,myints+pos);
//    std::sort (myvector.begin(), myvector.end());
//    std::vector<double>::iterator it;
//    it = std::unique (myvector.begin(), myvector.end());
//    myvector.resize( std::distance(myvector.begin(),it) );
//    if (myvector.size()< (pos*3/4)){
//        stat.Rpears = 0;
//        stat.rmse = 1000;
//        stat.devZero = 1000;
//    }
//    svm1.close();
//    cout << "all Rpears, rmse & devzero : " << stat.Rpears <<" " <<stat.rmse  <<" " << stat.devZero << endl;
//    return stat;

//}

//statsdrug classification(string filename, model *models,bool verbose,const int nbword){

//    WORD *words;DOC *doc =(DOC*)NULL;
//    statsdrug stat;
//    stat.devZero=0;
//    char *comment;
//    long max_words_doc=1000000,queryid=0,slackid=0, wnum=0;
//    double dist=0,doc_label,costfactor;
//    string line;
//    //    MODEL * models;
//    //    models= read_model((char*)modelname.c_str());
//    ifstream svm1((char*)filename.c_str() , ios::in );
//    double TP=0, TN=0,FP=0,FN=0;
////    double precision=0,recall=0, fmesure=0;
//    double myints[2000];
//    for (int i=0; i < 2000; i++) myints[i]=0.;
//    unsigned int pos=0;
//    while(getline(svm1, line))
//    {
//        words = (WORD *)my_malloc(sizeof(WORD)*(nbword));
//        parse_document((char*)line.c_str(),words,&doc_label,&queryid,&slackid,&costfactor,&wnum,
//                       max_words_doc,&comment);

//        doc = create_example(doc_label,queryid,slackid,costfactor,create_svector(words,comment,1.0));
//        dist=classify_example(models,doc);
//        myints[pos]=dist;
//        pos++;
//        if (verbose) cout << "set: "<< comment <<"  ||  real: "<< line[0] << "  ||  pred:"<< dist << endl;
//        free_example (doc,1);
//        //        for(int i=1;i<models->sv_num;i++) {
//        //            free_example(models->supvec[i],1);
//        //        }
//        //        cout << line[0] << " ";
//        if (line[0] == '1' && dist > 0 ) TP++;
//        if (line[0] == '-' && dist < 0 ) TN++;
//        if (line[0] == '1' && dist < 0 ) {
//            FN++;
//            stat.devZero -= dist;
//        }
//        if (line[0] == '-' && dist > 0 ){
//            FP++;
//            stat.devZero += dist;
//        }
//        //        cout << dist << endl;
//        free(words);
//    }

//    stat.pres = TP / (TP + FP);
//    stat.recall = TP / (TP + FN);
//    stat.fmeasure = (2*stat.pres*stat.recall)/(stat.pres+stat.recall);
//    stat.MCC = ((TP*TN)-(FP*FN))/ sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
//    stat.accuracy = (TP+TN)/ (TP+TN+FN+FP);

//    std::vector<double> myvector (myints,myints+pos);
//    std::sort (myvector.begin(), myvector.end());
//    std::vector<double>::iterator it;
//    it = std::unique (myvector.begin(), myvector.end());
//    myvector.resize( std::distance(myvector.begin(),it) );
////    cout << " >>>>>>>>>>>>> ici la taille : " << myvector.size() << endl;
//    if (myvector.size()< (pos*3/4)){
//        stat.pres = 0.5;
//        stat.recall = 0.5;
//        stat.fmeasure = 0.5;
//        stat.MCC = 0.5;
//        stat.accuracy = 0.5;
//        stat.devZero = 100;
//    }

//    svm1.close();
//    //    free_model(models,1);

//    return stat;

//}


string model_to_char(MODEL *model)
{

    long j,i,sv_num;
    SVECTOR *v;
ostringstream oss;
    if(verbosity>=1) {
        printf("Writing model ..."); fflush(stdout);
    }
    oss << "SVM-light Version " << VERSION<< "|_|";
    oss << model->kernel_parm.kernel_type << " # kernel type|_|";

    oss << model->kernel_parm.poly_degree << " # kernel parameter -d|_|";

    oss <<   model->kernel_parm.rbf_gamma<<" # kernel parameter -g|_|";

    oss << model->kernel_parm.coef_lin<< " # kernel parameter -s|_|";
    oss << model->kernel_parm.coef_const<< " # kernel parameter -r|_|";
    oss << model->kernel_parm.custom<< " # kernel parameter -u|_|";
    oss << model->totwords<< " # highest feature index|_|";
    oss << model->totdoc<< " # number of training documents|_|";

    /*

    sprintf(modelchar,"SVM-light Version %s|_|",VERSION);
    sprintf(modelchar,"%s%ld # kernel type|_|",modelchar,
            model->kernel_parm.kernel_type);
    sprintf(modelchar,"%s%ld # kernel parameter -d|_|",modelchar,
            model->kernel_parm.poly_degree);
    sprintf(modelchar,"%s%.8g # kernel parameter -g|_|",modelchar,
            model->kernel_parm.rbf_gamma);
    sprintf(modelchar,"%s%.8g # kernel parameter -s|_|",modelchar,
            model->kernel_parm.coef_lin);
    sprintf(modelchar,"%s%.8g # kernel parameter -r|_|",modelchar,
            model->kernel_parm.coef_const);
    sprintf(modelchar,"%s%s# kernel parameter -u|_|",modelchar,model->kernel_parm.custom);
    sprintf(modelchar,"%s%ld # highest feature index|_|",modelchar,model->totwords);
    sprintf(modelchar,"%s%ld # number of training documents|_|",modelchar,model->totdoc);
*/
    sv_num=1;
    for(i=1;i<model->sv_num;i++) {
        for(v=model->supvec[i]->fvec;v;v=v->next)
            sv_num++;
    }

    oss << sv_num << " # number of support vectors plus 1 |_|";
    oss << model->b <<  " # threshold b, each following line is a SV (starting with alpha*y)|_|";

//    sprintf(modelchar,"%s%.8g # threshold b, each following line is a SV (starting with alpha*y)|_|",modelchar,model->b);

    for(i=1;i<model->sv_num;i++) {
        for(v=model->supvec[i]->fvec;v;v=v->next) {
            oss << model->alpha[i]*v->factor << " ";
//            sprintf(modelchar,"%s%.32g ",modelchar,model->alpha[i]*v->factor);
            for (j=0; (v->words[j]).wnum; j++) {
                oss << (long)(v->words[j]).wnum <<":"<< (double)(v->words[j]).weight << " ";
//                sprintf(modelchar,"%s%ld:%.8g ",modelchar,
//                        (long)(v->words[j]).wnum,
//                        (double)(v->words[j]).weight);
            }
            oss << "#"<< v->userdefined <<"|_|";
//            sprintf(modelchar,"%s#%s|_|",modelchar,v->userdefined);
            /* NOTE: this could be made more efficient by summing the
       alpha's of identical vectors before writing them to the
       file. */
        }
    }
    if(verbosity>=1) {
        printf("done\n");
    }
    //cout << oss.str()<<endl;
    return oss.str();
}






