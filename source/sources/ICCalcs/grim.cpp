#include <iostream>
#include "headers/ICCalcs/grim.h"
#include "headers/ICMole/complex.h"
#include "headers/ICMole/atom.h"
#include "headers/ICMole/rigidbody.h"
#include "headers/ICMole/vertex.h"

using namespace std;
using namespace ICMole;


Grim::Grim(InterResults &refcomplex,  InterResults &compcomplex):
    reference(refcomplex),

    comparison(compcomplex),

    graphMatch(cliqueData),
    wGrScore(true),
    distclose(0.75),
    distsame(0.5),
    distfar(1),
    diststd(1)

{
    cliqueData.graph.reserveVertex(refcomplex.Ints.numAtom()*compcomplex.Ints.numAtom());
    cliqueData.graph.reserveEdges(refcomplex.Ints.numAtom()*compcomplex.Ints.numAtom()*3);
}






/**
 * @brief Tell whether two atoms are identical in graph theory
 * @param Ref:  Reference atom to check
 * @param Comp : Comparison atom to check
 * @return  Return true if the reference and the comparison atoms have the same name and the same residu name. Othewise, return false.
 */
bool Grim::isPair(const Atom&  Ref, const Atom& Comp) const
{
    if (Ref.getResidu()->getName() != Comp.getResidu()->getName()) return false;
    if (Ref.getName() == "CA"  && Comp.getName() == "CA") { return true;}
    if (Ref.getName() == "CZ"  && Comp.getName() == "CZ") { return true;}
    if (Ref.getName() == "OD1" && Comp.getName() == "OD1"){ return true;}
    if (Ref.getName() == "NZ"  && Comp.getName() == "NZ") { return true;}
    if (Ref.getName() == "O"   && Comp.getName() == "O")  { return true;}
    if (Ref.getName() == "N"   && Comp.getName() == "N")  { return true;}
    if (Ref.getName() == "Zn"  && Comp.getName() == "Zn") { return true;}
    return false;
}




/**
  \fn void Grim::calcsCliques(const unsigned int inter_match,const unsigned int min_cli_size,const bool searchAllcliques, const bool scoretype, const double d1, const double d2, const double d3)
  \brief Detect maximal cliques over the given complexes
  \param inter_match : 0=> Merged interactions. 1=> InterLig only. 2=> InterCent only. 3=> InterProt only. Default 0.
  \param min_cli_size : Minimal size of all clique. If used for 3D alignment, use at least 3. (Default 3)
  \param searchAllcliques : Boolean defining whether the function should find all cliques or only maximal one. (Default false=> Maximal one only).
  \param scoretype : Use the scoring function to order cliques (default false).
  \param d1 : Distance threshold to tell whether two pairs with the same positionning can be linked together in the product graph (Default 0.5 Angstroem)
  \param d2 : Distance threshold to tell whether two pairs (InterLig<->InterCent or InterCent-InterProt) can be linked together in the product graph (Default 0.75 Angstroem)
  \param d3 : Distance threshold to tell whether two pairs (InterLig<->InterProt only) can be linked together in the product graph (Default 1 Angstroem)
  */
void Grim::calcsCliques(const unsigned int& interMatch)
{
    const size_t NrefAtm=reference.Ints.numAtom();
    const size_t NcompAtm=comparison.Ints.numAtom();
    const Molecule& refMole = reference.Ints;
    const Molecule& compMole = comparison.Ints;
//    if (refMole.SizeRes() > 0 || compMole.SizeRes() > 0 )
//    {
//        cout << "Grim Failed on " << refMole.getName() << " / " << compMole.getName() << "alignment"<< endl;

//    }
//    else{
        unsigned short ptype;
        double weight,sumweight=0;
        for (size_t i=0; i< NrefAtm;++i)
        {
            const Atom &atomR = refMole.getAtom(i);

            const std::string position=atomR.getResiduName().substr(2,1);
            //cout << atomR.getIdentifier()<< " " << position<<endl;
            if (position=="P"
                    && (interMatch==0
                        || interMatch==4
                        || interMatch==5
                        ||interMatch==6)) ptype=7;
            else if
                    (position=="L"
                     && (interMatch==0
                         || interMatch==2
                         || interMatch==3
                         || interMatch==6)) ptype=1;
            else if
                    (position=="C"
                     && (interMatch==0
                         || interMatch==1
                         || interMatch==3
                         || interMatch==5)) ptype=3;
            else throw MoleExcept(1111111,"","Wrong interMatch");

            for (size_t j=0; j< NcompAtm;++j)
            {
                const Atom &atomC = compMole.getAtom(j);
                if (!isPair(atomR,atomC))continue;
                weight=1;
                if (atomR.getName() == "CA")  weight=0.299;
                else if (atomR.getName() == "CZ")  weight=0.99;
                else if (atomR.getName() == "O")   weight=0.930;
                else if (atomR.getName() == "OD1") weight=0.993;
                else if (atomR.getName() == "N")   weight=0.834;
                else if (atomR.getName() == "NZ")  weight=0.966;
                else if (atomR.getName() == "ZN")  weight=0.985;
                Pair pair(&atomR,&atomC,&cliqueData.graph.addVertex(),weight,ptype);
                // And keeping the pair into the listPair
                cliqueData.listPairs.push_back(pair);
                // Adding this pair weight to the sumweigh
                sumweight+= pair.weight;


            }
        }
        const size_t ncli=cliqueData.listPairs.size();
        for (unsigned int i=0  ; i < ncli; ++i){
            Pair &Prr = cliqueData.listPairs.at(i);
            Prr.vertex->reserve(ncli);
        }
        double dist=0;
        float    tmaxdist=0;
        unsigned int nlink=0,nlinkref=0;
        vector<double> WList;double Wscore, ClScore=0;
        vector<Coords> CRef,CComp;
        std::multimap<double,Results >   Scored_cliques;


        // STEP2 : Generating edges :
        // All against all pair comparison to create edges
        for (unsigned int i=0  ; i < ncli; ++i){
            const Pair &Prr = cliqueData.listPairs.at(i);
            nlinkref=0;
            for (unsigned int j=i+1; j < ncli; ++j){
                const Pair &Prc = cliqueData.listPairs.at(j);
//                cout << "case : " << Prr.PType+Prc.PType;
//                cout << "  tmxdist : " << tmaxdist ;
                switch(Prr.PType+Prc.PType) // Finding which kind of pairs we are looking at => Giving the threshold
                {
                case 2: tmaxdist=distsame; break; // LIG LIG pairs
                case 6: tmaxdist=distsame; break; // CENTER CENTER PAIRS
                case 14:tmaxdist=distsame; break; // PROT PROT PAIRS
                case 4: tmaxdist=distclose; break;// LIG CENTER PAIRS
                case 10:tmaxdist=distclose; break;// CENTER PROT
                case 8: tmaxdist=distfar; break;  // LIG PROT PAIRS

                default :tmaxdist=diststd;
                }
//                cout << " apres : " << tmaxdist << endl;
//                tmaxdist = 1;

                nlink=0;
                if (Prr.ref->hasBondWith(*Prc.ref)) nlink++;
                if (Prr.comp->hasBondWith(*Prc.comp)) nlink++;
                if (nlink == 1)continue;
                // Calculating distance deviation between d(Pair1.Ref-Pair2.Ref) and d(Pair1.Comp-Pair2.Comp)
                dist = fabs(Prr.comp->calcFixpos(*Prc.comp)-Prr.ref->calcFixpos(*Prc.ref));
                // Distance must be less than above threshold and Pair1.Ref != Pair2.Ref and Pair1.Comp != Pair2.Comp
                // So we can add an edges between the two vertexes
                if (Prr.comp != Prc.comp &&
                        Prr.ref != Prc.ref &&
                        dist < tmaxdist){cliqueData.graph.addEdge(*Prr.vertex,*Prc.vertex);}

            }}


        graphMatch.createMatrix();

        // cout << "PAIRS : "<< cliqueData.graph.numVertex()<<endl;
        // cout << "LINKS : "<< cliqueData.graph.numEdges()<<endl;

        graphMatch.runSearchAllMaxCliques(false);

        RigidBody rotator;
        const size_t ncompAtom= comparison.Ints.numAtom();
        const size_t nCliques= cliqueData.cliques.size();

        for (size_t posCli=0; posCli < nCliques; ++posCli)
        {
            const vector<int> &clique=cliqueData.cliques.at(posCli);
            if (clique.size() < 3)continue;

            WList.clear();
            Wscore=0;
            CRef.clear();
            CComp.clear();
            Results res(clique);
            for (size_t nP=0; nP < clique.size(); ++nP) {
                // Retaining the related pair :
                const Pair& pair = cliqueData.listPairs.at(clique.at(nP));
                // Adding weight of the pair to the score :
                Wscore+= pair.weight;
                WList.push_back(pair.weight);
                // Recording coordinates for alignment :
                CRef.push_back(pair.ref->fixpos);
                CComp.push_back(pair.comp->fixpos);

                if (pair.PType == 7)res.NProt  +=1;
                else if (pair.PType == 3)res.NCenter+=1;
                else                   res.NLig   +=1;
                if (pair.ref->getName() != "CA"&& pair.ref->getName() != "CZ") res.NPol++;
            }
            rotator.setWeigthList(WList);
            rotator.loadCoordsToRigid(CRef);
            rotator.loadCoordsToMobile(CComp);
            // So we can get the rmsd of the alignment :
            res.rmsd = rotator.calcRotation();
            res.Tc = Wscore/sumweight;
            res.score=0;
            res.matrix=rotator.getRotMatrix();
            res.TransMobil=rotator.getTransMobil();
            res.TransRigid=rotator.getTransRigid();
            comparison.Ints.rotateMolecule(rotator.getRotMatrix(),rotator.getTransMobil(),rotator.getTransRigid());

            for (size_t nAtom =0; nAtom <ncompAtom;++nAtom)
            {
                const Atom& compAtm = comparison.Ints.getAtom(nAtom);
                const double dist =compAtm.fixpos.calcDist(compAtm.rotpos);
                res.devAl=dist*dist;

            }
            res.devAl=sqrt(res.devAl)/ncompAtom;


            // By default, no scoring function is used, so we order by decreasing Tc :
            if (!wGrScore) scores.insert(std::pair<double,Results>( res.Tc,res));
            else
            {
                // Otherwise we use the scoring function :
                const int diff = (comparison.Ints.numAtom()-ncompAtom);
                ClScore = 0.5006
                        +static_cast<double>(res.NLig)   *0.0151
                        +static_cast<double>(res.NCenter)*0.0039
                        +static_cast<double>(res.NProt)  *0.0143
                        +static_cast<double>(res.Tc)     *0.2098
                        -static_cast<double>(res.rmsd)   *0.0720
                        -static_cast<double>(abs(diff))  *0.0003;
                res.score=ClScore;
                //  cout << res.score <<endl;
                scores.insert(std::pair<double,Results>(ClScore,res));
            }


            //

        }
//    }



}





void Grim::filterWithScore(const double& minscore)
{
    multimap<double,Results>::iterator todel;
    for (multimap<double,Results>::iterator it=scores.begin(); it != scores.end();++it)
    {
        if ((*it).first >= minscore) continue;
        todel=it;--it;
        scores.erase(it);
    }
}
void Grim::filterWithRMSD(const double& minscore)
{
    multimap<double,Results>::iterator todel;
    for (multimap<double,Results>::iterator it=scores.begin(); it != scores.end();++it)
    {
        const Results& result=(*it).second;
        if (result.rmsd >= minscore) continue;
        todel=it;--it;
        scores.erase(it);
    }
}
void Grim::filterWithDevAl(const double& minscore)
{
    multimap<double,Results>::iterator todel;
    for (multimap<double,Results>::iterator it=scores.begin(); it != scores.end();++it)
    {
        const Results& result=(*it).second;
        if (result.devAl >= minscore) continue;
        todel=it;--it;
        scores.erase(it);
    }
}
void Grim::filterWithNProt(const double& minscore)
{
    multimap<double,Results>::iterator todel;
    for (multimap<double,Results>::iterator it=scores.begin(); it != scores.end();++it)
    {
        const Results& result=(*it).second;
        if (result.NProt >= minscore) continue;
        todel=it;--it;
        scores.erase(it);
    }
}
void Grim::filterWithNLig(const double& minscore)
{
    multimap<double,Results>::iterator todel;
    for (multimap<double,Results>::iterator it=scores.begin(); it != scores.end();++it)
    {
        const Results& result=(*it).second;
        if (result.NLig >= minscore) continue;
        todel=it;--it;
        scores.erase(it);
    }
}
void Grim::filterWithNCent(const double& minscore)
{
    multimap<double,Results>::iterator todel;
    for (multimap<double,Results>::iterator it=scores.begin(); it != scores.end();++it)
    {
        const Results& result=(*it).second;
        if (result.NCenter >= minscore) continue;
        todel=it;--it;
        scores.erase(it);
    }
}
void Grim::filterWithNPol(const double& minscore)
{
    multimap<double,Results>::iterator todel;
    for (multimap<double,Results>::iterator it=scores.begin(); it != scores.end();++it)
    {
        const Results& result=(*it).second;
        if (result.NPol >= minscore) continue;
        todel=it;--it;
        scores.erase(it);
    }
}





