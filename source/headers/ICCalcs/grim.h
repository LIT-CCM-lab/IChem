#ifndef GRIM_H
#define GRIM_H
#include "headers/ICCalcs/interaction.h"

namespace ICMole
{

struct Results
{
    double score;
    double rmsd;
    double Tc;
    double devAl;
    unsigned int NProt, NLig,NCenter,NPol;
    const std::vector<int> &clique;
    std::vector<double> matrix;
    Coords TransRigid;
    Coords TransMobil;
    Results(const std::vector<int>& cli):score(0),rmsd(0),Tc(0),devAl(0),
        NProt(0),NLig(0),NCenter(0),NPol(0),clique(cli){}
};

class Grim
{
  InterResults& reference;
  InterResults& comparison;
  CliqueData   cliqueData;
  GraphMatch   graphMatch;
  bool         wGrScore;
  double       distclose;
  double       distsame;
  double       distfar;
  double       diststd;
  std::multimap<double, Results> scores;


  bool isPair(const Atom &Ref, const Atom &Comp) const;
public:
  Grim( InterResults& refcomplex,  InterResults& compcomplex);

  const InterResults& getRefInterResult()const {return reference;}
  const InterResults& getCompInterResult()const {return comparison;}

  void setMinCliqueSize(const size_t& size ){graphMatch.setMinSizeClique(size); }

void calcsCliques(const unsigned int& interMatch);
void setDistClose(const double& distance){distclose=distance;}
void setdistfar(const double& distance){distfar=distance;}
void setDistSame (const double& distance){distsame=distance;}
void setDiststd  (const double& distance){diststd=distance;}
const CliqueData& getResults() const { return  cliqueData;}
void filterWithScore(const double& minscore);
void filterWithRMSD(const double& minscore);
void filterWithDevAl(const double& minscore);
void filterWithNProt(const double& minscore);
void filterWithNLig(const double& minscore);
void filterWithNCent(const double& minscore);
void filterWithNPol(const double& minscore);
std::multimap<double, Results>::const_reverse_iterator firstScore() const {return scores.rbegin();}
std::multimap<double, Results>::const_reverse_iterator lastScore() const {return scores.rend();}
};

}


#endif // GRIM_H
