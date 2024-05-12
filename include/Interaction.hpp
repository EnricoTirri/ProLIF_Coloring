#ifndef INTERACTION_INTERACTION_H
#define INTERACTION_INTERACTION_H

#include <string>
#include <MoleculeMesh.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

class Interaction {
protected:
    RDKit::ROMol *matchMol;

    std::vector<RDKit::MatchVectType> *findMatch(const RDKit::ROMol *molecule) {
        auto *res = new std::vector<RDKit::MatchVectType>();
        RDKit::SubstructMatch(*molecule, *matchMol, *res);
        return res;
    }

    explicit Interaction(const std::string &smart) {
        matchMol = RDKit::SmartsToMol(smart);
    }

public:
    virtual bool getInteraction(const RDKit::ROMol *molecule, MoleculeMesh &mask) = 0;
};

#endif //INTERACTION_INTERACTION_H
