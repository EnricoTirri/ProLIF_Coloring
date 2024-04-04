//
// Created by Enrico on 29/03/2024.
//

#ifndef INTERACTION_INTERACTION_H
#define INTERACTION_INTERACTION_H

#include <string>
#include <InteractionMask.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

class Interaction{
protected:
    RDKit::ROMol* mol;

    RDKit::MatchVectType * findMatch(RDKit::ROMol* molecule){
        auto *res = new RDKit::MatchVectType();
        RDKit::SubstructMatch(*molecule,*mol,*res);
        return res;
    }

    Interaction(const std::string &smart){
        mol = RDKit::SmartsToMol(smart);
    }

public:
    virtual std::vector<InteractionMask> getInteractions(RDKit::ROMol *molecule) = 0;
};

#endif //INTERACTION_INTERACTION_H
