#ifndef PROLIF_COLORING_INTERACTION
#define PROLIF_COLORING_INTERACTION

#include <string>
#include <Mesh.hpp>
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
    virtual bool getInteraction(const RDKit::ROMol *molecule,
                                MoleculeMesh &interactionMask,
                                MoleculeMesh &subtractionMask) = 0;
};

#endif //PROLIF_COLORING_INTERACTION
