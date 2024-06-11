#ifndef PROLIF_COLORING_INTERACTION
#define PROLIF_COLORING_INTERACTION

#include <string>
#include <Mesh.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

/**
 * This is the abstract class that defines the methods needed to calculate an interaction by an input molecule
 */
class Interaction {
protected:
    /**
     * The continuous molecule definition of the match pattern required by interaction
     */
    RDKit::ROMol *matchMol;

    /**
     * This function return all the matches between input molecule and match-pattern molecule
     * @param molecule The input molecule
     * @return All matches between input molecule and match-pattern molecule
     */
    std::vector<RDKit::MatchVectType> *findMatch(const RDKit::ROMol *molecule) {
        auto *res = new std::vector<RDKit::MatchVectType>();
        RDKit::SubstructMatch(*molecule, *matchMol, *res);
        return res;
    }

    /**
     * This constructor initialize the match-pattern molecule from the input SMART match definition
     * @param smart The input SMART definition for interaction match-pattern
     */
    explicit Interaction(const std::string &smart) {
        matchMol = RDKit::SmartsToMol(smart);
    }

public:
    /**
     * This function calculate the discrete space the interaction is acting on
     * @param molecule The reference input continuous molecule
     * @param interactionMask The output discrete space definition of interaction acting space
     * @param subtractionMask A input discrete space definition of a subtraction mask for the interaction space
     * @return False if no interaction has been found, True otherwise
     */
    virtual bool getInteraction(const RDKit::ROMol *molecule,
                                MoleculeMesh &interactionMask,
                                MoleculeMesh &subtractionMask) = 0;
};

#endif //PROLIF_COLORING_INTERACTION
