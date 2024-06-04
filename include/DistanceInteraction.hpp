#ifndef INTERACTION_DISTANCE_H
#define INTERACTION_DISTANCE_H

#include <Interaction.hpp>
#include <GraphMol/GraphMol.h>
#include <vector>


class DistanceInteraction : public Interaction {
private:
    const double distance;

public:
    DistanceInteraction(const std::string &smart, double distance) : Interaction(smart), distance(distance) {};

    bool getInteraction(const RDKit::ROMol *molecule, MoleculeMesh &interactionMask,
                        MoleculeMesh &subtractionMask) override;
};

#endif //INTERACTION_DISTANCE_H
