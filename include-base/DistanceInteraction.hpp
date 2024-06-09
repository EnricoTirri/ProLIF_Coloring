#ifndef PROLIF_COLORING_DISTANCE_INTERACTION
#define PROLIF_COLORING_DISTANCE_INTERACTION

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

#endif //PROLIF_COLORING_DISTANCE_INTERACTION
