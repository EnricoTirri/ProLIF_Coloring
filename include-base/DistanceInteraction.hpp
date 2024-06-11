#ifndef PROLIF_COLORING_DISTANCE_INTERACTION
#define PROLIF_COLORING_DISTANCE_INTERACTION

#include <Interaction.hpp>
#include <GraphMol/GraphMol.h>

/**
 * This class extends Interaction class
 * Defines a type of interaction in which the space the interaction is acting on is defined by:
 *      - distance(point_in_space, centroid_of_interaction) <= reference_distance
 */
class DistanceInteraction : public Interaction {
private:
    /**
     * The reference distance for the interaction
     */
    const double distance;

public:
    /**
     * This constructor extends the Interaction class constructor by receiving also the reference distance value
     * @param smart The input SMART definition for interaction match-pattern
     * @param distance The reference distance for the interaction
     */
    DistanceInteraction(const std::string &smart, double distance) : Interaction(smart), distance(distance) {};

    /**
     * This function overrides the Interaction class one
     */
    bool getInteraction(const RDKit::ROMol *molecule, MoleculeMesh &interactionMask,
                        MoleculeMesh &subtractionMask) override;
};

#endif //PROLIF_COLORING_DISTANCE_INTERACTION
