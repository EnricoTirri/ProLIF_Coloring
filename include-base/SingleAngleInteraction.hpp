
#ifndef PROLIF_COLORING_SINGLEANGLE_INTERACTION
#define PROLIF_COLORING_SINGLEANGLE_INTERACTION

#include <GraphMol/GraphMol.h>
#include <Interaction.hpp>

/**
 * This class extends Interaction class
 * Defines a type of interaction in which the space the interaction is acting on is defined by:
 *      - distance(point_in_space, centroid_of_interaction) <= reference_distance
 *      - Be p1,p2 two centroid of interaction on the molecule,
 *           l1 point in space,
 *           p2l1 the vector from p2 to l1,
 *           p2p1 the vector from p2 to p1,
 *        reference_min_angle <= angle(p2p1,p2l1) <= reference_max_angle)
 */
class SingleAngleInteraction : public Interaction {
private:
    /**
     * Reference min/max angles and distance
     */
    const double min_angle, max_angle, distance;

    /**
     * Variable 0/1 that defines if p1 or p2 has to be used as centroid for distance calculation
     */
    int cp;
public:

    /**
     * This constructor extends the Interaction class constructor by receiving also
     * the reference min/max angle and distance value
     * @param smart The input SMART definition for interaction match-pattern
     * @param angle The pair <reference_min_angle, reference_max_angle>
     * @param distance The reference distance for the interaction
     * @param centerPoint The switch-variable for distance centroid
     */
    SingleAngleInteraction(const std::string &smart,
                           std::pair<double, double> angle,
                           double distance,
                           int centerPoint) : Interaction(smart),
                                              min_angle(angle.first),
                                              max_angle(angle.second),
                                              distance(distance) {
        if (centerPoint == 0 || centerPoint == 1) cp = centerPoint;
        //Default centroid is p1
        else cp = 0;
    };

    /**
     * This function overrides the Interaction class one
     */
    bool getInteraction(const RDKit::ROMol *molecule,
                        MoleculeMesh &interactionMask, MoleculeMesh &subtractionMask) override;
};

#endif //PROLIF_COLORING_SINGLEANGLE_INTERACTION
