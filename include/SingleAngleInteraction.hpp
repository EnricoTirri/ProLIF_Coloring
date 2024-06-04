
#ifndef INTERACTION_ANGLE_H
#define INTERACTION_ANGLE_H

#include <GraphMol/GraphMol.h>
#include <Interaction.hpp>

class SingleAngleInteraction : public Interaction {
private:
    const double min_angle, max_angle, distance;
    int cp;
public:
    SingleAngleInteraction(const std::string &smart,
                           std::pair<double, double> angle,
                           const int centerPoint,
                           double distance) : Interaction(smart),
                                              min_angle(angle.first),
                                              max_angle(angle.second),
                                              distance(distance) {
        if (centerPoint == 0 || centerPoint == 1) cp = centerPoint;
        else cp = 0;
    };

    bool getInteraction(const RDKit::ROMol *molecule,
                        MoleculeMesh &interactionMask, MoleculeMesh &subtractionMask) override;
};

#endif //INTERACTION_ANGLE_H
