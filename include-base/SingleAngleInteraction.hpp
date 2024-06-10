
#ifndef PROLIF_COLORING_SINGLEANGLE_INTERACTION
#define PROLIF_COLORING_SINGLEANGLE_INTERACTION

#include <GraphMol/GraphMol.h>
#include <Interaction.hpp>

class SingleAngleInteraction : public Interaction {
private:
    const double min_angle, max_angle, distance;
    int cp;
public:
    SingleAngleInteraction(const std::string &smart,
                           std::pair<double, double> angle,
                           double distance,
                           const int centerPoint) : Interaction(smart),
                                              min_angle(angle.first),
                                              max_angle(angle.second),
                                              distance(distance) {
        if (centerPoint == 0 || centerPoint == 1) cp = centerPoint;
        else cp = 0;
    };

    bool getInteraction(const RDKit::ROMol *molecule,
                        MoleculeMesh &interactionMask, MoleculeMesh &subtractionMask) override;
};

#endif //PROLIF_COLORING_SINGLEANGLE_INTERACTION
