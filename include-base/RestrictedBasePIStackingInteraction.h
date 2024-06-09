
#ifndef PROLIF_COLORING_RBSP_INTERACTION
#define PROLIF_COLORING_RBSP_INTERACTION

#include <GraphMol/GraphMol.h>
#include <Interaction.hpp>

class RestrictedBasePIStackingInteraction : public Interaction {
private:
    const double distance,
            min_angle_ring, max_angle_ring,
            min_angle_cent, max_angle_cent,
            intersect_radius;
    bool intersect;
    int cp;
public:
    RestrictedBasePIStackingInteraction(
            const double distance,
            const std::pair<double, double> angle,
            const std::pair<double, double> normal_to_centroid_angle,
            const std::string &pi_ring,
            bool intersect = false,
            double intersect_radius = 1.5) : Interaction(pi_ring),
                                             distance(distance),
                                             min_angle_ring(angle.first),
                                             max_angle_ring(angle.second),
                                             min_angle_cent(normal_to_centroid_angle.first),
                                             max_angle_cent(normal_to_centroid_angle.second),
                                             intersect_radius(intersect_radius),
                                             intersect(intersect) {};

    bool getInteraction(const RDKit::ROMol *molecule,
                        MoleculeMesh &interactionMask, MoleculeMesh &subtractionMask) override;
};

#endif //PROLIF_COLORING_RBSP_INTERACTION
