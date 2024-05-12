//
// Created by Enrico on 29/03/2024.
//

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

    bool getInteraction(const RDKit::ROMol *molecule, MoleculeMesh &mask) override {

        // Get molecule conformer and retrive matches of smart into given molecule
        RDKit::Conformer conformer = molecule->getConformer();
        std::vector<RDKit::MatchVectType> *matches = Interaction::findMatch(molecule);

        if (matches->empty()) return false;

        // Discretize mask radius and calculate mask dimension
        int scaledMaskRadius = static_cast<int>(ceil(distance * GRAIN));
        int maskDim = 2 * scaledMaskRadius;

        // Generate a pattern-mesh
        MoleculeMesh bubble(maskDim, maskDim, maskDim);

        // Over all size of pattern-mesh assign if (point-distance <= #distance) from the center of mesh
        double scaledDistance = distance * GRAIN;
        double ds = scaledDistance * scaledDistance;
        for (int z = 0; z < maskDim; ++z) {
            int dz = z - scaledMaskRadius;
            int z_res = dz * dz;
            for (int y = 0; y < maskDim; ++y) {
                int dy = y - scaledMaskRadius;
                int y_res = dy * dy;
                for (int x = 0; x < maskDim; ++x) {
                    int dx = x - scaledMaskRadius;
                    int x_res = dx * dx;
                    if (x_res + y_res + z_res <= ds)
                        bubble.at(x, y, z) = true;
                }
            }
        }

        // For each interaction-centroid apply the pattern-mesh centered at centroid onto the support-mesh
        auto paddingDisplacement = static_cast<double>(mask.internalDisplacement - scaledMaskRadius);
        bool found = false;
        for (RDKit::MatchVectType match: *matches) {
            if (!match.empty()) {
                found = true;
                // Get interaction match centroid position
                auto atomId = match.at(0).second;
                RDGeom::Point3D pos = conformer.getAtomPos(atomId);

                // Find the zero-point displacement of pattern from the zero-point of support-mask
                double px = (pos.x - mask.globalDisplacement.x) * GRAIN + paddingDisplacement;
                double py = (pos.y - mask.globalDisplacement.y) * GRAIN + paddingDisplacement;
                double pz = (pos.z - mask.globalDisplacement.z) * GRAIN + paddingDisplacement;

                // Discretize the displacement
                int displ_x = static_cast<int>(round(px));
                int displ_y = static_cast<int>(round(py));
                int displ_z = static_cast<int>(round(pz));

                // Apply pattern at displacement onto support-mesh
                mask.add(bubble, displ_x, displ_y, displ_z);
            }
        }

        return found;
    }
};

#endif //INTERACTION_DISTANCE_H
