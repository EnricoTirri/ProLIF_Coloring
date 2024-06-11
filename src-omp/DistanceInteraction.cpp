
#include "DistanceInteraction.hpp"
#include <vector>

bool DistanceInteraction::getInteraction(const RDKit::ROMol *molecule, MoleculeMesh &interactionMask,
                                          MoleculeMesh &subtractionMask) {
    // Get molecule conformer and retrive matches of smart into given molecule
    RDKit::Conformer conformer = molecule->getConformer();
    std::vector<RDKit::MatchVectType> *matches = Interaction::findMatch(molecule);

    if (matches->empty()) return false;

    // Discretize mask radius and calculate mask dimension
    int scaledMaskRadius = static_cast<int>(ceil(distance * GRAIN));
    int maskDim = 2 * scaledMaskRadius;

    // Generate a pattern-mesh
    MoleculeMesh bubble(maskDim, maskDim, maskDim);

    bool found = false;
    // Over all size of pattern-mesh assign if (point-distance <= #distance) from the center of mesh
#pragma omp parallel
    {
        double scaledDistance = distance * GRAIN;
        double ds = scaledDistance * scaledDistance;

        int layerS = maskDim * maskDim;
        int maskS = layerS * maskDim;

#pragma omp for simd
        for (int k = 0; k < maskS; ++k) {
            int z = k / layerS;
            int y = (k % layerS) / maskDim;
            int x = (k % layerS) % maskDim;
            int dz = z - scaledMaskRadius;
            int z_res = dz * dz;
            int dy = y - scaledMaskRadius;
            int y_res = dy * dy;
            int dx = x - scaledMaskRadius;
            int x_res = dx * dx;
            if (x_res + y_res + z_res <= ds)
                bubble.at(x, y, z) = true;
        }

        // For each interaction-centroid apply the pattern-mesh centered at centroid onto the support-mesh
        auto paddingDisplacement = static_cast<double>(interactionMask.internalDisplacement - scaledMaskRadius);

#pragma omp for simd
        for (int i = 0; i < matches->size(); ++i) {
            RDKit::MatchVectType match = (*matches)[i];
            if (!match.empty()) {
                found = true;
                // Get interaction match centroid position
                auto atomId = match.at(0).second;
                RDGeom::Point3D pos = conformer.getAtomPos(atomId);

                // Find the zero-point displacement of pattern from the zero-point of support-mask
                double px = (pos.x - interactionMask.globalDisplacement.x) * GRAIN + paddingDisplacement;
                double py = (pos.y - interactionMask.globalDisplacement.y) * GRAIN + paddingDisplacement;
                double pz = (pos.z - interactionMask.globalDisplacement.z) * GRAIN + paddingDisplacement;

                // Discretize the displacement
                int displ_x = static_cast<int>(round(px));
                int displ_y = static_cast<int>(round(py));
                int displ_z = static_cast<int>(round(pz));

                // Apply pattern at displacement onto support-mesh
                interactionMask.add(bubble, displ_x, displ_y, displ_z);
            }
        }
    }

    interactionMask.sub(subtractionMask, 0, 0, 0);

    return found;
}