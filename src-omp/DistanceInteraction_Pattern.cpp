
#include "DistanceInteraction.hpp"
#include "omp.h"

#if USEPATTERN

#define V1 0

#if V1
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

    //interactionMask.sub(subtractionMask, 0, 0, 0);

    return found;
}
#endif

#if !V1

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

    double scaledDistance, ds;
    int layerS, maskS;
    double paddingDisplacement;
    RDGeom::Point3D pos;


    bool found = false;
    // Over all size of pattern-mesh assign if (point-distance <= #distance) from the center of mesh
#pragma omp parallel default(shared)
    {

#pragma omp single
        {
            scaledDistance = distance * GRAIN;
            ds = scaledDistance * scaledDistance;
        }

#pragma omp single
        {
            layerS = maskDim * maskDim;
            maskS = layerS * maskDim;
        }

#pragma omp single
        {
            paddingDisplacement = static_cast<double>(interactionMask.internalDisplacement - scaledMaskRadius);
        }

#pragma omp barrier

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

        for (int i = 0; i < matches->size(); ++i) {
            RDKit::MatchVectType match = (*matches)[i];
            if (!match.empty()) {

#pragma omp single
                { found = true; }

#pragma omp single
                {
                    int atomId = match.at(0).second;
                    pos = conformer.getAtomPos(atomId);
                }

#pragma omp barrier


                double px = (pos.x - interactionMask.globalDisplacement.x) * GRAIN + paddingDisplacement;
                int displ_x = static_cast<int>(round(px));

                double py = (pos.y - interactionMask.globalDisplacement.y) * GRAIN + paddingDisplacement;
                int displ_y = static_cast<int>(round(py));

                double pz = (pos.z - interactionMask.globalDisplacement.z) * GRAIN + paddingDisplacement;
                int displ_z = static_cast<int>(round(pz));

                int sx = displ_x < 0 ? 0 : displ_x;
                int sy = displ_y < 0 ? 0 : displ_y;
                int sz = displ_z < 0 ? 0 : displ_z;

                int ex = maskDim + displ_x < interactionMask.dim_x ? maskDim + displ_x : interactionMask.dim_x;
                int ey = maskDim + displ_y < interactionMask.dim_y ? maskDim + displ_y : interactionMask.dim_y;
                int ez = maskDim + displ_z < interactionMask.dim_z ? maskDim + displ_z : interactionMask.dim_z;

                int xR = ex - sx;
                int yR = ey - sy;
                int zR = ez - sz;
                int R = xR * yR * zR;
                int L = xR * yR;

#pragma omp barrier

#pragma omp for simd
                for (int k = 0; k < R; ++k) {
                    int z = k / L + sz;
                    int y = (k % L) / xR + sy;
                    int x = (k % L) % xR + sx;

                    int az = z - displ_z;
                    int ay = y - displ_y;
                    int ax = x - displ_x;
                    if(bubble.at(ax,ay,az))
                        interactionMask.at(x, y, z) = true;
                }
            }
        }
    }

    //interactionMask.sub(subtractionMask, 0, 0, 0);

    return found;
}

#endif

#endif