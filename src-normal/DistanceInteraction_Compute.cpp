
#include <GraphMol/FileParsers/FileParsers.h>
#include "DistanceInteraction.hpp"
#include "Discretizer.hpp"


#if !USEPATTERN

bool DistanceInteraction::getInteraction(const RDKit::ROMol *molecule, MoleculeMesh &interactionMask,
                                         MoleculeMesh &subtractionMask) {
    // Get molecule conformer and retrive matches of smart into given molecule
    RDKit::Conformer conformer = molecule->getConformer();
    std::vector<RDKit::MatchVectType> *matches = Interaction::findMatch(molecule);

    if (matches->empty()) return false;

    // Discretize mask radius and calculate mask dimension
    int scaledMaskRadius = static_cast<int>(ceil(distance * GRAIN));
    int maskDim = 2 * scaledMaskRadius;

    // For each interaction-centroid apply the pattern-mesh centered at centroid onto the support-mesh
    auto paddingDisplacement = static_cast<double>(interactionMask.internalDisplacement - scaledMaskRadius);
    bool found = false;
    for (RDKit::MatchVectType match: *matches) {
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

            // calculate operative window
            int sx = displ_x;
            if (sx < 0) {
                sx = 0;
            }

            int sy = displ_y;
            if (sy < 0) {
                sy = 0;
            }

            int sz = displ_z;
            if (sz < 0) {
                sz = 0;
            }

            int ex = interactionMask.dim_x;
            if (maskDim + displ_x < ex) {
                ex = maskDim + displ_x;
            }

            int ey = interactionMask.dim_y;
            if (maskDim + displ_y < ey) {
                ey = maskDim + displ_y;
            }

            int ez = interactionMask.dim_z;
            if (maskDim + displ_z < ez) {
                ez = maskDim + displ_z;
            }

            // Over all size of pattern-mesh assign if (point-distance <= #distance) from the center of mesh
            double scaledDistance = distance * GRAIN;
            double ds = scaledDistance * scaledDistance;

            /* Execute operation over operative window */
            for (int z = sz; z < ez; z++) {
                int cz = z - displ_z - scaledMaskRadius;
                int z_res = cz * cz;
                for (int y = sy; y < ey; y++) {
                    int cy = y - displ_y - scaledMaskRadius;
                    int y_res = cy * cy;
                    for (int x = sx; x < ex; x++) {
                        int cx = x - displ_x - scaledMaskRadius;
                        int x_res = cx * cx;
                        if(x_res + y_res + z_res <= ds)
                            interactionMask.at(x,y,z) = true;//!subtractionMask.at(x,y,z);
                    }
                }
            }
        }
    }

    //interactionMask.sub(subtractionMask, 0,0,0);

    return found;
}

#endif