
#include "DistanceInteraction.hpp"

#if !USEPATTERN

#if MP

bool DistanceInteraction::getInteraction(const RDKit::ROMol *molecule, MoleculeMesh &interactionMask,
                                         MoleculeMesh &subtractionMask) {
    // Get molecule conformer and retrive matches of smart into given molecule
    RDKit::Conformer conformer = molecule->getConformer();
    std::vector<RDKit::MatchVectType> *matches = Interaction::findMatch(molecule);

    if (matches->empty()) return false;

    bool found = false;

    int i;

    int scaledMaskRadius;
    int maskDim;

    double scaledDistance;
    double ds;
    double paddingDisplacement;

    RDGeom::Point3D pos;


    int displ_x, displ_y, displ_z;
    int cx, cy, cz;

    // calculate operative window
    int sx, sy, sz, ex, ey, ez;

    int xR;

    int R, L;

    int atomId;


#pragma omp parallel default(shared) private(i)
    {


#pragma omp single
        {
            // Discretize mask radius and calculate mask dimension
            scaledMaskRadius = static_cast<int>(ceil(distance * GRAIN));
            maskDim = 2 * scaledMaskRadius;
            paddingDisplacement = static_cast<double>(interactionMask.internalDisplacement - scaledMaskRadius);
        }


        // Over all size of pattern-mesh assign if (point-distance <= #distance) from the center of mesh
#pragma omp single
        {
            scaledDistance = distance * GRAIN;
            ds = scaledDistance * scaledDistance;
        }

#pragma omp barrier

        for (i = 0; i < matches->size(); ++i) {
            RDKit::MatchVectType match = (*matches)[i];
            if (!match.empty()) {


#pragma omp single
                { found = true; }

#pragma omp single
                {
                    atomId = match.at(0).second;
                    pos = conformer.getAtomPos(atomId);
                }

#pragma omp barrier

#pragma omp single
                {
                    double px = (pos.x - interactionMask.globalDisplacement.x) * GRAIN + paddingDisplacement;
                    displ_x = static_cast<int>(round(px));
                }

#pragma omp single
                {
                    double py = (pos.y - interactionMask.globalDisplacement.y) * GRAIN + paddingDisplacement;
                    displ_y = static_cast<int>(round(py));
                }

#pragma omp single
                {
                    double pz = (pos.z - interactionMask.globalDisplacement.z) * GRAIN + paddingDisplacement;
                    displ_z = static_cast<int>(round(pz));
                }

#pragma omp barrier


#pragma omp single
                { sx = displ_x < 0 ? 0 : displ_x; }

#pragma omp single
                { sy = displ_y < 0 ? 0 : displ_y; }

#pragma omp single
                { sz = displ_z < 0 ? 0 : displ_z; }

#pragma omp single
                { ex = maskDim + displ_x < interactionMask.dim_x ? maskDim + displ_x : interactionMask.dim_x; }

#pragma omp single
                { ey = maskDim + displ_y < interactionMask.dim_y ? maskDim + displ_y : interactionMask.dim_y; }

#pragma omp single
                { ez = maskDim + displ_z < interactionMask.dim_z ? maskDim + displ_z : interactionMask.dim_z; }

#pragma omp barrier

#pragma omp single
                {
                    xR = ex - sx;
                    int yR = ey - sy;
                    int zR = ez - sz;
                    R = xR * yR * zR;
                    L = xR * yR;
                }

#pragma omp single
                {
                    cz = displ_z + scaledMaskRadius;
                    cy = displ_y + scaledMaskRadius;
                    cx = displ_x + scaledMaskRadius;
                }

#pragma omp barrier

#pragma omp for simd
                for (int k = 0; k < R; ++k) {
                    int z = k / L + sz;
                    int y = (k % L) / xR + sy;
                    int x = (k % L) % xR + sx;

                    int dz = z - cz;
                    int z_res = dz * dz;
                    int dy = y - cy;
                    int y_res = dy * dy;
                    int dx = x - cx;
                    int x_res = dx * dx;
                    if (x_res + y_res + z_res <= ds)
                        interactionMask.at(x, y, z) = true;
                }
            }
        }
    }

    return found;
}

#endif

#if !MP

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


#pragma omp parallel for shared(found)
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
                            interactionMask.at(x,y,z) = true;
                    }
                }
            }
        }
    }

    return found;
}

#endif

#endif