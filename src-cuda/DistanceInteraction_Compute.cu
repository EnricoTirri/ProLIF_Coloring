

#include <GraphMol/FileParsers/FileParsers.h>
#include "DistanceInteraction.hpp"
#include "Discretizer.hpp"

#if !USEPATTERN

__global__
void addBubble(MoleculeMesh::data_t *mask, int displ_x, int displ_y, int displ_z, int dim_x, int dim_y, int dim_z,
               double dist, int maskRad) {

    auto maskDim = 2 * maskRad;

    int k = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);

    int sx = displ_x < 0 ? 0 : displ_x;
    int sy = displ_y < 0 ? 0 : displ_y;
    int sz = displ_z < 0 ? 0 : displ_z;
    int ex = maskDim + displ_x < dim_x ? maskDim + displ_x : dim_x;
    int ey = maskDim + displ_y < dim_y ? maskDim + displ_y : dim_y;
    int ez = maskDim + displ_z < dim_z ? maskDim + displ_z : dim_z;
    int xR = ex - sx;
    int yR = ey - sy;
    int zR = ez - sz;
    int R = xR * yR * zR;
    int L = xR * yR;


    int cz = displ_z + maskRad;
    int cy = displ_y + maskRad;
    int cx = displ_x + maskRad;

    double ds = dist * dist;

    if (k < R) {
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
            mask[dim_x * (z * dim_y + y) + x] = true;
    }
}


bool DistanceInteraction::getInteraction(const RDKit::ROMol *molecule, MoleculeMesh &interactionMask,
                                         MoleculeMesh &/*subtractionMask*/) {

    MoleculeMesh::data_t *interaction_data = nullptr;
    bool ris = false;

    try {
        // Get molecule conformer and retrive matches of smart into given molecule
        RDKit::Conformer conformer = molecule->getConformer();
        std::vector<RDKit::MatchVectType> *matches = Interaction::findMatch(molecule);

        if (matches->empty()) return false;

        // Discretize mask radius and calculate mask dimension
        // Over all size of pattern-mesh assign if (point-distance <= #distance) from the center of mesh
        double scaledDistance = distance * GRAIN;
        int scaledMaskRadius = static_cast<int>(ceil(scaledDistance));
        int maskDim = 2 * scaledMaskRadius;

        int bubbleDim = maskDim * maskDim * maskDim;

        unsigned int numBlocks = (bubbleDim + BLOCK_SIZE) / BLOCK_SIZE;
        cudaMalloc((void **) &interaction_data, sizeof(MoleculeMesh::data_t) * interactionMask.getDataSize());

        auto paddingDisplacement = static_cast<double>(interactionMask.internalDisplacement - scaledMaskRadius);

        for (RDKit::MatchVectType match: *matches) {
            if (!match.empty()) {

                ris = true;
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


                addBubble<<<numBlocks, BLOCK_SIZE, 0>>>(interaction_data,
                                                     displ_x, displ_y, displ_z,
                                                     interactionMask.dim_x, interactionMask.dim_y, interactionMask.dim_z,
                                                     scaledDistance, scaledMaskRadius);
            }
        }
        cudaMemcpy(interactionMask.getData(), interaction_data,
                     sizeof(MoleculeMesh::data_t) * interactionMask.getDataSize(), cudaMemcpyDeviceToHost);

    } catch (...) {
        ris = false;
    }

    if (interaction_data != nullptr) cudaFree(interaction_data);

    return ris;
}

#endif