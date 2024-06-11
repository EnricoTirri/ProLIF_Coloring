
#include <GraphMol/FileParsers/FileParsers.h>
#include "DistanceInteraction.hpp"
#include "Discretizer.hpp"

__global__
void buildBubble_ker(MoleculeMesh::data_t *bubble, const double inter_d, const int maskEdge) {
    int thr_id = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    const int layerDim = maskEdge * maskEdge;

    if (thr_id < layerDim * maskEdge) {
        const int z_cord = thr_id / layerDim;
        const int y_cord = (thr_id % layerDim) / maskEdge;
        const int x_cord = (thr_id % layerDim) % maskEdge;

        const int maskRadius = maskEdge / 2;

        // Over all size of pattern-mesh assign if (point-distance <= #distance) from the center of mesh
        double ds = inter_d * inter_d;
        int dz = z_cord - maskRadius;
        int z_res = dz * dz;
        int dy = y_cord - maskRadius;
        int y_res = dy * dy;
        int dx = x_cord - maskRadius;
        int x_res = dx * dx;

        bubble[thr_id] = (x_res + y_res + z_res <= ds);
    }
}


bool DistanceInteraction::getInteraction(const RDKit::ROMol *molecule, MoleculeMesh &interactionMask,
                                         MoleculeMesh &subtractionMask) {
    cudaError_t err;
    MoleculeMesh::data_t *bubble_data = nullptr;
    MoleculeMesh::data_t *interaction_data = nullptr;
    MoleculeMesh::data_t *subtraction_data = nullptr;
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

        err = cudaMalloc((void **) &bubble_data, sizeof(MoleculeMesh::data_t) * (bubbleDim));
        if (err != cudaSuccess) throw;

        unsigned int numBlocks = (interactionMask.getDataSize() + BLOCK_SIZE) / BLOCK_SIZE;


        buildBubble_ker<<<numBlocks, BLOCK_SIZE>>>(bubble_data, scaledDistance, maskDim);
        err = cudaGetLastError();
        if (err != cudaSuccess) throw;


        err = cudaMalloc((void **) &interaction_data, sizeof(MoleculeMesh::data_t) * interactionMask.getDataSize());
        if (err != cudaSuccess) throw;

        // For each interaction-centroid apply the pattern-mesh centered at centroid onto the support-mesh
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

                MoleculeMesh::addMeshes(interaction_data, bubble_data,
                                        displ_x, displ_y, displ_z,
                                        interactionMask.dim_x, interactionMask.dim_y, interactionMask.dim_z,
                                        maskDim, maskDim, maskDim);
                err = cudaGetLastError();
                if (err != cudaSuccess) throw;
            }
        }

        if (subtractionMask.getDataSize()!=0) {
            err = cudaMalloc((void **) &subtraction_data, sizeof(MoleculeMesh::data_t) * subtractionMask.getDataSize());
            if (err != cudaSuccess) throw;

            err = cudaMemcpy(subtraction_data, subtractionMask.getData(),
                             sizeof(MoleculeMesh::data_t) * subtractionMask.getDataSize(), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) throw;

            MoleculeMesh::subMeshes(interaction_data, subtraction_data,
                                    0, 0, 0,
                                    interactionMask.dim_x, interactionMask.dim_y, interactionMask.dim_z,
                                    subtractionMask.dim_x, subtractionMask.dim_y, subtractionMask.dim_z);
            err = cudaGetLastError();
            if (err != cudaSuccess) throw;
        }


        err = cudaMemcpy(interactionMask.getData(), interaction_data,
                         sizeof(MoleculeMesh::data_t) * interactionMask.getDataSize(), cudaMemcpyDeviceToHost);
        if(err != cudaSuccess) throw;

    } catch (...) {
        ris = false;
        std::cout << cudaGetErrorString(err) << std::endl;
    }

    if(bubble_data != nullptr) cudaFree(bubble_data);
    if(interaction_data != nullptr) cudaFree(interaction_data);
    if(subtraction_data != nullptr) cudaFree(subtraction_data);

    return ris;
}