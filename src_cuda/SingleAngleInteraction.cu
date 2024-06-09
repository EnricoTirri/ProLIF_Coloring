
#include "SingleAngleInteraction.hpp"
#include <cuda/std/cmath>

__device__
double lengthSq(double x, double y, double z) {
    return x * x + y * y + z * z;
}

__device__
double dotProduct(double p1_x, double p1_y, double p1_z,
                  double p2_x, double p2_y, double p2_z) {
    return p1_x * (p2_x) + p1_y * (p2_y) + p1_z * (p2_z);
}

__device__
double angleTo(double p1_x, double p1_y, double p1_z,
               double p2_x, double p2_y, double p2_z) {
    double lsq = lengthSq(p1_x, p1_y, p1_z) * lengthSq(p2_x, p2_y, p2_z);
    double dotProd = dotProduct(p1_x, p1_y, p1_z, p2_x, p2_y, p2_z);
    dotProd /= cuda::std::sqrt(lsq);

    // watch for roundoff error:
    if (dotProd <= -1.0) {
        return M_PI;
    }
    if (dotProd >= 1.0) {
        return 0.0;
    }

    return cuda::std::acos(dotProd);
}

__global__
void buildBubbleSlice_ker(MoleculeMesh::data_t *bubble,
                          const double inter_d, const double min_angle, const double max_angle,
                          const double center_x, const double center_y, const double center_z,
                          const double p1_x, const double p1_y, const double p1_z,
                          const double p2_x, const double p2_y, const double p2_z,
                          const int maskEdge) {

    int thr_id = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);
    const int layerDim = maskEdge * maskEdge;

    if (thr_id < layerDim * maskEdge) {

        const int z_cord = thr_id / layerDim;
        const int y_cord = (thr_id % layerDim) / maskEdge;
        const int x_cord = (thr_id % layerDim) % maskEdge;

        const int maskRadius = maskEdge / 2;

        // Over all size of pattern-mesh assign if (point-distance <= #distance) from the center of mesh

        // Calculate vector p2 --> p1
        double p2p1_x = p1_x - p2_x;
        double p2p1_y = p1_y - p2_y;
        double p2p1_z = p1_z - p2_z;

        double p2p1_l = cuda::std::sqrt(lengthSq(p2p1_x, p2p1_y, p2p1_z));
        p2p1_x /= p2p1_l;
        p2p1_y /= p2p1_l;
        p2p1_z /= p2p1_l;

        /*
         * Over all size of pattern-mesh assign if:
         *      - (point-distance <= #distance) from the center of mesh
         *      - angle between l1 <-- p2 --> p1 is (#min <= #angle <= #max)
         */
        double ds = inter_d * inter_d;
        int dz = z_cord - maskRadius;
        int z_res = dz * dz;
        int dy = y_cord - maskRadius;
        int y_res = dy * dy;
        int dx = x_cord - maskRadius;
        int x_res = dx * dx;

        // Position of l1 is calculated taking account of pattern-center position

        double l1_z = z_cord + center_z;
        double l1_y = y_cord + center_y;
        double l1_x = x_cord + center_x;

        // Calculate vector p2 --> l1
        double p2l1_x = l1_x - p2_x;
        double p2l1_y = l1_y - p2_y;
        double p2l1_z = l1_z - p2_z;

        double p2l1_l = cuda::std::sqrt(lengthSq(p2l1_x, p2l1_y, p2l1_z));
        p2l1_x /= p2l1_l;
        p2l1_y /= p2l1_l;
        p2l1_z /= p2l1_l;

        // Calculate angle l1 <-- p2 --> p1
        double angle = angleTo(p2p1_x, p2p1_y, p2p1_z, p2l1_x, p2l1_y, p2l1_z);

        bubble[thr_id] = (x_res + y_res + z_res <= ds && angle >= min_angle && angle <= max_angle);
    }
}

bool SingleAngleInteraction::getInteraction(const RDKit::ROMol *molecule,
                                            MoleculeMesh &interactionMask, MoleculeMesh &subtractionMask) {
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

        // Calculate mask size and centering coordinates
        double scaledDistance = distance * GRAIN;
        auto scaledMaskCenter = static_cast<int>(ceil(scaledDistance));
        int maskDim = 2 * scaledMaskCenter;


        int bubbleDim = maskDim * maskDim * maskDim;

        unsigned int numBlocks = (interactionMask.getDataSize() + BLOCK_SIZE) / BLOCK_SIZE;

        err = cudaMalloc((void **) &bubble_data, sizeof(MoleculeMesh::data_t) * (bubbleDim));
        if (err != cudaSuccess) throw;

        err = cudaMalloc((void **) &interaction_data, sizeof(MoleculeMesh::data_t) * interactionMask.getDataSize());
        if (err != cudaSuccess) throw;

        for (RDKit::MatchVectType match: *matches) {
            if (match.size() >= 2) {
                ris = true;

                // Get molecule match and its centroids position
                auto p1Id = match.at(0).second;
                auto p1 = conformer.getAtomPos(p1Id);
                auto p2Id = match.at(1).second;
                auto p2 = conformer.getAtomPos(p2Id);

                RDGeom::Point3D center;
                if (cp) center = p1;
                else center = p2;

                buildBubbleSlice_ker<<<numBlocks, BLOCK_SIZE>>>(bubble_data,
                                                                scaledDistance, min_angle, max_angle,
                                                                center.x, center.y, center.z,
                                                                p1.x, p1.y, p1.z,
                                                                p2.x, p2.y, p2.z,
                                                                maskDim);
                err = cudaGetLastError();
                if (err != cudaSuccess) throw;

                // Find the zero-point displacement of pattern from the zero-point of support-mask
                auto paddingDisplacement = static_cast<double>(interactionMask.internalDisplacement - scaledMaskCenter);
                double px = (center.x - interactionMask.globalDisplacement.x) * GRAIN + paddingDisplacement;
                double py = (center.y - interactionMask.globalDisplacement.y) * GRAIN + paddingDisplacement;
                double pz = (center.z - interactionMask.globalDisplacement.z) * GRAIN + paddingDisplacement;

                // Discretize the displacement
                int displ_x = static_cast<int>(round(px));
                int displ_y = static_cast<int>(round(py));
                int displ_z = static_cast<int>(round(pz));

                // Apply pattern at displacement onto support-mesh
                MoleculeMesh::addMeshes(interaction_data, bubble_data,
                                        displ_x, displ_y, displ_z,
                                        interactionMask.dim_x, interactionMask.dim_y, interactionMask.dim_z,
                                        maskDim, maskDim, maskDim);
                err = cudaGetLastError();
                if (err != cudaSuccess) throw;
            }
        }

        if (!subtractionMask.getDataSize()!=0) {
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
        if (err != cudaSuccess) throw;

    } catch (...) {
        ris = false;
        std::cout << cudaGetErrorString(err) << std::endl;
    }

    if (bubble_data != nullptr) cudaFree(bubble_data);
    if (interaction_data != nullptr) cudaFree(interaction_data);
    if (subtraction_data != nullptr) cudaFree(subtraction_data);

    return ris;
}