
#include "Mesh.hpp"

__global__
void addMask_ker(MoleculeMesh::data_t *data, const MoleculeMesh::data_t *to_add,
                 const int displ_x, const int displ_y, const int displ_z,
                 const int data_dim_x, const int data_dim_y, const int data_dim_z,
                 const int add_dim_x, const int add_dim_y, const int add_dim_z) {

    const int thr_id = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);


    const int data_layer_size = data_dim_x * data_dim_y;
    const int data_z = thr_id / data_layer_size;

    if (data_z < data_dim_z) {
        const int data_y = (thr_id % data_layer_size) / data_dim_y;
        const int data_x = (thr_id % data_layer_size) % data_dim_y;

        const int add_z = data_z - displ_z;
        const int add_y = data_y - displ_y;
        const int add_x = data_x - displ_x;

        if (add_z >= 0 && add_y >= 0 && add_x >= 0 && add_z < add_dim_z && add_y < add_dim_y && add_x < add_dim_x) {
            const int data_id = data_dim_x * (data_z * data_dim_y + data_y) + data_x;
            const int add_id = add_dim_x * (add_z * add_dim_y + add_y) + add_x;
            data[data_id] = data[data_id] || to_add[add_id];
        }
    }
}

void MoleculeMesh::addMeshes(MoleculeMesh::data_t *data, MoleculeMesh::data_t *to_add,
                             const int displ_x, const int displ_y, const int displ_z,
                             const int data_dim_x, const int data_dim_y, const int data_dim_z,
                             const int add_dim_x, const int add_dim_y, const int add_dim_z) {

    int dataDim = data_dim_x * data_dim_y * data_dim_z;
    unsigned int numBlocks = (dataDim + BLOCK_SIZE) / (BLOCK_SIZE);

    addMask_ker<<<numBlocks, BLOCK_SIZE>>>(data, to_add, displ_x, displ_y, displ_z,
                                           data_dim_x, data_dim_y, data_dim_z,
                                           add_dim_x, add_dim_y, add_dim_z);
}

__global__
void subMask_ker(MoleculeMesh::data_t *data, const MoleculeMesh::data_t *to_subtract,
                 const int displ_x, const int displ_y, const int displ_z,
                 const int data_dim_x, const int data_dim_y, const int data_dim_z,
                 const int sub_dim_x, const int sub_dim_y, const int sub_dim_z) {

    int thr_id = static_cast<int>(blockIdx.x * blockDim.x + threadIdx.x);

    const int data_layer_size = data_dim_x * data_dim_y;
    const int data_z = thr_id / data_layer_size;

    if (data_z < data_dim_z) {
        const int data_y = (thr_id % data_layer_size) / data_dim_y;
        const int data_x = (thr_id % data_layer_size) % data_dim_y;

        const int sub_z = data_z + displ_z;
        const int sub_y = data_y + displ_y;
        const int sub_x = data_x + displ_x;

        if (sub_z >= 0 && sub_y >= 0 && sub_x >= 0 && sub_z < sub_dim_z && sub_y < sub_dim_y && sub_x < sub_dim_x) {
            data[thr_id] = data[thr_id] && !to_subtract[sub_dim_x * (sub_z * sub_dim_y + sub_y) + sub_x];
        }
    }
}

void MoleculeMesh::subMeshes(MoleculeMesh::data_t *data, MoleculeMesh::data_t *to_subtract,
                             const int displ_x, const int displ_y, const int displ_z,
                             const int data_dim_x, const int data_dim_y, const int data_dim_z,
                             const int sub_dim_x, const int sub_dim_y, const int sub_dim_z) {

    int dataDim = data_dim_x * data_dim_y * data_dim_z;
    unsigned int numBlocks = (dataDim + BLOCK_SIZE) / (BLOCK_SIZE);

    subMask_ker<<<numBlocks, BLOCK_SIZE>>>(data, to_subtract, displ_x, displ_y, displ_z,
                                           data_dim_x, data_dim_y, data_dim_z,
                                           sub_dim_x, sub_dim_y, sub_dim_z);
}
