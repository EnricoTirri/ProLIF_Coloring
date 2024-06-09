
#include "Mesh.hpp"

void MoleculeMesh::subMeshes(MoleculeMesh::data_t *data, MoleculeMesh::data_t *to_subtract,
                             const int displ_x, const int displ_y, const int displ_z,
                             const int data_dim_x, const int data_dim_y, const int data_dim_z,
                             const int sub_dim_x, const int sub_dim_y, const int sub_dim_z) {

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

    int ex = data_dim_x;
    if (sub_dim_x + displ_x < data_dim_x) {
        ex = sub_dim_x + displ_x;
    }

    int ey = data_dim_y;
    if (sub_dim_y + displ_y < data_dim_y) {
        ey = sub_dim_y + displ_y;
    }

    int ez = data_dim_z;
    if (sub_dim_z + displ_z < data_dim_z) {
        ez = sub_dim_z + displ_z;
    }

    /* Execute operation over operative window */
    for (int z = sz; z < ez; z++) {
        int az = z - displ_z;
        for (int y = sy; y < ey; y++) {
            int ay = y - displ_y;
            for (int x = sx; x < ex; x++) {
                int ax = x - displ_x;
                MoleculeMesh::ref(data, x, y, z, data_dim_x, data_dim_y, data_dim_z) =
                        MoleculeMesh::ref(data, x, y, z, data_dim_x, data_dim_y, data_dim_z)
                        && !MoleculeMesh::ref(to_subtract, ax, ay, az, sub_dim_x, sub_dim_y, sub_dim_z);
            }
        }
    }
}

void MoleculeMesh::addMeshes(MoleculeMesh::data_t *data, MoleculeMesh::data_t *to_add,
                             const int displ_x, const int displ_y, const int displ_z,
                             const int data_dim_x, const int data_dim_y, const int data_dim_z,
                             const int add_dim_x, const int add_dim_y, const int add_dim_z) {

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

    int ex = data_dim_x;
    if (add_dim_x + displ_x < data_dim_x) {
        ex = add_dim_x + displ_x;
    }

    int ey = data_dim_y;
    if (add_dim_y + displ_y < data_dim_y) {
        ey = add_dim_y + displ_y;
    }

    int ez = data_dim_z;
    if (add_dim_z + displ_z < add_dim_z) {
        ez = add_dim_z + displ_z;
    }

    /* Execute operation over operative window */
    for (int z = sz; z < ez; z++) {
        int az = z - displ_z;
        for (int y = sy; y < ey; y++) {
            int ay = y - displ_y;
            for (int x = sx; x < ex; x++) {
                int ax = x - displ_x;
                MoleculeMesh::ref(data, x, y, z, data_dim_x, data_dim_y, data_dim_z) =
                        MoleculeMesh::ref(data, x, y, z, data_dim_x, data_dim_y, data_dim_z)
                        || MoleculeMesh::ref(to_add, ax, ay, az, add_dim_x, add_dim_y, add_dim_z);
            }
        }
    }
}
