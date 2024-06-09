#ifndef PROLIF_COLORING_MESH
#define PROLIF_COLORING_MESH

#include <vector>
#include "Geometry/point.h"

#define GRAIN 2

class MoleculeMesh {
public:
    typedef int data_t;

private:
    std::vector<data_t> voxels;

public:
    const int dim_x, dim_y, dim_z;
    RDGeom::Point3D globalDisplacement;
    int internalDisplacement;

    MoleculeMesh(int p_dim_x, int p_dim_y, int p_dim_z, const RDGeom::Point3D &globalDisplacement,
                 int internalDisplacement) :
            dim_x(p_dim_x),
            dim_y(p_dim_y),
            dim_z(p_dim_z),
            globalDisplacement(globalDisplacement),
            internalDisplacement(internalDisplacement) {
        voxels = std::vector<data_t>(dim_x * dim_y * dim_z);
    }

    MoleculeMesh(int p_dim_x, int p_dim_y, int p_dim_z) :
            MoleculeMesh(p_dim_x, p_dim_y, p_dim_z, {0, 0, 0}, 0) {}

    inline size_t getDataSize() {
        return voxels.size();
    }

    inline data_t *getData() {
        return voxels.data();
    }

    inline data_t &at(int x, int y, int z) {
        return MoleculeMesh::ref(voxels.data(), x, y, z, dim_x, dim_y, dim_z);
    }

    inline static data_t &ref(MoleculeMesh::data_t *data, int x, int y, int z,
                              const int dim_x, const int dim_y, const int /*dim_z*/) {
        return data[dim_x * (z * dim_y + y) + x];
    }

    inline void add(MoleculeMesh &addend, int displ_x, int displ_y, int displ_z) {
        MoleculeMesh::addMeshes(getData(), addend.getData(),
                                displ_x, displ_y, displ_z,
                                dim_x, dim_y, dim_z,
                                addend.dim_x, addend.dim_y, addend.dim_z);
    }

    inline void sub(MoleculeMesh &addend, int displ_x, int displ_y, int displ_z) {
        MoleculeMesh::subMeshes(getData(), addend.getData(),
                                displ_x, displ_y, displ_z,
                                dim_x, dim_y, dim_z,
                                addend.dim_x, addend.dim_y, addend.dim_z);
    }

    static void addMeshes(MoleculeMesh::data_t *data, MoleculeMesh::data_t *to_add,
                          int displ_x, int displ_y, int displ_z,
                          int data_dim_x, int data_dim_y, int data_dim_z,
                          int add_dim_x, int add_dim_y, int add_dim_z);

    static void subMeshes(MoleculeMesh::data_t *data, MoleculeMesh::data_t *to_subtract,
                          int displ_x, int displ_y, int displ_z,
                          int data_dim_x, int data_dim_y, int data_dim_z,
                          int sub_dim_x, int sub_dim_y, int sub_dim_z);
};

#endif //PROLIF_COLORING_MESH
