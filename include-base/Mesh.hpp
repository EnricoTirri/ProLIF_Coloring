#ifndef PROLIF_COLORING_MESH
#define PROLIF_COLORING_MESH

#include <vector>
#include "Geometry/point.h"

#ifndef GRAIN
#define GRAIN 3
#endif

/**
 * This class defines the model for the discrete molecule
 */
class MoleculeMesh {
public:
    /**
     * The type of data the discrete space is based on
     */
    typedef int data_t;

private:
    /**
     * The data structure that contains the discrete space description
     */
    std::vector<data_t> voxels;

public:
    /**
     * The 3D sizes of the discrete space
     */
    const int dim_x, dim_y, dim_z;

    /**
     * The displacement this discrete space have in relation to a "global" one
     * (this is useful if we have to manage independently multiple discrete spaces all related to each other)
     */
    RDGeom::Point3D globalDisplacement;

    /**
     * The displacement data have internally in this discrete space
     * (this is useful if we have to manage paddings due to transformation from continuous space to discrete)
     */
    int internalDisplacement;

    /**
     * This constructor initialize the discrete space
     * @param p_dim_x X dimension of the space
     * @param p_dim_y Y dimension of the space
     * @param p_dim_z Z dimension of the space
     * @param globalDisplacement Global displacement of the space
     * @param internalDisplacement Internal displacement of data
     */
    MoleculeMesh(int p_dim_x, int p_dim_y, int p_dim_z, const RDGeom::Point3D &globalDisplacement,
                 int internalDisplacement) :
            dim_x(p_dim_x),
            dim_y(p_dim_y),
            dim_z(p_dim_z),
            globalDisplacement(globalDisplacement),
            internalDisplacement(internalDisplacement) {
        voxels = std::vector<data_t>(dim_x * dim_y * dim_z);
    }

    /**
     * This constructor initialize a discrete space with no Global nor Internal displacement
     * @param p_dim_x X dimension of the space
     * @param p_dim_y Y dimension of the space
     * @param p_dim_z Z dimension of the space
     */
    MoleculeMesh(int p_dim_x, int p_dim_y, int p_dim_z) :
            MoleculeMesh(p_dim_x, p_dim_y, p_dim_z, {0, 0, 0}, 0) {}

    /**
     * This function returns the number of data the space contains
     * @return
     */
    inline size_t getDataSize() {
        return voxels.size();
    }

    /**
     * This function returns the data of the space
     * @return
     */
    inline data_t *getData() {
        return voxels.data();
    }

    /**
     * This function returns the data at a specific discrete position of the space
     * @param x X discrete coordinates
     * @param y Y discrete coordinates
     * @param z Z discrete coordinates
     * @return The data at (X,Y,Z) discrete position in space
     */
    inline data_t &at(int x, int y, int z) {
        return MoleculeMesh::ref(voxels.data(), x, y, z, dim_x, dim_y, dim_z);
    }

    /**
     * This function defines how space data structure is managed in relation of spatial access
     * @param data The space data structure
     * @param x X discrete coordinates
     * @param y Y discrete coordinates
     * @param z Z discrete coordinates
     * @param dim_x The X dimension of space data structure
     * @param dim_y The Y dimension of space data structure
     * @param dim_z the Z dimension of space data structure
     * @return The data at (X,Y,Z) discrete position in input space data structure
     */
    inline static data_t &ref(MoleculeMesh::data_t *data, int x, int y, int z,
                              const int dim_x, const int dim_y, const int /*dim_z*/) {
        return data[dim_x * (z * dim_y + y) + x];
    }

    /**
     * This function allow to integrate a discrete space performing a boolean addition to the class managed one
     * @param addend The discrete space we want to integrate
     * @param displ_x The X displacement we want the input space to be placed
     * @param displ_y The Y displacement we want the input space to be placed
     * @param displ_z The Z displacement we want the input space to be placed
     */
    inline void add(MoleculeMesh &addend, int displ_x, int displ_y, int displ_z) {
        MoleculeMesh::addMeshes(getData(), addend.getData(),
                                displ_x, displ_y, displ_z,
                                dim_x, dim_y, dim_z,
                                addend.dim_x, addend.dim_y, addend.dim_z);
    }

    /**
     * This function allow to integrate a discrete space performing a boolean subtraction to the class managed one
     * @param addend The discrete space we want to integrate
     * @param displ_x The X displacement we want the input space to be placed
     * @param displ_y The Y displacement we want the input space to be placed
     * @param displ_z The Z displacement we want the input space to be placed
     */
    inline void sub(MoleculeMesh &addend, int displ_x, int displ_y, int displ_z) {
        MoleculeMesh::subMeshes(getData(), addend.getData(),
                                displ_x, displ_y, displ_z,
                                dim_x, dim_y, dim_z,
                                addend.dim_x, addend.dim_y, addend.dim_z);
    }

    /**
     * This function defines how the logical addition between two discrete spaces has to performed
     * @param data The base data on which the function integrate addend
     * @param to_add The addend data
     * @param displ_x The X displacement we want addend to be placed
     * @param displ_y The Y displacement we want addend to be placed
     * @param displ_z The Z displacement we want addend to be placed
     * @param data_dim_x The base data X dimension
     * @param data_dim_y The base data Y dimension
     * @param data_dim_z The base data Z dimension
     * @param add_dim_x The addend data X dimension
     * @param add_dim_y The addend data Y dimension
     * @param add_dim_z The addend data Z dimension
     */
    static void addMeshes(MoleculeMesh::data_t *data, MoleculeMesh::data_t *to_add,
                          int displ_x, int displ_y, int displ_z,
                          int data_dim_x, int data_dim_y, int data_dim_z,
                          int add_dim_x, int add_dim_y, int add_dim_z);

    /**
     * This function defines how the logical subtraction between two discrete spaces has to performed
     * @param data The base data on which the function integrate addend
     * @param to_add The addend data
     * @param displ_x The X displacement we want addend to be placed
     * @param displ_y The Y displacement we want addend to be placed
     * @param displ_z The Z displacement we want addend to be placed
     * @param data_dim_x The base data X dimension
     * @param data_dim_y The base data Y dimension
     * @param data_dim_z The base data Z dimension
     * @param sub_dim_x The addend data X dimension
     * @param sub_dim_y The addend data Y dimension
     * @param sub_dim_z The addend data Z dimension
     */
    static void subMeshes(MoleculeMesh::data_t *data, MoleculeMesh::data_t *to_subtract,
                          int displ_x, int displ_y, int displ_z,
                          int data_dim_x, int data_dim_y, int data_dim_z,
                          int sub_dim_x, int sub_dim_y, int sub_dim_z);
};

#endif //PROLIF_COLORING_MESH
