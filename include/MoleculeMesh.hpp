//
// Created by Enrico on 14/04/2024.
//

#ifndef INTERACTION_MOLECULEVOXEL_HPP
#define INTERACTION_MOLECULEVOXEL_HPP

#include <vector>
#include "Geometry/point.h"

#define GRAIN 2

class MoleculeMesh {
private:
    typedef int data_t;
    std::vector<data_t> voxels;

public:
    const int dim_x, dim_y, dim_z;
    RDGeom::Point3D globalDisplacement;
    int internalDisplacement;

    MoleculeMesh(int p_dim_x, int p_dim_y, int p_dim_z, const RDGeom::Point3D &globalDisplacement, int internalDisplacement) :
            dim_x(p_dim_x),
            dim_y(p_dim_y),
            dim_z(p_dim_z),
            globalDisplacement(globalDisplacement),
            internalDisplacement(internalDisplacement){
        voxels = std::vector<data_t>(dim_x * dim_y * dim_z);
    }

    MoleculeMesh(int p_dim_x, int p_dim_y, int p_dim_z) :
            MoleculeMesh(p_dim_x, p_dim_y, p_dim_z, {0, 0, 0}, 0) {}


    inline data_t &at(int x, int y, int z) {
        return voxels[dim_x * (z * dim_y + y) + x];
    };


    inline void add(MoleculeMesh &addend, int displ_x, int displ_y, int displ_z) {
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

        int ex = dim_x;
        if(addend.dim_x + displ_x < dim_x){
            ex = addend.dim_x + displ_x;
        }

        int ey = dim_y;
        if(addend.dim_y + displ_y < dim_y){
            ey = addend.dim_y + displ_y;
        }

        int ez = dim_z;
        if(addend.dim_z + displ_z < dim_z){
            ez = addend.dim_z + displ_z;
        }

        /* Execute operation over operative window */
        for(int z = sz; z < ez; z++){
            int az = z - displ_z;
            for(int y = sy; y < ey; y++){
                int ay = y - displ_y;
                for(int x = sx; x < ex; x++){
                    int ax = x - displ_x;
                    at(x,y,z) = at(x,y,z) || addend.at(ax,ay,az);
                }
            }
        }
    }

    inline void sub(MoleculeMesh &addend, int displ_x, int displ_y, int displ_z) {
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

        int ex = dim_x;
        if(addend.dim_x + displ_x < dim_x){
            ex = addend.dim_x + displ_x;
        }

        int ey = dim_y;
        if(addend.dim_y + displ_y < dim_y){
            ey = addend.dim_y + displ_y;
        }

        int ez = dim_z;
        if(addend.dim_z + displ_z < dim_z){
            ez = addend.dim_z + displ_z;
        }

        /* Execute operation over operative window */
        for(int z = sz; z < ez; z++){
            int az = z - displ_z;
            for(int y = sy; y < ey; y++){
                int ay = y - displ_y;
                for(int x = sx; x < ex; x++){
                    int ax = x - displ_x;
                    at(x,y,z) = at(x,y,z) && !addend.at(ax,ay,az);
                }
            }
        }
    }
};

#endif //INTERACTION_MOLECULEVOXEL_HPP
