//
// Created by Enrico on 29/03/2024.
//
#ifndef INTERACTION_MASK_H
#define INTERACTION_MASK_H

#include <vector>
#include <GraphMol/GraphMol.h>

class InteractionMask {
private:

    int size_x, size_y, size_z;
    int size_layer;
    std::vector<bool> mask;
    RDGeom::Point3D globalPos;

public:
    typedef typename std::tuple<int, int, int> mask_size;

    InteractionMask(int size_x, int size_y, int size_z) : size_x(size_x), size_y(size_y), size_z(size_z) {
        size_layer = size_x * size_y;
        mask = std::vector<bool>(size_z * size_y * size_x);
        globalPos = RDGeom::Point3D(0, 0, 0);
    }

    InteractionMask(int size_x, int size_y, int size_z, const RDGeom::Point3D &pos) : size_x(size_x), size_y(size_y),
                                                                                      size_z(size_z) {
        size_layer = size_x * size_y;
        mask = std::vector<bool>(size_z * size_y * size_x);
        globalPos = RDGeom::Point3D(pos.x, pos.y, pos.z);
    }

    mask_size size() {
        return {size_x, size_y, size_z};
    }

    inline std::_Bit_reference operator[](int pos) { return mask[pos]; }

    inline std::_Bit_reference at(int x, int y, int z) {
        return mask[(z * size_layer) + (y * size_y) + x];
    }

    RDGeom::Point3D &getGlobalPos() {
        return globalPos;
    }

    void setGlobalPos(const RDGeom::Point3D &newPos) {
        globalPos.x = newPos.x;
        globalPos.y = newPos.y;
        globalPos.z = newPos.z;
    }

    bool matchSize(InteractionMask &match) {
        return size() == match.size();
    }

    void operator+=(InteractionMask &mask2) {
        if (!matchSize(mask2)) return;

        for (int i = 0; i < size_z * size_layer; ++i) {
            mask[i] = mask[i] && mask2[i];
        }
    }

    void placeMask(RDKit::RWMol *molecule) {
        RDKit::Conformer &conformer = molecule->getConformer();

        for (int i = 0; i < size_z; i++) {
            for (int j = 0; j < size_y; j++) {
                for (int k = 0; k < size_x; k++) {
                    if (at(k, j, i)) {
                        auto px = static_cast<double>(k);
                        auto py = static_cast<double>(j);
                        auto pz = static_cast<double>(i);
                        unsigned int autoId = molecule->addAtom();
                        auto *pos = new RDGeom::Point3D(px + globalPos.x, py + globalPos.y, pz + globalPos.z);
                        conformer.setAtomPos(autoId, *pos);
                    }
                }
            }
        }
    }
};

#endif //INTERACTION_MASK_H
