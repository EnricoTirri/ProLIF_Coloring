//
// Created by Enrico on 29/03/2024.
//

#ifndef INTERACTION_DISTANCE_H
#define INTERACTION_DISTANCE_H

#include <Interaction.hpp>
#include <cstdlib>
#include <GraphMol/GraphMol.h>
#include <InteractionMask.hpp>
#include <vector>


class DistanceInteraction : public Interaction {
private:
    const double distance;

public:
    DistanceInteraction(const std::string &smart, double distance) : Interaction(smart), distance(distance) {};

    std::vector<InteractionMask> getInteractions(RDKit::ROMol *molecule) override {
        std::vector<InteractionMask> ris;

        //get molecule conformer and retrive matches of smart into given molecule
        RDKit::Conformer conformer = molecule->getConformer();
        RDKit::MatchVectType *match = Interaction::findMatch(molecule);

        if(!match->empty()) {
            //get molecule match and his position
            auto atomId = match->at(0).second;
            auto pos = conformer.getAtomPos(atomId);

            //calculate mask size and centering coordinates
            auto maskCenter = static_cast<int>(ceil(distance));
            auto maskDim = maskCenter * 2 + 1;
            pos.x -= maskCenter;
            pos.y -= maskCenter;
            pos.z -= maskCenter;

            //build mask
            InteractionMask mask(maskDim, maskDim, maskDim, pos);

            double ds = distance * distance;
            for (int z = 0; z < maskDim; ++z) {
                int dz = z - maskCenter;
                int z_res = dz * dz;
                for (int y = 0; y < maskDim; ++y) {
                    int dy = y - maskCenter;
                    int y_res = dy * dy;
                    for (int x = 0; x < maskDim; ++x) {
                        int dx = x - maskCenter;
                        int x_res = dx * dx;
                        if (x_res + y_res + z_res <= ds)
                            mask.at(x, y, z) = true;
                    }
                }
            }

            ris.push_back(mask);
        }

        return ris;
    }
};

#endif //INTERACTION_DISTANCE_H
