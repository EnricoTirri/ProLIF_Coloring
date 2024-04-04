//
// Created by Enrico on 29/03/2024.
//

#ifndef INTERACTION_ANGLE_H
#define INTERACTION_ANGLE_H

#include <cstdlib>
#include <InteractionMask.hpp>
#include <GraphMol/GraphMol.h>
#include <Interaction.hpp>

class SingleAngleInteraction : public Interaction {
private:
    const double min_angle, max_angle, distance;
public:
    SingleAngleInteraction(const std::string &smart,
                           std::pair<double, double> angle,
                           double distance
                           //TODO add distance-to point choice
                           )
                           : Interaction(smart),
                            min_angle(angle.first),
                            max_angle(angle.second),
                            distance(distance) {};

    std::vector<InteractionMask> getInteractions(RDKit::ROMol *molecule) override {
        std::vector<InteractionMask> ris;

        //get molecule conformer and retrive matches of smart into given molecule
        RDKit::Conformer conformer = molecule->getConformer();
        RDKit::MatchVectType *match = Interaction::findMatch(molecule);

        if(match->size()>=2){
            //get molecule match and his position
            auto p1Id = match->at(0).second;
            auto p1 = conformer.getAtomPos(p1Id);
            auto p2Id = match->at(1).second;
            auto p2 = conformer.getAtomPos(p2Id);


            //calculate mask size and centering coordinates
            auto maskCenter = static_cast<int>(ceil(distance));
            auto maskDim = maskCenter * 2 + 1;

            //TODO distance is calculated always to p1, this should be changed on implementation choice
            //Since distance build a sphere around the atom, mask is centered by that atom
            RDGeom::Point3D maskPos(p1.x-maskCenter, p1.y-maskCenter, p1.z-maskCenter);

            //build mask
            InteractionMask mask(maskDim, maskDim, maskDim, maskPos);

            RDGeom::Point3D p2p1 = p2.directionVector(p1);

            double ds = distance * distance;
            for (int z = 0; z < maskDim; ++z) {
                double pz = z + maskPos.z;
                int dz = z - maskCenter;
                int z_res = dz * dz;
                for (int y = 0; y < maskDim; ++y) {
                    double py = y + maskPos.y;
                    int dy = y - maskCenter;
                    int y_res = dy * dy;
                    for (int x = 0; x < maskDim; ++x) {
                        double px = x + maskPos.x;
                        int dx = x - maskCenter;
                        int x_res = dx * dx;

                        RDGeom::Point3D l1(px, py, pz);
                        RDGeom::Point3D p2l1 = p2.directionVector(l1);
                        double angle = p2p1.angleTo(p2l1);
                        if (x_res + y_res + z_res <= ds && angle >= min_angle && angle <= max_angle)
                            mask.at(x, y, z) = true;
                    }
                }
            }

            ris.push_back(mask);
        }

        return ris;
    }
};

#endif //INTERACTION_ANGLE_H
