
#include "SingleAngleInteraction.hpp"

bool SingleAngleInteraction::getInteraction(const RDKit::ROMol *molecule,
                                            MoleculeMesh &interactionMask, MoleculeMesh &subtractionMask) {

    // Get molecule conformer and retrive matches of smart into given molecule
    RDKit::Conformer conformer = molecule->getConformer();
    std::vector<RDKit::MatchVectType> *matches = Interaction::findMatch(molecule);

    if (matches->empty()) return false;


    bool found = false;
    for (RDKit::MatchVectType match: *matches) {
        if (match.size() >= 2) {
            found = true;

            // Get molecule match and its centroids position
            auto p1Id = match.at(0).second;
            auto p1 = conformer.getAtomPos(p1Id);
            auto p2Id = match.at(1).second;
            auto p2 = conformer.getAtomPos(p2Id);

            RDGeom::Point3D center;
                if (cp) center = p1;
                else center = p2;

            // Calculate mask size and centering coordinates
            auto scaledMaskCenter = static_cast<int>(ceil(distance * GRAIN));
            auto maskDim = 2 * scaledMaskCenter;

            // Generate pattern-mesh
            MoleculeMesh bubble(maskDim, maskDim, maskDim);

            // Calculate vector p2 --> p1
            RDGeom::Point3D p2p1 = p2.directionVector(p1);

            /*
             * Over all size of pattern-mesh assign if:
             *      - (point-distance <= #distance) from the center of mesh
             *      - angle between l1 <-- p2 --> p1 is (#min <= #angle <= #max)
             */
            double scaledDistance = distance * GRAIN;
            double ds = scaledDistance * scaledDistance;
            for (int z = 0; z < maskDim; ++z) {
                double pz = z + center.z;
                int dz = z - scaledMaskCenter;
                int z_res = dz * dz;
                for (int y = 0; y < maskDim; ++y) {
                    double py = y + center.y;
                    int dy = y - scaledMaskCenter;
                    int y_res = dy * dy;
                    for (int x = 0; x < maskDim; ++x) {
                        double px = x + center.x;
                        int dx = x - scaledMaskCenter;
                        int x_res = dx * dx;

                        // Position of l1 is calculated taking account of pattern-center position
                        RDGeom::Point3D l1(px, py, pz);

                        // Calculate vector p2 --> l1
                        RDGeom::Point3D p2l1 = p2.directionVector(l1);

                        // Calculate angle l1 <-- p2 --> p1
                        double angle = p2p1.angleTo(p2l1);

                        if (x_res + y_res + z_res <= ds && angle >= min_angle && angle <= max_angle)
                            bubble.at(x, y, z) = true;
                    }
                }
            }

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
            interactionMask.add(bubble, displ_x, displ_y, displ_z);
        }
    }

    interactionMask.sub(subtractionMask, 0, 0, 0);

    return found;
}