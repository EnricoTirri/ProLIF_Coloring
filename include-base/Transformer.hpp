#ifndef PROLIF_COLORING_DISCRETIZER
#define PROLIF_COLORING_DISCRETIZER

#include "GraphMol/RWMol.h"
#include "Mesh.hpp"

/**
 * This class allow the transformation from the continuous space of RDKit-molecule
 * definition to the discrete space of MoleculeMesh class and vice-versa
 */
class Transformer {

    /**
     * Standard atom radius used in transformation
     */
    static constexpr double atomRadius = 1.1;
    /**
     * Minimum padding applied to MoleculeMesh boundary
     */
    static constexpr int minPadding = 2;

public:
    /**
     * This function allow the transformation from the RDKit-molecule to MoleculeMesh
     * @param molecule The RDKit-molecule to get discrete definition
     * @param padding The padding to add to discrete definition
     * @return The discrete definition of the input molecule
     */
    static MoleculeMesh *discretize(const RDKit::ROMol &molecule, int padding = minPadding) {

        /* Padding less than the minimum one is not allowed */
        if (padding < minPadding) padding = minPadding;

        /* Retrieve all atom position of input molecule */
        std::vector<RDGeom::Point3D> atoms = molecule.getConformer().getPositions();

        /* Find min and max of the span of molecule */
        double t_max_x = 0, t_max_y = 0, t_max_z = 0, t_min_x = atoms[0].x, t_min_y = atoms[0].y, t_min_z = atoms[0].z;

        for (const auto &pos: atoms) {
            if (t_max_x < pos.x)
                t_max_x = pos.x;
            if (t_max_y < pos.y)
                t_max_y = pos.y;
            if (t_max_z < pos.z)
                t_max_z = pos.z;

            if (t_min_x > pos.x)
                t_min_x = pos.x;
            if (t_min_y > pos.y)
                t_min_y = pos.y;
            if (t_min_z > pos.z)
                t_min_z = pos.z;
        }

        /* Discretize min and max */
        int max_x = static_cast<int>(ceil(t_max_x) * GRAIN);
        int max_y = static_cast<int>(ceil(t_max_y) * GRAIN);
        int max_z = static_cast<int>(ceil(t_max_z) * GRAIN);

        int min_x = static_cast<int>(floor(t_min_x) * GRAIN);
        int min_y = static_cast<int>(floor(t_min_y) * GRAIN);
        int min_z = static_cast<int>(floor(t_min_z) * GRAIN);

        /* Calculate span of molecule */
        int span_x = max_x - min_x;
        int span_y = max_y - min_y;
        int span_z = max_z - min_z;

        /* Size of mesh is molecule_span + border_padding, all scaled to a granularity factor */
        int scaledPadding = padding * GRAIN;

        int size_x = span_x + scaledPadding * 2;
        int size_y = span_y + scaledPadding * 2;
        int size_z = span_z + scaledPadding * 2;

        /* Generate support-mesh */
        auto mesh = new MoleculeMesh(size_x, size_y, size_z,
                                     {floor(t_min_x),
                                      floor(t_min_y),
                                      floor(t_min_z)},
                                     scaledPadding);

        /* Calculate discrete atom size */
        double scaledAtomRadius = atomRadius * GRAIN;
        int scaledAtomPadding = static_cast<int>(ceil(scaledAtomRadius));

        /* For each atom of molecule */
        for (auto &pos: atoms) {
            /* Calculate atom position on support-mesh reference system */
            double px = pos.x * GRAIN - min_x + scaledPadding;
            double py = pos.y * GRAIN - min_y + scaledPadding;
            double pz = pos.z * GRAIN - min_z + scaledPadding;

            /* Calculate operative ranges of atom */
            std::pair<int, int> range_x = {
                    static_cast<int>(floor(px)) - scaledAtomPadding,
                    static_cast<int>(ceil(px)) + scaledAtomPadding
            };
            std::pair<int, int> range_y = {
                    static_cast<int>(floor(py)) - scaledAtomPadding,
                    static_cast<int>(ceil(py)) + scaledAtomPadding
            };
            std::pair<int, int> range_z = {
                    static_cast<int>(floor(pz)) - scaledAtomPadding,
                    static_cast<int>(ceil(pz)) + scaledAtomPadding
            };

            /* Over operative ranges */
            double ds = scaledAtomRadius * scaledAtomRadius;
            for (int z = range_z.first; z < range_z.second; ++z) {
                double dz = z - pz;
                double z_res = dz * dz;
                for (int y = range_y.first; y < range_y.second; ++y) {
                    double dy = y - py;
                    double y_res = dy * dy;
                    for (int x = range_x.first; x < range_x.second; ++x) {
                        double dx = x - px;
                        double x_res = dx * dx;
                        /* Find if the point (x,y,z) has distance <= #atomRadius from atom-position (#pos) */
                        if (x_res + y_res + z_res <= ds)
                            mesh->at(x, y, z) = 1;
                    }
                }
            }
        }

        return mesh;
    }


    /**
     * This function allow the transformation from the MoleculeMesh to RDKit-molecule
     * @param mesh The RDKit-molecule to get continuous definition
     * @return The continuous definition of the input molecule
     */
    static RDKit::RWMol *sintetize(MoleculeMesh &mesh) {
        /* Generate a new molecule and assign a conformer for atoms position */
        auto *molecule = new RDKit::RWMol();
        auto *conformer = new RDKit::Conformer();
        molecule->addConformer(conformer);

        /* Over the entire size of mesh */
        for (int i = 0; i < mesh.dim_z; i++) {
            auto pz = static_cast<double>(i - mesh.internalDisplacement) / GRAIN + mesh.globalDisplacement.z;
            for (int j = 0; j < mesh.dim_y; j++) {
                auto py = static_cast<double>(j - mesh.internalDisplacement) / GRAIN + mesh.globalDisplacement.y;
                for (int k = 0; k < mesh.dim_x; k++) {
                    if (mesh.at(k, j, i)) {
                        auto px =
                                static_cast<double>(k - mesh.internalDisplacement) / GRAIN + mesh.globalDisplacement.x;

                        /* Get atom position back from mesh reference system */
                        auto *pos = new RDGeom::Point3D(px, py, pz);

                        /* Add atom to the molecule and assign the calculated position */
                        unsigned int autoId = molecule->addAtom();
                        molecule->getAtomWithIdx(autoId)->setAtomicNum(1);
                        conformer->setAtomPos(autoId, *pos);
                    }
                }
            }
        }
        return molecule;
    }
};

#endif //PROLIF_COLORING_DISCRETIZER
