#include <iostream>
#include <filesystem>
#include "GraphMol/FileParsers/FileParsers.h"
#include "MoleculeMesh.hpp"
#include "Discretizer.hpp"
#include "InteractionCollection.hpp"

int main(int argc, char *argv[]) {
    /* Get molecule file path */
    if (argc != 2) {
        std::cout << "Usage:\tProLIF_coloring <molecule_path>" << std::endl;
        return 1;
    }
    std::string molPath = argv[1];

    /* Setup directory for output files */
    std::filesystem::create_directory("./outs/");

    /* Read molecule file and generate molecule mesh */
    RDKit::ROMol *molecule = RDKit::PDBFileToMol(molPath, true, false);

    MoleculeMesh *moleculeMesh = Discretizer::discretize(*molecule, 5);

    /* Generate discrete molecule from molecule-mesh */
    RDKit::RWMol *discrMolecule = Discretizer::sintetize(*moleculeMesh);

    /* Save discrete molecule */
    std::string discrMoleculePath = "./outs/Molecule.pdb";
    RDKit::MolToPDBFile(*discrMolecule, discrMoleculePath);
    std::cout << "Generated discrete molecule file -> " << discrMoleculePath << std::endl;

    /* Retrive interaction list */
    std::vector<std::pair<std::string, Interaction *>> interactions = InteractionCollection::buildList();

    /* Iterate over interaction list */
    for (const std::pair<std::string, Interaction *> &interaction: interactions) {
        Interaction *inter = interaction.second;
        std::string desc = interaction.first;

        /* Generate a support-mesh for interaction as large as molecule one */
        MoleculeMesh interactionMesh(moleculeMesh->dim_x, moleculeMesh->dim_y, moleculeMesh->dim_z,
                                     moleculeMesh->globalDisplacement, moleculeMesh->internalDisplacement);

        /* If interaction mesh generation has succeeded */
        if (inter->getInteraction(molecule, interactionMesh)) {
            /* Subtract molecule-mesh from interaction-mesh */
            //interactionMesh.sub(*moleculeMesh, 0, 0, 0);

            /* Generate discrete molecule from interaction-mesh */
            RDKit::RWMol *discrInteraction = Discretizer::sintetize(interactionMesh);

            /* Save discrete interaction */
            std::string interactionPath = "./outs/" + desc + ".pdb";
            RDKit::MolToPDBFile(*discrInteraction, interactionPath);
            std::cout << "Generated discrete interaction file -> " << interactionPath << std::endl;
        }
    }

    return EXIT_SUCCESS;
}