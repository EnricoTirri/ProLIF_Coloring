#include <iostream>
#include <filesystem>
#include "GraphMol/FileParsers/FileParsers.h"
#include "Mesh.hpp"
#include "Transformer.hpp"
#include "InteractionCollection.hpp"

int main(int argc, char *argv[]) {
    /* Get molecule file path */
    if (argc != 2) {
        std::cout << "Usage:\tProLIF_coloring <molecule_path>" << std::endl;
        return 1;
    }
    std::string molPath = argv[1];

    timespec startTime, endTime;

    /* Setup directory for output files */
    std::filesystem::create_directory("./outs/");

    /* Read molecule file and generate molecule mesh */
    RDKit::ROMol *molecule = RDKit::PDBFileToMol(molPath, true, false);

    std::cout << "Discretizing molecule" << std::endl;
    clock_gettime(CLOCK_MONOTONIC, &startTime);
    MoleculeMesh *moleculeMesh = Transformer::discretize(*molecule, 5);
    clock_gettime(CLOCK_MONOTONIC, &endTime);

    auto elapsed = static_cast<double>((endTime.tv_sec - startTime.tv_sec));
    elapsed += static_cast<double>((endTime.tv_nsec - startTime.tv_nsec)) / 1000000000.0;
    std::cout << "\t-> elapsed time : " << elapsed << std::endl;


    /* Generate discrete molecule from molecule-mesh */
    RDKit::RWMol *discrMolecule = Transformer::sintetize(*moleculeMesh);

    /* Save discrete molecule */
    std::string discrMoleculePath = "./outs/Molecule.pdb";

    std::cout << "\t-> saving discrete molecule file -> ";
    RDKit::MolToPDBFile(*discrMolecule, discrMoleculePath);
    std::cout << discrMoleculePath << std::endl;

    /* Retrive interaction list */
    std::vector<std::pair<std::string, Interaction *>> interactions = InteractionCollection::buildList();

    /* Iterate over interaction list */
    for (const std::pair<std::string, Interaction *> &interaction: interactions) {
        Interaction *inter = interaction.second;
        std::string desc = interaction.first;

        /* Generate a support-mesh for interaction as large as molecule one */
        MoleculeMesh interactionMesh(moleculeMesh->dim_x, moleculeMesh->dim_y, moleculeMesh->dim_z,
                                     moleculeMesh->globalDisplacement, moleculeMesh->internalDisplacement);

        std::cout << "Calculating interaction: " << desc << std::endl;
        clock_gettime(CLOCK_MONOTONIC, &startTime);
        bool succeed = inter->getInteraction(molecule, interactionMesh, *moleculeMesh);
        clock_gettime(CLOCK_MONOTONIC, &endTime);

        /* If interaction mesh generation has succeeded */
        if (succeed) {
            auto elapsed = static_cast<double>((endTime.tv_sec - startTime.tv_sec));
            elapsed += static_cast<double>((endTime.tv_nsec - startTime.tv_nsec)) / 1000000000.0;
            std::cout << "\t-> elapsed time : " << elapsed << std::endl;

            /* Generate discrete molecule from interaction-mesh */
            RDKit::RWMol *discrInteraction = Transformer::sintetize(interactionMesh);

            /* Save discrete interaction */
            std::string interactionPath = "./outs/" + desc + ".pdb";

            std::cout << "\t-> saving discrete interaction file -> ";
            RDKit::MolToPDBFile(*discrInteraction, interactionPath);
            std::cout << interactionPath << std::endl;

        }else{
            std::cout << "\t-> no interaction found" << std::endl;
        }
    }

    return EXIT_SUCCESS;
}