
#include <iostream>
#include "GraphMol/FileParsers/FileParsers.h"
#include "Mesh.hpp"
#include "Discretizer.hpp"
#include "InteractionCollection.hpp"

#if USEOMP
#include "omp.h"
#endif

#define ITERATIONS 10

int main(int argc, char *argv[]) {
#if USEOMP
    omp_set_num_threads(OMP_NUM_THREADS);
#endif

    /* Get molecule file path */
    if (argc != 2) {
        std::cout << "#Usage:\tProLIF_coloring <molecule_path>" << std::endl;
        return 1;
    }
    std::string molPath = argv[1];

    timespec startTime, endTime;

    /* Read molecule file and generate molecule mesh */
    RDKit::ROMol *molecule = RDKit::PDBFileToMol(molPath, true, false);

    std::cout << "#Discretizing molecule" << std::endl;
    MoleculeMesh *moleculeMesh = Discretizer::discretize(*molecule, 5);
    std::cout << "#-----------------------------------" << std::endl;

    Interaction *testInteraction = new HydrophobicInteraction();

    std::string method = "Serial";
    std::string type = "Compute";

    if(USECUDA) {
        std::cout << "#TEST ON CUDA IMPLEMENTATION" << std::endl;
        method = "Cuda";
    } else if(USEOMP) {
        std::cout << "#TEST ON OMP IMPLEMENTATION" << std::endl;
        method = "OpenMP";
    } else{
        std::cout << "#TEST ON SERIAL IMPLEMENTATION" << std::endl;
    }
    if(USEPATTERN){
        std::cout << "#Method type: Pattern placing" << std::endl;
        type = "Pattern";
    }else{
        std::cout << "#METHOD TYPE: Pure compute" << std::endl;
    }

    std::cout << "#MESH GRAINING: " << GRAIN << std::endl;
    std::cout << "#ITERATIONS: " << ITERATIONS << std::endl;
    std::cout << "#-----------------------------------" << std::endl;

    for(int i=0; i< ITERATIONS; ++i) {
        /* Generate a support-mesh for interaction as large as molecule one */
        MoleculeMesh interactionMesh(moleculeMesh->dim_x, moleculeMesh->dim_y, moleculeMesh->dim_z,
                                     moleculeMesh->globalDisplacement, moleculeMesh->internalDisplacement);

        std::cout << "#Iteration: " << i << std::endl;
        clock_gettime(CLOCK_MONOTONIC, &startTime);
        bool succeed = testInteraction->getInteraction(molecule, interactionMesh, *moleculeMesh);
        clock_gettime(CLOCK_MONOTONIC, &endTime);

        /* If interaction mesh generation has succeeded */
        if (succeed) {
            auto elapsed = static_cast<double>((endTime.tv_sec - startTime.tv_sec));
            elapsed += static_cast<double>((endTime.tv_nsec - startTime.tv_nsec)) / 1000000000.0;
            std::cout << "#Elapsed: " << elapsed << std::endl;
        } else {
            std::cout << "#\t-> no interaction found" << std::endl;
            break;
        }
        std::cout << "#-----------------------------------" << std::endl;
    }


    return EXIT_SUCCESS;
}