#include <cstdlib>
#include "GraphMol/FileParsers/FileParsers.h"
#include <GraphMol/Substruct/SubstructMatch.h>
#include <InteractionMask.hpp>
#include <filesystem>
#include <InteractionCollection.hpp>

int main(int argc, char *argv[]) {
    if(argc != 2){
        std::cout << "Usage:\tProLIF_coloring <molecule_path>" << std::endl;
        return 1;
    }

    std::filesystem::create_directory("./outs/");

    std::string molPath = argv[1];
    RDKit::ROMol* molecule = RDKit::PDBFileToMol(molPath, true, false);
    RDKit::Conformer conformer = molecule->getConformer();

    std::vector<std::pair<std::string,Interaction *>> interactions = InteractionCollection::buildList();

    for(const std::pair<std::string,Interaction *>& interaction : interactions){
        auto *bubble = new RDKit::RWMol();
        auto *bconf = new RDKit::Conformer();
        bubble->addConformer(bconf);

        Interaction * inter = interaction.second;
        std::string desc = interaction.first;
        std::vector<InteractionMask> masks = inter->getInteractions(molecule);
        if(!masks.empty()) {
            masks[0].placeMask(bubble);
            RDKit::MolToPDBFile(*bubble, "./outs/" + desc + ".pdb");
        }
    }

    return EXIT_SUCCESS;
}