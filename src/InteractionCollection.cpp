//
// Created by Enrico on 04/04/2024.
//
#include <InteractionCollection.hpp>

std::vector<std::pair<std::string, Interaction*>> InteractionCollection::buildList(){
    std::vector<std::pair<std::string, Interaction*>> interactionsList;
    interactionsList.push_back({"Hydrophobic", new HydrophobicInteraction()});
    interactionsList.push_back({"HBAcceptor", new HBAcceptorInteraction()});
    interactionsList.push_back({"HBDonor", new HBDonorInteraction()});
    interactionsList.push_back({"Cationic", new AnionicInteraction()});
    interactionsList.push_back({"Anionic", new CationicInteraction()});
    interactionsList.push_back({"MetalAcceptor", new MetalAcceptorInteraction()});
    interactionsList.push_back({"MetalDonor", new MetalDonorInteraction()});
    return interactionsList;
}