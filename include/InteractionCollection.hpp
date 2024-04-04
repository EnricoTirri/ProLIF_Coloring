//
// Created by Enrico on 04/04/2024.
//

#ifndef INTERACTION_INTERACTIONCOLLECTION_HPP
#define INTERACTION_INTERACTIONCOLLECTION_HPP

#include <vector>
#include <string>
#include <Interaction.hpp>
#include <HydrophobicInteraction.hpp>
#include <HBInteraction.hpp>
#include <IonicInteraction.hpp>
#include <MetalInteraction.hpp>

class InteractionCollection{
public:
    static std::vector<std::pair<std::string, Interaction *>> buildList();
};
#endif //INTERACTION_INTERACTIONCOLLECTION_HPP
