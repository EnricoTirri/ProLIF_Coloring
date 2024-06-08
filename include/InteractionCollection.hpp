#ifndef PROLIF_COLORING_INTERACTION_COLLECTION
#define PROLIF_COLORING_INTERACTION_COLLECTION

#include <vector>
#include <string>
#include <Interaction.hpp>
#include <HydrophobicInteraction.hpp>
#include <HBInteraction.hpp>
#include <IonicInteraction.hpp>
#include <MetalInteraction.hpp>

class InteractionCollection {
public:
    static std::vector<std::pair<std::string, Interaction *>> buildList();
};

#endif //PROLIF_COLORING_INTERACTION_COLLECTION
