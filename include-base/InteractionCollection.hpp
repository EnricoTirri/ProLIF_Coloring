#ifndef PROLIF_COLORING_INTERACTION_COLLECTION
#define PROLIF_COLORING_INTERACTION_COLLECTION

#include <vector>
#include <string>
#include "Interaction.hpp"

/**
 * This abstract class defines the method necessary to build a collection of Interaction types
 */
class InteractionCollection {
public:
    /**
     * This function defines a list of Interaction type and return a list-map that associate each interaction to an id
     * @return A list-map: Interaction-ID <--> Interaction type
     */
    static std::vector<std::pair<std::string, Interaction *>> buildList();
};

#endif //PROLIF_COLORING_INTERACTION_COLLECTION
