#ifndef PROLIF_COLORING_IONIC_INTERACTION
#define PROLIF_COLORING_IONIC_INTERACTION

#include <DistanceInteraction.hpp>

class CationicInteraction : public DistanceInteraction {
public:
    CationicInteraction() : DistanceInteraction(
            "[-{1-},$(O=[C,S,P]-[O-])]",
            4.5
    ) {}
};

class AnionicInteraction : public DistanceInteraction {
public:
    AnionicInteraction() : DistanceInteraction(
            "[+{1-},$([NX3&!$([NX3]-O)]-[C]=[NX3+])]",
            4.5
    ) {}
};

#endif //PROLIF_COLORING_IONIC_INTERACTION
