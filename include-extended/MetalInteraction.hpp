#ifndef PROLIF_COLORING_METAL_INTERACTION
#define PROLIF_COLORING_METAL_INTERACTION

#include "DistanceInteraction.hpp"

class MetalAcceptorInteraction : public DistanceInteraction {
public:
    MetalAcceptorInteraction() : DistanceInteraction(
            "[Ca,Cd,Co,Cu,Fe,Mg,Mn,Ni,Zn]",
            2.8
    ) {}
};

class MetalDonorInteraction : public DistanceInteraction {
public:
    MetalDonorInteraction() : DistanceInteraction(
            "[O,#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4]),-{1-};!+{1-}]",
            2.8
    ) {}
};

#endif //PROLIF_COLORING_METAL_INTERACTION
