#ifndef PROLIF_COLORING_HYDRO_INTERACTION
#define PROLIF_COLORING_HYDRO_INTERACTION

#include <DistanceInteraction.hpp>

class HydrophobicInteraction : public DistanceInteraction {
public:
    HydrophobicInteraction() : DistanceInteraction(
            "[c,s,Br,I,S&H0&v2,$([D3,D4;#6])&!$([#6]~[#7,#8,#9])&!$([#6X4H0]);+0]", 4.5) {}
};

#endif //PROLIF_COLORING_HYDRO_INTERACTION
