#ifndef PROLIF_COLORING_HB_INTERACTION
#define PROLIF_COLORING_HB_INTERACTION

#include "SingleAngleInteraction.hpp"

class HBDonorInteraction : public SingleAngleInteraction {
public:
    HBDonorInteraction() : SingleAngleInteraction(
            "[#7&!$([nX3])&!$([NX3]-*=[O,N,P,S])&!$([NX3]-[a])&!$([Nv4&+1]),O&!$([OX2](C)C=O)&!$(O(~a)~a)&!$(O=N-*)&!$([O-]-N=O),o+0,F&$(F-[#6])&!$(F-[#6][F,Cl,Br,I])]",
            {M_PI * 130 / 180, M_PI},
            3.5,
            0
    ) {}
};

class HBAcceptorInteraction : public SingleAngleInteraction {
public:
    HBAcceptorInteraction() : SingleAngleInteraction(
            "[$([O,S;+0]),$([N;v3,v4&+1]),n+0]-[H]",
            {M_PI * 130 / 180, M_PI},
            3.5,
            0
    ) {}
};

#endif //PROLIF_COLORING_HB_INTERACTION
