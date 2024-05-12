#ifndef INTERACTION_METALINTERACTION_HPP
#define INTERACTION_METALINTERACTION_HPP

#include <DistanceInteraction.hpp>

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

#endif //INTERACTION_METALINTERACTION_HPP
