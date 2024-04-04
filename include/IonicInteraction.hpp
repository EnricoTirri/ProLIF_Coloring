//
// Created by Enrico on 03/04/2024.
//

#ifndef INTERACTION_IONICINTERACTION_HPP
#define INTERACTION_IONICINTERACTION_HPP

#include <DistanceInteraction.hpp>

class CationicInteraction : public DistanceInteraction{
public:
    CationicInteraction() : DistanceInteraction(
            "[-{1-},$(O=[C,S,P]-[O-])]",
            4.5
            ){}
};

class AnionicInteraction : public DistanceInteraction{
public:
    AnionicInteraction() : DistanceInteraction(
            "[+{1-},$([NX3&!$([NX3]-O)]-[C]=[NX3+])]",
            4.5
    ){}
};

#endif //INTERACTION_IONICINTERACTION_HPP
