#pragma once
#include <soillayer.h>
#include <vector>
class Soilprofile {
   public:
    std::vector<Soillayer> mSoillayers;
    double waterHeight=0;

    Soilprofile(std::vector<Soillayer> Soillayers,double waterHeight =0);

    void addSoilLayer(Soillayer soillayer);
    void replaceSoilprofile(std::vector<Soillayer> soillayers);
    double getEffectiveSoilePressure(double depth);
};
