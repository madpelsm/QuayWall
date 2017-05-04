#pragma once
#include <soillayer.h>
#include <vector>
class Soilprofile {
   public:
    std::vector<Soillayer> mSoillayers;
    double waterHeight;

    Soilprofile(std::vector<Soillayer> Soillayers);

    void addSoilLayer(Soillayer soillayer);
    void replaceSoilprofile(std::vector<Soillayer> soillayers);
    double getEffectiveSoilePressure(double depth);
};
