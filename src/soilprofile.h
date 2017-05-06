#pragma once
#include <math.h>
#include <soillayer.h>
#include <cmath>
#include <vector>
class Soilprofile {
   public:
    std::vector<Soillayer> mSoillayers;
    double waterHeight = 0;

    Soilprofile(std::vector<Soillayer> Soillayers, double waterHeight = 0);
    Soilprofile();
    ~Soilprofile();

    void addSoilLayer(Soillayer soillayer);
    void replaceSoilprofile(std::vector<Soillayer> soillayers);
    double getEffectiveSoilePressure(double depth);
    double getQa(double lambda_a, Soillayer& soil, double upper, double lower,
                 double depth, double alpha = 0, double epsilon = 0);
    double getPoE(double gamma, double thickness, double q_0);
    double getLambda_a(double phi, double alpha, double psi, double epsilon);
};
