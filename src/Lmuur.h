#pragma once
#include <forcevector.h>
#include <soilprofile.h>
#include <algorithm>
#include <fstream>
#include <glm/glm.hpp>
#include <iostream>
#include <string>
#include <vector>

class Lmuur {
   public:
    // intermediate results
    double kastnerH = 0;
    // end of intermediate result storage
    double gamma_water = 9.81, mq = 30;
    double mHm, mHv, mBl, mBm, mBr, mBz, mxI, mxII, myI, myII, mx, my, mAI,
        mAII, gamma, mA, mtoe, mxtoe, mytoe, mAtoe;
    double mForce, mSoilHeightDifference;
    ForceVector mOwnWeight;  // real weight, no buoyance effect
    ForceVector mBuoyantForce;
    ForceVector mlhsWaterPressure, mrhsWaterPressure;
    ForceVector mResultingForce;
    std::vector<ForceVector> mSoilWedgeWeight;  // effective weight
    std::vector<ForceVector> mActiveSoilPressure;
    std::vector<ForceVector> mBoussinesqResultant;
    std::vector<ForceVector> mPassiveSoilPressure;
    Soilprofile rightProfile;
    Soilprofile leftProfile;

    Lmuur(double mHm, double mHv, double mBl, double mBm, double mBr,
          double gamma, double toe = 0);
    ~Lmuur();
    void calculateProperties();
    void calculateActiveSoilPressures(Soilprofile& soilprofile, double side,
                                      double footwidth);
    void calculateSoilWedgeWeight(Soilprofile& soilprofile, double base,
                                  double H, double offSet, double side);
    // side means: if it is on the left side of the wall -1, else 1
    void calculateWaterPressures();
    void calculateBoussinesqLoads();  // for general purposes a class for this
                                      // load should handle this itself, but
                                      // since it is extremely simple in our
                                      // case. it is faster to implement the
                                      // result right away.
    void calculateBuoyancy();
    void setSolidHeightDifference(double height);
    void calculateActiveSoilPressureLeft();
    void calculateResultingForce();
    void calculateAll();
    void writeToCSV();
    void addSoilprofiles(Soilprofile Right, Soilprofile Left);
    double squareSurface(double height, double width);
    double triangleSurface(double height, double base);
};
