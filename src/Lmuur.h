#pragma once
#include <forcevector.h>
#include <vector>
#include <soilprofile.h>
#include <algorithm>

class Lmuur {
   public:
    double mHm, mHv, mBl, mBm, mBr, mBz, mxI, mxII, myI, myII, mx, my, mAI,
        mAII, gamma, mA, mtoe, mxtoe, mytoe, mAtoe;
    double mForce;
    ForceVector mForceVector;
    std::vector<ForceVector> activeSoilPressure;
    Soilprofile rightProfile;
    Soilprofile leftProfile;

    Lmuur(double mHm, double mHv, double mBl, double mBm, double mBr,
          double gamma, double toe = 0);
    ~Lmuur();
    void calculateProperties();
    void calculateActiveSoilPressures();
};
