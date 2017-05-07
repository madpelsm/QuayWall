#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
// weights are in kN/m^3

class Soillayer {
   public:
    double mUpperbounds, mLowerBounds, mDryWeight, mPorosity, mphi;
    double mPhiA, mPhiB, mPhiC;
    double mEffectiveWeight, mWetWeight;
    // calculating safetyweight is less cluttery than the phi's
    // input in degrees, internal format is radians
    Soillayer(double mUpperbounds, double mLowerBounds, double mDryWeight,
              double mPorosity, double mPhi);
    ~Soillayer();

    void calculateSafetyValues();
    void calculateWeights();
};
