#include <soillayer.h>
#include <iostream>
Soillayer::Soillayer(double mUpperbounds, double mLowerBounds,
                     double mDryWeight, double mPorosity, double mPhi)
    : mUpperbounds(mUpperbounds),
      mLowerBounds(mLowerBounds),
      mDryWeight(mDryWeight),
      mPorosity(mPorosity),
      mphi(mPhi) {
    // input in degree, internal radians
    mphi = mphi * M_PI / 180.0;
    calculateSafetyValues();
    calculateWeights();
}

void Soillayer::calculateWeights() {
    mWetWeight = mDryWeight + mPorosity * 9.81;
    mEffectiveWeight = mWetWeight - 9.81;
}

void Soillayer::calculateSafetyValues() {
    // input is degree, internal is radians
    mPhiA = atan(tan(mphi) / 1.25);
    mPhiB = M_PI * mphi / 180.0;
    mPhiC = mPhiA;
}

Soillayer::~Soillayer() {}
