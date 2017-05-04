#include <soillayer.h>
#include <iostream>
Soillayer::Soillayer(double mUpperbounds, double mLowerBounds,
                     double mDryWeight, double mPorosity, double mPhi)
    : mUpperbounds(mUpperbounds),
      mLowerBounds(mLowerBounds),
      mDryWeight(mDryWeight),
      mPorosity(mPorosity),
      mphi(mPhi) {
    calculateSafetyValues();
    calculateWeights();
}

void Soillayer::calculateWeights() {
    mWetWeight = mDryWeight + mPorosity * 9.81;
    mEffectiveWeight = mWetWeight - 9.81;
}

void Soillayer::calculateSafetyValues() {
    mPhiA = 180.0 / M_PI * atan(tan(mphi * M_PI / 180.0) / 1.25);
    mPhiB = mphi;
    mPhiC = mPhiA;
}
