#pragma once

class Lmuur {
   public:
    double mHm, mHv, mBl, mBm, mBr, mBz, mxI, mxII, myI, myII, mx, my, mAI,
        mAII, gamma, mA, mtoe, mxtoe, mytoe, mAtoe;
    double mForce;
    Lmuur(double mHm, double mHv, double mBl, double mBm, double mBr,
          double gamma, double toe=0);
    void calculateProperties();
};
