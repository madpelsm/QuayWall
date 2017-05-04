#include <Lmuur.h>

Lmuur::Lmuur(double mHm, double mHv, double mBl, double mBm, double mBr,
             double gamma, double toe)
    : mHm(mHm),
      mHv(mHv),
      mBl(mBl),
      mBm(mBm),
      mBr(mBr),
      gamma(gamma),
      mtoe(toe) {
    calculateProperties();
}

void Lmuur::calculateProperties() {
    // calculate foot width
    mBz = mBl + mBm + mBr;
    // calculate areas of the two rectangular parts
    mAI = mHm * mBm;
    mAII = mHv * (mBz);
    mAtoe = mtoe * mtoe;
    // calculate points of gravity
    mxI = mBm / 2.0;
    myI = mHm / 2.0;
    mxII = mBz / 2.0 - mBl;
    myII = mHv / 2.0 + mHm;
    mxtoe = -mBl + mtoe / 2.0;
    mytoe = mHm + mHv + mtoe / 2.0;
    // general p.o.g.
    mA = mAI + mAII + mAtoe;
    mx = (mxI*mAI + mxII*mAII + mxtoe*mAtoe) / mA;
    my = (myI*mAI + myII*mAII + mytoe*mAtoe) / mA;
    // calculate weight force
    mForce = gamma * mA;
}
