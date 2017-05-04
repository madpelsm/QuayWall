#include <soilprofile.h>

Soilprofile::Soilprofile(std::vector<Soillayer> Soillayers)
    : mSoillayers(Soillayers) {  // added soillayers
}

void Soilprofile::addSoilLayer(Soillayer soillayer) {
    mSoillayers.push_back(soillayer);
}
void Soilprofile::replaceSoilprofile(std::vector<Soillayer> soillayers) {
    mSoillayers = soillayers;
}

double Soilprofile::getEffectiveSoilePressure(double depth) {
    double effectivePressure = 0;
    double lowbounds = 0, upbounds = 0;
    bool tooDeep = false;
    for (std::size_t i = 0; i < mSoillayers.size() && !tooDeep; ++i) {
        lowbounds = mSoillayers[i].mLowerBounds;
        upbounds = mSoillayers[i].mUpperbounds;
        if (upbounds > depth) tooDeep = true;

        if (mSoillayers[i].mLowerBounds > depth && !tooDeep) {
            // set the lowbounds to this depth
            lowbounds = depth;
        }
        if (!tooDeep) {
            if (waterHeight > lowbounds) {
                effectivePressure +=
                    mSoillayers[i].mDryWeight * (lowbounds - upbounds);
            } else if (waterHeight < upbounds) {
                effectivePressure +=
                    (lowbounds - upbounds) * mSoillayers[i].mEffectiveWeight;
            } else if (waterHeight < lowbounds && waterHeight > upbounds) {
                effectivePressure +=
                    (lowbounds - waterHeight) * mSoillayers[i].mEffectiveWeight;
                effectivePressure +=
                    (waterHeight - upbounds) * mSoillayers[i].mDryWeight;
            }
        }
    }
    return effectivePressure;
}
