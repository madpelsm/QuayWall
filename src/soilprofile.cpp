#include <soilprofile.h>
#include <iostream>

Soilprofile::Soilprofile(std::vector<Soillayer> Soillayers, double waterHeight)
    : mSoillayers(Soillayers), waterHeight(waterHeight) {  // added soillayers
}
Soilprofile::Soilprofile() {}

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
            } else if (waterHeight <= lowbounds && waterHeight >= upbounds) {
                effectivePressure +=
                    (lowbounds - waterHeight) * mSoillayers[i].mEffectiveWeight;
                effectivePressure +=
                    (waterHeight - upbounds) * mSoillayers[i].mDryWeight;
            }
        }
    }
    return effectivePressure;
}

double Soilprofile::getQa(double lambda_a, Soillayer &soil, double upper,
                          double lower, double depth, double alpha,
                          double epsilon) {
    double gamma = 0;
    double Qa_1 = 0;
    double Qa_2 = 0;
    if (waterHeight < upper || waterHeight > lower) {
        if (waterHeight > lower) {
            gamma = soil.mDryWeight;
        } else {
            gamma = soil.mEffectiveWeight;
        }
        Qa_1 = 0.5 * lambda_a *
               (gamma * (lower - upper) * (lower - upper) +
                2 * getEffectiveSoilePressure(depth) * cos(alpha) /
                    (cos(alpha - epsilon)) * (lower - upper));

    } else {
        Qa_1 =
            0.5 * lambda_a *
            (soil.mDryWeight * (waterHeight - upper) * (waterHeight - upper) +
             2 * getEffectiveSoilePressure(depth) * (waterHeight - upper));
        Qa_2 = 0.5 * lambda_a *
               (soil.mEffectiveWeight * (lower - waterHeight) *
                    (lower - waterHeight) +
                2 * getEffectiveSoilePressure(depth + (waterHeight - upper)));
    }
    return Qa_1;
}

double Soilprofile::getLambda_a(double phi, double alpha, double psi,
                                double epsilon) {
    return cos(phi - alpha) * cos(phi - alpha) /
           (cos(alpha) * cos(alpha) * cos(alpha + psi) *
            pow((1 + std::sqrt(sin(phi + psi) * sin(phi - epsilon) /
                               (cos(alpha + psi) * cos(epsilon - alpha)))),
                2));
}

double Soilprofile::getPoE(double gamma, double thickness, double q_0) {
    return (1.0 / 3.0 * gamma * thickness + 0.5 * q_0) /
           (0.5 * gamma * thickness + q_0) * thickness;
}

Soilprofile::~Soilprofile() {}
