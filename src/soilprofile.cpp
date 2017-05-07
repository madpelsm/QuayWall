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
            if (waterHeight >= lowbounds) {
                effectivePressure +=
                    mSoillayers[i].mDryWeight * (lowbounds - upbounds);
            } else if (waterHeight <= upbounds) {
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

ForceVector Soilprofile::getQa(double lambda_a, Soillayer &soil, double upper,
                               double lower, double depth, double alpha,
                               double epsilon) {
    // levert een kracht Qa voor een horizontaal maaiveld (epsilon =0)
    double gamma = 0;
    double Qa_1 = 0;
    double Qa_2 = 0;
    double Qa = 0;
    double Qa_1PoE = 0, Qa_2PoE = 0;
    double resultingPoE = 0;
    if (waterHeight < upper || waterHeight > lower) {
        if (waterHeight > lower) {
            gamma = soil.mDryWeight;
        } else {
            gamma = soil.mEffectiveWeight;
        }
        Qa_1 = 0.5 * lambda_a *
               (gamma * (lower - upper) * (lower - upper) +
                2 * getEffectiveSoilePressure(depth) * cos(alpha) /
                    (cos(alpha - epsilon)) * ((lower - upper)));
        Qa_1PoE =
            getPoE(gamma, (lower - upper), getEffectiveSoilePressure(depth));

    } else {
        Qa_1 =
            0.5 * lambda_a *
            (soil.mDryWeight * (waterHeight - upper) * (waterHeight - upper) +
             2 * getEffectiveSoilePressure(depth) * (waterHeight - upper));
        Qa_2 = 0.5 * lambda_a *
               (soil.mEffectiveWeight * (lower - waterHeight) *
                    (lower - waterHeight) +
                2 * getEffectiveSoilePressure(depth + (waterHeight - upper)));
        Qa_1PoE = getPoE(soil.mDryWeight, (waterHeight - upper),
                         getEffectiveSoilePressure(depth));
        Qa_2PoE =
            getPoE(soil.mEffectiveWeight, (lower - waterHeight),
                   getEffectiveSoilePressure(depth + (waterHeight - upper)));
    }
    resultingPoE = (Qa_1PoE * Qa_1 + Qa_2PoE * Qa_2) / (Qa_1 + Qa_2);
    Qa = Qa_1 + Qa_2;
    ForceVector Qactief =
        ForceVector(glm::vec2(Qa, 0), glm::vec2(0, resultingPoE));
    // correction to be made when putting it on the construction, this is the
    // magnitude in x,positive, with the right PoE_y, but not the right
    // direction also the right x is dependent on the surface
    return Qactief;
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
    // want horizontaal maaiveld, dus d=h laagdikte, y willen we weten in
    // yCoords tov de bovenkant van de beschouwde grondlaag.
    // want yt = y=cos(alpha), maar we willen yt*cos(alpha) zie cursus
    return (1.0 / 3.0 * gamma * thickness + 0.5 * q_0) /
           (0.5 * gamma * thickness + q_0) * thickness;
}

Soilprofile::~Soilprofile() {}
