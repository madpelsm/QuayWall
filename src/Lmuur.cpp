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
    mx = (mxI * mAI + mxII * mAII + mxtoe * mAtoe) / mA;
    my = (myI * mAI + myII * mAII + mytoe * mAtoe) / mA;
    // calculate weight force
    mForce = gamma * mA;
    mOwnWeight = ForceVector(glm::vec2(0, mForce), glm::vec2(mx, my));
}

void Lmuur::calculateActiveSoilPressures(Soilprofile& soilprofile, double side,
                                         double footwidth) {
    double correctionHeight = 0;
    double forceDirectionCorrection = 1;
    double wallBorder = mBm;
    if (side != 1) {
        forceDirectionCorrection = -1;
        correctionHeight = mSoilHeightDifference;
        wallBorder = 0;
    }

    double phi_d = M_PI / 2.0;
    for (size_t i = 0; i < soilprofile.mSoillayers.size(); ++i) {
        phi_d = std::min(phi_d, soilprofile.mSoillayers[i].mPhiA);
    }

    ForceVector Qa;
    double psi = 0, phi = 0, alpha = 0;
    double yTemp =
        (mHm - correctionHeight) - footwidth * tan(M_PI / 4.0 + phi_d);
    for (size_t i = 0; i < soilprofile.mSoillayers.size(); ++i) {
        // dependencies for lambda_a
        psi = phi_d * (2.0 / 3.0);
        alpha = M_PI / 4.0 - phi / 2.0;

        // dependencies for Qa

        // h = diepte tov maaiveld (pos = naar beneden)
        double h = soilprofile.mSoillayers[i].mUpperbounds;
        double lambda_a = 0;
        double lower = std::min(soilprofile.mSoillayers[i].mLowerBounds,
                                mHm - correctionHeight);
        if (h < yTemp) {
            // yTemp zal dus sowieso positief zijn, want h is positief
            if (lower < yTemp) {
                // glijlijn is nog niet bereikt
                alpha = 0;
                // helling van oppervlak is helling van muur
                // (alpha=0,psi=2/3phi)
                //
                lambda_a = soilprofile.getLambda_a(phi_d, 0, psi, 0);
                Qa = soilprofile.getQa(lambda_a, soilprofile.mSoillayers[i], h,
                                       lower, h, alpha);
                // correctie voor de juiste richting van de vector te bekomen
                // magnitude zit in x.
                Qa.mForce.y = Qa.mForce.x * sin(alpha + psi);
                Qa.mForce.x =
                    -forceDirectionCorrection * Qa.mForce.x * cos(alpha + psi);
                // correctie voor aangrijpingspunt. Assenstelsel linker
                // verticale aan de muur is de y-as. indien de kracht aangrijpt
                // rechts van de wand is de xCo dus mBm.
                // de yCoordinaat werd berekent relatief tov de bovenkant van de
                // grondlaag
                Qa.mPoE.x = wallBorder;
                Qa.mPoE.y += h;
                mActiveSoilPressure.push_back(Qa);
            } else {
                // glijlijn ligt boven de onderkant van de beschouwde
                // grondlaag
                // bovenkant van de grondlaag ligt nog bove de glijlijn
                lower = yTemp;
                alpha = 0;
                lambda_a = soilprofile.getLambda_a(phi_d, alpha, psi, 0);
                Qa = soilprofile.getQa(lambda_a, soilprofile.mSoillayers[i], h,
                                       lower, h, alpha);
                Qa.mForce.y = Qa.mForce.x * sin(alpha + psi);
                Qa.mForce.x =
                    -forceDirectionCorrection * Qa.mForce.x * cos(alpha + psi);
                Qa.mPoE.x = wallBorder;
                Qa.mPoE.y += h;
                mActiveSoilPressure.push_back(Qa);
                // Deel van de grond op grond actie
                lower = std::min(soilprofile.mSoillayers[i].mLowerBounds,
                                 mHm - correctionHeight);
                alpha = M_PI / 4.0 - phi_d / 2.0;
                lambda_a = soilprofile.getLambda_a(phi_d, alpha, phi_d, 0);
                Qa = soilprofile.getQa(lambda_a, soilprofile.mSoillayers[i],
                                       yTemp, lower, yTemp);
                Qa.mForce.y = Qa.mForce.x * sin(alpha + psi);
                Qa.mForce.x =
                    -forceDirectionCorrection * Qa.mForce.x * cos(alpha + psi);
                // Qa.mPoE.y is nu nog de ycoordinaat tov de bovenkant van de
                // laag
                Qa.mPoE.x = wallBorder +
                            forceDirectionCorrection * Qa.mPoE.y /
                                tan(M_PI / 4.0 + phi_d / 2.0);
                // naar globaal assenstelsel
                Qa.mPoE.y += yTemp;
                mActiveSoilPressure.push_back(Qa);
            }
        } else {
            // nu zitten we dus volledig onder de grondwig zijn top
            // yTemp kan dus negatief zijn
            h = soilprofile.mSoillayers[i].mUpperbounds;
            alpha = M_PI / 4.0 - phi_d / 2.0;
            psi = phi_d;
            lambda_a = soilprofile.getLambda_a(phi_d, alpha, phi_d, 0);
            Qa = soilprofile.getQa(lambda_a, soilprofile.mSoillayers[i], h,
                                   lower, h, alpha);
            Qa.mForce.y = Qa.mForce.x * sin(alpha + psi);
            Qa.mForce.x =
                -forceDirectionCorrection * Qa.mForce.x * cos(alpha + psi);
            Qa.mPoE.x = wallBorder +
                        forceDirectionCorrection * Qa.mPoE.y /
                            tan(M_PI / 4.0 + phi_d / 2.0);
            Qa.mPoE.y += h;
        }
        // nu zijn we tot op de bovenkant van de voet gekomen
        // enkel nog Qa op de wand van die voet, hoogte mHv
        for (size_t i = 0; i < soilprofile.mSoillayers.size(); ++i) {
            h = std::max(mHm - correctionHeight,
                         soilprofile.mSoillayers[i].mUpperbounds);
            if (soilprofile.mSoillayers[i].mUpperbounds <
                (mHm - correctionHeight + mHv)) {
                lower = std::min(soilprofile.mSoillayers[i].mLowerBounds,
                                 mHm + mHv - correctionHeight);
                if (lower > mHm - correctionHeight) {
                    lambda_a = soilprofile.getLambda_a(phi_d, 0, psi, 0);
                    Qa = soilprofile.getQa(lambda_a, soilprofile.mSoillayers[i],
                                           h, lower, h, 0);
                    Qa.mForce.y = Qa.mForce.x * sin(alpha + psi);
                    Qa.mForce.x = -forceDirectionCorrection * Qa.mForce.x *
                                  cos(alpha + psi);
                    Qa.mPoE.x =
                        wallBorder + footwidth * forceDirectionCorrection;
                    Qa.mPoE.y += h;
                }
            }
        }
    }
}

void Lmuur::calculateSoilWedgeWeight(Soilprofile& soilprofile, double base,
                                     double H, double offSet, double side) {
    // voor de rechterkant: offset =mBm, base= mBr, H= mHm
    // voor de linkerkant: offset =0, base = mBl, H=mHm-grondhoogteverschil
    double phi_d = M_PI / 2.0;
    for (size_t i = 0; i < soilprofile.mSoillayers.size(); ++i) {
        phi_d = std::min(phi_d, soilprofile.mSoillayers[i].mPhiA);
    }
    double tanphi_d = tan(M_PI / 4.0 + phi_d / 2.0);
    double triangleHeight = base * tanphi_d;
    double up = 0, low = 0;
    for (size_t i = 0; i < soilprofile.mSoillayers.size(); ++i) {
        if (soilprofile.mSoillayers[i].mUpperbounds < H &&
            soilprofile.mSoillayers[i].mLowerBounds > H - triangleHeight) {
            up = H - std::max(soilprofile.mSoillayers[i].mUpperbounds,
                              H - triangleHeight);
            low = H - std::min(soilprofile.mSoillayers[i].mLowerBounds, H);
            double triangleSurf = 0, squareSurf = 0;
            if (low < up) {
                // valid area
                // now check for water
                double gamma = 0;
                if (H - soilprofile.waterHeight >= up ||
                    H - soilprofile.waterHeight <= low) {
                    // water is out the layer
                    gamma = (H - soilprofile.waterHeight >= up)
                                ? soilprofile.mSoillayers[i].mEffectiveWeight
                                : soilprofile.mSoillayers[i].mDryWeight;
                    triangleSurf =
                        triangleSurface(up - low, (up - low) / tanphi_d);
                    squareSurf = squareSurface((triangleHeight - up) / tanphi_d,
                                               up - low);
                    double sx = (triangleHeight - up) / (2.0 * tanphi_d);
                    double sy = H - (up - low) / 2.0;
                    double tx = (triangleHeight - up) / (tanphi_d) +
                                (2.0 / 3.0) * (up - low) / tanphi_d;
                    double ty = (2.0 / 3.0) * (up - low) + (H - up);
                    ForceVector Fweight = ForceVector(
                        glm::vec2(0, (triangleSurf + squareSurf) * gamma),
                        glm::vec2(
                            side * (offSet +
                                    (sx * squareSurf + tx * triangleSurf) /
                                        (squareSurf + triangleSurf)),
                            (sy * squareSurf + ty * triangleSurf) /
                                (squareSurf + triangleSurf)));
                    mSoilWedgeWeight.push_back(Fweight);
                } else {
                    // water goes through the layer
                    double gammaEf =
                        soilprofile.mSoillayers[i].mEffectiveWeight;
                    double gammaD = soilprofile.mSoillayers[i].mDryWeight;
                    double waterLevel = H - soilprofile.waterHeight;
                    triangleSurf = triangleSurface(
                        up - waterLevel, (up - waterLevel) / tanphi_d);
                    squareSurf = squareSurface((triangleHeight - up) / tanphi_d,
                                               up - waterLevel);
                    double triangleSurf2 = triangleSurface(
                        waterLevel - low, (waterLevel - low) / tanphi_d);
                    double squareSurf2 =
                        squareSurface((triangleHeight - waterLevel) / tanphi_d,
                                      waterLevel - low);
                    double sx = (triangleHeight - up) / (2.0 * tanphi_d);
                    double sy = H - (up - waterLevel) / 2.0;
                    double tx = (triangleHeight - up) / (tanphi_d) +
                                (2.0 / 3.0) * (up - waterLevel) / tanphi_d;
                    double ty = (2.0 / 3.0) * (up - waterLevel) + (H - up);
                    double sx2 =
                        (triangleHeight - waterLevel) / (2.0 * tanphi_d);
                    double sy2 = H - (waterLevel - low) / 2.0;
                    double tx2 = (triangleHeight - waterLevel) / (tanphi_d) +
                                 (2.0 / 3.0) * (waterLevel - low) / tanphi_d;
                    double ty2 =
                        (2.0 / 3.0) * (waterLevel - low) + (H - waterLevel);
                    double resultingPoEx =
                        (gammaD * (triangleSurf * tx + squareSurf * sx) +
                         gammaEf * (triangleSurf2 * tx2 + squareSurf2 * sx2)) /
                        (gammaD * (triangleSurf + squareSurf) +
                         gammaEf * (triangleSurf2 + squareSurf2));
                    double resultingPoEy =
                        (gammaD * (triangleSurf * ty + squareSurf * sy) +
                         gammaEf * (triangleSurf2 * ty2 + squareSurf2 * sy2)) /
                        (gammaD * (triangleSurf + squareSurf) +
                         gammaEf * (triangleSurf2 + squareSurf2));
                    mSoilWedgeWeight.push_back(ForceVector(
                        glm::vec2(0,
                                  gammaD * (squareSurf + triangleSurf) +
                                      gammaEf * (squareSurf2 + triangleSurf2)),
                        glm::vec2(side * (offSet + resultingPoEx),
                                  resultingPoEy)));
                }
            }
        }
    }
}

void Lmuur::calculateWaterPressures() {
    double waterHeightLeft = mSoilHeightDifference + leftProfile.waterHeight;
    double waterHeightRight = rightProfile.waterHeight;
    double h_a = waterHeightLeft - waterHeightRight;
    double d = mHm + mHv - h_a;
    double kastnerHcorr = d + h_a / (1 + cbrt(1 + h_a / d));
    ForceVector rhsWaterpressure = ForceVector(
        glm::vec2(-gamma_water * kastnerHcorr * (mHm + mHv) * 0.5, 0),
        glm::vec2(mBm, (2.0 / 3.0) * (mHm + mHv)));
    ForceVector lhsWaterpressure =
        ForceVector(glm::vec2(gamma_water * kastnerHcorr * d * 0.5, 0),
                    glm::vec2(0, h_a + (2.0 / 3.0) * d));
}

void Lmuur::calculateBoussinesqLoads() {
    mBoussinesqResultant.push_back(
        ForceVector(glm::vec2(-0.5 * mq * mHm, 0), glm::vec2(mBm, mHm * 0.5)));
    mBoussinesqResultant.push_back(ForceVector(
        glm::vec2(0, 0.5 * mq * mBr), glm::vec2(mBm + mBr * 0.5, mHm)));
    mBoussinesqResultant.push_back(
        ForceVector(glm::vec2(-0.5 * mq * mHv, mq * mHv / M_PI),
                    glm::vec2(mBm + mBr, mHm + mHv * 0.5)));
}

void Lmuur::calculateBuoyancy() {
    double wetHeight = mHm - rightProfile.waterHeight;
    double newWetArea = mAI * (wetHeight / mHm);
    double y1 = rightProfile.waterHeight + wetHeight * 0.5;
    mBuoyantForce = ForceVector(
        glm::vec2(0, -((newWetArea + mAII) * gamma_water)),
        glm::vec2((newWetArea * mxI + mAII * mxII) / (newWetArea + mAII),
                  (newWetArea * y1 + mAII * myII) / (newWetArea + mAII)));
}
void Lmuur::calculateAll() {
    calculateProperties();
    calculateActiveSoilPressures(rightProfile, 1, mBr);
    calculateActiveSoilPressures(leftProfile, -1, mBl);
    calculateSoilWedgeWeight(rightProfile, mBr, mHm, mBm, 1);
    calculateSoilWedgeWeight(leftProfile, mBl, mHm - mSoilHeightDifference, 0,
                             -1);
    calculateWaterPressures();
    calculateBuoyancy();
    calculateBoussinesqLoads();
}
void Lmuur::calculateActiveSoilPressureLeft() {}
double Lmuur::squareSurface(double height, double width) {
    return height * width;
}

double Lmuur::triangleSurface(double height, double base) {
    return height * base * 0.5;
}

void Lmuur::addSoilprofiles(Soilprofile Right, Soilprofile Left) {
    leftProfile = Left;
    rightProfile = Right;
}

Lmuur::~Lmuur() {}
