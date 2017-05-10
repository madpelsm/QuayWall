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
    std::cout << "Wall properties calculated!" << std::endl;
}

void Lmuur::calculateActiveSoilPressures(Soilprofile& soilprofile, double side,
                                         double footwidth, int safetyCase) {
    double correctionHeight = 0;
    double forceDirectionCorrection = 1;
    double wallBorder = mBm;
    double toeHeight = 0;
    if (side != 1) {
        forceDirectionCorrection = -1;
        correctionHeight = mSoilHeightDifference;
        wallBorder = 0;
        toeHeight = mtoe;
    }
    if (safetyCase > soilprofile.mSoillayers[0].mSafetyPhi.size()) {
        abort();
    }

    double phi_d = M_PI / 2.0;
    for (size_t i = 0; i < soilprofile.mSoillayers.size(); ++i) {
        phi_d =
            std::min(phi_d, soilprofile.mSoillayers[i].mSafetyPhi[safetyCase]);
    }

    ForceVector Qa;
    double psi = 0, phi = 0, alpha = 0;
    double lambda_a = 0;
    double yTemp =
        (mHm - correctionHeight) - footwidth * tan(M_PI / 4.0 + phi_d);
    for (size_t i = 0; i < soilprofile.mSoillayers.size(); ++i) {
        // dependencies for lambda_a
        psi = phi_d * (2.0 / 3.0);
        alpha = M_PI / 4.0 - phi_d / 2.0;
        // dependencies for Qa

        // h = diepte tov maaiveld (pos = naar beneden)
        double h = soilprofile.mSoillayers[i].mUpperbounds;
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
                psi = (2.0 / 3.0) *
                      soilprofile.mSoillayers[i].mSafetyPhi[safetyCase];
                lambda_a = soilprofile.getLambda_a(
                    soilprofile.mSoillayers[i].mSafetyPhi[safetyCase], 0, psi,
                    0);
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
                Qa.mPoE.y += h + correctionHeight;
                mActiveSoilPressure.push_back(Qa);
            } else {
                // glijlijn ligt boven de onderkant van de beschouwde
                // grondlaag
                // bovenkant van de grondlaag ligt nog boven de glijlijn
                lower = yTemp;
                alpha = 0;
                psi = (2.0 / 3.0) *
                      soilprofile.mSoillayers[i].mSafetyPhi[safetyCase];
                lambda_a = soilprofile.getLambda_a(
                    soilprofile.mSoillayers[i].mSafetyPhi[safetyCase], 0, psi,
                    0);
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
                psi = soilprofile.mSoillayers[i].mSafetyPhi[safetyCase];
                lambda_a = soilprofile.getLambda_a(
                    soilprofile.mSoillayers[i].mSafetyPhi[safetyCase], 0, psi,
                    0);
                Qa = soilprofile.getQa(lambda_a, soilprofile.mSoillayers[i],
                                       yTemp, lower, yTemp);
                Qa.mForce.y = Qa.mForce.x * sin(alpha + psi);
                Qa.mForce.x =
                    -forceDirectionCorrection * Qa.mForce.x * cos(alpha + psi);
                // Qa.mPoE.y is nu nog de ycoordinaat tov de bovenkant van de
                // laag
                Qa.mPoE.x =
                    wallBorder +
                    forceDirectionCorrection * Qa.mPoE.y /
                        tan(M_PI / 4.0 +
                            soilprofile.mSoillayers[i].mSafetyPhi[safetyCase] /
                                2.0);
                // naar globaal assenstelsel
                Qa.mPoE.y += yTemp + correctionHeight;
                mActiveSoilPressure.push_back(Qa);
            }
        } else {
            // nu zitten we dus volledig onder de grondwig zijn top
            // yTemp kan dus negatief zijn
            h = soilprofile.mSoillayers[i].mUpperbounds;
            if (h < lower) {
                alpha = M_PI / 4.0 - phi_d / 2.0;
                psi = soilprofile.mSoillayers[i].mSafetyPhi[safetyCase];
                lambda_a = soilprofile.getLambda_a(
                    soilprofile.mSoillayers[i].mSafetyPhi[safetyCase], 0, psi,
                    0);
                Qa = soilprofile.getQa(lambda_a, soilprofile.mSoillayers[i], h,
                                       lower, h, alpha);
                Qa.mForce.y =
                    Qa.mForce.x * abs(sin(alpha + psi)) + correctionHeight;
                Qa.mForce.x = -forceDirectionCorrection * Qa.mForce.x *
                              abs(cos(alpha + psi));
                Qa.mPoE.x =
                    wallBorder +
                    forceDirectionCorrection * Qa.mPoE.y /
                        tan(M_PI / 4.0 +
                            soilprofile.mSoillayers[i].mSafetyPhi[safetyCase] /
                                2.0);
                Qa.mPoE.y += h + correctionHeight;
                mActiveSoilPressure.push_back(Qa);
            }
        }
        // nu zijn we tot op de bovenkant van de voet gekomen
        // enkel nog Qa op de wand van die voet, hoogte mHv
    }
    for (size_t i = 0; i < soilprofile.mSoillayers.size(); ++i) {
        double h = std::max(soilprofile.mSoillayers[i].mUpperbounds,
                            mHm - correctionHeight);
        double lower = std::min(soilprofile.mSoillayers[i].mLowerBounds,
                                mHm - correctionHeight);
        if (soilprofile.mSoillayers[i].mUpperbounds <
                (mHm - correctionHeight + mHv + toeHeight) &&
            soilprofile.mSoillayers[i].mLowerBounds > mHm - correctionHeight) {
            lower = std::min(soilprofile.mSoillayers[i].mLowerBounds,
                             mHm + mHv - correctionHeight + toeHeight);
            if (lower > mHm - correctionHeight) {
                alpha = 0;
                psi = (2.0 / 3.0) *
                      soilprofile.mSoillayers[i].mSafetyPhi[safetyCase];
                lambda_a = soilprofile.getLambda_a(
                    soilprofile.mSoillayers[i].mSafetyPhi[safetyCase], 0, psi,
                    0);
                Qa = soilprofile.getQa(lambda_a, soilprofile.mSoillayers[i], h,
                                       lower, h, 0);
                Qa.mForce.y = Qa.mForce.x * sin(alpha + psi);
                Qa.mForce.x =
                    -forceDirectionCorrection * Qa.mForce.x * cos(alpha + psi);
                Qa.mPoE.x = wallBorder + footwidth * forceDirectionCorrection;
                Qa.mPoE.y += h + correctionHeight;
                mActiveSoilPressure.push_back(Qa);
            }
        }
    }
}

void Lmuur::calculateSoilWedgeWeight(Soilprofile& soilprofile, double base,
                                     double H, double offSet, double side,
                                     int safetyCase) {
    // voor de rechterkant: offset =mBm, base= mBr, H= mHm
    // voor de linkerkant: offset =0, base = mBl, H=mHm-grondhoogteverschil
    double phi_d = M_PI / 2.0;
    for (size_t i = 0; i < soilprofile.mSoillayers.size(); ++i) {
        phi_d =
            std::min(phi_d, soilprofile.mSoillayers[i].mSafetyPhi[safetyCase]);
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

void Lmuur::calculatePassiveSoilPressure(Soilprofile& soilprofile, double side,
                                         double footwidth, int safetyCase) {
    double correctionHeight = 0;
    double forceDirectionCorrection = 1;
    double wallBorder = mBm;
    if (side != 1) {
        // if side !=1 thus -1, we are on the negative side of the x axis. i.e.
        // on the left side of the wall
        forceDirectionCorrection = -1;
        correctionHeight = mSoilHeightDifference;
        wallBorder = 0;
    }

    double phi_d = M_PI / 2.0;
    for (size_t i = 0; i < soilprofile.mSoillayers.size(); ++i) {
        phi_d =
            std::min(phi_d, soilprofile.mSoillayers[i].mSafetyPhi[safetyCase]);
    }
    phi_d *= -1;
    double psi = (3.0 / 2.0) * phi_d;
    ForceVector Qp;
    double lambda_p = 1, alpha = 0;
    for (size_t i = 0; i < soilprofile.mSoillayers.size(); ++i) {
        double h = std::max(soilprofile.mSoillayers[i].mUpperbounds,
                            mHm - correctionHeight);
        double lower = std::min(soilprofile.mSoillayers[i].mLowerBounds,
                                mHm - correctionHeight);
        if (soilprofile.mSoillayers[i].mUpperbounds <
                (mHm - correctionHeight + mHv + mtoe) &&
            soilprofile.mSoillayers[i].mLowerBounds > mHm - correctionHeight) {
            lower = std::min(soilprofile.mSoillayers[i].mLowerBounds,
                             mHm + mHv + mtoe - correctionHeight);
            if (lower > mHm - correctionHeight) {
                alpha = 0;
                psi = -(2.0 / 3.0) *
                      soilprofile.mSoillayers[i].mSafetyPhi[safetyCase];
                Qp = soilprofile.getQa(lambda_p, soilprofile.mSoillayers[i], h,
                                       lower, h, 0);
                Qp.mForce.y = Qp.mForce.x * sin(alpha + psi);
                Qp.mForce.x =
                    -forceDirectionCorrection * Qp.mForce.x * cos(alpha + psi);
                Qp.mPoE.x = wallBorder + footwidth * forceDirectionCorrection;
                Qp.mPoE.y += h + correctionHeight;
                mPassiveSoilPressure.push_back(Qp);
            }
        }
    }
    std::cout << "Passive soil pressure calculation completed!" << std::endl;
}

void Lmuur::calculateWaterPressures() {
    double waterHeightLeft = mSoilHeightDifference + leftProfile.waterHeight;
    double waterHeightRight = rightProfile.waterHeight;
    double d = mHm + mHv - waterHeightLeft;
    double kastnerHcorr =
        d + (waterHeightLeft / (1 + cbrt((1 + waterHeightLeft / d))));
    kastnerH = kastnerHcorr;
    mrhsWaterPressure = ForceVector(
        glm::vec2(-gamma_water * kastnerHcorr * (mHm + mHv) * 0.5, 0),
        glm::vec2(mBm, (2.0 / 3.0) * (mHm + mHv)));
    mlhsWaterPressure =
        ForceVector(glm::vec2(gamma_water * kastnerHcorr * (mtoe + d) * 0.5, 0),
                    glm::vec2(0, waterHeightLeft + (2.0 / 3.0) * (d + mtoe)));
    std::cout << "Water pressures calculatio completed!" << std::endl;
}

void Lmuur::calculateBoussinesqLoads() {
    mBoussinesqResultant.push_back(
        ForceVector(glm::vec2(-0.5 * mq * mHm, 0), glm::vec2(mBm, mHm * 0.5)));
    mBoussinesqResultant.push_back(ForceVector(
        glm::vec2(0, 0.5 * mq * mBr), glm::vec2(mBm + mBr * 0.5, mHm)));
    mBoussinesqResultant.push_back(
        ForceVector(glm::vec2(-0.5 * mq * mHv, mq * mHv / M_PI),
                    glm::vec2(mBm + mBr, mHm + mHv * 0.5)));
    std::cout << "Boussinesq calculation completed!" << std::endl;
}

void Lmuur::calculateBuoyancy() {
    double wetHeight = mHm - rightProfile.waterHeight;
    double newWetArea = mAI * (wetHeight / mHm);
    double y1 = rightProfile.waterHeight + wetHeight * 0.5;
    mBuoyantForce =
        ForceVector(glm::vec2(0, -((newWetArea + mAII + mAtoe) * gamma_water)),
                    glm::vec2((newWetArea * mxI + mAII * mxII + mxtoe * mAtoe) /
                                  (newWetArea + mAII + mAtoe),
                              (newWetArea * y1 + mAII * myII + mAtoe * mytoe) /
                                  (newWetArea + mAII + mAtoe)));
    std::cout << "Buoyancy calculation completed!" << std::endl;
}

void Lmuur::setSolidHeightDifference(double height) {
    mSoilHeightDifference = height;
}

bool Lmuur::isBeneficial(ForceVector& _forcevector, int _failuremode,
                         int _safetycase) {
    bool b = false;
    if (_failuremode == 0) {
        // evenwichtsdraagvermogen
        b = _forcevector.mForce.y < 0;
    } else if (_failuremode == 1) {
        if (_forcevector.mForce.x >= 0 && _forcevector.mForce.y >= 0) {
            b = true;
        } else if (_forcevector.mForce.x > 0 && _forcevector.mForce.y < 0) {
            b = _forcevector.mForce.x >
                (-1.0) * _forcevector.mForce.y *
                    std::tan(getPhiAtConstructionFoot(_safetycase));
        } else if (_forcevector.mForce.x <= 0 && _forcevector.mForce.y > 0) {
            b = (-1.0) * _forcevector.mForce.x <
                _forcevector.mForce.y *
                    std::tan(getPhiAtConstructionFoot(_safetycase));
        }
    }
    return b;
}

ForceVector Lmuur::calculateResultingForce(int failureMode, int safetyCase) {
    // Calculate the resulting force
    double sumFx = 0, sumFy = 0, poex = 0, poey = 0, mOx = 0, mOy = 0;
    // mOx = sum(Fx*y), mOy = sum(Fy*x)
    double safetyCoef = 1;
    double ssoil = mSafetyGSup[safetyCase], sactive = mSafetyGSup[safetyCase],
           sbous = mSafetyQ[safetyCase], spassive = mSafetyGInf[safetyCase],
           sownweight = mSafetyGSup[safetyCase],
           sbuoy = mSafetyGInf[safetyCase], swater = mSafetyGSup[safetyCase];
    if (failureMode == 0) {
        // evenwichtsdraagvermogen
        ssoil = mSafetyGSup[safetyCase], sactive = mSafetyGSup[safetyCase],
        sbous = mSafetyQ[safetyCase], spassive = mSafetyGInf[safetyCase],
        sownweight = mSafetyGSup[safetyCase], sbuoy = mSafetyGInf[safetyCase],
        swater = mSafetyGSup[safetyCase];
    }
    if (failureMode == 1) {
        // schuiven
        ssoil = mSafetyGInf[safetyCase], sactive = mSafetyGSup[safetyCase],
        sbous = mSafetyQ[safetyCase], spassive = mSafetyGInf[safetyCase],
        sownweight = mSafetyGInf[safetyCase], sbuoy = mSafetyGInf[safetyCase],
        swater = mSafetyGSup[safetyCase];
    } else {
        // alles op 1 zetten
        ssoil = 1.0, sactive = 1.0, sbous = 1.0, spassive = 1.0,
        sownweight = 1.0, sbuoy = 1.0, swater = 1.0;
    }

    for (size_t i = 0; i < mSoilWedgeWeight.size(); ++i) {
        sumFx += ssoil * mSoilWedgeWeight[i].mForce.x;
        sumFy += ssoil * mSoilWedgeWeight[i].mForce.y;
        mOy +=
            ssoil * mSoilWedgeWeight[i].mForce.y * mSoilWedgeWeight[i].mPoE.x;
        mOx +=
            ssoil * mSoilWedgeWeight[i].mForce.x * mSoilWedgeWeight[i].mPoE.y;
    }
    for (size_t i = 0; i < mActiveSoilPressure.size(); ++i) {
        sumFx += sactive * mActiveSoilPressure[i].mForce.x;
        sumFy += sactive * mActiveSoilPressure[i].mForce.y;
        mOy += sactive * mActiveSoilPressure[i].mForce.y *
               mActiveSoilPressure[i].mPoE.x;
        mOx += sactive * mActiveSoilPressure[i].mForce.x *
               mActiveSoilPressure[i].mPoE.y;
    }
    for (size_t i = 0; i < mBoussinesqResultant.size(); ++i) {
        sumFx += sbous * mBoussinesqResultant[i].mForce.x;
        sumFy += sbous * mBoussinesqResultant[i].mForce.y;
        mOy += sbous * mBoussinesqResultant[i].mForce.y *
               mBoussinesqResultant[i].mPoE.x;
        mOx += sbous * mBoussinesqResultant[i].mForce.x *
               mBoussinesqResultant[i].mPoE.y;
    }
    // Passive
    for (size_t i = 0; i < mPassiveSoilPressure.size(); ++i) {
        sumFx += spassive * mPassiveSoilPressure[i].mForce.x;
        sumFy += spassive * mPassiveSoilPressure[i].mForce.y;
        mOy += spassive * mPassiveSoilPressure[i].mForce.y *
               mPassiveSoilPressure[i].mPoE.x;
        mOx += spassive * mPassiveSoilPressure[i].mForce.x *
               mPassiveSoilPressure[i].mPoE.y;
    }
    // own weight
    sumFx += sownweight * mOwnWeight.mForce.x;
    sumFy += sownweight * mOwnWeight.mForce.y;
    mOy += sownweight * mOwnWeight.mForce.y * mOwnWeight.mPoE.x;
    mOx += sownweight * mOwnWeight.mForce.x * mOwnWeight.mPoE.y;
    // buoyant
    sumFx += sbuoy * mBuoyantForce.mForce.x;
    sumFy += sbuoy * mBuoyantForce.mForce.y;
    mOy += sbuoy * mBuoyantForce.mForce.y * mBuoyantForce.mPoE.x;
    mOx += sbuoy * mBuoyantForce.mForce.x * mBuoyantForce.mPoE.y;
    // rhsWaterPressure
    sumFx += swater * mrhsWaterPressure.mForce.x;
    sumFy += swater * mrhsWaterPressure.mForce.y;
    mOy += swater * mrhsWaterPressure.mForce.y * mrhsWaterPressure.mPoE.x;
    mOx += swater * mrhsWaterPressure.mForce.x * mrhsWaterPressure.mPoE.y;
    // lhsWaterPressure
    sumFx += swater * mlhsWaterPressure.mForce.x;
    sumFy += swater * mlhsWaterPressure.mForce.y;
    mOy += swater * mlhsWaterPressure.mForce.y * mlhsWaterPressure.mPoE.x;
    mOx += swater * mlhsWaterPressure.mForce.x * mlhsWaterPressure.mPoE.y;

    // mOx = sum(Fx*y), mOy = sum(Fy*x)

    std::cout << "resulting force calculation completed." << std::endl;
    return ForceVector(glm::vec2(sumFx, sumFy),
                       glm::vec2(mOy / sumFy, mOx / sumFx));
}

void Lmuur::calculateTiltMomentAtFoot(int safetyCase) {
    double distanceX = 0;
    double distanceY = 0;
    double momenti = 0;
    for (size_t i = 0; i < mSoilWedgeWeight.size(); ++i) {
        distanceX = mSoilWedgeWeight[i].mPoE.x - mBl;
        distanceY = mSoilWedgeWeight[i].mPoE.y - mHv - mHm - mtoe;
        momenti = distanceX * mSoilWedgeWeight[i].mForce.y -
                  distanceY * mSoilWedgeWeight[i].mForce.x;
        if (momenti >= 0) {
            momentST += mSafetyGInf[safetyCase] * momenti;
        } else {
            momentDST += mSafetyGSup[safetyCase] * momenti;
        }
    }
    for (size_t i = 0; i < mActiveSoilPressure.size(); ++i) {
        distanceX = mActiveSoilPressure[i].mPoE.x - mBl;
        distanceY = mActiveSoilPressure[i].mPoE.y - mHv - mHm - mtoe;
        momenti = distanceX * mActiveSoilPressure[i].mForce.y -
                  distanceY * mActiveSoilPressure[i].mForce.x;
        if (momenti >= 0) {
            momentST += mSafetyGInf[safetyCase] * momenti;
        } else {
            momentDST += mSafetyGSup[safetyCase] * momenti;
        }
    }
    for (size_t i = 0; i < mPassiveSoilPressure.size(); ++i) {
        distanceX = mPassiveSoilPressure[i].mPoE.x - mBl;
        distanceY = mPassiveSoilPressure[i].mPoE.y - mHv - mHm - mtoe;
        momenti = distanceX * mPassiveSoilPressure[i].mForce.y -
                  distanceY * mPassiveSoilPressure[i].mForce.x;
        if (momenti >= 0) {
            momentST += mSafetyGInf[safetyCase] * momenti;
        } else {
            momentDST += mSafetyGSup[safetyCase] * momenti;
        }
    }
    for (size_t i = 0; i < mBoussinesqResultant.size(); ++i) {
        distanceX = mBoussinesqResultant[i].mPoE.x - mBl;
        distanceY = mBoussinesqResultant[i].mPoE.y - mHv - mHm - mtoe;
        momenti = distanceX * mBoussinesqResultant[i].mForce.y -
                  distanceY * mBoussinesqResultant[i].mForce.x;
        if (momenti >= 0) {
            momentST += mSafetyQ[3] * momenti;
        } else {
            momentDST += mSafetyQ[3] * momenti;
        }
    }
    distanceX = mOwnWeight.mPoE.x - mBl;
    distanceY = mOwnWeight.mPoE.y - mHv - mHm - mtoe;
    momenti = distanceX * mOwnWeight.mForce.y - distanceY * mOwnWeight.mForce.x;
    if (momenti >= 0) {
        momentST += mSafetyGInf[safetyCase] * momenti;
    } else {
        momentDST += mSafetyGSup[safetyCase] * momenti;
    }
    distanceX = mBuoyantForce.mPoE.x - mBl;
    distanceY = mBuoyantForce.mPoE.y - mHv - mHm - mtoe;
    momenti =
        distanceX * mBuoyantForce.mForce.y - distanceY * mBuoyantForce.mForce.x;
    if (momenti >= 0) {
        momentST += mSafetyGInf[safetyCase] * momenti;
    } else {
        momentDST += mSafetyGSup[safetyCase] * momenti;
    }
    distanceX = mlhsWaterPressure.mPoE.x - mBl;
    distanceY = mlhsWaterPressure.mPoE.y - mHv - mHm - mtoe;
    momenti = distanceX * mlhsWaterPressure.mForce.y -
              distanceY * mlhsWaterPressure.mForce.x;
    if (momenti >= 0) {
        momentST += mSafetyGInf[safetyCase] * momenti;
    } else {
        momentDST += mSafetyGSup[safetyCase] * momenti;
    }
    distanceX = mrhsWaterPressure.mPoE.x - mBl;
    distanceY = mrhsWaterPressure.mPoE.y - mHv - mHm - mtoe;
    momenti = distanceX * mrhsWaterPressure.mForce.y -
              distanceY * mrhsWaterPressure.mForce.x;
    if (momenti >= 0) {
        momentST += mSafetyGInf[safetyCase] * momenti;
    } else {
        momentDST += mSafetyGSup[safetyCase] * momenti;
    }
    std::cout << "Tilt calculation completed!" << std::endl;
}

void Lmuur::clearForces() {
    mSoilWedgeWeight.clear();
    mActiveSoilPressure.clear();
    mBoussinesqResultant.clear();
    mPassiveSoilPressure.clear();
}

void Lmuur::calculateAll(int safetyCase) {
    std::cout << "starting complete calculation" << std::endl;

    clearForces();
    calculateProperties();
    calculateActiveSoilPressures(rightProfile, 1, mBr, safetyCase);
    calculateActiveSoilPressures(leftProfile, -1, mBl, safetyCase);
    calculatePassiveSoilPressure(leftProfile, -1, mBl, safetyCase);
    calculateSoilWedgeWeight(rightProfile, mBr, mHm, mBm, 1, safetyCase);
    calculateSoilWedgeWeight(leftProfile, mBl, mHm - mSoilHeightDifference, 0,
                             -1, safetyCase);
    calculateWaterPressures();
    calculateBuoyancy();
    calculateBoussinesqLoads();

    mResultingR_d = calculateResultingForce(0, safetyCase);
    calculateExcentricity();

    mResultingR_dSchuiven = calculateResultingForce(1, safetyCase);

    makeUnityChecks(safetyCase);
}

void Lmuur::calculateExcentricity() {
    mExcentricity = mResultingR_d.mPoE.x - mxII;
}

double Lmuur::squareSurface(double height, double width) {
    return height * width;
}

double Lmuur::triangleSurface(double height, double base) {
    return height * base * 0.5;
}

void Lmuur::writeToCSV(std::string file_name) {
    std::ofstream file;
    file.open(file_name);
    if (file.is_open()) {
        file << "L vormige kaaimuur\n";
        file << "Eigenschappen L-wand,[m],,[m]\n";
        file << "Hoogte vanaf de voet," << mHm << ",Breedte verticale wand,"
             << mBm << "\n";
        file << "Breedte linker voet," << mBl << ",Breedte rechter voet," << mBr
             << "\n";
        file << "Breedte verticale wand," << mBm << ",Dikte funderingszool,"
             << mHv << "\n";
        file << "Afmeting dextrateen:," << mtoe << "\n";
        file << "Opmerkingen\n";
        file << "Het assenstelsel gaat met positieve richting naar "
                "beneden\nDeze y as valt samen met de linkerzijde van de "
                "verticale "
                "wand\nDe x as "
                "gaat positief naar rechts\n samenvallend met het "
                "horizontale "
                "maaiveld.\n\n";
        file << "Aangrijpende krachten. \n";
        file << "Kracht,Fx[kN],Fy[kN],x[m],y[m]\n";
        file << "Eigengewicht.\n";
        file << "Eigen," << mOwnWeight.mForce.x << "," << mOwnWeight.mForce.y
             << "," << mOwnWeight.mPoE.x << "," << mOwnWeight.mPoE.y << "\n";

        file << "Opwaartse stuwkracht van het beton.\n";
        file << "Archimedes," << mBuoyantForce.mForce.x << ","
             << mBuoyantForce.mForce.y << "," << mBuoyantForce.mPoE.x << ","
             << mBuoyantForce.mPoE.y << "\n";

        file << "Waterdrukken mbv Kastner.\n";
        file << "Rechterkant," << mrhsWaterPressure.mForce.x << ","
             << mrhsWaterPressure.mForce.y << "," << mrhsWaterPressure.mPoE.x
             << "," << mrhsWaterPressure.mPoE.y << "\n";
        file << "Linkerkant," << mlhsWaterPressure.mForce.x << ","
             << mlhsWaterPressure.mForce.y << "," << mlhsWaterPressure.mPoE.x
             << "," << mlhsWaterPressure.mPoE.y << "\n";
        file << "Actieve gronddrukken.\n";
        for (size_t i = 0; i < mActiveSoilPressure.size(); ++i) {
            file << "Qa" << i << "," << mActiveSoilPressure[i].mForce.x << ","
                 << mActiveSoilPressure[i].mForce.y << ","
                 << mActiveSoilPressure[i].mPoE.x << ","
                 << mActiveSoilPressure[i].mPoE.y << "\n";
        }
        file << "Passieve gronddrukken.\n";
        for (size_t i = 0; i < mPassiveSoilPressure.size(); ++i) {
            file << "Qa" << i << "," << mPassiveSoilPressure[i].mForce.x << ","
                 << mPassiveSoilPressure[i].mForce.y << ","
                 << mPassiveSoilPressure[i].mPoE.x << ","
                 << mPassiveSoilPressure[i].mPoE.y << "\n";
        }
        file << "GrondGewicht\n";
        for (size_t i = 0; i < mSoilWedgeWeight.size(); ++i) {
            file << "G" << i << "," << mSoilWedgeWeight[i].mForce.x << ","
                 << mSoilWedgeWeight[i].mForce.y << ","
                 << mSoilWedgeWeight[i].mPoE.x << ","
                 << mSoilWedgeWeight[i].mPoE.y << "\n";
        }
        file << "Boussinesq.\n";
        for (size_t i = 0; i < mBoussinesqResultant.size(); ++i) {
            file << "F" << i << "," << mBoussinesqResultant[i].mForce.x << ","
                 << mBoussinesqResultant[i].mForce.y << ","
                 << mBoussinesqResultant[i].mPoE.x << ","
                 << mBoussinesqResultant[i].mPoE.y << "\n";
        }
        file << "\n";
        file << "Resultante\n";
        file << "F_d," << mResultingR_d.mForce.x << ","
             << mResultingR_d.mForce.y << "," << mResultingR_d.mPoE.x << ","
             << mResultingR_d.mPoE.y << "\n";
        file << "\nUNITY CHECKS\n";
        file << "type,R_d,E_d,veiligheid\n";
        file << "Evewichtsdraagvermogen," << R_d << ","
             << mResultingR_d.mForce.y * 100 << ","
             << R_d / (100 * mResultingR_d.mForce.y) << ",excentriciteit:,"
             << mExcentricity << "\n";
        file << "Schuiven," << RH_d << "," << mResultingR_dSchuiven.mForce.x
             << "," << abs(RH_d / mResultingR_dSchuiven.mForce.x) << "\n";
        file << "Kantelen," << momentST << "," << momentDST << ","
             << (momentDST != 0 ? (std::abs(momentST / momentDST)) : 0)
             << ",moment in kNm,\n";
        if (!file) {
            std::cout << "Write failed." << std::endl;
        }
        file.close();
        if (file) {
            std::cout << "Write complete!" << std::endl;
        } else {
            std::cout << "Write failed." << std::endl;
        }
    } else {
        std::cout << "Write failed." << std::endl;
    }
}

void Lmuur::makeUnityChecks(int safetyCase) {
    double phi_d = getPhiAtConstructionFoot(safetyCase);
    R_d =
        calculateR_d(phi_d, leftProfile, mHm + mHv - mSoilHeightDifference, 1);
    RH_d = mResultingR_dSchuiven.mForce.y * tan(phi_d);

    calculateTiltMomentAtFoot(safetyCase);
}

double Lmuur::getPhiAtConstructionFoot(int safetyCase) {
    // since left and right have to be compatible left or right doesnt
    // matter
    double phi_d = rightProfile.mSoillayers[0].mSafetyPhi[0];
    for (size_t i = 0; i < rightProfile.mSoillayers.size(); ++i) {
        if (rightProfile.mSoillayers[i].mUpperbounds <= (mHm + mHv) &&
            rightProfile.mSoillayers[i].mLowerBounds > (mHm + mHv)) {
            phi_d = rightProfile.mSoillayers[i].mSafetyPhi[safetyCase];
        }
    }
    return phi_d;
}

double Lmuur::calculateR_d(double phi_d, Soilprofile& soilprofile, double depth,
                           double effectiveCohesion_safetyF) {
    // Safetyfactors on cohesion
    double c_d = mCohesion / effectiveCohesion_safetyF;

    double B = mBz - 2 * mExcentricity;
    double L = 100;
    // invloedsfactoren
    double N_q = exp(M_PI * tan(phi_d)) * tan(M_PI / 4.0 + phi_d / 2.0) *
                 tan(M_PI / 4.0 + phi_d / 2.0);
    double N_c = 0;
    if (phi_d != 0) {
        N_c = (N_q - 1) / tan(phi_d);
    } else if (phi_d == 0) {
        N_c = M_PI + 2;
    }
    double N_g = 2 * (N_q - 1) * tan(phi_d);
    // korrespanning aan aanzet
    double p_t = soilprofile.getEffectiveSoilePressure(depth);
    // dieptefactore d = 1
    double d_c = 1;
    double d_q = 1;
    // vormfactoren
    double s_q = 1 + B / L * sin(phi_d);
    double s_c = (s_q * N_q - 1) / (N_q - 1);
    double s_g = 1 - 0.3 * B / L;
    // effectief volumegewicht van steungevende grond die er naast ligt
    double g_k = 10000;
    for (size_t i = 0; i < soilprofile.mSoillayers.size(); ++i) {
        if (depth > soilprofile.mSoillayers[i].mUpperbounds) {
            g_k = std::min(g_k, soilprofile.mSoillayers[i].mEffectiveWeight);
        }
    }
    if (soilprofile.mSoillayers.size() != 0) {
        g_k = std::min(soilprofile.mSoillayers[0].mEffectiveWeight, g_k);
    }
    // hellingsfactoren
    double i_q =
        pow((1 -
             0.7 * abs(mResultingR_d.mForce.x) /
                 (abs(mResultingR_d.mForce.y) + B * L / (tan(phi_d) * c_d))),
            3);
    double i_c = (i_q * N_q - 1) / (N_q - 1);
    double i_g =
        pow((1 -
             abs(mResultingR_d.mForce.x) /
                 (abs(mResultingR_d.mForce.y) + B * L / (tan(phi_d) * c_d))),
            3);
    double R_dToReturn =
        B * L * (d_q * s_q * N_q * p_t * i_q + d_c * s_c * N_c * i_c * c_d +
                 s_g * N_g * g_k * i_g * B / 2.0);
    // std::cout << "q_u = d_q*s_q*N_q*p_t+d_c*s_c*N_c*c+s_g*N_g*g_k*B/2\n"
    //           << d_q << "*" << s_q << "*" << N_q << "*" << p_t << "+" <<
    //           d_c
    //           << "*" << s_c << "*" << N_c << "*" << c << "/" <<
    //           effectiveCohesion_safetyF  << "+" << s_g << "*"
    //           << N_g << "*" << g_k << "*" << B << "/" << 2.0 << "=" <<
    //           _qu
    //           << "\n"
    //           << std::endl;
    return R_dToReturn;
}

void Lmuur::addSoilprofiles(Soilprofile Right, Soilprofile Left) {
    leftProfile = Left;
    rightProfile = Right;
    std::cout << "soil profiles added" << std::endl;
}

Lmuur::~Lmuur() {}
