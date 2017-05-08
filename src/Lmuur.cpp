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
    double lambda_a = 0;
    double yTemp =
        (mHm - correctionHeight) - footwidth * tan(M_PI / 4.0 + phi_d);
    for (size_t i = 0; i < soilprofile.mSoillayers.size(); ++i) {
        // dependencies for lambda_a
        psi = phi_d * (2.0 / 3.0);
        alpha = M_PI / 4.0 - phi / 2.0;
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
            if (h < lower) {
                alpha = M_PI / 4.0 - phi_d / 2.0;
                psi = phi_d;
                lambda_a = soilprofile.getLambda_a(phi_d, alpha, phi_d, 0);
                Qa = soilprofile.getQa(lambda_a, soilprofile.mSoillayers[i], h,
                                       lower, h, alpha);
                Qa.mForce.y = Qa.mForce.x * abs(sin(alpha + psi));
                Qa.mForce.x = -forceDirectionCorrection * Qa.mForce.x *
                              abs(cos(alpha + psi));
                Qa.mPoE.x = wallBorder +
                            forceDirectionCorrection * Qa.mPoE.y /
                                tan(M_PI / 4.0 + phi_d / 2.0);
                Qa.mPoE.y += h;
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
                (mHm - correctionHeight + mHv) &&
            soilprofile.mSoillayers[i].mLowerBounds > mHm - correctionHeight) {
            lower = std::min(soilprofile.mSoillayers[i].mLowerBounds,
                             mHm + mHv - correctionHeight);
            if (lower > mHm - correctionHeight) {
                alpha = 0;
                lambda_a = soilprofile.getLambda_a(phi_d, 0, psi, 0);
                Qa = soilprofile.getQa(lambda_a, soilprofile.mSoillayers[i], h,
                                       lower, h, 0);
                Qa.mForce.y = Qa.mForce.x * sin(alpha + psi);
                Qa.mForce.x =
                    -forceDirectionCorrection * Qa.mForce.x * cos(alpha + psi);
                Qa.mPoE.x = wallBorder + footwidth * forceDirectionCorrection;
                Qa.mPoE.y += h;
                mActiveSoilPressure.push_back(Qa);
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
    double d = mHm + mHv - waterHeightLeft;
    double kastnerHcorr =
        d + (waterHeightLeft / (1 + cbrt((1 + waterHeightLeft / d))));
    kastnerH = kastnerHcorr;
    mrhsWaterPressure = ForceVector(
        glm::vec2(-gamma_water * kastnerHcorr * (mHm + mHv) * 0.5, 0),
        glm::vec2(mBm, (2.0 / 3.0) * (mHm + mHv)));
    mlhsWaterPressure =
        ForceVector(glm::vec2(gamma_water * kastnerHcorr * d * 0.5, 0),
                    glm::vec2(0, waterHeightLeft + (2.0 / 3.0) * d));
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

void Lmuur::setSolidHeightDifference(double height) {
    mSoilHeightDifference = height;
}

void Lmuur::calculateResultingForce() {
    // Calculate the resulting force
    double sumFx = 0, sumFy = 0, poex = 0, poey = 0, mOx = 0, mOy = 0;
    // mOx = sum(Fx*y), mOy = sum(Fy*x)
    for (size_t i = 0; i < mSoilWedgeWeight.size(); ++i) {
        sumFx += mSoilWedgeWeight[i].mForce.x;
        sumFy += mSoilWedgeWeight[i].mForce.y;
        mOy += mSoilWedgeWeight[i].mForce.y * mSoilWedgeWeight[i].mPoE.x;
        mOx += mSoilWedgeWeight[i].mForce.x * mSoilWedgeWeight[i].mPoE.y;
    }
    for (size_t i = 0; i < mActiveSoilPressure.size(); ++i) {
        sumFx += mActiveSoilPressure[i].mForce.x;
        sumFy += mActiveSoilPressure[i].mForce.y;
        mOy += mActiveSoilPressure[i].mForce.y * mActiveSoilPressure[i].mPoE.x;
        mOx += mActiveSoilPressure[i].mForce.x * mActiveSoilPressure[i].mPoE.y;
    }
    for (size_t i = 0; i < mBoussinesqResultant.size(); ++i) {
        sumFx += mBoussinesqResultant[i].mForce.x;
        sumFy += mBoussinesqResultant[i].mForce.y;
        mOy +=
            mBoussinesqResultant[i].mForce.y * mBoussinesqResultant[i].mPoE.x;
        mOx +=
            mBoussinesqResultant[i].mForce.x * mBoussinesqResultant[i].mPoE.y;
    }
    // Passive not yet implemented
    for (size_t i = 0; i < mPassiveSoilPressure.size() && false; ++i) {
        sumFx += mPassiveSoilPressure[i].mForce.x;
        sumFy += mPassiveSoilPressure[i].mForce.y;
        mOy +=
            mPassiveSoilPressure[i].mForce.y * mPassiveSoilPressure[i].mPoE.x;
        mOx +=
            mPassiveSoilPressure[i].mForce.x * mPassiveSoilPressure[i].mPoE.y;
    }
    // own weight
    sumFx += mOwnWeight.mForce.x;
    sumFy += mOwnWeight.mForce.y;
    mOy += mOwnWeight.mForce.y * mOwnWeight.mPoE.x;
    mOx += mOwnWeight.mForce.x * mOwnWeight.mPoE.y;
    // buoyant
    sumFx += mBuoyantForce.mForce.x;
    sumFy += mBuoyantForce.mForce.y;
    mOy += mBuoyantForce.mForce.y * mBuoyantForce.mPoE.x;
    mOx += mBuoyantForce.mForce.x * mBuoyantForce.mPoE.y;
    // rhsWaterPressure
    sumFx += mrhsWaterPressure.mForce.x;
    sumFy += mrhsWaterPressure.mForce.y;
    mOy += mrhsWaterPressure.mForce.y * mrhsWaterPressure.mPoE.x;
    mOx += mrhsWaterPressure.mForce.x * mrhsWaterPressure.mPoE.y;
    // lhsWaterPressure
    sumFx += mlhsWaterPressure.mForce.x;
    sumFy += mlhsWaterPressure.mForce.y;
    mOy += mlhsWaterPressure.mForce.y * mlhsWaterPressure.mPoE.x;
    mOx += mlhsWaterPressure.mForce.x * mlhsWaterPressure.mPoE.y;

    // mOx = sum(Fx*y), mOy = sum(Fy*x)

    mResultingForce = ForceVector(glm::vec2(sumFx, sumFy),
                                  glm::vec2(mOy / sumFy, mOx / sumFx));
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
    calculateResultingForce();

    calculateExcentricity();

    makeUnityChecks();
}

void Lmuur::calculateExcentricity() {
    mExcentricity = mResultingForce.mPoE.x - mxII;
}

void Lmuur::calculateActiveSoilPressureLeft() {}
double Lmuur::squareSurface(double height, double width) {
    return height * width;
}

double Lmuur::triangleSurface(double height, double base) {
    return height * base * 0.5;
}

void Lmuur::writeToCSV() {
    std::ofstream file;
    file.open("output.csv");
    file << "L vormige kaaimuur\n";
    file << "Eigenschappen L-wand,[m],,[m]\n";
    file << "Hoogte vanaf de voet," << mHm << ",Breedte verticale wand," << mBm
         << "\n";
    file << "Breedte linker voet," << mBl << ",Breedte rechter voet," << mBr
         << "\n";
    file << "Breedte verticale wand," << mBm << ",Dikte funderingszool," << mHv
         << "\n";
    file << "Opmerkingen\n";
    file << "Het assenstelsel gaat met positieve richting naar "
            "beneden\nDeze y as valt samen met de linkerzijde van de verticale "
            "wand\nDe x as "
            "gaat positief naar rechts\n samenvallend met het horizontale "
            "maaiveld.\n\n";
    file << "Aangrijpende krachten. \n";
    file << "Kracht,Fx[kN],Fy[kN],x[m],y[m]\n";
    file << "Eigengewicht.\n";
    file << "Eigen," << mOwnWeight.mForce.x << "," << mOwnWeight.mForce.y << ","
         << mOwnWeight.mPoE.x << "," << mOwnWeight.mPoE.y << "\n";

    file << "Opwaartse stuwkracht van het beton.\n";
    file << "Archimedes," << mBuoyantForce.mForce.x << ","
         << mBuoyantForce.mForce.y << "," << mBuoyantForce.mPoE.x << ","
         << mBuoyantForce.mPoE.y << "\n";

    file << "Waterdrukken mbv Kastner.\n";
    file << "Rechterkant," << mrhsWaterPressure.mForce.x << ","
         << mrhsWaterPressure.mForce.y << "," << mrhsWaterPressure.mPoE.x << ","
         << mrhsWaterPressure.mPoE.y << "\n";
    file << "Linkerkant," << mlhsWaterPressure.mForce.x << ","
         << mlhsWaterPressure.mForce.y << "," << mlhsWaterPressure.mPoE.x << ","
         << mlhsWaterPressure.mPoE.y << "\n";
    file << "Actieve gronddrukken.\n";
    for (size_t i = 0; i < mActiveSoilPressure.size(); ++i) {
        file << "Qa" << i << "," << mActiveSoilPressure[i].mForce.x << ","
             << mActiveSoilPressure[i].mForce.y << ","
             << mActiveSoilPressure[i].mPoE.x << ","
             << mActiveSoilPressure[i].mPoE.y << "\n";
    }
    file << "GrondGewicht\n";
    for (size_t i = 0; i < mSoilWedgeWeight.size(); ++i) {
        file << "G" << i << "," << mSoilWedgeWeight[i].mForce.x << ","
             << mSoilWedgeWeight[i].mForce.y << ","
             << mSoilWedgeWeight[i].mPoE.x << "," << mSoilWedgeWeight[i].mPoE.y
             << "\n";
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
    file << "Fk," << mResultingForce.mForce.x << "," << mResultingForce.mForce.y
         << "," << mResultingForce.mPoE.x << "," << mResultingForce.mPoE.y
         << "\n";
    file << "UNITY CHECKS\n";
    file << "type,R_d,E_d,veiligheid\n";
    file << "Evewichtsdraagvermogen," << R_d << ","
         << mResultingForce.mForce.y * 100 << ","
         << R_d / (100 * mResultingForce.mForce.y) << ",excentriciteit:,"
         << mExcentricity << "\n";
    file << "Schuiven," << RH_d << "," << mResultingForce.mForce.x << ","
         << abs(RH_d / mResultingForce.mForce.x) << "\n";

    file.close();
    std::cout << "Write complete!" << std::endl;
}

void Lmuur::makeUnityChecks() {
    R_d = calculateR_d(leftProfile.mSoillayers[0].mPhiA, leftProfile,
                       mHm + mHv - mSoilHeightDifference, 1);
    RH_d = mResultingForce.mForce.y * tan(leftProfile.mSoillayers[0].mPhiA);
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
             0.7 * abs(mResultingForce.mForce.x) /
                 (abs(mResultingForce.mForce.y) + B * L / (tan(phi_d) * c_d))),
            3);
    double i_c = (i_q * N_q - 1) / (N_q - 1);
    double i_g =
        pow((1 -
             abs(mResultingForce.mForce.x) /
                 (abs(mResultingForce.mForce.y) + B * L / (tan(phi_d) * c_d))),
            3);
    double R_dToReturn =
        B * L * (d_q * s_q * N_q * p_t * i_q + d_c * s_c * N_c * i_c * c_d +
                 s_g * N_g * g_k * i_g * B / 2.0);
    // std::cout << "q_u = d_q*s_q*N_q*p_t+d_c*s_c*N_c*c+s_g*N_g*g_k*B/2\n"
    //           << d_q << "*" << s_q << "*" << N_q << "*" << p_t << "+" << d_c
    //           << "*" << s_c << "*" << N_c << "*" << c << "/" <<
    //           effectiveCohesion_safetyF  << "+" << s_g << "*"
    //           << N_g << "*" << g_k << "*" << B << "/" << 2.0 << "=" << _qu
    //           << "\n"
    //           << std::endl;
    return R_dToReturn;
}

void Lmuur::addSoilprofiles(Soilprofile Right, Soilprofile Left) {
    leftProfile = Left;
    rightProfile = Right;
}

Lmuur::~Lmuur() {}
