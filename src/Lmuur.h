#pragma once
#include <forcevector.h>
#include <soilprofile.h>
#include <algorithm>
#include <fstream>
#include <glm/glm.hpp>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

class Lmuur {
   public:
    // intermediate results
    double kastnerH = 0;
    // end of intermediate result storage
    // begin unityCheck quantities
    double R_d = 0, RH_d = 0, momentST = 0, momentDST = 0;
    // end unity check quantities
    double gamma_water = 9.81, mq = 30;
    double mExcentricity = 0;
    double mCohesion = 2;  // kNm^2/m
    double mHm, mHv, mBl, mBm, mBr, mBz, mxI, mxII, myI, myII, mx, my, mAI,
        mAII, gamma, mA, mtoe, mxtoe, mytoe, mAtoe;
    double mForce, mSoilHeightDifference;
    // FORCES
    ForceVector mOwnWeight;  // real weight, no buoyance effect
    ForceVector mBuoyantForce;
    ForceVector mlhsWaterPressure, mrhsWaterPressure;
    ForceVector mResultingForce;
    std::vector<ForceVector> mSoilWedgeWeight;  // effective weight
    std::vector<ForceVector> mActiveSoilPressure;
    std::vector<ForceVector> mBoussinesqResultant;
    std::vector<ForceVector> mPassiveSoilPressure;

    // SOILPROFILES
    Soilprofile rightProfile;
    Soilprofile leftProfile;

    Lmuur(double mHm, double mHv, double mBl, double mBm, double mBr,
          double gamma, double toe = 0);
    ~Lmuur();
    void calculateProperties();
    void calculateActiveSoilPressures(Soilprofile& soilprofile, double side,
                                      double footwidth, int safetyCase = 0);
    // if 0, mPhiA, 1 -> mPhiB, 2-> mPhiC
    void calculateSoilWedgeWeight(Soilprofile& soilprofile, double base,
                                  double H, double offSet, double side);
    // side means: if it is on the left side of the wall -1, else 1
    void calculatePassiveSoilPressure(Soilprofile& soilprofile, double side,
                                      double footwidth, int safetyCase = 0);
    void calculateWaterPressures();
    void calculateBoussinesqLoads();  // for general purposes a class for this
                                      // load should handle this itself, but
                                      // since it is extremely simple in our
                                      // case. it is faster to implement the
                                      // result right away.
    void calculateBuoyancy();
    void setSolidHeightDifference(double height);
    void calculateResultingForce();
    void calculateTiltMomentAtFoot(int safeteCase = 0);
    void calculateAll(int safetyCase = 0);
    void calculateExcentricity();
    void writeToCSV(std::string file_name);
    void makeUnityChecks();
    double calculateR_d(double phi_d, Soilprofile& soilprofile, double depth,
                        double effectiveCohesion_safetyF);
    void addSoilprofiles(Soilprofile Right, Soilprofile Left);
    double squareSurface(double height, double width);
    double triangleSurface(double height, double base);
};
