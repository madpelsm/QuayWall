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
    double b_a = 0;
    // end of intermediate result storage
    // begin unityCheck quantities
    double R_d = 0, RH_d = 0, momentST = 0, momentDST = 0;
    // end unity check quantities
    bool mexcentricitycalculated = false;
    double gamma_water = 9.81, mq = 30;
    std::vector<double> mSoilWedgeHeight;
    double mExcentricity = 0;
    double mCohesion = 2;  // kNm^2/m
    double mHm, mHv, mBl, mBm, mBr, mBz, mxI, mxII, myI, myII, mx, my, mAI,
        mAII, gamma, mA, mtoe, mxtoe, mytoe, mAtoe, mtoewidth;
    double mForce, mSoilHeightDifference;
    // FORCES
    ForceVector mOwnWeight;  // real weight, no buoyance effect
    ForceVector mBuoyantForce;
    ForceVector mlhsWaterPressure, mrhsWaterPressure;
    ForceVector mResultingR_d;
    ForceVector mResultingR_dSchuiven;
    std::vector<ForceVector> mSoilWedgeWeight;  // effective weight
    std::vector<ForceVector> mActiveSoilPressure;
    std::vector<ForceVector> mBoussinesqResultant;
    std::vector<ForceVector> mPassiveSoilPressure;

    // SOILPROFILES
    Soilprofile rightProfile;
    Soilprofile leftProfile;

    // SAFETY FACTORS
    std::vector<double> mSafetyGSup = {1.1, 1.35, 1};
    std::vector<double> mSafetyGInf = {0.9, 1, 1};
    std::vector<double> mSafetyQ = {1.5, 1.5, 1.3, 0};
    std::vector<double> mSafetyC = {1.25,1,1.25};

    Lmuur(double mHm, double mHv, double mBl, double mBm, double mBr,
          double gamma, double toe = 0, double toewidth = 0);
    ~Lmuur();
    void clearForces();
    void calculateProperties();
    void calculateActiveSoilPressures(Soilprofile& soilprofile, double side,
                                      double footwidth, int safetyCase = 0);
    // if 0, mPhiA, 1 -> mPhiB, 2-> mPhiC
    void calculateSoilWedgeWeight(Soilprofile& soilprofile, double base,
                                  double H, double offSet, double side,
                                  int safetyCase = 0);
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
    bool isBeneficial(ForceVector& _forcevector, int _failuremode,
                      int _safetycase);
    // failure mode 0: grondbreuk(evenwichtsdraagvermogen)
    // failure mode 1: schuiven
    ForceVector calculateResultingForce(int failureMode, int safetyCase = 0);
    // directions: y = 1, x=2 (negative for oposite direction)

    void calculateTiltMomentAtFoot(int safeteCase = 0);
    void calculateAll(int safetyCase = 0);
    void calculateExcentricity();
    void writeToCSV(std::string file_name);
    void makeUnityChecks(int safetyCase = 0);
    double getPhiAtConstructionFoot(int safetyCase = 0);
    double calculateR_d(double phi_d, Soilprofile& soilprofile, double depth,
                        double effectiveCohesion_safetyF);
    void addSoilprofiles(Soilprofile Right, Soilprofile Left);
    double squareSurface(double height, double width);
    double triangleSurface(double height, double base);
};
