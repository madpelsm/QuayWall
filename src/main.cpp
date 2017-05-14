#include <Lmuur.h>
#include <soillayer.h>
#include <soilprofile.h>
#include <iostream>
#include <vector>

int main() {
    // Lmuur(double mHm, double mHv, double mBl, double mBm, double mBr, double
    // gamma, double toe);
    // mHm, mHv, mBl,mBm,mBr,gamma,toe
    double verticalHeight = 26, footheight = 2, leftWidth = 2.5,
           middleWidth = 2, rightWidth = 16.5, gamma = 25, toe = 0.5,
           toeWidth = 1;
    Lmuur Lm = Lmuur(verticalHeight, footheight, leftWidth, middleWidth,
                     rightWidth, gamma, toe, toeWidth);
    Soillayer R1 = Soillayer(0.0, 3.0, 16.0, 0.39, 15.0);
    Soillayer R2 = Soillayer(3.0, 28.0, 18.0, 0.35, 25.0);
    Soillayer R3 = Soillayer(28.0, 100.0, 18.0, 0.35, 35.0);
    Soillayer L1 = Soillayer(0.0, 12.0, 18.0, 0.35, 25.0);
    Soillayer L2 = Soillayer(12.0, 100.0, 18.0, 0.35, 35.0);

    std::vector<Soillayer> rhs{R1, R2, R3};
    std::vector<Soillayer> lhs{L1, L2};

    Soilprofile Right = Soilprofile(rhs, 3.0);
    Soilprofile Left = Soilprofile(lhs, -12.0);

    Lm.addSoilprofiles(Right, Left);

    Lm.setSolidHeightDifference(21);

    Lm.calculateAll(0);
    Lm.writeToCSV("CaseAfin.csv");

    Lm.calculateAll(1);
    Lm.writeToCSV("CaseBfin.csv");

    Lm.calculateAll(2);
    Lm.writeToCSV("CaseCfin.csv");
    std::cout << "kaaimuur met hoogte vanaf de voet: " << Lm.mHm
              << "[m] en dikte " << Lm.mBm << "[m]"
              << "\neen funderingszooldikte van " << Lm.mHv << "[m]"
              << "\nlinkervoetbreedte: " << Lm.mBl << "[m]"
              << "\nDe rechter voet is " << Lm.mBr << "[m]" << std::endl;
    std::cout << "Force with magintude (Fx,Fy): (" << Lm.mResultingR_d.mForce.x
              << "," << Lm.mResultingR_d.mForce.y
              << ")\n and point of engagement (x,y): ("
              << Lm.mResultingR_d.mPoE.x << "," << Lm.mResultingR_d.mPoE.y
              << ")" << std::endl;

    return 0;
}
