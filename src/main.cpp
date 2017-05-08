#include <Lmuur.h>
#include <soillayer.h>
#include <soilprofile.h>
#include <iostream>
#include <vector>

int main() {
    // Lmuur(double mHm, double mHv, double mBl, double mBm, double mBr, double
    // gamma, double toe);
    // mHm, mHv, mBl,mBm,mBr,gamma,toe
    double verticalHeight = 25, footheight = 3, leftWidth = 2, middleWidth = 2,
           rightWidth = 20, gamma = 25, toe = 0;
    Lmuur Lm = Lmuur(verticalHeight, footheight, leftWidth, middleWidth,
                     rightWidth, gamma, toe);
    std::cout << Lm.mOwnWeight.mForce.y << "," << Lm.mOwnWeight.mPoE.x << ","
              << Lm.mOwnWeight.mPoE.y << std::endl;
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

    Lm.calculateAll();

    std::cout << "kaaimuur met hoogte vanaf de voet: " << Lm.mHm
              << "[m] en dikte " << Lm.mBm << "[m]"
              << "\neen funderingszooldikte van " << Lm.mHv << "[m]"
              << "\nlinkervoetbreedte: " << Lm.mBl << "[m]"
              << "\nDe rechter voet is " << Lm.mBr << "[m]" << std::endl;
    std::cout << "Force with magintude (Fx,Fy): ("
              << Lm.mResultingForce.mForce.x << ","
              << Lm.mResultingForce.mForce.y
              << ")\n and point of engagement (x,y): ("
              << Lm.mResultingForce.mPoE.x << "," << Lm.mResultingForce.mPoE.y
              << ")" << std::endl;
    Lm.writeToCSV();

    std::cout << Lm.kastnerH << std::endl;

    return 0;
}
