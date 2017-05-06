#include <iostream>

#include <Lmuur.h>
#include <soillayer.h>
#include <soilprofile.h>
#include <vector>
int main() {
    // Lmuur(double mHm, double mHv, double mBl, double mBm, double mBr, double
    // gamma, double toe);
    Lmuur L1 = Lmuur(21, 3, 2, 2, 28, 25);
    std::cout << L1.ForceVector.mForce.y << "," << L1.ForceVector.mPoE.x << ","
              << L1.ForceVector.mPoE.y << std::endl;
    Soillayer Layer1 = Soillayer(0, 3, 16, 0.39, 15);
    std::vector<Soillayer> v1;
    v1.push_back(Layer1);
    std::cout << Layer1.mEffectiveWeight << "," << Layer1.mWetWeight
              << std::endl;
    Soilprofile p1 = Soilprofile(v1, 0);
    std::cout << p1.getEffectiveSoilePressure(2) << std::endl;

    return 0;
}
