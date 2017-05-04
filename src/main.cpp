#include <iostream>

#include <Lmuur.h>
#include <soillayer.h>
int main() {
    // Lmuur(double mHm, double mHv, double mBl, double mBm, double mBr, double
    // gamma, double toe);
    Lmuur L1 = Lmuur(21, 3, 2, 2, 28, 25);
    std::cout << L1.mForce << "," << L1.mx << "," << L1.my << std::endl;
    Soillayer Layer1 = Soillayer(0,3,16,0.39,15);
    std::cout<<Layer1.mPhiA<<std::endl;
    return 0;
}
