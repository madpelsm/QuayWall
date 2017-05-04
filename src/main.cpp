#include <Lmuur.h>
#include <iostream>
int main() {
    // Lmuur(double mHm, double mHv, double mBl, double mBm, double mBr, double
    // gamma, double toe);
    Lmuur L1 = Lmuur(21, 3, 2, 2, 28, 25);
    std::cout << L1.mForce << "," << L1.mx << "," << L1.my << std::endl;
    return 0;
}
