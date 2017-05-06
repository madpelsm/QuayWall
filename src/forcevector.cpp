#include <forcevector.h>

ForceVector::ForceVector(glm::vec2 Force, glm::vec2 PoE) : mForce(Force), mPoE(PoE) {}
ForceVector::ForceVector() {
    mForce = glm::vec2(0, 0);
    mPoE = glm::vec2(0, 0);
}
ForceVector::~ForceVector() {}
