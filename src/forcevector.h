#pragma once
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <sstream>
#include <string>
class ForceVector {
   public:
    glm::vec2 mForce;
    glm::vec2 mPoE;  // point of engagement

    ForceVector();
    ForceVector(glm::vec2 Force, glm::vec2 PoE);
    std::string toString();
    ~ForceVector();
};
