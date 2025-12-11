// File: Transfer/src/Entities/PhysicsStructures.cpp

#include "Entities/PhysicsStructures.h"

std::ostream& operator<<(std::ostream& os, const Vector2D& Vec) {
    // 3. Insert the desired formatting into the output stream (os)
    // os << "Point(" << p.label << ") - X: " << p.x << ", Y: " << p.y;
    os << "{ "<< Vec.x_val << ", " << Vec.y_val << " }";
    
    // 4. Return the stream reference to allow chaining
    return os;
}