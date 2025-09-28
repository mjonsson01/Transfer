//#include "UtilityHelpers.h"
//
//Color GetColorForMass(double mass) {
//    if (mass <= BodyMassOptions::Light) {
//        return ColorLibrary::White;
//    }
//    else if (mass <= BodyMassOptions::Medium) {
//        return ColorLibrary::Green;
//    }
//    else {
//        return ColorLibrary::Blue;
//    }
//};
//
//void updatePhysicsFrame(std::vector<GravitationalBody>& bodies) {
//    for (auto& body : bodies) {
//        body.x += body.netVelocity.v_x * TIME_STEP;
//        body.y += body.netVelocity.v_y * TIME_STEP;
//    }
//}