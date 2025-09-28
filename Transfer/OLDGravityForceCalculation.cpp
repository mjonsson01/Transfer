//#include "GravityForceCalculation.h"
//// if km scale masses differently
//// mass1 and mass2 in kg, radius1 and radius2 in kilometers, distance in meters, where distance is the distance between the centers of the two bodies
//double calculateGravitationalForce(double mass1, double mass2, int radius1, int radius2, int distance)
//{
//	if (distance <= radius1 + radius2) {
//		return 0; // Avoid division by zero
//	}
//	return (GRAVITATIONAL_CONSTANT * mass1 * mass2) / (distance * distance);
//};
//
//// Calculate effects on bodies as a result of gravitational forces
//
//int distanceBetweenBodies(const GravitationalBody& body1, const GravitationalBody& body2) {
//	int dx = body2.x - body1.x;
//	int dy = body2.y - body1.y;
//	return static_cast<int>(sqrt(dx * dx + dy * dy)); // return rounded int distance
//}
//
//int signedDisplacementBetweenBodies(const GravitationalBody& body1, const GravitationalBody& body2, char axis) {
//	if (axis == 'x') {
//		return body2.x - body1.x;
//	}
//	else if (axis == 'y') {
//		return body2.y - body1.y;
//	}
//	else {
//		std::cerr << "Error: Invalid axis specified. Use 'x' or 'y'." << std::endl;
//		return 0;
//	}
//}
//
//void applyGravitationalForceEffects(GravitationalBody& body1, GravitationalBody& body2)
//{
//	// Calculate distance between bodies
//	int displacement_X = signedDisplacementBetweenBodies(body1, body2, 'x');
//	int displacement_Y = signedDisplacementBetweenBodies(body1, body2, 'y');
//	//double measuredForce = calculateGravitationalForce(body1.mass, body2.mass, body1.radius, body2.radius, distance);
//	// Will need to call all forces acting on a body and sum them to get net force.
//	// Apply the force. We will assume that we only have two bodies for now. We will limit the calculation to 62.5 Hz as well.
//	// F = m * a  -> a = F / m. We will apply the force over a time interval of 0.016 seconds (62.5 Hz). So the impulse applied is F*dt = m*dv -> dv = F*dt/m
//
//	// in order to correctly scale these, we must have the components scaled correctly x/sqrt(x^2+y^2) = cos(theta), y/sqrt(x^2+y^2) = sin(theta). So F cos(theta) = F * x/sqrt(x^2+y^2) and F sin(theta) = F * y/sqrt(x^2+y^2)
//
//	int distance = distanceBetweenBodies(body1, body2);
//	double measuredForce = calculateGravitationalForce(body1.mass, body2.mass, body1.radius, body2.radius, distance);
//	double force_x = measuredForce* displacement_X / sqrt(distance);
//	double force_y = measuredForce* displacement_Y / sqrt(distance);
//	// apply to body 1
//	body1.netForce.f_x += force_x;
//	body1.netForce.f_y += force_y;
//	body1.netVelocity.v_x += ( force_x * TIME_STEP) / body1.mass;
//	body1.netVelocity.v_y += ( force_y * TIME_STEP) / body1.mass;
//
//	// apply to body 2
//	body2.netForce.f_x -= force_x; // equal and opposite force
//	body2.netForce.f_y -= force_y;
//	body2.netVelocity.v_x -= (force_x* TIME_STEP) / body2.mass;
//	body2.netVelocity.v_y -= (force_y* TIME_STEP) / body2.mass;	
//}
//
//void applyGravityToSystem(std::vector<GravitationalBody>& bodies)
//{
//	for (size_t i = 0; i < bodies.size(); ++i) {
//		for (size_t j = i + 1; j < bodies.size(); ++j) {
//			applyGravitationalForceEffects(bodies[i], bodies[j]);
//		}
//	}
//}