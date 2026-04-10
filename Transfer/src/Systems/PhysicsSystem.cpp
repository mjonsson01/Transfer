// File: Transfer/src/Systems/PhysicsSystem.cpp

#include "PhysicsSystem.h"

PhysicsSystem::PhysicsSystem() {
  // Initialize physics system variables if needed
}

PhysicsSystem::~PhysicsSystem() {}

// --------- ANONYMOUS NAMESPACE HELPER STRUCTS --------- //
struct CollisionInfo {
  double distance;
  Vector2D unitNormalVector;       // unit normal (from bodyA to bodyB)
  Vector2D relativeVelocityVector; // vB - vA
  double normalSpeed;              // signed speed along normal vector
  double absNormalSpeed; // abs value of signed speed along normal vector
};

static inline CollisionInfo getCollisionInfo(const GravitationalBody& a,
                                             const GravitationalBody& b) {
  Vector2D r_vector = b.position - a.position;
  double distance = r_vector.magnitude();
  Vector2D unit_normal_vector =
      (distance > 1e-8) ? (r_vector / distance) : Vector2D(1.0, 0.0);

  Vector2D relative_velocity_vector = b.velocity - a.velocity;
  double normal_speed = relative_velocity_vector.dot(unit_normal_vector);
  double abs_normal_speed = std::abs(normal_speed);
  return {distance, unit_normal_vector, relative_velocity_vector, normal_speed,
          abs_normal_speed};
}

struct GravitationalBodyPair {
  GravitationalBody* heavy;
  GravitationalBody* light;
  double ratio; // heavy/light (>= 1)
  bool equal;
};

static inline GravitationalBodyPair pickMassPair(GravitationalBody& a,
                                                 GravitationalBody& b) {
  if (a.mass == b.mass)
    return {&a, &b, 1.0, true};
  if (a.mass > b.mass)
    return {&a, &b, a.mass / b.mass, false};
  return {&b, &a, b.mass / a.mass, false};
}

// --------- SYSTEM-LEVEL METHOD --------- //

void PhysicsSystem::UpdateSystemFrame(GameState& gameState, UIState& UIState) {

  // Check for dirty input that requires immediate updates to the Physics System
  updateGravBodyInstantiations(gameState, UIState);

  // Check if the Physics System is paused for early exit
  if (UIState.getInputState().isPhysicsPaused)
    return;

  // ------------------------ MAIN PHYSICS LOOP ------------------------ //
  // Integrate first half step
  integrateForwardsPhase1(gameState);

  // Update all the gravity
  updateGravityForSystem(gameState);

  // Integrate second half step
  integrateForwardsPhase2(gameState);

  // Check total energy
  // calculateTotalEnergy(gameState);

  // Handle all the collisions
  handleCollisions(gameState);

  // Advance Collision artifact lifetime
  auto& macro_bodies = gameState.getMacroBodiesMutable();
  for (auto& body : macro_bodies) {
    if (body.isTransient) {
      if (body.lifetime > 0) {
        --body.lifetime; // decrement lifetime
      }

      if (body.lifetime == 0) {
        body.isMarkedForDeletion = true; // flag for removal
      }
    }
  }
  // Clean up the marked for deletion particles and macro bodies
  cleanupParticles(gameState);
  cleanupMacroBodies(gameState);
}

// --------- CLEANUP METHOD --------- //

void PhysicsSystem::CleanUp() {
  // Any necessary cleanup code for the physics system
}

// --------- BODY INSTANTIATION METHOD --------- //

void PhysicsSystem::updateGravBodyInstantiations(GameState& gameState,
                                                 UIState& UIState) {
  InputState& input_state = UIState.getMutableInputState();
  if (input_state.instantiateDirty) {
    if (input_state.isCreatingMacro &&
        macroLimiter.canSpawn(1 / UIState.getFPS())) {
      createMacroBody(gameState, input_state);
      input_state.resetTransientFlags();
    } else if (input_state.isCreatingParticle &&
               particleLimiter.canSpawn(1 / UIState.getFPS())) {
      createParticle(gameState, input_state);
      input_state.resetTransientFlags();
    } else if (input_state.isCreatingParticleCluster &&
               clusterLimiter.canSpawn(1 / UIState.getFPS())) {
      createParticleCluster(gameState, input_state);
      input_state.resetTransientFlags();
    }
  } else {
    macroLimiter.reset();
    particleLimiter.reset();
    clusterLimiter.reset();
  }
  // Check if all Gravitational Bodies are supposed to be wiped
  if (input_state.clearAll) {
    gameState.getMacroBodiesMutable().clear();
    gameState.getParticlesMutable().clear();
    input_state.clearAll = false;
  }
}

// --------- GRAVITY METHODS --------- //

void PhysicsSystem::updateGravityForSystem(GameState& gameState) {
  // Get Particles and Macro Bodies
  auto& particles = gameState.getParticlesMutable();
  int num_particles = particles.size();
  auto& macro_bodies = gameState.getMacroBodiesMutable();
  int num_macro_bodies = macro_bodies.size();

  // Update gravity for macro-macro interactions
  for (size_t i = 0; i < num_macro_bodies; i++) {
    for (size_t j = i + 1; j < num_macro_bodies; j++) {
      calculateGravity(macro_bodies[i], macro_bodies[j]);
    }
  }
  // Update gravity for particle-macro interactions
  for (size_t i = 0; i < num_particles; i++) {
    for (size_t j = 0; j < num_macro_bodies; j++) {
      calculateGravity(particles[i], macro_bodies[j]);
    }
  }
}

void PhysicsSystem::calculateGravity(GravitationalBody& body1,
                                     GravitationalBody& body2) {
  if (body1.mass == 0 || body2.mass == 0)
    return;
  static const double G = GRAVITATIONAL_CONSTANT;

  // Define softening epsilon
  double epsilon = (body1.radius + body2.radius) / 2.0; // Simple average radius
  double epsilon_squared = epsilon * epsilon;

  // 1. Calculate direction Vector
  Vector2D direction_vector = body2.position - body1.position;

  // 2. Calculate Distance Squared
  double r_squared = direction_vector.square_magnitude();

  // 3. Calculate the Denominator Term (r^2 + epsilon^2)^(3/2)
  double denominator_1 = sqrt(r_squared + epsilon_squared);
  double denominator_3 = denominator_1 * denominator_1 * denominator_1;

  // 4. Calculate Coefficient (F = direction vector * G * m1 * m2 / Denominator)
  double coefficient = (G * body1.mass * body2.mass) / denominator_3;

  // Force vector is C * direction_vector (r)
  Vector2D force = direction_vector * coefficient;

  // 6. Apply Forces (Newton's Third Law)
  body1.netForce += force;
  body2.netForce -= force;
}

// --------- INTEGRATION METHODS --------- //

void PhysicsSystem::integrateForwardsPhase1(GameState& gameState) {
  // Get particles and macro bodies
  auto& particles = gameState.getParticlesMutable();
  auto& macro_bodies = gameState.getMacroBodiesMutable();

  // Phase 1 integrate particles
  for (auto& particle : particles) {
    if (abs(particle.mass) <= EPSILON)
      continue;
    Vector2D new_acceleration = particle.isStatic
                                    ? Vector2D(0.0, 0.0)
                                    : particle.netForce / particle.mass;
    // Calculate current acceleration (a_t)
    // check for 0 mass?

    // Step 1. Kick 1: Update velocity by half the acceleration over one PHYSICS
    // TIME STEP (dt)
    particle.velocity += new_acceleration * (PHYSICS_TIME_STEP / 2.0);

    // Step 2. Drift: Update position using the half-step velocity
    particle.previousPosition = particle.position;
    particle.position += particle.velocity * PHYSICS_TIME_STEP;

    // Reset netForce for the next step's force accumulation
    particle.netForce = Vector2D(0.0, 0.0);
  }

  // Phase 1 integrate macro bodies
  for (auto& body : macro_bodies) {
    if (body.isStatic || abs(body.mass) <= EPSILON)
      continue;

    // Calculate current acceleration (a_t)
    //
    // std::cout<<body.mass<<std::endl;
    Vector2D current_acceleration = body.netForce / body.mass;
    // std::cout<<current_acceleration<<std::endl;
    // Step 1. Kick 1: Update velocity by half the acceleration over one PHYSICS
    // TIME STEP (dt)
    body.velocity += current_acceleration * (PHYSICS_TIME_STEP / 2.0);

    // Step 2. Drift: Update position using the half-step velocity
    body.previousPosition = body.position;
    body.position += body.velocity * PHYSICS_TIME_STEP;

    // Reset netForce for the next step's force accumulation
    body.netForce = Vector2D(0.0, 0.0);
  }
}

void PhysicsSystem::integrateForwardsPhase2(GameState& gameState) {
  // Get particles and macro bodies
  auto& particles = gameState.getParticlesMutable();
  auto& macro_bodies = gameState.getMacroBodiesMutable();

  // Phase 2 integrate particles
  for (auto& particle : particles) {
    if (particle.isStatic || abs(particle.mass) <= EPSILON)
      continue;

    // Calculate new acceleration (a_t + dt) using the newly calculated netForce
    Vector2D new_acceleration = particle.netForce / particle.mass;

    // Step 3. Kick 2: Update velocity by half the acceleration over one PHYSICS
    // TIME STEP (dt)
    particle.velocity += new_acceleration * (PHYSICS_TIME_STEP / 2.0);
  }

  // Phase 2 integrate macro bodies
  for (auto& body : macro_bodies) {
    if (abs(body.mass) <= EPSILON)
      continue;

    // Calculate new acceleration (a_t + dt) using the newly calculated netForce
    Vector2D new_acceleration =
        body.isStatic ? Vector2D(0.0, 0.0) : body.netForce / body.mass;

    // Step 3. Kick 2: Update velocity by half the acceleration over one PHYSICS
    // TIME STEP (dt)
    body.velocity += new_acceleration * (PHYSICS_TIME_STEP / 2.0);
  }
}

// --------- COLLISION METHODS --------- //

void PhysicsSystem::handleCollisions(GameState& gameState) {
  // Get bodies and deal with body-body collisions
  auto& macro_bodies = gameState.getMacroBodiesMutable();
  size_t num_macro_bodies = macro_bodies.size();

  // Loop through all macro body collision pairs
  for (size_t i = 0; i < num_macro_bodies; ++i) {
    auto& body_a = macro_bodies[i];
    if (body_a.isMarkedForDeletion)
      continue; // ignore bodies that will be deleted this frame

    for (size_t j = i + 1; j < num_macro_bodies; ++j) {
      auto& body_b = macro_bodies[j];
      if (body_b.isMarkedForDeletion)
        continue; // ignore bodies that will be deleted this frame

      CollisionInfo info = getCollisionInfo(body_a, body_b);

      if (info.distance >= body_a.radius + body_b.radius)
        continue;

      GravitationalBodyPair pair = pickMassPair(body_a, body_b);
      if (body_a.isBounce || body_b.isBounce) {
        handleElasticCollisions(*pair.light, *pair.heavy);
      } else if (info.absNormalSpeed < MIN_SHATTER_SPEED) {
        if (pair.ratio >= MIN_BODY_BODY_ACCRETION_THRESHOLD_RATIO) {
          handleAccretion(*pair.light, *pair.heavy);
        } else {
          handleElasticCollisions(*pair.light, *pair.heavy);
        }
      } else {

        // as long as the heavier/lighter is within 40%, blow em both to shit.
        // will tweak as needed and place in engine constants when im happy with
        // it note that ratio alwasys >= 1.0 since I force the heavier to be in
        // numerator
        if (pair.ratio <= 1.4) {
          handleDynamicExplosionCollision(body_a, body_b, gameState);
          // Leave these here as baseline.
          // substituteWithParticles(*pair.light, gameState);
          // substituteWithParticles(*pair.heavy, gameState);
        }
        // otherwise just kill the light one
        else {
          substituteWithParticles(*pair.light, gameState);
          // handleDynamicExplosionCollision(body_a, body_b, gameState);
        }
      }
    }
  }

  auto& particles = gameState.getParticlesMutable();
  for (auto& particle : particles) {
    if (particle.isMarkedForDeletion)
      continue;

    for (auto& body : macro_bodies) {
      if (body.isMarkedForDeletion)
        continue;

      auto info = getCollisionInfo(particle, body);
      if (body.isBounce) {
        if (info.distance < particle.radius + body.radius) {
          handleElasticCollisions(particle, body);
          break;
        }
      } else {
        if (info.distance >= particle.radius + body.radius)
          continue;

        if (info.absNormalSpeed > MAX_ACCRETION_COLLISION_SPEED) {
          // handleElasticCollisions(particle, body);
          handleAccretion(particle, body);
        } else {
          handleAccretion(particle, body);
          break; // Don't accrete same particle into multiple bodies this frame
        }
      }
    }
  }
}

void PhysicsSystem::handleElasticCollisions(GravitationalBody& smallerBody,
                                            GravitationalBody& largerBody) {
  Vector2D r_vector = largerBody.position - smallerBody.position;
  double distance = r_vector.magnitude();

  if (distance == 0)
    return;
  // Solution 1
  // if (smallerBody.isBounce || largerBody.isBounce)
  // {
  //     GravitationalBody &bounceBody = smallerBody.isBounce ? smallerBody :
  //     largerBody; GravitationalBody &other = smallerBody.isBounce ?
  //     largerBody : smallerBody;

  //     Vector2D normal = (other.position - bounceBody.position).normalize();
  //     double v_dot_n = other.velocity.dot(normal);

  //     if (v_dot_n < 0) // only reflect if moving inward
  //     {
  //         other.velocity -= normal * (2.0 * v_dot_n);
  //     }

  //     return;
  // }
  // if ((smallerBody.mass == 0 || largerBody.mass == 0))
  //     return;

  // Standalone Soln2
  // if (!smallerBody.isBounce && !largerBody.isBounce && (smallerBody.mass == 0
  // || largerBody.mass == 0))
  //     return;

  // Handle bounce bodies first
  if (smallerBody.isBounce || largerBody.isBounce) {
    GravitationalBody& bounceBody =
        smallerBody.isBounce ? smallerBody : largerBody;
    GravitationalBody& other = smallerBody.isBounce ? largerBody : smallerBody;

    Vector2D normal = (other.position - bounceBody.position).normalize();
    double v_dot_n = other.velocity.dot(normal);

    if (v_dot_n < 0) // only reflect if moving inward
    {
      other.velocity -= normal * (2.0 * v_dot_n);
    }

    return;
  }

  // Skip massless bodies only if neither is a bounce body
  if (smallerBody.mass == 0 || largerBody.mass == 0)
    return;
  Vector2D normal_vector = r_vector / distance; // collision normal
  Vector2D v_smaller = smallerBody.velocity;
  Vector2D v_larger = largerBody.velocity;
  double m_smaller = smallerBody.mass;
  double m_larger = largerBody.mass;
  double v_smaller_n = v_smaller.dot(normal_vector);
  double v_larger_n = v_larger.dot(normal_vector);

  // 1D elastic collision formula (normal direction only)
  double v_smaller_n_new =
      (v_smaller_n * (m_smaller - m_larger) + 2 * m_larger * v_larger_n) /
      (m_smaller + m_larger);
  double v_larger_n_new =
      (v_larger_n * (m_larger - m_smaller) + 2 * m_smaller * v_smaller_n) /
      (m_smaller + m_larger);

  // Change in normal components
  Vector2D v_smaller_change =
      normal_vector * (v_smaller_n_new - v_smaller_n) * ELASTIC_LOSS_FACTOR;
  Vector2D v_larger_change =
      normal_vector * (v_larger_n_new - v_larger_n) * ELASTIC_LOSS_FACTOR;

  // Apply
  // smallerBody.velocity = v_smaller + v_smaller_change;
  // largerBody.velocity = v_larger + v_larger_change;
  if (!smallerBody.isStatic)
    smallerBody.velocity = v_smaller + v_smaller_change;

  if (!largerBody.isStatic)
    largerBody.velocity = v_larger + v_larger_change;

  double overlap = smallerBody.radius + largerBody.radius - distance;

  // Positional correction to disperse overlapping bodies
  if (abs(overlap) > 0) {
    // Small slop to avoid micro jitter
    const double slop = 0.01;
    const double percent = 0.8; // correction strength

    double correction_mag = std::max(overlap - slop, 0.0) * percent;
    Vector2D correction = normal_vector * correction_mag;

    double totalMass = smallerBody.mass + largerBody.mass;

    // smallerBody.position -= correction * (largerBody.mass / totalMass);
    // largerBody.position += correction * (smallerBody.mass / totalMass);
    if (!smallerBody.isStatic && !largerBody.isStatic) {
      smallerBody.position -= correction * (largerBody.mass / totalMass);
      largerBody.position += correction * (smallerBody.mass / totalMass);
    } else if (smallerBody.isStatic && !largerBody.isStatic) {
      largerBody.position += normal_vector * correction_mag;
    } else if (!smallerBody.isStatic && largerBody.isStatic) {
      smallerBody.position -= normal_vector * correction_mag;
    }
    // both static → do nothing

    // Prevent fake velocity from Verlet-style position jumps
    smallerBody.previousPosition = smallerBody.position;
    largerBody.previousPosition = largerBody.position;
  }
}

void PhysicsSystem::handleDynamicExplosionCollision(
    GravitationalBody& macroBody1, GravitationalBody& macroBody2,
    GameState& gameState) {
  // Buffers to prevent iterator invalidation
  std::vector<GravitationalBody> pending_particles;
  std::vector<GravitationalBody> pending_macros;

  // instantiate collision proxy for both bodies. Consider a shrinking proxy?Add
  // later?

  GravitationalBody proxyForMB1;
  populateProxyFromMacroBody();

  // 2. Instantiate the "Invisible Bouncy Circle" (Proxy Body)
  // GravitationalBody proxy;
  // proxy.position = center_of_explosion;
  // proxy.previousPosition = center_of_explosion;
  // proxy.velocity = blended_vel;
  // proxy.radius = falloff_radius; // Size of the shockwave
  // proxy.mass = 0.0;              // Massless so it doesn't affect gravity
  // proxy.isBounce = true;         // Makes fragments reflect
  // proxy.isStatic = true;         // Makes the object fixed.
  // proxy.isCollidable = true;
  // proxy.isMacroGhost = true; // proxy is macro particle but does not collide
  // with macro bodies. proxy.isTransient = true; proxy.lifetime = TARGET_FPS;
  // // sloppy solution for now but probably fine. NEED TO FIX double avgSpeed =
  // (body1.velocity.magnitude() + body2.velocity.magnitude()) * 0.5; double
  // avgMass = (body1.mass + body2.mass) / 2; double sizeFactor = (body1.radius
  // + body2.radius) * 0.5; lifetime in frames proxy.lifetime =
  // static_cast<uint32_t>(TARGET_FPS * (5 * sizeFactor / (avgSpeed + 1e-5)));
  // proxy.lifetime = std::min(static_cast<int>(proxy.lifetime), TARGET_FPS *
  // 3); // max 3 seconds proxy.lifetime = TARGET_FPS * (1 + 5 * avgMass /
  // MAX_MASS); Note: You'll need to add a 'lifeSpan' or similar to your struct
  // and decrement it in your main loop to eventually delete this.
  pending_macros.push_back(proxy);

  // 3. Helper to process each body's fragments
  auto processBody = [&](GravitationalBody& source, const Vector2D& otherVel) {
    int R_int = static_cast<int>(source.radius);
    double sourceMass = source.mass;

    // Count pixels first for mass distribution
    int pixel_count = 0;
    for (int i = -R_int; i <= R_int; ++i) {
      for (int j = -R_int; j <= R_int; ++j) {
        if ((i * i) + (j * j) <= (source.radius * source.radius))
          pixel_count++;
      }
    }
    if (pixel_count == 0)
      return;
    double p_mass = sourceMass / pixel_count;

    for (int i = -R_int; i <= R_int; ++i) {
      for (int j = -R_int; j <= R_int; ++j) {
        if ((i * i) + (j * j) <= (source.radius * source.radius)) {
          // new solution
          // Vector2D pos = source.position + Vector2D{(double)i, (double)j};

          // Vector2D r = pos - center_of_explosion;
          // double dist = r.magnitude();
          // Vector2D n = (dist > 1e-8) ? (r / dist) : Vector2D(1, 0);
          // Vector2D t = Vector2D(-n.yVal, n.xVal);

          // // Calculate Falloff (1.0 at center, ~0.0 at edges)
          // double d = dist / (falloff_radius + 1e-8);
          // double falloff = 1.0 / (1.0 + d * d);

          // // Velocity Calculation
          // // Blends between the original body velocity and the collision
          // momentum Vector2D base_vel = lerp(source.velocity, blended_vel,
          // falloff);

          // double factor = source.velocity.magnitude();
          // constexpr double EXPLOSION_STRENGTH = 0.5; // 0.0 = no
          // explosion, 1.0 = current strength double blast = factor * falloff *
          // EXPLOSION_STRENGTH; double shear = factor * falloff * 1.5 *
          // EXPLOSION_STRENGTH; double jitter = randomDouble(-1.0, 1.0) *
          // falloff * EXPLOSION_STRENGTH; constexpr double RADIAL_SCALE = 0.6;
          // constexpr double TANGENTIAL_SCALE = 0.3;

          // GravitationalBody p;
          // p.mass = p_mass;
          // p.radius = 1.0;
          // p.position = pos;
          // p.isFragment = true;
          // p.velocity = base_vel + n * (blast * RADIAL_SCALE + jitter) + t *
          // (shear * TANGENTIAL_SCALE * randomDouble(-1.0, 1.0));
          // p.previousPosition = p.position - p.velocity * PHYSICS_TIME_STEP;
          // pending_particles.push_back(p);

          // old solution
          Vector2D pos = source.position + Vector2D{(double)i, (double)j};

          Vector2D r = pos - center_of_explosion;
          double dist = r.magnitude();
          Vector2D n = (dist > 1e-8) ? (r / dist) : Vector2D(1, 0);
          Vector2D t = Vector2D(-n.yVal, n.xVal);

          // Calculate Falloff (1.0 at center, ~0.0 at edges)
          double d = dist / (falloff_radius + 1e-8);
          double falloff = 1.0 / (1.0 + d * d);

          // Velocity Calculation
          // Blends between the original body velocity and the collision
          // momentum
          Vector2D base_vel = lerp(source.velocity, blended_vel, falloff);

          double factor = source.velocity.magnitude();
          double blast = factor * falloff * 0.25;
          double shear =
              factor * falloff * 1.5; // Slight reduction in shear for stability
          double jitter = randomDouble(-1.0, 1.0) * falloff;

          GravitationalBody p;
          p.mass = p_mass;
          p.radius = 1.0;
          p.position = pos;
          p.isFragment = true;
          p.velocity = base_vel + n * (blast + jitter) +
                       t * (shear * randomDouble(-1.0, 1.0));
          p.previousPosition = p.position - p.velocity * PHYSICS_TIME_STEP;

          pending_particles.push_back(p);
        }
      }
    }
  };

  processBody(body1, body2.velocity);
  processBody(body2, body1.velocity);

  // 4. Final Cleanup and Vector Updates
  body1.isMarkedForDeletion = true;
  body2.isMarkedForDeletion = true;

  // Transfer pending items to the actual GameState
  auto& particles = gameState.getParticlesMutable();
  particles.insert(particles.end(), pending_particles.begin(),
                   pending_particles.end());

  auto& macros = gameState.getMacroBodiesMutable();
  macros.insert(macros.end(), pending_macros.begin(), pending_macros.end());
}

void PhysicsSystem::handleAccretion(GravitationalBody& particle,
                                    GravitationalBody& body) {
  // Accrete mass
  body.mass += particle.mass;
  // Add to radius
  body.radius =
      body.radius * pow((body.mass + particle.mass) / body.mass, 1.0 / 3.0);
  particle.isMarkedForDeletion = true;
}

void PhysicsSystem::createMacroBody(GameState& gameState,
                                    InputState& inputState) {
  GravitationalBody body;
  body.mass = inputState.selectedMass;
  body.radius = inputState.selectedRadius;
  body.position = inputState.mouseCurrPosition;
  body.previousPosition = body.position;
  body.isPlanet = true;
  body.isStatic = inputState.isCreatingStatic;

  // --- NEW LOGIC: Nudge fragments out of the new body's radius ---
  auto& particles = gameState.getParticlesMutable();
  double nudge_factor = 1.01; // Nudge fragments out by 1% more than the radius

  for (auto& particle : particles) {
    Vector2D displacement = particle.position - body.position;
    double dist = displacement.magnitude();
    double min_dist = body.radius + particle.radius;

    if (dist < min_dist) {
      // Calculate the required outward displacement
      // We want the new distance to be min_dist * nudge_factor
      double required_separation = (min_dist * nudge_factor) - dist;

      // Normalize the displacement vector (safety check for dist=0, though
      // unlikely here)
      Vector2D direction =
          (dist == 0) ? Vector2D(1.0, 0.0) : displacement / dist;

      // Apply the positional nudge
      particle.position += direction * required_separation;

      // Crucial for stability in Verlet integration:
      // Set prevPosition to the new position to prevent the next frame's
      // velocity calculation from being massive and inaccurate.
      particle.previousPosition = particle.position;
    }
  }

  gameState.getMacroBodiesMutable().push_back(body);
}

void PhysicsSystem::createParticle(GameState& gameState,
                                   InputState& inputState) {
  GravitationalBody body;
  body.mass = inputState.selectedMass;
  body.radius = inputState.selectedRadius;
  body.position = inputState.mouseCurrPosition;
  body.previousPosition = body.position;
  body.isDust = inputState.isCreatingDust;
  body.isStatic = inputState.isCreatingStatic;
  inputState.isCreatingDust = false;
  gameState.getParticlesMutable().push_back(body);
  inputState.instantiateDirty = false;
}

void PhysicsSystem::createParticleCluster(GameState& gameState,
                                          InputState& inputState) {
  GravitationalBody body;
  body.mass = inputState.selectedMass;
  body.radius = inputState.selectedRadius;
  body.position = inputState.mouseCurrPosition;
  body.previousPosition = body.position;
  body.isPlanet = true;
  body.isStatic = inputState.isCreatingStatic;

  substituteWithParticles(body, gameState);
}
// --------- TOTAL ENERGY CALCULATION METHOD --------- //

void PhysicsSystem::calculateTotalEnergy(GameState& gameState) {
  // Calculate total energy logic here
  auto& bodies = gameState.getMacroBodies();
  int num_bodies = bodies.size();
  double totalEnergy = 0.0;
  for (int i = 0; i < num_bodies; ++i) {
    totalEnergy += bodies[i].mass * bodies[i].velocity.square_magnitude() / 2.0;
  }

  for (size_t i = 0; i < num_bodies; i++) {
    for (size_t j = i + 1; j < num_bodies; j++) {
      double epsilon =
          (bodies[i].radius + bodies[j].radius) / 2.0; // Simple average radius
      double epsilonSq = epsilon * epsilon;

      // 1. Calculate Distance Vector
      Vector2D distance = bodies[i].position - bodies[j].position;

      // 2. Calculate Distance Squared (r^2)
      double rSq = distance.square_magnitude();

      // 3. Calculate the Denominator Term (r^2 + epsilon^2)^(3/2)
      // The term inside the parenthesis: rSq + epsilonSq
      // The final term in the denominator: pow(rSq + epsilonSq, 1.5)
      double denominator = sqrt(rSq + epsilonSq);

      totalEnergy -= GRAVITATIONAL_CONSTANT * bodies[i].mass * bodies[j].mass /
                     denominator;
    }
  }
  // std::cout<<"Total E: "<< totalEnergy << std::endl;
}

// --------- PARTICLE SUBSTITUTION METHOD --------- //

void PhysicsSystem::substituteWithParticles(GravitationalBody& originalBody,
                                            GameState& gameState) {
  // A simplified value for the collision radius (R)
  const double R = originalBody.radius;
  const Vector2D center = originalBody.position;
  const double originalMass = originalBody.mass;
  const Vector2D originalVelocity = originalBody.velocity;

  std::vector<Vector2D> pixelPositions;

  // 1. RASTERIZATION (Find all pixel positions)
  // Assuming a 1-to-1 mapping where 1 unit = 1 pixel for simplicity.
  int R_int = static_cast<int>(R);
  for (int i = -R_int; i <= R_int; ++i) {
    for (int j = -R_int; j <= R_int; ++j) {
      // Check if (i, j) is inside the circle
      if ((i * i) + (j * j) <= (R * R)) {
        // Store the particle's center position
        pixelPositions.push_back(center + Vector2D{(double)i, (double)j});
      }
    }
  }

  // Check if we found any pixels (safety)
  if (pixelPositions.empty()) {
    // std::cout<<"no pixels"<<std::endl;
    return;
  }

  // 2. MASS CALCULATION
  const double particleMass = originalMass / pixelPositions.size();

  // 3. PARTICLE INSTANTIATION
  for (const auto& pos : pixelPositions) {
    GravitationalBody p;
    p.mass = particleMass;
    p.radius = 1.0;
    p.position = pos;
    p.previousPosition = pos; // Initialize for Verlet integration
    p.isFragment = true;

    // Start with the original body's momentum
    p.velocity = originalVelocity;

    // Add a random outward "explosion" vector
    // This is necessary to make the cloud expand.
    Vector2D explosionVector = (pos - center); // Vector from center to particle
    // Normalize and scale by a collision factor (e.g., 0.1)
    explosionVector =
        explosionVector.normalize() * (0.0 * originalVelocity.magnitude());

    p.velocity += explosionVector;
    p.previousPosition = p.position - p.velocity * PHYSICS_TIME_STEP;

    gameState.getParticlesMutable().push_back(p);
  }

  // 4. CLEANUP
  originalBody.isMarkedForDeletion =
      true; // Destroy the original macroscopic body
}

// --------- CLEANUP GRAVITATIONAL BODIES METHODS --------- //

void PhysicsSystem::cleanupParticles(GameState& gameState) {
  auto& particles = gameState.getParticlesMutable();

  // 1. Use std::remove_if to move all elements marked for deletion
  //    to the end of the vector. It returns an iterator to the new end.
  auto new_end = std::remove_if(particles.begin(), particles.end(),
                                [](const GravitationalBody& p) {
                                  // The predicate returns true for elements to
                                  // be 'removed' (moved to end)
                                  return p.isMarkedForDeletion;
                                });

  // 2. Use vector::erase to destroy the elements in the range [new_end,
  // particles.end())
  //    This efficiently shrinks the vector to the correct size.
  particles.erase(new_end, particles.end());
}

void PhysicsSystem::cleanupMacroBodies(GameState& gameState) {
  auto& particles = gameState.getMacroBodiesMutable();

  // 1. Use std::remove_if to move all elements marked for deletion
  //    to the end of the vector. It returns an iterator to the new end.
  auto new_end = std::remove_if(particles.begin(), particles.end(),
                                [](const GravitationalBody& b) {
                                  // The predicate returns true for elements to
                                  // be 'removed' (moved to end)
                                  return b.isMarkedForDeletion;
                                });

  // 2. Use vector::erase to destroy the elements in the range [new_end,
  // particles.end())
  //    This efficiently shrinks the vector to the correct size.
  particles.erase(new_end, particles.end());
}
