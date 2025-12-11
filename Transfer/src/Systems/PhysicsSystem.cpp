// File: Transfer/src/Systems/PhysicsSystem.cpp

#include "PhysicsSystem.h"

PhysicsSystem::PhysicsSystem()
{
    // Initialize physics system variables if needed
    // physicsData.grid.resize(physicsData.gridWidth * physicsData.gridHeight);
}

PhysicsSystem::~PhysicsSystem()
{
    // Cleanup if necessary
}
void PhysicsSystem::CleanUp()
{
    // Any necessary cleanup code for the physics system
}



void PhysicsSystem::UpdateSystemFrame(GameState& state, UIState& UIState)
{
    // Apply to particles
    auto& inputState = UIState.getMutableInputState();
    if (inputState.dirty)
    {
        if (inputState.isCreatingCluster)
        {   
            createParticleCluster(state.physicsData, inputState);
            // state.physicsData.particles.push_back(inputState.newCluster);
            inputState.dirty = false;
        }
    }

    applyGravity(state.physicsData);
    predictPositions(state.physicsData);
    buildSpatialGrid(state.physicsData);
    solveXPBDConstraints(state.physicsData);
    for (auto& p : state.physicsData.particles)
    {
        if (p.inv_mass == 0.0) continue;
        p.velocity = (p.position - p.prevPosition) / PHYSICS_TIME_STEP;
    }
    // if (state.physicsData.particles.size() > 0)
    // {
    //     std::cout<<"particle 0 position: " << state.physicsData.particles[0].position<< std::endl;
    // }
    // std::cout<<"PhysicsSystem contains "<< state.physicsData.particles.size() << " particles." << std::endl;

}

void PhysicsSystem::applyGravity(PhysicsData& physicsData)
{
    const double G = 1e-4; // tune for your simulation scale

    int N = physicsData.particles.size();
    for (int i = 0; i < N; ++i)
    {
        auto& pi = physicsData.particles[i];
        if (pi.inv_mass == 0.0) continue; // skip static particles

        Vector2D accel = {0.0, 0.0};

        for (int j = 0; j < N; ++j)
        {
            if (i == j) continue;

            auto& pj = physicsData.particles[j];

            Vector2D dir = pj.position - pi.position;
            double r2 = dir.square_magnitude() + 1.0; // avoid division by zero
            Vector2D dirNorm = dir / sqrt(r2);

            double mass_j = pj.inv_mass > 0 ? 1.0 / pj.inv_mass : 0.0;

            accel += dirNorm * (G * mass_j / r2);
        }
        // Verlet displacement due to acceleration
        pi.position += accel * PHYSICS_TIME_STEP * PHYSICS_TIME_STEP;
    }

}

// Verlet Integration to predict the positions of the particle without constraints. These are later corrected
void PhysicsSystem::predictPositions(PhysicsData& physicsData)
{
    for (auto& p : physicsData.particles)
    {
        if (p.inv_mass == 0.0) continue; // skip inifinite mass (static bodies)

        Vector2D temp = p.position;                   // store current position
        p.position += (p.position - p.prevPosition);  // displacement from last frame -> used to approximate the velocity 
        p.prevPosition = temp;                        // save previous for next step
    }
}

void PhysicsSystem::buildSpatialGrid(PhysicsData& physicsData)
{
    for (auto& cell : physicsData.grid) cell.particleIndices.clear(); // wipe previous fram particle indices

    for (int i = 0; i < physicsData.particles.size(); ++i)
    {
        const auto& p = physicsData.particles[i];
        int cellX = int(p.position.x_val / physicsData.cellSize);
        int cellY = int(p.position.y_val / physicsData.cellSize);
        int idx = cellY * physicsData.gridWidth + cellX;
        if (idx >= 0 && idx < physicsData.grid.size())
            physicsData.grid[idx].particleIndices.push_back(i);
    }
}

void PhysicsSystem::solveXPBDConstraints(PhysicsData& physicsData)
{
    for (int iter = 0; iter < SOLVER_ITERATIONS; ++iter)
    {
        for (auto& c : physicsData.constraints)
        {
            if (c.isBroken) continue;  // skip broken constraints

            XPBDSoftParticle& p1 = physicsData.particles[c.idx1];
            XPBDSoftParticle& p2 = physicsData.particles[c.idx2];

            Vector2D delta = p2.position - p1.position;
            double dist = delta.magnitude();

            if (dist == 0.0) continue; // avoid divide by zero

            // compliance term from XPBD
            double w1 = p1.inv_mass;
            double w2 = p2.inv_mass;
            double wSum = w1 + w2;
            if (wSum == 0.0) continue;

            double alpha = c.compliance / (PHYSICS_TIME_STEP * PHYSICS_TIME_STEP);  // XPBD softening

            // compute constraint correction
            double C = dist - c.restLength;
            double lambda = -C / (wSum + alpha);

            Vector2D correction = delta / dist * lambda;

            // apply corrections
            if (w1 > 0.0) p1.position += correction * w1;
            if (w2 > 0.0) p2.position -= correction * w2;

            // store lambda for next frame (optional)
            c.lambda = lambda;

            // mark as broken if stretched too far
            if (dist > c.restLength * 1.5)  // tweak threshold
                c.isBroken = true;
        }
    }
}
void PhysicsSystem::createParticleCluster(PhysicsData& physicsData, InputState& inputState)
{
    Vector2D center = inputState.mouseCurrPosition;
    double clusterRadius = inputState.selectedRadius;
    double totalMass = inputState.selectedMass;
    double compliance = 0.01;
    int numParticles = 300;

    double particleMass = totalMass / numParticles;
    double inv_mass = 1.0 / particleMass;

    std::vector<int> newParticleIndices;
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> jitter(-PARTICLE_RADIUS * 0.25, PARTICLE_RADIUS * 0.25);

    int gridSize = static_cast<int>(std::sqrt(numParticles));
    double spacing = 2.1 * PARTICLE_RADIUS;

    for (int i = 0; i < gridSize; ++i)
    {
        for (int j = 0; j < gridSize; ++j)
        {
            Vector2D offset = {
                (i - gridSize / 2) * spacing + jitter(rng),
                (j - gridSize / 2) * spacing + jitter(rng)
            };

            if (offset.magnitude() > clusterRadius) continue; // stay inside circle

            XPBDSoftParticle p;
            p.position = center + offset;
            p.prevPosition = p.position;
            p.velocity = {0.0, 0.0};
            p.inv_mass = inv_mass;
            p.mass = particleMass;
            p.radius = PARTICLE_RADIUS;

            physicsData.particles.push_back(p);
            newParticleIndices.push_back(int(physicsData.particles.size() - 1));
        }
    }

    // Create constraints between nearby particles
    for (size_t i = 0; i < newParticleIndices.size(); ++i)
    {
        int idx1 = newParticleIndices[i];
        auto& p1 = physicsData.particles[idx1];

        for (size_t j = i + 1; j < newParticleIndices.size(); ++j)
        {
            int idx2 = newParticleIndices[j];
            auto& p2 = physicsData.particles[idx2];

            Vector2D delta = p2.position - p1.position;
            double dist = delta.magnitude();

            if (dist <= PARTICLE_RADIUS * 2.5)
            {
                XPBDConstraints c;
                c.idx1 = idx1;
                c.idx2 = idx2;
                c.restLength = dist;
                c.compliance = compliance;
                c.lambda = 0.0;
                c.isBroken = false;

                physicsData.constraints.push_back(c);
            }
        }
    }
}
