// File: Transfer/src/Utilities/Physics/UniformParticleGrid.hpp

#pragma once

#include "Entities/Physics/GravitationalBody.hpp"

#include <cstdint>
#include <vector>

// Broad-phase spatial index for particle-particle collision candidates. Rebuilt from scratch
// every physics tick (particles move every tick, so nothing would be gained by persisting
// bucket contents across ticks). Cell size is derived from the largest particle radius present
// *this tick*, so it stays correct regardless of what generated the particles or how much their
// sizes vary -- it never needs to know about density factors, impact skew, or anything else
// upstream.
class UniformParticleGrid
{
  public:
    // Builds the grid from the current particle list. Must be called before queryCandidates().
    void build(const std::vector<GravitationalBody>& particles);

    // Fills outCandidates with indices (into the same particle vector passed to build()) of
    // particles that may be colliding with particles[index]. Each unordered pair is returned
    // at most once across a full pass over all indices for one build() -- callers don't need
    // their own de-duplication.
    void queryCandidates(size_t index, const std::vector<GravitationalBody>& particles,
                         std::vector<size_t>& outCandidates) const;

  private:
    struct Entry
    {
        int64_t cellKey;
        size_t particleIndex;
    };

    double cellSize = 1.0;
    std::vector<Entry> sortedEntries;
};