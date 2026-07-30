#include "Utilities/Physics/UniformParticleGrid.hpp"

#include <algorithm>
#include <cmath>

namespace
{
// Packs two cell coordinates into one sortable key. Each coordinate is offset by a large bias
// before packing so negative world positions produce non-negative packed components; floor()
// (not truncation) is used when deriving the coordinate itself so negative positions bucket
// the same consistent way positive ones do.
constexpr int64_t CELL_COORD_BIAS = 1'000'000;

int64_t packCell(int64_t cx, int64_t cy)
{
    return (cx + CELL_COORD_BIAS) * (2 * CELL_COORD_BIAS) + (cy + CELL_COORD_BIAS);
}
} // namespace

void UniformParticleGrid::build(const std::vector<GravitationalBody>& particles)
{
    double maxRadius = 0.5; // floor, so cellSize never collapses near-0 with an empty/tiny particle set
    for (const auto& p : particles)
    {
        maxRadius = std::max(maxRadius, p.radius);
    }
    cellSize = 2.0 * maxRadius;

    sortedEntries.clear();
    sortedEntries.reserve(particles.size());
    for (size_t i = 0; i < particles.size(); ++i)
    {
        int64_t cx = static_cast<int64_t>(std::floor(particles[i].position.xVal / cellSize));
        int64_t cy = static_cast<int64_t>(std::floor(particles[i].position.yVal / cellSize));
        sortedEntries.push_back({packCell(cx, cy), i});
    }

    std::sort(sortedEntries.begin(), sortedEntries.end(),
              [](const Entry& a, const Entry& b) { return a.cellKey < b.cellKey; });
}

void UniformParticleGrid::queryCandidates(size_t index, const std::vector<GravitationalBody>& particles,
                                          std::vector<size_t>& outCandidates) const
{
    outCandidates.clear();

    const Vector2D& pos = particles[index].position;
    int64_t cx = static_cast<int64_t>(std::floor(pos.xVal / cellSize));
    int64_t cy = static_cast<int64_t>(std::floor(pos.yVal / cellSize));

    // Forward half-stencil: self + 4 directional neighbors. This specific 5-cell shape,
    // combined with the self-cell index check below, guarantees every adjacent cell pair
    // (including diagonals) is visited from exactly one side -- no separate de-dup pass needed.
    static constexpr int64_t offsets[5][2] = {{0, 0}, {1, 0}, {0, 1}, {1, 1}, {-1, 1}};

    for (const auto& offset : offsets)
    {
        int64_t neighborKey = packCell(cx + offset[0], cy + offset[1]);

        auto rangeBegin = std::lower_bound(sortedEntries.begin(), sortedEntries.end(), neighborKey,
                                           [](const Entry& e, int64_t key) { return e.cellKey < key; });
        auto rangeEnd = std::upper_bound(sortedEntries.begin(), sortedEntries.end(), neighborKey,
                                         [](int64_t key, const Entry& e) { return key < e.cellKey; });

        bool isSelfCell = (offset[0] == 0 && offset[1] == 0);
        for (auto it = rangeBegin; it != rangeEnd; ++it)
        {
            if (isSelfCell && it->particleIndex <= index)
            {
                continue; // each same-cell pair reported once, from the lower index's query
            }
            outCandidates.push_back(it->particleIndex);
        }
    }
}