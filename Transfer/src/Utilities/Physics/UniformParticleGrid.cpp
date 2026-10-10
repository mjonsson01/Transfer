// File: Transfer/src/Utilities/Physics/UniformParticleGrid.cpp

#include "Utilities/Physics/UniformParticleGrid.hpp"

#include <algorithm>
#include <cmath>

namespace
{
// Packs two cell coordinates into one sortable key. Each coordinate is offset by a large bias
// before packing so negative world positions produce non-negative packed components; floor()
// (not truncation) is used when deriving the coordinate itself so negative positions bucket
// the same consistent way positive ones do.
// queryCandidates RELIES on this layout: sorting by key sorts by column (cx), then by row (cy), so the
// cells (cx, cy) and (cx, cy + 1) have consecutive keys and a column can be walked with a forward scan.
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
        int64_t cx = static_cast<int64_t>(std::floor(particles[i].position.x_val / cellSize));
        int64_t cy = static_cast<int64_t>(std::floor(particles[i].position.y_val / cellSize));
        sortedEntries.push_back({packCell(cx, cy), i});
    }

    std::sort(sortedEntries.begin(), sortedEntries.end(),
              [](const Entry& a, const Entry& b) { return a.cellKey < b.cellKey; });
}

void UniformParticleGrid::queryCandidates(size_t index, const std::vector<GravitationalBody>& particles,
                                          std::vector<size_t>& outCandidates) const
{
    outCandidates.clear();

    const DynamoEngine::Vector2D& pos = particles[index].position;
    int64_t cx = static_cast<int64_t>(std::floor(pos.x_val / cellSize));
    int64_t cy = static_cast<int64_t>(std::floor(pos.y_val / cellSize));

    // Forward half-stencil: the own cell, the next row's cell (0,1), and the three cells of the next column
    // (1,-1), (1,0), (1,1). For every pair of neighbouring cells exactly one sees the other in this shape, so
    // each pair is reported once. Each column's cells have consecutive keys (see packCell), so a column costs
    // ONE binary search plus a forward scan, instead of one search per cell.

    // Column cx: the own cell, then the next row's cell (cy + 1)
    int64_t own_key = packCell(cx, cy);
    size_t i = firstEntryAtOrAfter(own_key);
    while (i < sortedEntries.size() && sortedEntries[i].cellKey == own_key)
    {
        if (sortedEntries[i].particleIndex > index) // each same-cell pair reported once, from the lower index's query
        {
            outCandidates.push_back(sortedEntries[i].particleIndex);
        }
        ++i;
    }
    while (i < sortedEntries.size() && sortedEntries[i].cellKey == own_key + 1)
    {
        outCandidates.push_back(sortedEntries[i].particleIndex);
        ++i;
    }

    // Column cx + 1: the cells from row cy - 1 to cy + 1
    int64_t first_key = packCell(cx + 1, cy - 1);
    int64_t last_key = packCell(cx + 1, cy + 1);
    i = firstEntryAtOrAfter(first_key);
    while (i < sortedEntries.size() && sortedEntries[i].cellKey <= last_key)
    {
        outCandidates.push_back(sortedEntries[i].particleIndex);
        ++i;
    }
}

size_t UniformParticleGrid::firstEntryAtOrAfter(int64_t cellKey) const
{
    // Binary search: the first entry whose key is >= cellKey (sortedEntries.size() if there is none)
    auto it = std::lower_bound(sortedEntries.begin(), sortedEntries.end(), cellKey,
                               [](const Entry& e, int64_t key) { return e.cellKey < key; });
    return static_cast<size_t>(it - sortedEntries.begin());
}
