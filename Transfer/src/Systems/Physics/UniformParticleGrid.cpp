// File: Transfer/src/Systems/Physics/UniformParticleGrid.cpp

#include "Systems/Physics/UniformParticleGrid.hpp"

#include <algorithm>
#include <cmath>

namespace
{
// How many slots per particle. More slots = fewer cells sharing a slot, but a bigger table to clear every tick.
constexpr size_t SLOTS_PER_PARTICLE = 2;
constexpr size_t MIN_SLOT_COUNT = 64;

// The smallest power of two that is >= value (the slot count is a power of two so a slot is picked with "& mask")
size_t roundUpToPowerOfTwo(size_t value)
{
    size_t power = 1;
    while (power < value)
    {
        power *= 2;
    }
    return power;
}
} // namespace

int64_t UniformParticleGrid::cellColumn(double x) const
{
    // floor, not a plain cast: a cast rounds toward zero, which would put x = -0.5 and x = +0.5 in the same cell
    return static_cast<int64_t>(std::floor(x / m_cell_size));
}

int64_t UniformParticleGrid::cellRow(double y) const { return static_cast<int64_t>(std::floor(y / m_cell_size)); }

uint32_t UniformParticleGrid::slotOfCell(int64_t column, int64_t row) const
{
    // Multiply the column and row by two big odd constants (they scatter the bits upward) and mix them with XOR, then
    // fold the high half down (the multiplications leave the best-mixed bits up there). Neighboring cells land in
    // unrelated slots, and two occupied cells rarely share one.
    uint64_t hash = static_cast<uint64_t>(column) * 0x9E3779B97F4A7C15ULL;
    hash ^= static_cast<uint64_t>(row) * 0xC2B2AE3D27D4EB4FULL;
    hash ^= hash >> 32;
    return static_cast<uint32_t>(hash & m_slot_mask);
}

void UniformParticleGrid::chooseCellSize(const std::vector<GravitationalBody>& particles)
{
    // As wide as the biggest particle, so touching particles are always in the same or neighboring cells
    double largest_radius = 0.5; // floor, so the cell size never collapses to ~0 with no or tiny particles
    for (const GravitationalBody& particle : particles)
    {
        largest_radius = std::max(largest_radius, particle.radius);
    }
    m_cell_size = 2.0 * largest_radius;
}

void UniformParticleGrid::build(const std::vector<GravitationalBody>& particles)
{
    chooseCellSize(particles);

    const size_t particle_count = particles.size();
    const size_t slot_count = roundUpToPowerOfTwo(std::max(MIN_SLOT_COUNT, SLOTS_PER_PARTICLE * particle_count));
    m_slot_mask = static_cast<uint32_t>(slot_count - 1);

    // Each particle's cell and slot. m_slot_start counts the particles per slot for now (pass 1 of the counting sort):
    // slot s is counted in m_slot_start[s + 1], so that pass 2 turns the counts into starts in place.
    m_particle_column.resize(particle_count);
    m_particle_row.resize(particle_count);
    m_particle_slot.resize(particle_count);
    m_slot_start.assign(slot_count + 1, 0);
    for (size_t i = 0; i < particle_count; ++i)
    {
        int64_t column = cellColumn(particles[i].position.x_val);
        int64_t row = cellRow(particles[i].position.y_val);
        uint32_t slot = slotOfCell(column, row);
        m_particle_column[i] = column;
        m_particle_row[i] = row;
        m_particle_slot[i] = slot;
        m_slot_start[slot + 1] += 1;
    }

    // Pass 2: each slot's section starts where the previous one ended (a running total of the counts)
    for (size_t slot = 0; slot < slot_count; ++slot)
    {
        m_slot_start[slot + 1] += m_slot_start[slot];
    }

    // Pass 3: put each particle into the next free spot of its slot's section. Going through the particles in index
    // order keeps each section in index order.
    m_next_free.assign(m_slot_start.begin(), m_slot_start.end() - 1);
    m_sorted_entries.resize(particle_count);
    for (size_t i = 0; i < particle_count; ++i)
    {
        uint32_t slot = m_particle_slot[i];
        m_sorted_entries[m_next_free[slot]] = {m_particle_column[i], m_particle_row[i], static_cast<uint32_t>(i)};
        m_next_free[slot] += 1;
    }
}

void UniformParticleGrid::sortByCell(std::vector<GravitationalBody>& particles)
{
    const size_t particle_count = particles.size();
    if (particle_count < 2)
    {
        return;
    }
    chooseCellSize(particles);

    // 1. Each particle's cell, and the smallest box of cells that holds them all
    m_particle_column.resize(particle_count);
    m_particle_row.resize(particle_count);
    int64_t min_column = cellColumn(particles[0].position.x_val);
    int64_t max_column = min_column;
    int64_t min_row = cellRow(particles[0].position.y_val);
    int64_t max_row = min_row;
    for (size_t i = 0; i < particle_count; ++i)
    {
        m_particle_column[i] = cellColumn(particles[i].position.x_val);
        m_particle_row[i] = cellRow(particles[i].position.y_val);
        min_column = std::min(min_column, m_particle_column[i]);
        max_column = std::max(max_column, m_particle_column[i]);
        min_row = std::min(min_row, m_particle_row[i]);
        max_row = std::max(max_row, m_particle_row[i]);
    }

    // 2. The sort key: the cell's number when the box's cells are counted column by column, top to bottom.
    //    Sorting by it puts columns left to right, and each column top to bottom.
    const uint64_t box_height = static_cast<uint64_t>(max_row - min_row) + 1;
    m_sort_keys.resize(particle_count);
    for (size_t i = 0; i < particle_count; ++i)
    {
        uint64_t column_in_box = static_cast<uint64_t>(m_particle_column[i] - min_column);
        uint64_t row_in_box = static_cast<uint64_t>(m_particle_row[i] - min_row);
        m_sort_keys[i] = column_in_box * box_height + row_in_box;
    }
    const uint64_t largest_key = static_cast<uint64_t>(max_column - min_column) * box_height + (box_height - 1);

    // 3. Radix sort: the counting sort, run once per byte of the key, lowest byte first. Each pass is STABLE (it
    //    keeps the order the previous passes made among particles whose byte is equal), so after the pass for the
    //    highest byte the whole key is in order. Only as many passes as the largest key has bytes (a 2000-wide box of
    //    6-wide cells needs about 2), so it costs a few straight sweeps over the particles.
    m_sort_order.resize(particle_count);
    m_sort_order_next.resize(particle_count);
    for (size_t i = 0; i < particle_count; ++i)
    {
        m_sort_order[i] = static_cast<uint32_t>(i);
    }
    for (int shift = 0; shift < 64 && (largest_key >> shift) != 0; shift += 8)
    {
        // Count each byte value (byte b is counted in starts[b + 1], so the running total turns counts into starts)
        uint32_t starts[257] = {};
        for (size_t k = 0; k < particle_count; ++k)
        {
            uint32_t byte = static_cast<uint32_t>((m_sort_keys[m_sort_order[k]] >> shift) & 0xFF);
            starts[byte + 1] += 1;
        }
        for (int byte = 0; byte < 256; ++byte)
        {
            starts[byte + 1] += starts[byte];
        }
        // Place each particle (in the current order) at the next free spot of its byte's section
        for (size_t k = 0; k < particle_count; ++k)
        {
            uint32_t particle_index = m_sort_order[k];
            uint32_t byte = static_cast<uint32_t>((m_sort_keys[particle_index] >> shift) & 0xFF);
            m_sort_order_next[starts[byte]] = particle_index;
            starts[byte] += 1;
        }
        m_sort_order.swap(m_sort_order_next);
    }

    // 4. Copy the particles in the sorted order, then swap the lists: `particles` takes the reordered copy, and its
    //    old storage becomes the destination next time (swap only exchanges the vectors' insides: nothing is copied)
    m_reordered_particles.resize(particle_count);
    for (size_t k = 0; k < particle_count; ++k)
    {
        m_reordered_particles[k] = particles[m_sort_order[k]];
    }
    particles.swap(m_reordered_particles);
}

void UniformParticleGrid::addCellCandidates(int64_t column, int64_t row, size_t index, bool only_after_index,
                                            std::vector<size_t>& candidates_out) const
{
    uint32_t slot = slotOfCell(column, row);
    uint32_t section_end = m_slot_start[slot + 1];
    for (uint32_t k = m_slot_start[slot]; k < section_end; ++k)
    {
        const Entry& entry = m_sorted_entries[k];
        if (entry.column != column || entry.row != row)
        {
            continue; // another cell that shares this slot
        }
        if (only_after_index && entry.particle_index <= index)
        {
            continue; // this same-cell pair is reported by the other particle's query (or it's the particle itself)
        }
        candidates_out.push_back(entry.particle_index);
    }
}

void UniformParticleGrid::queryCandidates(size_t index, std::vector<size_t>& candidates_out) const
{
    candidates_out.clear();
    const int64_t column = m_particle_column[index];
    const int64_t row = m_particle_row[index];

    // Forward half-stencil: the own cell, the next row's cell (0, +1), and the three cells of the next column
    // (+1, -1), (+1, 0), (+1, +1). For any two neighboring cells, exactly one sees the other in this shape, so each
    // pair is reported once:
    //        column  column+1
    //   row-1   .       X
    //   row    own      X
    //   row+1   X       X
    addCellCandidates(column, row, index, true, candidates_out);
    addCellCandidates(column, row + 1, index, false, candidates_out);
    addCellCandidates(column + 1, row - 1, index, false, candidates_out);
    addCellCandidates(column + 1, row, index, false, candidates_out);
    addCellCandidates(column + 1, row + 1, index, false, candidates_out);
}
