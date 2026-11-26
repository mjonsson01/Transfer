// File: Transfer/src/Systems/Physics/UniformParticleGrid.hpp

#pragma once

#include "Systems/Physics/GravitationalBody.hpp"

#include <cstdint>
#include <vector>

// Broad phase for particle-particle collisions: which particles are close enough that they MIGHT touch.
//
// The world is cut into square cells, each as wide as the biggest particle, so two particles can only touch if they
// are in the same cell or in neighboring cells. The cells form an endless 2D grid (the world has no fixed size), too
// big to store, so each cell is HASHED (folded) into one of a fixed number of SLOTS, about two per particle. Usually a
// slot holds one occupied cell; now and then two cells share one, and queryCandidates skips the particles of the
// cell it didn't ask for.
//
// build() fills two flat arrays with a COUNTING SORT (three straight passes, no comparisons):
//     m_slot_start: one entry per slot, where that slot's section of m_sorted_entries begins
//     m_sorted_entries: one entry per particle, grouped by slot, so each slot's particles sit side by side
// Finding a cell's particles is then: hash the cell, read m_slot_start[slot] and m_slot_start[slot + 1].
//
// Rebuilt from scratch every physics tick (particles move every tick). The cell size follows the biggest particle
// present THIS tick, so it is right whatever made the particles.
class UniformParticleGrid
{
  public:
    // Builds the grid from the current particle list. Must be called before queryCandidates().
    void build(const std::vector<GravitationalBody>& particles);

    // Fills candidates_out with indices (into the particle vector passed to build()) of particles that may be
    // touching particle `index`, as it was at build(). Each pair of particles is returned at most once over a full
    // pass over every index for one build(), so callers don't need to remove duplicates. Candidates come out in a
    // fixed order (by cell, then by particle index), so the same particles always give the same candidates.
    void queryCandidates(size_t index, std::vector<size_t>& candidates_out) const;

    // Reorders `particles` by position: column by column (left to right), and top to bottom within a column. A query
    // reads its own column and the next one, so in this order the pair pass reads memory in two smooth streams
    // instead of jumping all over the list. Sorted with a RADIX SORT (the counting sort, one byte of the key at a
    // time), the way a GPU sorts. Call build() afterwards: the grid isn't left describing the new order.
    void sortByCell(std::vector<GravitationalBody>& particles);

  private:
    // One particle in the sorted list, with its cell, so a query can skip particles from another cell that happens
    // to share the slot without looking the particle itself up
    struct Entry
    {
        int64_t column;
        int64_t row;
        uint32_t particle_index;
    };

    // Sets m_cell_size from the biggest particle (build() and sortByCell() both start with it)
    void chooseCellSize(const std::vector<GravitationalBody>& particles);
    int64_t cellColumn(double x) const;
    int64_t cellRow(double y) const;
    // The slot a cell is folded into (0 .. slot count - 1)
    uint32_t slotOfCell(int64_t column, int64_t row) const;
    // Adds the particles in cell (column, row) to candidates_out. only_after_index: same-cell pairs are reported
    // only from the query of the particle with the lower index, so each pair comes out once.
    void addCellCandidates(int64_t column, int64_t row, size_t index, bool only_after_index,
                           std::vector<size_t>& candidates_out) const;

    double m_cell_size = 1.0;
    uint32_t m_slot_mask = 0; // slot count - 1 (the slot count is a power of two, so "& mask" picks the slot)
    std::vector<uint32_t> m_slot_start;  // slot count + 1 entries: the last one is the end of the list
    std::vector<Entry> m_sorted_entries; // one per particle, grouped by slot
    // Each particle's cell, by particle index, so a query knows where its own particle is without recomputing it
    std::vector<int64_t> m_particle_column;
    std::vector<int64_t> m_particle_row;
    // Work space for build(), kept between ticks so building doesn't allocate
    std::vector<uint32_t> m_particle_slot;
    std::vector<uint32_t> m_next_free;
    // sortByCell's work space: each particle's sort key, the order being sorted (and a second buffer for each
    // radix pass), and the destination of the reordered particles (swapped with the particle list at the end)
    std::vector<uint64_t> m_sort_keys;
    std::vector<uint32_t> m_sort_order;
    std::vector<uint32_t> m_sort_order_next;
    std::vector<GravitationalBody> m_reordered_particles;
};
