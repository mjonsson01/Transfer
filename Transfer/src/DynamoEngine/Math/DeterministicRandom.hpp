// File: Transfer/src/DynamoEngine/Math/DeterministicRandom.hpp

#pragma once

// Standard Library Imports
#include <cstdint>
#include <cstring>

namespace DynamoEngine
{
// Repeatable random numbers, the same on every machine and every run.
//
// Instead of a generator that hands out the "next" number (whose answers depend on how many numbers were drawn before,
// and in what order), a random number is HASHED from what it is for: a seed, the tick, which body, which fragment, and
// a "purpose" so two different decisions about the same fragment don't get the same number. Same inputs -> same
// number, whatever order things happen in, so it also stays the same when the work is later split across threads.
//     uint64_t key = combineHash(seed, tick);
//     key = combineHash(key, body_id);
//     key = combineHash(key, fragment_index);
//     double speed_scale = hashToRange(combineHash(key, PURPOSE_SPEED), 0.8, 1.1);
// (std::mt19937 and std::uniform_real_distribution aren't used: the distributions are allowed to differ between the
// standard libraries of Apple and Microsoft.)

// Scrambles 64 bits so that inputs differing in a single bit give unrelated outputs (SplitMix64's finalizer, a
// well-tested mixing function). Only integer multiply, xor and shift: identical everywhere.
inline uint64_t mixBits(uint64_t value)
{
    value += 0x9E3779B97F4A7C15ULL;
    value = (value ^ (value >> 30)) * 0xBF58476D1CE4E5B9ULL;
    value = (value ^ (value >> 27)) * 0x94D049BB133111EBULL;
    return value ^ (value >> 31);
}

// Adds one more input to a hash (order matters: combineHash(combineHash(h, a), b) != combineHash(combineHash(h, b), a))
inline uint64_t combineHash(uint64_t hash, uint64_t value) { return mixBits(hash ^ mixBits(value)); }

// A double's exact bits, so positions and masses can be hashed
inline uint64_t doubleBits(double value)
{
    uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

// A hash as a number in [0, 1): its top 53 bits (a double holds 53 bits exactly) divided by 2^53
inline double hashToUnitInterval(uint64_t hash)
{
    return static_cast<double>(hash >> 11) * (1.0 / 9007199254740992.0); // 2^53 = 9007199254740992
}

// A hash as a number in [min, max)
inline double hashToRange(uint64_t hash, double min, double max) { return min + (max - min) * hashToUnitInterval(hash); }
} // namespace DynamoEngine
