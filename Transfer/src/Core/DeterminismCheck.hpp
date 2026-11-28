// File: Transfer/src/Core/DeterminismCheck.hpp

#pragma once

// The cross-machine repeatability check (Notes/PHYSICS_REQUIREMENTS.md, requirement A1).
//
// Run the game with --determinism-check: instead of opening a window it plays a fixed set of scenes (debris piles,
// planets merging and shattering, negative mass, the ship flying through debris) with no player input, and prints a
// checksum of every body's exact bits along the way, ending with one OVERALL line. The same lines are written to
// determinism_check.txt in the current folder.
//
// Every build on every machine must print exactly the same lines: compare the file from the Mac with the one from
// Windows. Any difference means some physics code gives different floating-point results on the two machines.
//
// Returns the process exit code (0).
int runDeterminismCheck();
