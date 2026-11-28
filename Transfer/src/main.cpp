// File: Transfer/src/main.cpp

// Custom Imports
#include "Core/DeterminismCheck.hpp"
#include "Core/Game.hpp"

// Standard Library Imports
#include <cstring>

int main(int argc, char* argv[])
{
    // TransferGame --determinism-check: play the fixed check scenes and print their checksums instead of starting the
    // game (see Core/DeterminismCheck.hpp)
    for (int i = 1; i < argc; ++i)
    {
        if (std::strcmp(argv[i], "--determinism-check") == 0)
        {
            return runDeterminismCheck();
        }
    }

    Game gameInstance;
    gameInstance.StartGame();
    return 0;
}
