// File: Transfer/src/Systems/AudioSystem.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>

// Custom Imports
#include "Core/GameState.h"
#include "Core/UIState.h"
#include "Utilities/SystemPathUtility.h"

//Standard Library Imports
#include <queue>
#include <string>
#include <filesystem>
#include <iostream>
class AudioSystem
{
    public:
        // Constructor and Destructor
        AudioSystem();
        ~AudioSystem();

    public:
        // Main method to process audio each frame
        void ProcessSystemAudioFrame(GameState& gameState, UIState& UIState);

        
        bool LoadTrack(const std::string& musicFilePath); // Loads a new track into the audio system, replacing any existing track
        void AddAllMusicToPlaylist(const std::string& folderPath);

        // Clean up helper
        void CleanUp();

    private: 
        // Add any private helper methods or member variables for audio management here
        SDL_AudioSpec wavSpec = {};
        Uint8* wavBuffer = nullptr;
        Uint32 wavLength = 0;
        SDL_AudioStream* stream = nullptr;
        SDL_AudioDeviceID device = 0;
        bool isQueued = false; // true when wavBuffer has been queued to the stream/device
        std::queue<std::string> musicPlaylist; // Queue to manage songs to be played
};