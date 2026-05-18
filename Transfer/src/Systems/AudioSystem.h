// File: Transfer/src/Systems/AudioSystem.h

#pragma once

// SDL3 Imports
#include <SDL3/SDL.h>

// Custom Imports
#include "Core/GameState.h"
#include "Core/UIState.h"
#include "Entities/Sound/MusicModeEnum.h"
#include "Entities/Sound/SoundEffect.h"
#include "Utilities/System/SystemPathUtility.h"

// Standard Library Imports
#include <filesystem>
#include <iostream>
#include <queue>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>
class AudioSystem
{
  public:
    // Constructor and Destructor
    AudioSystem();
    ~AudioSystem();

  public:
    // Main method to process audio each frame
    void ProcessSystemAudioFrame(GameState& gameState, UIState& UIState);

    void loadAndPlayTrack(const std::string& path);
    void addAllSoundEffectsToLibrary(const std::string& folderPath);
    void loadMusicLibrary(const std::string& folderPath);
    void onTrackFinished();
    void playSoundEffect(const std::string& name);
    void transitionMusicMode(MusicMode newMode);
    void stopCurrentMusic();
    void prepareShufflePlaylist();
    void playNextShuffleTrack();
    void cleanupFinishedSFX();

    // Clean up helper
    void CleanUp();

  private:
    std::mt19937 rng{std::random_device{}()};
    void processMusic();
    bool musicEnabled = false;                    // Will get populated by var in game state
    bool soundEffectsEnabled = false;             // Will get populated by var in game state
    MusicMode currentMusicMode = MusicMode::NONE; // owned by audio system
    SDL_AudioSpec wavSpec = {};
    SDL_AudioSpec outputSpec = {};
    Uint8* wavBuffer = nullptr;
    Uint32 wavLength = 0;
    SDL_AudioStream* musicStream = nullptr;
    SDL_AudioDeviceID device = 0;
    std::vector<std::string> shuffleTracks;
    std::string titleTrack;
    float masterVolume = 1.0f;
    float musicVolume = 1.0f;
    float soundEffectsVolume = 0.03125f;
    size_t currentShuffleIndex = 0;
    std::vector<SDL_AudioStream*> activeSFXStreams;
    std::unordered_map<std::string, SoundEffect> soundEffects;
};