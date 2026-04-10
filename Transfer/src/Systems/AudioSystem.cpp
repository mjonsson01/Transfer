// File: Transfer/src/Systems/AudioSystem.cpp

#include "Systems/AudioSystem.h"

AudioSystem::AudioSystem() {
    // Initialize audio system variables if needed
    SDL_InitSubSystem(SDL_INIT_AUDIO);

    device = SDL_OpenAudioDevice(SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK, NULL);
    if (!device) {
        SDL_Log("Failed to open audio device: %s", SDL_GetError());
    }
    AddAllMusicToPlaylist(Utilities::GetResourcePath("Music"));
    SDL_PauseAudioDevice(device);
}

AudioSystem::~AudioSystem() {}

void AudioSystem::CleanUp() {
    if (device)
        SDL_PauseAudioDevice(device);
    if (stream) {
        SDL_DestroyAudioStream(stream);
    }
    SDL_Delay(10);
    if (device) {
        SDL_CloseAudioDevice(device);
    }

    if (wavBuffer) {
        SDL_free(wavBuffer);
    }
    stream = nullptr;
    device = 0;
    wavBuffer = nullptr;
    wavLength = 0;
    SDL_QuitSubSystem(SDL_INIT_AUDIO);
}

void AudioSystem::ProcessSystemAudioFrame(GameState& gameState,
                                          UIState& UIState) {
    if (gameState.getIsShuttingDownAudioSystem())
        return;
    if (!device)
        return;

    if (gameState.getPlayMusic()) {
        if (!isQueued) {
            if (!musicPlaylist.empty()) {
                std::string nextTrack = musicPlaylist.front();
                musicPlaylist.pop();
                musicPlaylist.emplace(nextTrack);

                // Always try to load the next track
                if (LoadTrack(nextTrack)) {
                    if (wavBuffer && wavLength > 0 && stream) {
                        SDL_PutAudioStreamData(stream, wavBuffer, wavLength);
                        SDL_FlushAudioStream(stream);
                        isQueued = true;
                        SDL_ResumeAudioDevice(device);
                    }
                }
            }
        } else {
            // If the current track finished
            if (stream && SDL_GetAudioStreamAvailable(stream) == 0) {
                isQueued = false;
            }
        }
    } else {
        // Pause music
        if (stream)
            SDL_ClearAudioStream(stream);
        SDL_PauseAudioDevice(device);
        isQueued = false;
    }
}

bool AudioSystem::LoadTrack(const std::string& musicFilePath) {
    // Free previous buffer
    if (wavBuffer) {
        SDL_free(wavBuffer);
        wavBuffer = nullptr;
        wavLength = 0;
    }

    if (!SDL_LoadWAV(musicFilePath.c_str(), &wavSpec, &wavBuffer, &wavLength)) {
        SDL_Log("Failed to load WAV %s: %s", musicFilePath.c_str(),
                SDL_GetError());
        return false;
    }

    // Destroy old stream safely
    if (stream) {
        SDL_DestroyAudioStream(stream);
        stream = nullptr;
    }

    stream = SDL_CreateAudioStream(&wavSpec, &wavSpec);
    if (!stream) {
        SDL_Log("Failed to create audio stream: %s", SDL_GetError());
        SDL_free(wavBuffer);
        wavBuffer = nullptr;
        wavLength = 0;
        return false;
    }

    if (device) {
        SDL_BindAudioStream(device, stream);
    } else {
        SDL_Log("Audio device invalid!");
        SDL_DestroyAudioStream(stream);
        stream = nullptr;
        SDL_free(wavBuffer);
        wavBuffer = nullptr;
        wavLength = 0;
        return false;
    }

    return true;
}

void AudioSystem::AddAllMusicToPlaylist(const std::string& folderPath) {
    if (!std::filesystem::exists(folderPath) ||
        !std::filesystem::is_directory(folderPath)) {
        SDL_Log("Music folder not found: %s", folderPath.c_str());
        return;
    }

    for (const auto& entry : std::filesystem::directory_iterator(folderPath)) {
        if (entry.is_regular_file()) {
            std::string path = entry.path().string();

            // Only accept .wav files
            if (entry.path().extension() == ".wav") {
                musicPlaylist.push(path);
                // SDL_Log("Added to playlist: %s", path.c_str());
            }
        }
    }
}
