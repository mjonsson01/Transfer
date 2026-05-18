// File: Transfer/src/Systems/AudioSystem.cpp

#include "Systems/AudioSystem.h"

AudioSystem::AudioSystem()
{
    // Initialize audio system variables if needed
    SDL_InitSubSystem(SDL_INIT_AUDIO);

    device = SDL_OpenAudioDevice(SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK, NULL);
    if (!device)
    {
        // SDL_Log("Failed to open audio device: %s", SDL_GetError());
    }
    else
    {
        // SDL_Log("Audio device opened successfully.");
    }
    SDL_GetAudioDeviceFormat(device, &outputSpec, nullptr);
    // SDL_Log("Output Frequency: %d", outputSpec.freq);
    // SDL_Log("Output Channels: %d", outputSpec.channels);

    loadMusicLibrary(Utilities::GetResourcePath("Music"));
    addAllSoundEffectsToLibrary(Utilities::GetResourcePath("SoundEffects"));
}

AudioSystem::~AudioSystem() {}

void AudioSystem::CleanUp()
{
    // SDL_Log("Cleaning up Audio System...");
    if (device)
        SDL_PauseAudioDevice(device);
    if (musicStream)
    {
        SDL_DestroyAudioStream(musicStream);
    }
    SDL_Delay(10);
    if (device)
    {
        SDL_CloseAudioDevice(device);
    }

    if (wavBuffer)
    {
        SDL_free(wavBuffer);
    }
    musicStream = nullptr;
    device = 0;
    wavBuffer = nullptr;
    wavLength = 0;
    for (SDL_AudioStream* sfxStream : activeSFXStreams)
    {
        SDL_DestroyAudioStream(sfxStream);
    }

    activeSFXStreams.clear();

    for (auto& [name, sfx] : soundEffects)
    {
        if (sfx.buffer)
        {
            SDL_free(sfx.buffer);
            sfx.buffer = nullptr;
        }
    }
    SDL_QuitSubSystem(SDL_INIT_AUDIO);
}

void AudioSystem::ProcessSystemAudioFrame(GameState& gameState, UIState& UIState)
{
    musicEnabled = UIState.getPlayMusic();
    soundEffectsEnabled = UIState.getPlaySoundEffects();
    MusicMode requested_music_mode = UIState.getRequestedMusicMode();
    if (requested_music_mode != currentMusicMode)
    {
        transitionMusicMode(requested_music_mode);
    }

    if (gameState.getIsShuttingDownAudioSystem())
        return;

    if (!device)
        return;
    // syncSettings(gameState);
    cleanupFinishedSFX();

    processMusic();
    while (UIState.HasPendingSoundEffects())
    {
        playSoundEffect(UIState.PopNextSoundEffect());
    }
}
void AudioSystem::cleanupFinishedSFX()
{
    for (auto it = activeSFXStreams.begin(); it != activeSFXStreams.end();)
    {
        SDL_AudioStream* stream = *it;

        if (SDL_GetAudioStreamAvailable(stream) == 0)
        {
            SDL_DestroyAudioStream(stream);
            it = activeSFXStreams.erase(it);
        }
        else
        {
            ++it;
        }
    }
}
void AudioSystem::transitionMusicMode(MusicMode newMode)
{
    stopCurrentMusic();

    currentMusicMode = newMode;

    switch (newMode)
    {
    case MusicMode::TITLE_THEME:
        loadAndPlayTrack(titleTrack);
        break;

    case MusicMode::MAIN_SHUFFLE:
        prepareShufflePlaylist();
        playNextShuffleTrack();
        break;

    default:
        break;
    }
}
void AudioSystem::prepareShufflePlaylist()
{
    std::shuffle(shuffleTracks.begin(), shuffleTracks.end(), rng);

    currentShuffleIndex = 0;
}
void AudioSystem::playNextShuffleTrack()
{
    if (shuffleTracks.empty())
        return;

    if (currentShuffleIndex >= shuffleTracks.size())
    {
        prepareShufflePlaylist();
    }

    loadAndPlayTrack(shuffleTracks[currentShuffleIndex]);

    currentShuffleIndex++;
}
void AudioSystem::processMusic()
{

    if (!musicEnabled)
        return;

    if (!musicStream)
        return;

    if (SDL_GetAudioStreamAvailable(musicStream) == 0)
    {
        onTrackFinished();
    }
}

void AudioSystem::onTrackFinished()
{
    switch (currentMusicMode)
    {
    case MusicMode::MAIN_SHUFFLE:
        playNextShuffleTrack();
        break;

    case MusicMode::TITLE_THEME:
        loadAndPlayTrack(titleTrack);
        break;

    default:
        break;
    }
}
void AudioSystem::loadMusicLibrary(const std::string& folderPath)
{
    if (!std::filesystem::exists(folderPath))
        return;

    for (const auto& entry : std::filesystem::directory_iterator(folderPath))
    {
        if (!entry.is_regular_file())
            continue;

        if (entry.path().extension() != ".wav")
            continue;

        std::string path = entry.path().string();
        std::string stem = entry.path().stem().string();

        if (stem == "TitleTheme")
        {
            titleTrack = path;
        }
        else
        {
            shuffleTracks.push_back(path);
        }
    }
}

void AudioSystem::stopCurrentMusic()
{
    if (musicStream)
    {
        SDL_DestroyAudioStream(musicStream);
        musicStream = nullptr;
    }

    if (wavBuffer)
    {
        SDL_free(wavBuffer);
        wavBuffer = nullptr;
    }

    wavLength = 0;
}
void AudioSystem::addAllSoundEffectsToLibrary(const std::string& folderPath)
{
    if (!std::filesystem::exists(folderPath) || !std::filesystem::is_directory(folderPath))
    {
        // SDL_Log("SoundEffects folder not found: %s", folderPath.c_str());
        return;
    }

    for (const auto& entry : std::filesystem::directory_iterator(folderPath))
    {
        if (entry.is_regular_file())
        {
            std::string path = entry.path().string();

            // Only accept .wav files
            if (entry.path().extension() == ".wav")
            {
                SoundEffect sfx;

                SDL_LoadWAV(path.c_str(), &sfx.spec, &sfx.buffer, &sfx.length);

                soundEffects[entry.path().stem().string()] = sfx;
                // SDL_Log("added to playlist: %s", path.c_str());
            }
        }
    }
}

void AudioSystem::loadAndPlayTrack(const std::string& path)
{
    stopCurrentMusic();

    if (!SDL_LoadWAV(path.c_str(), &wavSpec, &wavBuffer, &wavLength))
    {
        // SDL_Log("Failed to load music track: %s", path.c_str());
    }

    musicStream = SDL_CreateAudioStream(&wavSpec, &outputSpec);

    if (!musicStream)
    {
        // SDL_Log("Failed to create music musicStream: %s", SDL_GetError());
    }

    if (!SDL_BindAudioStream(device, musicStream))
    {
        // SDL_Log("Failed to bind music musicStream: %s", SDL_GetError());

        SDL_DestroyAudioStream(musicStream);
        musicStream = nullptr;
    }
    SDL_SetAudioStreamGain(musicStream, masterVolume * musicVolume);

    if (!SDL_PutAudioStreamData(musicStream, wavBuffer, wavLength))
    {
        // SDL_Log("Failed to queue music data: %s", SDL_GetError());
    }

    SDL_FlushAudioStream(musicStream);
}

void AudioSystem::playSoundEffect(const std::string& name)
{
    if (!soundEffectsEnabled)
        return;

    auto it = soundEffects.find(name);

    if (it == soundEffects.end())
    {
        // SDL_Log("Sound effect not found: %s", name.c_str());
        return;
    }

    SoundEffect& sfx = it->second;

    SDL_AudioStream* sfxStream = SDL_CreateAudioStream(&sfx.spec, &outputSpec);

    if (!sfxStream)
    {
        // SDL_Log("Failed to create SFX musicStream: %s", SDL_GetError());
        return;
    }

    if (!SDL_BindAudioStream(device, sfxStream))
    {
        // SDL_Log("Failed to bind SFX musicStream: %s", SDL_GetError());
        SDL_DestroyAudioStream(sfxStream);
        return;
    }
    SDL_SetAudioStreamGain(sfxStream, masterVolume * soundEffectsVolume);

    SDL_PutAudioStreamData(sfxStream, sfx.buffer, sfx.length);
    SDL_FlushAudioStream(sfxStream);

    activeSFXStreams.push_back(sfxStream);
}