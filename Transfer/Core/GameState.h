#pragma once


struct InputState
{
	// Body instantiation control vars
	double selectedRadius = 0.0;
	double selectedMass = 0.0;
	bool bodySelectionValidity = false;
};
class GameState
{
	public:
		// Constructor and Destructor
		GameState();
		~GameState();


	public: 
		// Getters and Setters for Game State
		bool IsPlaying() const { return isPlaying; }
		void SetPlaying(bool playing) { isPlaying = playing; }

	private:
		bool isPlaying;
};