//#include "EventListenerHelpers.h"
//
//void HandleBodyPropertyKeyEvent(const SDL_Event& event, double& selectedRadius, double& selectedMass)
//{
//	if (event.type == SDL_EVENT_KEY_DOWN)
//		switch (event.key.scancode) {
//			case SDL_SCANCODE_1:
//				selectedRadius = BodyRadiusOptions::Small;
//				break;
//			case SDL_SCANCODE_2:
//				selectedRadius = BodyRadiusOptions::Medium;
//				break;
//			case SDL_SCANCODE_3:
//				selectedRadius = BodyRadiusOptions::Large;
//				break;
//			case SDL_SCANCODE_Q:
//				selectedMass = BodyMassOptions::Light;
//				//std::cout << "Light mass selected and value is " << selectedMass << std::endl;
//				break;
//			case SDL_SCANCODE_W:
//				selectedMass = BodyMassOptions::Medium;
//				//std::cout << "Medium mass selected and value is " << selectedMass << std::endl;
//				break;
//			case SDL_SCANCODE_E:
//				selectedMass = BodyMassOptions::Heavy;
//				//std::cout << "Heavy mass selected and value is " << selectedMass << std::endl;
//				break;
//			default:
//				break;
//		}
//}
