#include "SDL3/SDL.h"


class MassSlider : public UIElement
{
};

class PauseMenu : public UIElement
{
};

// class VelocityVectorToggle : public UIElement
// {
// };

// class GravityToggle : public UIElement
// {
// };

class FPSCounter : public UIElement
{
};


class UIElement{
    public:
        UIElement();
        ~UIElement();

        virtual void renderElement(SDL_Renderer* renderer);
    private:
        int posX;
        int posY;
        int width;
        int height;
}