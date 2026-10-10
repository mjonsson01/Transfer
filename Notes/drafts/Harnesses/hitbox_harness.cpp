// Build (from the repo root, macOS):
//   clang++ -std=c++20 -O1 -w -ITransfer/src -IThirdParty/SDL3/include -IThirdParty/SDL3_ttf/include Notes/drafts/Harnesses/hitbox_harness.cpp $(find Transfer/src -name '*.cpp' ! -name main.cpp) -LThirdParty/SDL3/lib -LThirdParty/SDL3_ttf/lib -lSDL3 -lSDL3_ttf -Wl,-rpath,$PWD/ThirdParty/SDL3/lib -Wl,-rpath,$PWD/ThirdParty/SDL3_ttf/lib -o /tmp/hitbox_harness && /tmp/hitbox_harness
// The hitbox must face the same way as the ship (getPointingVector) at every rotation, and sit where the sprite is.
#include "Player/Starship.hpp"
#include "Core/UIState.hpp"
#include <cstdio>
#include <vector>
int main()
{
    Starship ship; UIState ui;
    DEPRECATED_InputState& in = ui.getMutableDEPRECATED_InputState();
    std::vector<DynamoEngine::Vector2D> p = ship.collisionPolygon();
    printf("rotation 0 corners:"); for (auto& c : p) printf(" (%.2f, %.2f)", c.x_val, c.y_val); printf("\n");
    double worst = 0;
    for (int step = 0; step < 130; ++step) // 130 x 0.05 rad = a full turn and a bit
    {
        p = ship.collisionPolygon();
        DynamoEngine::Vector2D center = (p[0] + p[3]) / 2.0;                      // wingtips' midpoint, on the centre line
        DynamoEngine::Vector2D nose = (p[1] + p[2]) / 2.0;
        DynamoEngine::Vector2D facing = (nose - center).normalize();
        DynamoEngine::Vector2D pointing = ship.getPointingVector();
        double err = (facing - pointing).magnitude();
        if (err > worst) worst = err;
        in.negativeRotation = true; ship.applyRotation(ui);
    }
    printf("worst |hitbox facing - getPointingVector| over a full turn: %.2e\n", worst);
}
