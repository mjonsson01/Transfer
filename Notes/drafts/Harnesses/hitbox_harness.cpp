// Build + run: see README.md in this folder.
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
