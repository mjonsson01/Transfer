// File: Transfer/src/Utilities/TransferInputs.h

#include "TransferInputs.h"

std::ostream& operator<<(std::ostream& os, const TransferInputs& input) {
    os << "Active Inputs: ";

    bool any = false;

    auto printIf = [&](bool condition, const char* name) {
        if (condition) {
            if (any)
                os << ", ";
            os << name;
            any = true;
        }
    };
    // Mouse
    os << "Mouse Current Position: " << input.mouseCurrPosition << ", ";
    os << "Mouse Start Position: " << input.mouseDragStartPosition << ", ";
    printIf(input.leftMousePressed, "LeftMouse");
    printIf(input.rightMousePressed, "RightMouse");
    printIf(input.middleMousePressed, "MiddleMouse");

    // Keyboard SpaceCraft Controls
    printIf(input.wPressed, "W");
    printIf(input.aPressed, "A");
    printIf(input.sPressed, "S");
    printIf(input.dPressed, "D");
    printIf(input.spacePressed, "Space");

    // Menu / Media
    printIf(input.escPressed, "Esc");
    printIf(input.zeroPressed, "0");

    // Modifiers
    // printIf(input.ctrlPressed, "Ctrl"); //Not needed yet
    printIf(input.shiftPressed, "Shift");
    printIf(input.altPressed, "Alt");

    if (!any)
        os << "None";

    return os;
}
