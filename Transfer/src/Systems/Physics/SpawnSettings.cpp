// File: Transfer/src/Core/SpawnSettings.cpp
#include "Core/SpawnSettings.hpp"

void SpawnSettings::setFlag(SpawnFlag flag, bool on)
{
    switch (flag)
    {
    case SpawnFlag::Static:
        is_force_static = on;
        break;
    case SpawnFlag::Shatterable:
        is_shatterable = on;
        break;
    case SpawnFlag::Accretable:
        is_accretable = on;
        break;
    case SpawnFlag::Bounce:
        is_bounce = on;
        break;
    case SpawnFlag::Ghost:
        is_ghost = on;
        break;
    }
    if (!on)
    {
        return; // switching off never switches anything else
    }

    // The settings this one makes meaningless (the rules are listed in SpawnSettings.hpp)
    switch (flag)
    {
    case SpawnFlag::Static:
        break;
    case SpawnFlag::Ghost:
        is_bounce = false;
        is_accretable = false;
        is_shatterable = false;
        break;
    case SpawnFlag::Bounce:
        is_ghost = false;
        is_accretable = false;
        is_shatterable = false;
        break;
    case SpawnFlag::Accretable:
    case SpawnFlag::Shatterable:
        is_ghost = false;
        is_bounce = false;
        break;
    }
}

bool SpawnSettings::flag(SpawnFlag flag) const
{
    switch (flag)
    {
    case SpawnFlag::Static:
        return is_force_static;
    case SpawnFlag::Shatterable:
        return is_shatterable;
    case SpawnFlag::Accretable:
        return is_accretable;
    case SpawnFlag::Bounce:
        return is_bounce;
    case SpawnFlag::Ghost:
        return is_ghost;
    }
    return false;
}
