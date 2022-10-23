#include "color.h"

#include <ostream>

std::ostream& operator << (std::ostream& os, const Color& color)
{
    os << "[" << color.c[0] << ", " << color.c[1] << ", " << color.c[2] << "]";

    return os;
}
