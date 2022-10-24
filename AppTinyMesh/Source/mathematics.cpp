#include "mathematics.h"

static unsigned long xorshift96_x=123456789, xorshift96_y=362436069, xorshift_96z=521288629;
unsigned int Math::xorshift96(void)
{
    unsigned long t;

    xorshift96_x ^= xorshift96_x << 16;
    xorshift96_x ^= xorshift96_x >> 5;
    xorshift96_x ^= xorshift96_x << 1;

    t = xorshift96_x;
    xorshift96_x = xorshift96_y;
    xorshift96_y = xorshift_96z;
    xorshift_96z = t ^ xorshift96_x ^ xorshift96_y;

    return xorshift_96z;
}
