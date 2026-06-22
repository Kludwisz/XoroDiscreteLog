#ifndef XOROSHIRO_H
#define XOROSHIRO_H

#include <cstdint>
#include <functional>

// Xoroshiro 128++ implementation taken from the Cubiomes library by Cubitect, see:
// https://github.com/Cubitect/cubiomes/blob/master/rng.h

/* 
MIT License

Copyright (c) 2020 Cubitect

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
///=============================================================================

struct Xoroshiro {
    uint64_t lo, hi;

    bool operator==(const Xoroshiro& other) const noexcept {
        return lo == other.lo && hi == other.hi;
    }
};

struct XoroshiroHash {
    std::size_t operator()(const Xoroshiro& x) const noexcept {
        std::size_t h1 = std::hash<uint64_t>{}(x.lo);
        std::size_t h2 = std::hash<uint64_t>{}(x.hi);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
};

static inline uint64_t rotl64(uint64_t x, uint8_t b)
{
    return (x << b) | (x >> (64-b));
}

static inline void xSetSeed(Xoroshiro *xr, uint64_t value)
{
    const uint64_t XL = 0x9e3779b97f4a7c15ULL;
    const uint64_t XH = 0x6a09e667f3bcc909ULL;
    const uint64_t A = 0xbf58476d1ce4e5b9ULL;
    const uint64_t B = 0x94d049bb133111ebULL;
    uint64_t l = value ^ XH;
    uint64_t h = l + XL;
    l = (l ^ (l >> 30)) * A;
    h = (h ^ (h >> 30)) * A;
    l = (l ^ (l >> 27)) * B;
    h = (h ^ (h >> 27)) * B;
    l = l ^ (l >> 31);
    h = h ^ (h >> 31);
    xr->lo = l;
    xr->hi = h;
}

static inline uint64_t xNextLong(Xoroshiro *xr)
{
    uint64_t l = xr->lo;
    uint64_t h = xr->hi;
    uint64_t n = rotl64(l + h, 17) + l;
    h ^= l;
    xr->lo = rotl64(l, 49) ^ h ^ (h << 21);
    xr->hi = rotl64(h, 28);
    return n;
}

static inline int xNextInt(Xoroshiro *xr, uint32_t n)
{
    uint64_t r = (xNextLong(xr) & 0xFFFFFFFF) * n;
    if ((uint32_t)r < n)
    {
        while ((uint32_t)r < (~n + 1) % n)
        {
            r = (xNextLong(xr) & 0xFFFFFFFF) * n;
        }
    }
    return r >> 32;
}

static inline int xNextIntFast(Xoroshiro *xr, uint32_t n)
{
    uint64_t r = (xNextLong(xr) & 0xFFFFFFFF) * n;
    return r >> 32;
}

static inline double xNextDouble(Xoroshiro *xr)
{
    return (xNextLong(xr) >> (64-53)) * 1.1102230246251565E-16;
}

static inline float xNextFloat(Xoroshiro *xr)
{
    return (xNextLong(xr) >> (64-24)) * 5.9604645E-8F;
}

static inline void xSkipN(Xoroshiro *xr, int count)
{
    while (count --> 0)
        xNextLong(xr);
}

static inline uint64_t xNextLongJ(Xoroshiro *xr)
{
    int32_t a = xNextLong(xr) >> 32;
    int32_t b = xNextLong(xr) >> 32;
    return ((uint64_t)a << 32) + b;
}

static inline int xNextIntJ(Xoroshiro *xr, uint32_t n)
{
    int bits, val;
    const int m = n - 1;

    if ((m & n) == 0) {
        uint64_t x = n * (xNextLong(xr) >> 33);
        return (int) ((int64_t) x >> 31);
    }

    do {
        bits = (xNextLong(xr) >> 33);
        val = bits % n;
    }
    while (bits - val + m < 0);
    return val;
}

#endif