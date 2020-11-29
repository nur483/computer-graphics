#include <nori/color.h>
using namespace nori;

// Linear interpolation
static float mix(const float a, const float b, const float t) {
    return (1 - t) * a + t * b;
}
static Color3f mix(const Color3f& a, const Color3f& b, const float t) {
    return (1 - t) * a + t * b;
}

// Bilinear interpolation
static Color3f bilinear(const Color3f &q11, const Color3f &q12, const Color3f &q21, const Color3f &q22, const float& tu, const float& tv) {
    return mix(mix(q11, q21, tu), mix(q12, q22, tu), tv);
}