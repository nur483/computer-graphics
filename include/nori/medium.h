#if !defined(__NORI_MEDIUM_H)
#define __NORI_MEDIUM_H

#include <nori/shape.h>
#include <openvdb/openvdb.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Data record for conveniently querying and sampling the medium
 */
struct MediumQueryRecord {
    /// Sampled interaction point
    Point3f p;

    /// Maximum length for free path
    float tMax;

    /// Whether the sample has an interaction, e.g. if t < tMax
    bool hasInteraction = false;

    explicit MediumQueryRecord(float tMax) : tMax(tMax) {}

};

class Medium : public NoriObject {

public:

    virtual Color3f sampleFreePath(const Ray3f& ray, Sampler *sampler, MediumQueryRecord& mRec) const = 0;

    virtual float Tr(const Ray3f &ray, Sampler *sampler, MediumQueryRecord& mRec) const = 0;

    std::string toString() const override = 0;

    bool rayIntersect(const Ray3f &ray, float &nearT, float &farT) const {
        return m_bbox.rayIntersect(ray, nearT, farT);
    }

    EClassType getClassType() const override {
        return EMedium;
    }

    BoundingBox3f getBbox() const {
        return m_bbox;
    }

private:
    BoundingBox3f m_bbox;

};

NORI_NAMESPACE_END
#endif //__
// NORI_MEDIUM_H
