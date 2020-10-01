#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class AvIntegrator : public Integrator {
public:
    AvIntegrator(const PropertyList &props) {
        m_length = props.getFloat("length");
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its)){
            return Color3f(1.0f);
        }

        /* Compute average visibility */
        Normal3f n = its.shFrame.n;
        Ray3f castedRay(its.p, Warp().sampleUniformHemisphere(sampler, n), Epsilon, m_length);

        if (scene->rayIntersect(castedRay, its)) {
            return Color3f(0.0f);
        } else {
            return Color3f(1.0);
        }
    }

    std::string toString() const {
        return "AvIntegrator[]";
    }

private:
    float m_length;

};

NORI_REGISTER_CLASS(AvIntegrator, "av");
NORI_NAMESPACE_END
