#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class DirectEmsIntegrator : public Integrator {
public:
    DirectEmsIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const override {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its)){
            return {0};
        }

        Color3f Lo(0);

        // Le(p, w0)
        if (its.mesh->isEmitter()) {
            EmitterQueryRecord lRec(ray.o, its.p, its.shFrame.n);
            auto Le = its.mesh->getEmitter()->eval(lRec);
            Lo += Le;
        }

        for (auto light : scene->getLights()) {
            EmitterQueryRecord lRec;
            lRec.ref = its.p;
            auto LeOverPdf = light->sample(lRec, sampler->next2D());

            // No direct ray to light source possible
            if (scene->rayIntersect(lRec.shadowRay)) {
                continue;
            }

            // Convert to local frame
            auto localRay = its.shFrame.toLocal(-ray.d); // wi
            auto localLRec = its.shFrame.toLocal(lRec.wi); // wo

            // Cosine value between shading normal and lRec
            auto cosTheta = Frame::cosTheta(localLRec);

            // Evaluate the BSDF for a pair of directions (lRec and ray)
            BSDFQueryRecord bsdfRec(localRay, localLRec, ESolidAngle);
            bsdfRec.uv = its.uv;
            auto fr = its.mesh->getBSDF()->eval(bsdfRec);

            auto F = fr * LeOverPdf * cosTheta;
            Lo += F;
        }

        return Lo;
    }

    std::string toString() const override {
        return "DirectEmsIntegrator[]";
    }

};

NORI_REGISTER_CLASS(DirectEmsIntegrator, "direct_ems");
NORI_NAMESPACE_END
