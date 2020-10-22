#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class DirectMatsIntegrator : public Integrator {
public:
    DirectMatsIntegrator(const PropertyList &props) {}

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

        // Sample from BSDF
        BSDFQueryRecord bRec(its.shFrame.toLocal(-ray.d));
        bRec.uv = its.uv;
        auto frCosThetaOverPdf = its.mesh->getBSDF()->sample(bRec, sampler->next2D());

        // Cast ray based on sample
        Ray3f wi(its.p, its.shFrame.toWorld(bRec.wo));

        // Check for intersection with emitter
        Intersection itsWi;
        if (!scene->rayIntersect(wi, itsWi) || !itsWi.mesh->isEmitter()) {
            return Lo;
        }

        // Eval emitter intersection
        EmitterQueryRecord lRec(wi.o, itsWi.p, itsWi.shFrame.n);
        auto Le = itsWi.mesh->getEmitter()->eval(lRec);


        Lo += Le * frCosThetaOverPdf;
        return Lo;
    }

    std::string toString() const override {
        return "DirectMatsIntegrator[]";
    }

};

NORI_REGISTER_CLASS(DirectMatsIntegrator, "direct_mats");
NORI_NAMESPACE_END
