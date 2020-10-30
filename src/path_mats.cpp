#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class PathMatsIntegrator : public Integrator {
public:
    PathMatsIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const override {

        Color3f t(1);
        Color3f Li(0);

        Intersection x0;
        Ray3f pathRay = ray;

        while (true) {

            if (!scene->rayIntersect(pathRay, x0)) {
                break;
            }

            auto Le = Color3f(0);
            if (x0.mesh->isEmitter()) {
                EmitterQueryRecord lRec(pathRay.o, x0.p, x0.shFrame.n);
                Le = x0.mesh->getEmitter()->eval(lRec);
            }
            Li += t * Le;

            // russian roulette with success probability p
            auto p = std::min(t.maxCoeff(), 1.f);
            if (sampler->next1D() > p)  {
                break;
            }
            t /= p;


            // Sample from BSDF
            BSDFQueryRecord bRec(x0.shFrame.toLocal(-pathRay.d));
            bRec.uv = x0.uv;
            auto frCosThetaOverPdf = x0.mesh->getBSDF()->sample(bRec, sampler->next2D());
            t *= frCosThetaOverPdf;

            pathRay = Ray3f(x0.p, x0.shFrame.toWorld(bRec.wo));
        }
        return Li;
    }

    std::string toString() const override {
        return "PathMatsIntegrator[]";
    }

};

NORI_REGISTER_CLASS(PathMatsIntegrator, "path_mats");
NORI_NAMESPACE_END
