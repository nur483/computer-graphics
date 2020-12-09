#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class PathMatsIntegrator : public Integrator {
public:
    PathMatsIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const override {

        Color3f t(1);
        Color3f tNew;
        Color3f Li(0);

        Intersection x0;
        Ray3f pathRay = ray;

        while (true) {

            bool intersectsScene = scene->rayIntersect(pathRay, x0);

            float tMax = pathRay.maxt;
            float nearT, farT;
            auto medium = scene->getMedium(pathRay, nearT, farT);
            if (medium != nullptr && intersectsScene) {
                tMax = (x0.p - pathRay.o).norm();
            }

            if (medium == nullptr && !intersectsScene) {
                break;
            }

            // Sample free path
            Color3f Tr;
            MediumQueryRecord mRec(tMax);
            if (medium != nullptr) {
                Tr = medium->sampleFreePath(pathRay, sampler, mRec);
            }

            // Volume interaction
            if (mRec.hasInteraction) {

                // Multiply with transmittance
                tNew = t * Tr;
                Li += tNew;

                // Isotropic Scattering
                auto direction = Warp::squareToUniformSphere(sampler->next2D());

                pathRay = Ray3f(mRec.p, direction);
            }

            // Surface interaction
            else if (intersectsScene) {
                Color3f Le(0);
                if (x0.mesh->isEmitter()) {
                    EmitterQueryRecord lRec(pathRay.o, x0.p, x0.shFrame.n);
                    Le = x0.mesh->getEmitter()->eval(lRec);
                }

                Li += t * Le;

                // Sample from BSDF
                BSDFQueryRecord bRec(x0.shFrame.toLocal(-pathRay.d));
                bRec.uv = x0.uv;
                auto frCosThetaOverPdf = x0.mesh->getBSDF()->sample(bRec, sampler->next2D());
                tNew = t * frCosThetaOverPdf;

                pathRay = Ray3f(x0.p, x0.shFrame.toWorld(bRec.wo));
            }
            else {
                break;
            }

            // russian roulette with success probability p
            auto p = std::min(t.maxCoeff(), .99f);
            if (sampler->next1D() > p)  {
                break;
            }
            t = tNew / p;
        }
        return Li;
    }

    std::string toString() const override {
        return "PathMatsIntegrator[]";
    }

};

NORI_REGISTER_CLASS(PathMatsIntegrator, "path_mats");
NORI_NAMESPACE_END
