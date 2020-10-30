#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class PathMisIntegrator : public Integrator {
public:
    PathMisIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const override {

        Color3f t(1);
        Color3f Li(0);

        Intersection x0;
        Ray3f pathRay = ray;

        auto wEm = 1.f;
        auto wMat = 1.f;

        while (true) {

            if (!scene->rayIntersect(pathRay, x0)) {
                break;
            }

            auto Le = Color3f(0);
            if (x0.mesh->isEmitter()) {
                EmitterQueryRecord lRec(pathRay.o, x0.p, x0.shFrame.n);
                Le = x0.mesh->getEmitter()->eval(lRec);
            }
            // Contrib from mats
            Li += wMat * t * Le;

            // russian roulette with termination probability (1 - p)
            auto p = std::min(t.maxCoeff(), 1.f);
            if (sampler->next1D() > p) {
                break;
            }
            t /= p;

            // Contribution from emitter sampling
            auto light = scene->getRandomEmitter(sampler->next1D());
            EmitterQueryRecord lRec(x0.p);
            auto LeOverPdf = light->sample(lRec, sampler->next2D()) * scene->getLights().size();
            if (!scene->rayIntersect(lRec.shadowRay)) {
                auto localRay = x0.shFrame.toLocal(-pathRay.d); // wi
                auto localLRec = x0.shFrame.toLocal(lRec.wi); // wo
                auto cosTheta = Frame::cosTheta(localLRec);
                BSDFQueryRecord bsdfRec(localRay, localLRec, ESolidAngle);
                bsdfRec.uv = x0.uv;
                auto fr = x0.mesh->getBSDF()->eval(bsdfRec);

                auto pdfEm = light->pdf(lRec);
                auto pdfMat = x0.mesh->getBSDF()->pdf(bsdfRec);
                if (pdfEm + pdfMat != 0) {
                    wEm = pdfEm / (pdfEm + pdfMat);
                }

                Li += wEm * t * (fr * LeOverPdf * cosTheta);
            }


            // Sample from BSDF
            BSDFQueryRecord bRec(x0.shFrame.toLocal(-pathRay.d));
            bRec.uv = x0.uv;
            auto frCosThetaOverPdf = x0.mesh->getBSDF()->sample(bRec, sampler->next2D());

            pathRay = Ray3f(x0.p, x0.shFrame.toWorld(bRec.wo));
            t *= frCosThetaOverPdf;

            // Compute new wMat
            if (bRec.measure == EDiscrete) {
                wMat = 1;
            }
            else {
                Intersection its;
                if (scene->rayIntersect(pathRay, its)) {
                    if (its.mesh->isEmitter()) {
                        EmitterQueryRecord itsERec(x0.p, its.p, its.shFrame.n);
                        auto pdfMat = x0.mesh->getBSDF()->pdf(bRec);
                        auto pdfEm = its.mesh->getEmitter()->pdf(itsERec);
                        if (pdfEm + pdfMat > 0) {
                            wMat = pdfMat / (pdfEm + pdfMat);
                        }
                    }
                }
            }
        }
        return Li;
    }

    std::string toString() const override {
        return "PathMisIntegrator[]";
    }

};

NORI_REGISTER_CLASS(PathMisIntegrator, "path_mis");
NORI_NAMESPACE_END
