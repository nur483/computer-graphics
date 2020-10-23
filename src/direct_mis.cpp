#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class DirectMisIntegrator : public Integrator {
public:
    DirectMisIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const override {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its)){
            return {0};
        }

        Color3f Le(0);

        // Le(p, w0)
        if (its.mesh->isEmitter()) {
            EmitterQueryRecord lRec(ray.o, its.p, its.shFrame.n);
            Le = its.mesh->getEmitter()->eval(lRec);
        }

        auto valueEm = Color3f(0);
        auto valueMat = Color3f(0);

        // Ems
        {
            for (auto light : scene->getLights()) {
                EmitterQueryRecord lRec;
                lRec.ref = its.p;
                auto LeOverPdf = light->sample(lRec, sampler->next2D());
                auto pdfEm = light->pdf(lRec);

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
                auto pdfMat = its.mesh->getBSDF()->pdf(bsdfRec);

                auto weight = 0.f;
                if (pdfEm + pdfMat != 0) {
                    weight = pdfEm / (pdfEm + pdfMat);
                }
                Color3f F = fr * LeOverPdf * cosTheta;
                valueEm += weight * F;
            }
        }


        // Mats
        {
            // Sample from BSDF
            BSDFQueryRecord bRec(its.shFrame.toLocal(-ray.d));
            bRec.uv = its.uv;
            auto frCosThetaOverPdf = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
            auto pdfMat = its.mesh->getBSDF()->pdf(bRec);

            // Cast ray based on sample
            Ray3f wi(its.p, its.shFrame.toWorld(bRec.wo));

            // Check for intersection with emitter
            Intersection itsWi;
            auto Le = Color3f(0);
            auto pdfEm = 0.f;
            if (scene->rayIntersect(wi, itsWi) && itsWi.mesh->isEmitter()) {
                // Eval emitter intersection
                EmitterQueryRecord lRec(wi.o, itsWi.p, itsWi.shFrame.n);
                Le = itsWi.mesh->getEmitter()->eval(lRec);
                pdfEm = itsWi.mesh->getEmitter()->pdf(lRec);
            }

            auto weight = 0.f;
            if (pdfEm + pdfMat != 0) {
                weight = pdfMat / (pdfEm + pdfMat);
            }
            valueMat = weight * Le * frCosThetaOverPdf;
        }

        return Le + valueEm + valueMat;
    }

    std::string toString() const override {
        return "DirectMisIntegrator[]";
    }

};

NORI_REGISTER_CLASS(DirectMisIntegrator, "direct_mis");
NORI_NAMESPACE_END
