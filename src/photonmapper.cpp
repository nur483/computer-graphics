/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/scene.h>
#include <nori/photon.h>

NORI_NAMESPACE_BEGIN

class PhotonMapper : public Integrator {
public:
    /// Photon map data structure
    typedef PointKDTree<Photon> PhotonMap;

    PhotonMapper(const PropertyList &props) {
        /* Lookup parameters */
        m_photonCount  = props.getInteger("photonCount", 1000000);
        m_photonRadius = props.getFloat("photonRadius", 0.0f /* Default: automatic */);
    }

    void preprocess(const Scene *scene) override {
        cout << "Gathering " << m_photonCount << " photons .. ";
        cout.flush();

        /* Create a sample generator for the preprocess step */
        Sampler *sampler = dynamic_cast<Sampler *>(
            NoriObjectFactory::createInstance("independent", PropertyList()));

        /* Allocate memory for the photon map */
        m_photonMap = std::unique_ptr<PhotonMap>(new PhotonMap());
        m_photonMap->reserve(m_photonCount);

		/* Estimate a default photon radius */
		if (m_photonRadius == 0)
			m_photonRadius = scene->getBoundingBox().getExtents().norm() / 500.0f;

        m_emittedCount = 0;

		int depositedPhotonsCount = 0;
		while (depositedPhotonsCount < m_photonCount) {
            ++m_emittedCount;

            Ray3f pathRay;
            Intersection xi;

            auto randomEmitter = scene->getRandomEmitter(sampler->next1D());
            Color3f W = randomEmitter->samplePhoton(pathRay, sampler->next2D(), sampler->next2D()) *
                     scene->getLights().size();

            while (true) {

                if (!scene->rayIntersect(pathRay, xi)) {
                    break;
                }

                if (xi.mesh->getBSDF()->isDiffuse()) {
                    m_photonMap->push_back(Photon(xi.p, -pathRay.d, W));
                    ++depositedPhotonsCount;
                }

                // russian roulette with success probability p
                auto p = std::min(W.maxCoeff(), .99f);
                if (sampler->next1D() > p) {
                    break;
                }
                W /= p;


                // Sample from BSDF
                BSDFQueryRecord bRec(xi.shFrame.toLocal(-pathRay.d));
                bRec.uv = xi.uv;
                auto bsdfCosThetaOverPdf = xi.mesh->getBSDF()->sample(bRec, sampler->next2D());
                W *= bsdfCosThetaOverPdf;

                pathRay = Ray3f(xi.p, xi.shFrame.toWorld(bRec.wo));
            }
        }

		/* Build the photon map */
        m_photonMap->build();
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &_ray) const override {

        Ray3f pathRay = _ray;
        Intersection xo;

        Color3f t(1);
        Color3f Li(0);

        while (true) {

            if (!scene->rayIntersect(pathRay, xo)) {
                break;
            }

            if (xo.mesh->isEmitter()) {
                EmitterQueryRecord eRec(pathRay.o, xo.p, xo.shFrame.n);
                Li += t * xo.mesh->getEmitter()->eval(eRec);
            }

            if (xo.mesh->getBSDF()->isDiffuse()) {
                Color3f photonDensityEstimation(0);
                std::vector<uint32_t> results;
                m_photonMap->search(xo.p, m_photonRadius,results);
                for (auto i : results) {
                    const Photon &photon = (*m_photonMap)[i];
                    BSDFQueryRecord bRec(xo.shFrame.toLocal(-pathRay.d), xo.shFrame.toLocal(photon.getDirection()), ESolidAngle);
                    bRec.uv = xo.uv;
                    auto fr = xo.mesh->getBSDF()->eval(bRec);
                    photonDensityEstimation += fr * photon.getPower();
                }
                photonDensityEstimation /= M_PI * pow(m_photonRadius, 2) * m_emittedCount;

                Li += t * photonDensityEstimation;
                break;
            }

            // russian roulette with success probability p
            auto p = std::min(t.maxCoeff(), .99f);
            if (sampler->next1D() > p) {
                break;
            }
            t /= p;


            // Sample from BSDF
            BSDFQueryRecord bRec(xo.shFrame.toLocal(-pathRay.d));
            bRec.uv = xo.uv;
            auto bsdfCosThetaOverPdf = xo.mesh->getBSDF()->sample(bRec, sampler->next2D());
            t *= bsdfCosThetaOverPdf;

            pathRay = Ray3f(xo.p, xo.shFrame.toWorld(bRec.wo));
        }
		return Li;
    }

    std::string toString() const override {
        return tfm::format(
            "PhotonMapper[\n"
            "  photonCount = %i,\n"
            "  photonRadius = %f\n"
            "]",
            m_photonCount,
            m_photonRadius
        );
    }
private:
    /* 
     * Important: m_photonCount is the total number of photons deposited in the photon map,
     * NOT the number of emitted photons. You will need to keep track of those yourself.
     */ 
    int m_photonCount;
    int m_emittedCount;
    float m_photonRadius;
    std::unique_ptr<PhotonMap> m_photonMap;
};

NORI_REGISTER_CLASS(PhotonMapper, "photonmapper");
NORI_NAMESPACE_END
