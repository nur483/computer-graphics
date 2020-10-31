/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Romain Prévost

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

#include <nori/emitter.h>
#include <nori/warp.h>
#include <nori/shape.h>

NORI_NAMESPACE_BEGIN

class AreaEmitter : public Emitter {
public:
    AreaEmitter(const PropertyList &props) {
        m_radiance = props.getColor("radiance");
    }

    std::string toString() const override {
        return tfm::format(
                "AreaLight[\n"
                "  radiance = %s,\n"
                "]",
                m_radiance.toString());
    }

    Color3f eval(const EmitterQueryRecord & lRec) const override {
        if(!m_shape)
            throw NoriException("There is no shape attached to this Area light!");

        if (lRec.n.dot(lRec.wi) < 0) {
            return m_radiance;
        }
        return 0;
    }

    Color3f sample(EmitterQueryRecord & lRec, const Point2f & sample) const override {
        if(!m_shape)
            throw NoriException("There is no shape attached to this Area light!");

        ShapeQueryRecord sRec(lRec.ref);
        m_shape->sampleSurface(sRec, sample);
        lRec.p = sRec.p;
        lRec.n = sRec.n;

        auto direction = (lRec.p - lRec.ref);
        lRec.wi = direction.normalized();
        lRec.shadowRay = Ray3f(lRec.ref, lRec.wi, Epsilon, direction.norm() - Epsilon);
        lRec.pdf = pdf(lRec);
        if (lRec.pdf > 0) {
            return eval(lRec) / lRec.pdf;
        }
        return 0;
    }

    float pdf(const EmitterQueryRecord &lRec) const override {
        if(!m_shape)
            throw NoriException("There is no shape attached to this Area light!");

        auto dotProduct = lRec.n.dot(lRec.wi);
        if (dotProduct < 0) {
            ShapeQueryRecord sRec(lRec.ref, lRec.p);
            auto pdf = m_shape->pdfSurface(sRec);
            auto sNorm = (lRec.ref - lRec.p).squaredNorm();
            return pdf * sNorm / -dotProduct;
        }
        return 0;
    }


    Color3f samplePhoton(Ray3f &ray, const Point2f &sample1, const Point2f &sample2) const override {
        ShapeQueryRecord sRec;
        m_shape->sampleSurface(sRec, sample1);

        auto pdf = sRec.pdf;
        auto direction = Frame(sRec.n).toWorld(Warp::squareToCosineHemisphere(sample2));
        ray = Ray3f(sRec.p, direction);

        if (pdf == 0) {
            return 0;
        }

        // Constant pdf => pdf = 1 / area
        auto area = 1 / pdf;

        // Get Le
        EmitterQueryRecord eRec(sRec.p + direction, sRec.p, sRec.n);
        auto Le = eval(eRec);

        return M_PI * area * Le;
    }


protected:
    Color3f m_radiance;
};

NORI_REGISTER_CLASS(AreaEmitter, "area")
NORI_NAMESPACE_END