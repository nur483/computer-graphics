#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/texture.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief The Disney BRDF model
 */
class Disney : public BSDF {

private:
    Texture<Color3f> * m_albedo;
    Color3f m_specular;
    float m_metallic;
    float m_roughness;
    float m_subsurface;
    float m_specularTint;

public:
    Disney(const PropertyList &propList) : m_albedo(nullptr) {
        PropertyList l;
        l.setColor("value", propList.getColor("albedo", 0));
        m_albedo = static_cast<Texture<Color3f> *>(NoriObjectFactory::createInstance("constant_color", l));

        m_specular = propList.getColor("specular", 0);
        m_metallic = propList.getFloat("metallic", 0);
        m_roughness = propList.getFloat("roughness", 0);
        m_subsurface = propList.getFloat("subsurface", 0);
        m_specularTint = propList.getFloat("specularTint", 0);

    }
    ~Disney() override {
        delete m_albedo;
    }

    // Linear interpolation
    static Color3f mix(const float & a, const float& b, const float& t) {
        return (1 - t) * a + t * b;
    }
    static Color3f mix(const Color3f& a, const Color3f& b, const float& t) {
        return (1 - t) * a + t * b;
    }

    static float SchlickFresnel(const float& u) {
        float m = clamp(1 - u, 0.f, 1.f);
        return pow(m, 5);
    }

    static float smithGGX(float cosTheta_v, float alpha) {
        float a = alpha * alpha;
        float b = cosTheta_v * cosTheta_v;
        return 1 / (cosTheta_v + sqrt(a + b - a * b));
    }


    /// Evaluate the BRDF model
    Color3f eval(const BSDFQueryRecord &bRec) const override {

        // Using naming scheme from disney paper
        auto v = bRec.wi;
        auto l = bRec.wo;
        auto h = (v + l).normalized();

        // Compute cosTheta
        auto cosTheta_v = Frame::cosTheta(v);
        auto cosTheta_l = Frame::cosTheta(l);
        auto cosTheta_h = Frame::cosTheta(h);
        auto cosTheta_d = l.dot(h);

        // Inside
        if (cosTheta_v < 0 || cosTheta_l < 0) {
            return 0;
        }

        // Define colors
        auto baseColor = m_albedo->eval(bRec.uv);
        auto luminance = baseColor.getLuminance();
        auto tintColor = luminance > 0 ? Color3f(baseColor / luminance) : Color3f(1);
        auto specularColor = mix(m_specular * .08f * mix(1, tintColor, m_specularTint), baseColor, m_metallic);


        // Diffuse (d)
        auto Fl = SchlickFresnel(cosTheta_l);
        auto Fv = SchlickFresnel(cosTheta_v);
        auto Fd90 = .5f + 2 * m_roughness * cosTheta_d * cosTheta_d;
        auto f_diffuse = mix(1, Fd90, Fl) * mix(1, Fd90, Fv) * baseColor;


        // Subsurface (ss)
        auto Fss90 = cosTheta_d * cosTheta_d * m_roughness;
        auto Fss =  mix(1, Fss90, Fl) * mix(1, Fss90, Fv);
        auto f_subsurface = 1.25 * (Fss * (1 / (cosTheta_l + cosTheta_v) - .5) + .5) * baseColor;


        // Specular (s)
        auto alpha = std::max(0.001f, m_roughness * m_roughness);
        auto Ds = Warp::squareToGTR2Pdf(cosTheta_h, alpha);
        auto Fh = SchlickFresnel(cosTheta_d);
        auto Gs = smithGGX(cosTheta_l, alpha) * smithGGX(cosTheta_v, alpha);
        auto f_specular = mix(specularColor, 1, Fh);


        // Combine everything
        return mix(f_diffuse, f_subsurface, m_subsurface) * (1 - m_metallic) * INV_PI
                + f_specular * Gs * Ds;
    }

    /// Compute the density of \ref sample() wrt. solid angles
    virtual float pdf(const BSDFQueryRecord &bRec) const override {
        auto v = bRec.wi;
        auto l = bRec.wo;
        auto h = (v + l).normalized();

        auto cosTheta_l = Frame::cosTheta(l);
        auto cosTheta_h = Frame::cosTheta(h);

        if (cosTheta_l <= 0) {
            return 0;
        }

        auto alpha = std::max(0.001f, m_roughness * m_roughness);
        auto Jh = 1 / (4 * h.dot(l));
        auto Ds = Warp::squareToGTR2Pdf(cosTheta_h, alpha);

        return (1 - m_metallic) * cosTheta_l * INV_PI
                + m_metallic * Ds * cosTheta_h * Jh;
    }

    /// Draw a a sample from the BRDF model
    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const override {
        auto v = bRec.wi;
        Vector3f l, h;

        if (Frame::cosTheta(v) <= 0) {
            return 0;
        }

        auto diffuse = (1 - m_metallic);

        if (sample.x() < diffuse) {
            l = Warp::squareToCosineHemisphere({sample.x() / diffuse, sample.y()});
        }
        else {
            auto alpha = std::max(0.001f, m_roughness * m_roughness);
            auto newSample = Point2f((sample.x() - diffuse) / (1 - diffuse), sample.y());
            h = Warp::squareToGTR2(newSample, alpha);
            l = (2 * v.dot(h) * h - v).normalized();
        }
        bRec.wo = l;

        auto cosTheta = Frame::cosTheta(l);
        if (cosTheta <= 0) {
            return 0;
        }
        return eval(bRec) * cosTheta / pdf(bRec);
    }

    // Mainly for photon mapping
    bool isDiffuse() const override {
        return true;
    }

    void activate() override {
        if(!m_albedo) {
            PropertyList l;
            l.setColor("value", Color3f(0.5f));
            m_albedo = static_cast<Texture<Color3f> *>(NoriObjectFactory::createInstance("constant_color", l));
            m_albedo->activate();
        }
    }

    /// Add texture for the albedo
    void addChild(NoriObject *obj) override {
        switch (obj->getClassType()) {
            case ETexture:
                if(obj->getIdName() == "albedo") {
                    if (m_albedo)
                        throw NoriException("There is already an albedo defined!");
                    m_albedo = static_cast<Texture<Color3f> *>(obj);
                }
                else {
                    throw NoriException("The name of this texture does not match any field!");
                }
                break;

            default:
                throw NoriException("Disney::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    }

    /// Return a human-readable summary
    std::string toString() const override {
        return tfm::format(
            "Disney[\n"
            "  albedo = %s\n"
            "  specular = %s\n"
            "  specularTint = %s\n"
            "  metallic = %s\n"
            "  roughness = %s\n"
            "  subsurface = %s\n"
            "]",
            m_albedo->toString(),
            m_specular.toString(),
            m_specularTint,
            m_metallic,
            m_roughness,
            m_subsurface
        );
    }

    EClassType getClassType() const override { return EBSDF; }
};

NORI_REGISTER_CLASS(Disney, "disney");
NORI_NAMESPACE_END
