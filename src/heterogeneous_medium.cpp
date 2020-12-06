#include <nori/medium.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class HeterogeneousMedium : public Medium {

private:
    openvdb::FloatGrid::Ptr m_density;
    float m_sigmaA; // Absorption coefficient
    float m_sigmaS; // Scattering coefficient
    float m_sigmaT; // Extinction coefficient

    float m_maxDensity;

public:

    explicit HeterogeneousMedium(const PropertyList &props) {
        auto filePath = props.getString("vdbPath");
        m_sigmaA = props.getFloat("sigmaA", 0);
        m_sigmaS = props.getFloat("sigmaS", 0);
        m_sigmaT = m_sigmaS + m_sigmaA;

        openvdb::initialize();
        openvdb::io::File file(filePath);
        file.open();

        for (auto nameIter = file.beginName(); nameIter != file.endName(); ++nameIter) {
            cout << "Found grid " << nameIter.gridName() << endl;
            if (nameIter.gridName() == "density") {
                auto baseGrid = file.readGrid(nameIter.gridName());
                m_density = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);
            }
        }
        openvdb::FloatGrid::Accessor accessor = m_density->getAccessor();
        auto bboxMin = m_density->metaValue<openvdb::Vec3I>("file_bbox_min");
        auto bboxMax = m_density->metaValue<openvdb::Vec3I>("file_bbox_max");
        getBbox().min << bboxMin.x(), bboxMin.y(), bboxMin.z();
        getBbox().max << bboxMax.x(), bboxMax.y(), bboxMax.z();

        file.close();

        m_maxDensity = 0;
        for (unsigned int x = bboxMin.x(); x < bboxMax.x(); ++x) {
            for (unsigned int y = bboxMin.y(); y < bboxMax.y(); ++y) {
                for (unsigned int z = bboxMin.z(); z < bboxMax.z(); ++z) {
                    auto density =  m_density->getAccessor().getValue(openvdb::Coord(x, y, z));
                    m_maxDensity = std::max(m_maxDensity, density);
                    if (density < 0) {
                        throw NoriException("A negative density value is not allowed.");
                    }
                }
            }
        }
        if (m_maxDensity == 0) {
            throw NoriException("The density grid need to have at least one positive value.");
        }
    }

    Color3f sampleFreePath(const Ray3f &ray, Sampler *sampler, MediumQueryRecord &mRec) const override {
        float t = Epsilon;

        while ((t += sampleDt(sampler)) < mRec.tMax) {
            if (evalDensity(ray(t)) / m_maxDensity > sampler->next1D()) {
                return m_sigmaS / m_sigmaT; // Real collision
            }
        }
        return 1; // No real collision
    }

    float Tr(const Ray3f &ray, Sampler *sampler, MediumQueryRecord &mRec) const override {
        float t = Epsilon;
        float tr = 1;

        while ((t += sampleDt(sampler)) < mRec.tMax) {
            tr *= 1 - evalDensity(ray(t)) / m_maxDensity;
        }
        return tr;
    }

    float evalDensity(const Point3f &p) const {
        if (!getBbox().contains(p)) {
            return 0;
        }
        openvdb::Coord xyz(p.x(), p.y(), p.z());
        return m_density->getAccessor().getValue(xyz);
    }

    float sampleDt(Sampler *sampler) const {
        return -log(1 - sampler->next1D()) / (m_maxDensity * m_sigmaT);
    }

    NoriObject::EClassType getClassType() const override {
        return EMedium;
    }

    std::string toString() const override {
        return "HeterogeneousMedium[]";
    }

};


NORI_REGISTER_CLASS(HeterogeneousMedium, "heterogeneous_medium")
NORI_NAMESPACE_END