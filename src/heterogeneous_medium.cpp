#include <nori/medium.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class HeterogeneousMedium : public Medium {

private:
    openvdb::FloatGrid::Ptr m_density;
    Color3f m_sigmaA; // Absorption coefficient
    Color3f m_sigmaS; // Scattering coefficient
    Color3f m_sigmaT; // Extinction coefficient

    float m_maxDensity;

    BoundingBox3f m_bbox;
    BoundingBox3i m_bboxVoxelGrid;

public:

    explicit HeterogeneousMedium(const PropertyList &props) {

        m_sigmaA = props.getColor("sigma_a", 1);
        m_sigmaS = props.getColor("sigma_s", 1);
        m_sigmaT = m_sigmaS + m_sigmaA;

        auto size = props.getVector3("size", Vector3f(0.4)).cwiseAbs();
        auto center = props.getPoint3("center", Vector3f(0.f));
        m_bbox = BoundingBox3f(center - size / 2, center + size / 2);


        auto filePath = props.getString("vdb_path");

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

        cout << "Reading meta data" << endl;
        for (auto iter = m_density->beginMeta(); iter != m_density->endMeta(); ++iter) {
            const std::string& name = iter->first;
            openvdb::Metadata::Ptr value = iter->second;
            std::string valueAsString = value->str();
            std::cout << name << " = " << valueAsString << " (" << value->typeName() << ")" << std::endl;
        }
        auto bboxMin = m_density->getMetadata<openvdb::Vec3IMetadata>("file_bbox_min")->value();
        auto bboxMax = m_density->getMetadata<openvdb::Vec3IMetadata>("file_bbox_max")->value();
        BoundingBox3i bbox(
            {bboxMin.x(), bboxMin.y(), bboxMin.z()},
            {bboxMax.x(), bboxMax.y(), bboxMax.z()}
        );
        m_bboxVoxelGrid = bbox;

        file.close();

        m_maxDensity = 0;
        for (int x = bboxMin.x(); x < bboxMax.x(); ++x) {
            for (int y = bboxMin.y(); y < bboxMax.y(); ++y) {
                for (int z = bboxMin.z(); z < bboxMax.z(); ++z) {
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

        float nearT, farT;
        rayIntersect(ray, nearT, farT);
        float tMax = std::min(mRec.tMax, farT);

        while ((t += sampleDt(sampler)) < tMax) {
            if (evalDensity(ray(t)) / m_maxDensity > sampler->next1D()) {
                mRec.hasInteraction = true;
                mRec.p = ray(t);
                return m_sigmaS / m_sigmaT; // Real collision
            }
        }
        return 1; // No real collision
    }

    Color3f Tr(const Ray3f &ray, Sampler *sampler, MediumQueryRecord &mRec) const override {
        float t = Epsilon;
        Color3f tr = 1;

        float nearT, farT;
        rayIntersect(ray, nearT, farT);
        float tMax = std::min(mRec.tMax, farT);

        while ((t += sampleDt(sampler)) < tMax) {
            tr *= 1 - evalDensity(ray(t)) / m_maxDensity;
        }
        return tr;
    }

    float evalDensity(const Point3f &p) const {
        if (!m_bbox.contains(p)) {
            //return 0;
        }
        Point3f pGrid = (p - m_bbox.min).cwiseQuotient(m_bbox.max - m_bbox.min);

        Vector3i gridSize = m_bboxVoxelGrid.max - m_bboxVoxelGrid.min;

        int x = m_bboxVoxelGrid.min.x() + int(round(pGrid.x() * gridSize.x()));
        int y = m_bboxVoxelGrid.min.y() + int(round(pGrid.y() * gridSize.y()));
        int z = m_bboxVoxelGrid.min.z() + int(round(pGrid.z() * gridSize.z()));

        openvdb::Coord xyz(x, y, z);
        return m_density->getAccessor().getValue(xyz);
    }

    float sampleDt(Sampler *sampler) const {
        return -log(1 - sampler->next1D()) / (m_maxDensity * m_sigmaT.maxCoeff());
    }

    bool rayIntersect(const Ray3f &ray, float &nearT, float &farT) const override {
        return m_bbox.rayIntersect(ray, nearT, farT);
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