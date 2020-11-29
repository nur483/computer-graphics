#include <nori/emitter.h>
#include <nori/warp.h>
#include <nori/shape.h>
#include <nori/bitmap.h>
#include "utils.cpp"

NORI_NAMESPACE_BEGIN

class EnvironmentEmitter : public Emitter {

private:
    Bitmap m_envMap;
    int m_rows;
    int m_cols;
    Eigen::VectorXf m_pdfTheta;
    Eigen::VectorXf m_cdfTheta;
    Eigen::MatrixXf m_conditionalPdfPhi;
    Eigen::MatrixXf m_conditionalCdfPhi;

public:
    explicit EnvironmentEmitter(const PropertyList &props) {
        auto envMapPath = props.getString("envMapPath");
        m_envMap = Bitmap(envMapPath);
        m_rows = m_envMap.rows();
        m_cols = m_envMap.cols();

        // Precompute (marginal) pdf and cdf
        preCompute();
    }

    void preCompute() {

        // Luminance * sinTheta
        Eigen::MatrixXf luminance(m_rows, m_cols);
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                auto theta = M_PI * i / (m_rows - 1);
                luminance(i, j) = m_envMap(i, j).getLuminance() * sin(theta);
            }
        }

        // Compute marginals for theta
        m_pdfTheta = Eigen::VectorXf(m_rows);
        for (int i = 0; i < m_rows; ++i) {
            m_pdfTheta(i) = luminance.row(i).sum();
        }
        m_pdfTheta /= m_pdfTheta.sum();
        computeCdf(m_pdfTheta, m_cdfTheta);

        // Compute conditional marginals for phi
        m_conditionalPdfPhi = Eigen::MatrixXf(m_rows, m_cols);
        m_conditionalCdfPhi = Eigen::MatrixXf(m_rows, m_cols + 1);
        for (int i = 0; i < m_rows; ++i) {
            m_conditionalPdfPhi.row(i) = luminance.row(i) / luminance.row(i).sum();
            Eigen::VectorXf cdf;
            computeCdf(m_conditionalPdfPhi.row(i).transpose(), cdf);
            m_conditionalCdfPhi.row(i) = cdf.transpose();
        }
    }

    static void computeCdf(const Eigen::VectorXf &pdf, Eigen::VectorXf &cdf) {
        int n = pdf.size();
        cdf = Eigen::VectorXf(n + 1);
        cdf(0) = 0;
        for (int i = 1; i < n; ++i) {
            cdf(i) = cdf(i - 1) + pdf(i - 1);
        }
        cdf(n) = 1;
    }

    static int sampleDiscrete(const Eigen::VectorXf &cdf, const double &sample) {
        int n = cdf.size();
        int a = 0;
        int b = n - 1;
        int index = -1;
        while (a <= b) {
            int m = (a + b + 1) / 2;
            if (cdf(m) <= sample) {
                index = m;
                a = m + 1;
            }
            else {
                b = m - 1;
            }
        }
        return index;
    }

    Color3f eval(const EmitterQueryRecord & eRec) const override {
        auto thetaPhi = sphericalCoordinates(eRec.wi);

        float u, v;
        std::tie(u, v) = get_uv(eRec);

        // Perform a bilinear interpolation
        int i1 = clamp(int(floor(u)), 0, m_rows - 1);
        int i2 = clamp(int(floor(u)) + 1, 0, m_rows - 1);
        int j1 = clamp(int(floor(v)), 0, m_cols - 1);
        int j2 = clamp(int(floor(v)) + 1, 0, m_cols - 1);

        auto q11 = m_envMap(i1, j1);
        auto q12 = m_envMap(i1, j2);
        auto q21 = m_envMap(i2, j1);
        auto q22 = m_envMap(i2, j2);

        auto tu = u - i1;
        auto tv = v - j1;

        return bilinear(q11, q12, q21, q22, tu, tv);
    }

    Color3f sample(EmitterQueryRecord& eRec, const Point2f& sample) const override {

        auto i = sampleDiscrete(m_cdfTheta, sample.x());
        auto j = sampleDiscrete(m_conditionalCdfPhi.row(i).transpose(), sample.y());

        auto theta = M_PI * i / (m_rows - 1);
        auto phi = 2 * M_PI * j / (m_cols - 1);

        eRec.wi = sphericalDirection(theta, phi);
        Ray3f ray(eRec.ref, eRec.wi, Epsilon, std::numeric_limits<double>::infinity());

        // Calculate self intersection, set maxt accordingly
        float u, v, t;
        m_shape->rayIntersect(0, eRec.shadowRay, u, v, t);
        eRec.shadowRay.maxt = t - Epsilon;

        auto pdfValue = pdf(eRec);
        if (pdfValue == 0) {
            return 0;
        }

        auto J = (m_cols - 1) * (m_rows - 1) / (2 * M_PI * M_PI * Frame::sinTheta(eRec.wi));
        return eval(eRec) / (J * pdfValue);
    }

    float pdf(const EmitterQueryRecord &eRec) const override {
        int i, j;
        std::tie(i, j) = get_ij(eRec);

        if (m_pdfTheta(i) == 0) {
            return 0;
        }
        return m_pdfTheta(i) * m_conditionalPdfPhi(i, j);
    }

    std::tuple<float, float> get_uv(const EmitterQueryRecord &eRec) const {
        auto thetaPhi = sphericalCoordinates(eRec.wi);
        auto theta = thetaPhi.x();
        auto phi = thetaPhi.y();

        auto u = theta * (m_rows - 1) * INV_PI;
        auto v = phi  * 0.5 * (m_cols - 1) * INV_PI;

        return {u, v};
    }

    std::tuple<int, int> get_ij(const EmitterQueryRecord &eRec) const {

        float u, v;
        std::tie(u, v) = get_uv(eRec);

        auto i = int(round(u));
        auto j = int(round(v));

        return {i, j};
    }

    std::string toString() const override {
        return tfm::format(
                "EnvironmentEmitter[]"
        );
    }
};

NORI_REGISTER_CLASS(EnvironmentEmitter, "environment")
NORI_NAMESPACE_END