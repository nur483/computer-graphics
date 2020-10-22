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

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

Vector3f Warp::sampleUniformHemisphere(Sampler *sampler, const Normal3f &pole) {
    // Naive implementation using rejection sampling
    Vector3f v;
    do {
        v.x() = 1.f - 2.f * sampler->next1D();
        v.y() = 1.f - 2.f * sampler->next1D();
        v.z() = 1.f - 2.f * sampler->next1D();
    } while (v.squaredNorm() > 1.f);

    if (v.dot(pole) < 0.f)
        v = -v;
    v /= v.norm();

    return v;
}

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    auto r = sqrt(sample.x());
    auto phi = 2 * M_PI * sample.y();
    return r * Point2f(cos(phi), sin(phi));
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    if (p.norm() > 1) {
        return 0;
    }
    return INV_PI;
}

Vector3f Warp::squareToUniformCylinder(const Point2f &sample) {
    auto phi = 2 * M_PI * sample.x();
    return {cos(phi), sin(phi), 2 * sample.y() - 1};
}

Vector3f Warp::squareToUniformSphereCap(const Point2f &sample, float cosThetaMax) {
    auto cylinder = squareToUniformCylinder(sample);
    auto z = .5f * (cylinder.z() + 1) * (1 - cosThetaMax) + cosThetaMax;
    auto r = sqrt(1 - pow(z, 2.f));
    return {r * cylinder.x(), r * cylinder.y(), z};
}

float Warp::squareToUniformSphereCapPdf(const Vector3f &v, float cosThetaMax) {
    if (abs(1 - v.norm()) > Epsilon || v.z() < cosThetaMax) { 
        return 0;
    }
    return INV_TWOPI / (1 - cosThetaMax);
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    auto cylinder = squareToUniformCylinder(sample);
    auto r = sqrt(1 - pow(cylinder.z(), 2.f));
    return {r * cylinder.x(), r * cylinder.y(), cylinder.z()};
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    if (abs(1 - v.norm()) > Epsilon) { // |v| == 1
        return 0;
    }
    return INV_FOURPI;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    return squareToUniformSphereCap(sample, 0);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    return squareToUniformSphereCapPdf(v, 0);
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    auto phi = 2 * M_PI * sample.x();
    auto theta =  acos(sqrt(sample.y()));
    return {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    return squareToUniformHemispherePdf(v);
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    auto phi = 2 * M_PI * sample.x();
    auto theta = atan(sqrt(-pow(alpha, 2.f) * log(1 - sample.y())));
    return {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    auto theta = acos(m.z());
    if (abs(1 - m.norm()) > Epsilon || m.z() < 0) { 
        return 0;
    }
    auto a2 = pow(alpha, 2);
    return exp(-pow(tan(theta), 2) / a2) / (M_PI * a2 * pow(cos(theta), 3));
}

Vector3f Warp::squareToUniformTriangle(const Point2f &sample) {
    float su1 = sqrtf(sample.x());
    float u = 1.f - su1, v = sample.y() * su1;
    return Vector3f(u,v,1.f-u-v);
}

NORI_NAMESPACE_END
