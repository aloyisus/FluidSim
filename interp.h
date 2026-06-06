#include <openvdb/tools/Interpolation.h>
#include <cmath>



using namespace openvdb;

namespace fluidsim::tools {

template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

struct CubicSampler
{
    static const char* name() { return "cubic"; }
    static int radius() { return 1; }
    static bool mipmap() { return true; }
    static bool consistent() { return false; }
    static bool staggered() { return false; }
    static size_t order() { return 3; }

    template<class TreeT>
    static bool sample(const TreeT& inTree, const Vec3R& inCoord,
                       typename TreeT::ValueType& result);

    template<class TreeT>
    static typename TreeT::ValueType sample(const TreeT& inTree, const Vec3R& inCoord);

    template<class ValueT, size_t N>
    static inline ValueT tricubicInterpolation(ValueT (&data)[N][N][N], const Vec3R& uvw);
};


// Thanks to James Bird here https://www.jb101.co.uk/2020/12/27/monotone-cubic-interpolation.html
// for making me aware of this implementation used in Field3D - https://github.com/imageworks/Field3D

/*
 * Copyright (c) 2009 Sony Pictures Imageworks Inc
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the
 * distribution.  Neither the name of Sony Pictures Imageworks nor the
 * names of its contributors may be used to endorse or promote
 * products derived from this software without specific prior written
 * permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */
template<class T>
T monotonicCubicInterpolant(const T &f1, const T &f2, const T &f3, const T &f4, double t)
{
    T d_k = T(.5) * (f3 - f1);
    T d_k1 = T(.5) * (f4 - f2);
    T delta_k = f3 - f2;

    // This bit comes from Fix monotonic cubic interpolation #69 https://github.com/imageworks/Field3D/issues/69
    if (delta_k == static_cast<T>(0) || (sign(d_k) != sign(delta_k) || sign(d_k1) != sign(delta_k))) {
        d_k = static_cast<T>(0);
        d_k1 = static_cast<T>(0);
    }

    T a0 = f2;
    T a1 = d_k;
    T a2 = (T(3) * delta_k) - (T(2) * d_k) - d_k1;
    T a3 = d_k + d_k1 - (T(2) * delta_k);

    T t1 = t;
    T t2 = t1 * t1;
    T t3 = t2 * t1;

    return a3 * t3 + a2 * t2 + a1 * t1 + a0;
}

template<class TreeT>
inline bool
CubicSampler::sample(const TreeT& inTree, const Vec3R& inCoord,
    typename TreeT::ValueType& result)
{
    using ValueT = typename TreeT::ValueType;

    const Vec3i inIdx = openvdb::tools::local_util::floorVec3(inCoord);
    const Vec3i inLoIdx = inIdx - Vec3i(1, 1, 1);

    // Fractions
    const Vec3R uvw = inCoord - inIdx;

    bool active = false;
    ValueT data[4][4][4];
    for (int dx = 0, ix = inLoIdx.x(); dx < 4; ++dx, ++ix) {
        for (int dy = 0, iy = inLoIdx.y(); dy < 4; ++dy, ++iy) {
            for (int dz = 0, iz = inLoIdx.z(); dz < 4; ++dz, ++iz) {
                if (inTree.probeValue(Coord(ix, iy, iz), data[dx][dy][dz])) active = true;
            }
        }
    }

    result = CubicSampler::tricubicInterpolation(data, uvw);

    return active;
}

template<class TreeT>
inline typename TreeT::ValueType
CubicSampler::sample(const TreeT& inTree, const Vec3R& inCoord)
{

    using ValueT = typename TreeT::ValueType;

    const Vec3i inIdx = openvdb::tools::local_util::floorVec3(inCoord);
    const Vec3i inLoIdx = inIdx - Vec3i(1, 1, 1);

    // Fractions
    const Vec3R uvw = inCoord - inIdx;

    ValueT data[4][4][4];
    for (int dx = 0, ix = inLoIdx.x(); dx < 4; ++dx, ++ix) {
        for (int dy = 0, iy = inLoIdx.y(); dy < 4; ++dy, ++iy) {
            for (int dz = 0, iz = inLoIdx.z(); dz < 4; ++dz, ++iz) {
                data[dx][dy][dz] = inTree.getValue(Coord(ix, iy, iz));
            }
        }
    }

    return CubicSampler::tricubicInterpolation(data, uvw);
}

template<class ValueT, size_t N>
inline ValueT
CubicSampler::tricubicInterpolation(ValueT (&data)[N][N][N], const Vec3R& uvw)
{
    ValueT vx[4];
    ValueT vy[4];
    ValueT vz[4];
    for (int dx = 0; dx < 4; ++dx) {
        for (int dy = 0; dy < 4; ++dy) {
            for (int dz = 0; dz < 4; ++dz) {
                vz[dz] = data[dx][dy][dz];
            }
            vy[dy] = monotonicCubicInterpolant(vz[0],vz[1],vz[2],vz[3],uvw[2]);
        }
        vx[dx] = monotonicCubicInterpolant(vy[0],vy[1],vy[2],vy[3],uvw[1]);
    }
    return monotonicCubicInterpolant(vx[0],vx[1],vx[2],vx[3],uvw[0]);
}

}


namespace openvdb::OPENVDB_VERSION_NAME::tools {
    template<>
    struct Sampler<3, false> : public fluidsim::tools::CubicSampler {};
}