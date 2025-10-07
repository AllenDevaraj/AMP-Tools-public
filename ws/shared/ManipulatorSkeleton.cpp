#include "ManipulatorSkeleton.h"
#include <cmath>
#include <vector>
#include <algorithm>

// Helpers
namespace {
inline double wrapPi(double a) { return std::atan2(std::sin(a), std::cos(a)); }
inline double clampUnit(double v){ return v < -1.0 ? -1.0 : (v > 1.0 ? 1.0 : v); }

inline Eigen::Vector2d fkEE(const std::vector<double>& L, const amp::ManipulatorState& q){
    Eigen::Vector2d p(0.0, 0.0);
    double c = 0.0;
    for (std::size_t i=0;i<L.size();++i){
        c += q(i);
        p.x() += L[i]*std::cos(c);
        p.y() += L[i]*std::sin(c);
    }
    return p;
}

inline void solve2R(double L1, double L2, double tx, double ty, int elbowSign,
                    double& th1, double& th2){
    const double r2 = tx*tx + ty*ty;
    double c2 = clampUnit((r2 - L1*L1 - L2*L2) / (2.0*L1*L2));
    double s2 = std::sqrt(std::max(0.0, 1.0 - c2*c2));
    if (elbowSign < 0) s2 = -s2;
    th2 = std::atan2(s2, c2);
    th1 = std::atan2(ty, tx) - std::atan2(L2*s2, L1 + L2*c2);
    th1 = wrapPi(th1); th2 = wrapPi(th2);
}
} // namespace

MyManipulator2D::MyManipulator2D()
    : LinkManipulator2D({1.0, 1.0}) {}

MyManipulator2D::MyManipulator2D(const std::vector<double>& link_lengths)
    : LinkManipulator2D(link_lengths) {}

Eigen::Vector2d
MyManipulator2D::getJointLocation(const amp::ManipulatorState& state, uint32_t joint_index) const {
    if (joint_index == 0) return {0.0, 0.0};
    const auto& L = getLinkLengths();
    const uint32_t J = std::min<uint32_t>(joint_index, static_cast<uint32_t>(L.size()));

    Eigen::Vector2d p(0.0, 0.0);
    double th_sum = 0.0;
    for (uint32_t k = 0; k < J; ++k) {
        th_sum += state(k);
        p.x() += L[k] * std::cos(th_sum);
        p.y() += L[k] * std::sin(th_sum);
    }
    return p;
}

amp::ManipulatorState
MyManipulator2D::getConfigurationFromIK(const Eigen::Vector2d& ee) const {
    const auto& L = getLinkLengths();
    const std::size_t N = nLinks();
    amp::ManipulatorState q; q.resize(N); q.setZero();

    const double x = ee.x(), y = ee.y();
    const double r = std::hypot(x, y);

    if (N == 2){
        double th1, th2;
        solve2R(L[0], L[1], x, y, +1, th1, th2);
        q(0) = th1; q(1) = th2;
        return q;
    }

    if (N == 3){
        const double L1=L[0], L2=L[1], L3=L[2];

        if (r < 1e-12){
            q.setZero();
            return q;
        }

        const double rmin = std::fabs(L1 - L2), rmax = L1 + L2;

        const double phi0   = std::atan2(y, x);
        const double rw_des = std::max(0.0, r - L3);
        const double rw_goal= std::min(std::max(rw_des, rmin), rmax);

        const double denom  = std::max(1e-12, 2.0*r*L3);
        const double cosDel = clampUnit((r*r + L3*L3 - rw_goal*rw_goal)/denom);
        const double delta  = std::acos(cosDel);

        double bestErr = 1e9;
        amp::ManipulatorState best = q;

        for (int sPhi : {+1, -1}){
            const double phi = phi0 + sPhi*delta;
            const double xw = x - L3*std::cos(phi);
            const double yw = y - L3*std::sin(phi);

            for (int elbowSign : {+1, -1}){ 
                double th1, th2;
                solve2R(L1, L2, xw, yw, elbowSign, th1, th2);
                const double th3 = wrapPi(phi - th1 - th2);

                amp::ManipulatorState cand(3); cand.setZero();
                cand(0)=th1; cand(1)=th2; cand(2)=th3;

                const double err = (fkEE(L, cand) - ee).norm();
                if (err < bestErr){ bestErr = err; best = cand; }
            }
        }
        q = best;
        return q;
    }

    // N link
    {
        const double totalLen = std::accumulate(L.begin(), L.end(), 0.0);
        std::vector<Eigen::Vector2d> p(N+1);
        p[0] = Eigen::Vector2d(0.0, 0.0);
        for (std::size_t i=0;i<N;++i) p[i+1] = p[i] + Eigen::Vector2d(L[i], 0.0);

        if (r >= totalLen - 1e-12){
            const Eigen::Vector2d dir = (r>1e-12) ? (ee / r) : Eigen::Vector2d(1.0, 0.0);
            p[0] = Eigen::Vector2d(0.0, 0.0);
            for (std::size_t i=0;i<N;++i) p[i+1] = p[i] + dir * L[i];
        } else {
            const Eigen::Vector2d base = p[0];
            const double tol = 1e-8;
            const int maxIter = 300;

            for (int it=0; it<maxIter; ++it){
                p[N] = ee;
                for (int i = static_cast<int>(N)-1; i >= 0; --i){
                    Eigen::Vector2d dir = p[i] - p[i+1];
                    double d = dir.norm();
                    if (d < 1e-12) dir = Eigen::Vector2d(1.0, 0.0), d = 1.0;
                    dir /= d;
                    p[i] = p[i+1] + dir * L[i];
                }
                p[0] = base;
                for (std::size_t i=0;i<N;++i){
                    Eigen::Vector2d dir = p[i+1] - p[i];
                    double d = dir.norm();
                    if (d < 1e-12) dir = Eigen::Vector2d(1.0, 0.0), d = 1.0;
                    dir /= d;
                    p[i+1] = p[i] + dir * L[i];
                }
                if ( (p[N] - ee).norm() <= tol ) break;
            }
        }

        double acc = 0.0;
        for (std::size_t i=0;i<N;++i){
            Eigen::Vector2d vi = p[i+1] - p[i];
            double abs_i = std::atan2(vi.y(), vi.x());
            double rel_i = wrapPi(abs_i - acc);
            q(i) = rel_i;
            acc = wrapPi(acc + rel_i);
        }
        return q;
    }
}
