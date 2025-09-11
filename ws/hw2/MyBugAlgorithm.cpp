#include "MyBugAlgorithm.h"
#include <algorithm>
#include <cmath>

#define _USE_MATH_DEFINES
#include <math.h>

// ===== Constants (do not edit order/signatures) =====
constexpr double MyBugAlgorithm::STEP_SIZE;
constexpr double MyBugAlgorithm::CONTACT_TOL;
constexpr double MyBugAlgorithm::GOAL_RADIUS;
constexpr double MyBugAlgorithm::HIT_RADIUS;
constexpr double MyBugAlgorithm::WALL_DIST;
constexpr double MyBugAlgorithm::SENSOR_RANGE;
constexpr double MyBugAlgorithm::IMPROVE_MIN;
constexpr double MyBugAlgorithm::LEAVE_EPS;
constexpr double MyBugAlgorithm::GOAL_BIAS_EPS;
constexpr double MyBugAlgorithm::MIN_LOOP_ARC;
constexpr int    MyBugAlgorithm::MAX_OUTER_ITERS;
constexpr int    MyBugAlgorithm::MAX_WALK_STEPS;

// ===== File-local helpers =====
namespace {
    using Vec = Eigen::Vector2d;
    constexpr double EPS = 1e-9;

    inline double dist(const Vec& a, const Vec& b) { return (a - b).norm(); }
    inline Vec    unit(const Vec& v) { double n = v.norm(); return (n > 0.0) ? (v / n) : Vec(1,0); }
    inline double angle(const Vec& v) { return std::atan2(v.y(), v.x()); }
    inline Vec    ang2vec(double th) { return Vec(std::cos(th), std::sin(th)); }
    inline double angDiff(double a, double b) {
        double d = std::fmod(a - b + M_PI, 2.0*M_PI);
        if (d < 0) d += 2.0*M_PI;
        return d - M_PI;
    }
    inline void append(std::vector<Vec>& poly, const Vec& q, double tol = 1e-9) {
        if (poly.empty() || (poly.back() - q).norm() > tol) poly.push_back(q);
    }

    // Segment/segment
    bool segIntersect(const Vec& p, const Vec& q, const Vec& a, const Vec& b) {
        const Vec r = q - p, s = b - a;
        const double denom = r.x()*s.y() - r.y()*s.x();
        const double num_t = (a.x()-p.x())*s.y() - (a.y()-p.y())*s.x();
        const double num_u = (a.x()-p.x())*r.y() - (a.y()-p.y())*r.x();
        if (std::abs(denom) < EPS) return false;
        const double t = num_t / denom, u = num_u / denom;
        return (t >= -EPS && t <= 1.0+EPS && u >= -EPS && u <= 1.0+EPS);
    }

    // CW verts
    std::vector<Vec> vertsCW(const amp::Obstacle2D& obs) {
        std::vector<Vec> V = obs.verticesCW();
        if (V.size() < 3) return V;
        double A = 0.0;
        for (size_t i=0;i<V.size();++i){
            const auto& p = V[i]; const auto& q = V[(i+1)%V.size()];
            A += p.x()*q.y() - q.x()*p.y();
        }
        if (A > 0) std::reverse(V.begin(), V.end()); // ensure CW
        return V;
    }

    // Point in polygon (ray cross)
    bool pointInObs(const Vec& P, const std::vector<Vec>& V) {
        bool inside = false; const size_t n = V.size();
        for (size_t i=0, j=n-1; i<n; j=i++){
            const Vec& A = V[j]; const Vec& B = V[i];
            const bool cond = ((A.y() > P.y()) != (B.y() > P.y())) &&
                              (P.x() < (B.x()-A.x())*(P.y()-A.y())/((B.y()-A.y())+0.0) + A.x());
            if (cond) inside = !inside;
        }
        return inside;
    }

    bool freePoint(const Vec& p, const std::vector<amp::Obstacle2D>& obs) {
        for (const auto& O: obs) {
            auto V = vertsCW(O);
            if (V.size() < 3) continue;
            if (pointInObs(p, V)) return false;
        }
        return true;
    }

    bool freeSeg(const Vec& a, const Vec& b, const std::vector<amp::Obstacle2D>& obs) {
        if (!freePoint(a, obs) || !freePoint(b, obs)) return false;
        for (const auto& O: obs) {
            auto V = vertsCW(O); const size_t n = V.size(); if (n < 3) continue;
            for (size_t i=0;i<n;++i) {
                if (segIntersect(a, b, V[i], V[(i+1)%n])) return false;
            }
        }
        return true;
    }
} // anon

// ===== “Range” sensor =====
MyBugAlgorithm::SensorHit
MyBugAlgorithm::checkSensor(const Vec& pos, const Vec& dir, double range,
                            const std::vector<amp::Obstacle2D>& obstacles)
{
    SensorHit out;
    const Vec u = unit(dir);
    const int N  = std::max(1, int(std::ceil(range / STEP_SIZE)));

    Vec a = pos;
    for (int i=0; i<N; ++i) {
        Vec b = pos + u * (STEP_SIZE * (i+1));
        if (!freeSeg(a, b, obstacles)) {
            // localize first contact on [a,b]
            double lo = 0.0, hi = STEP_SIZE;
            for (int it=0; it<40; ++it) {
                double mid = 0.5*(lo+hi);
                Vec m = a + u*mid;
                if (freeSeg(a, m, obstacles)) lo = mid; else hi = mid;
                if (hi - lo < CONTACT_TOL) break;
            }
            out.found    = true;
            out.distance = (STEP_SIZE * i) + lo;
            return out;
        }
        a = b;
    }

    // clear to 'range'
    if (freeSeg(pos, pos + u*range, obstacles)) {
        out.found    = false;
        out.distance = std::numeric_limits<double>::infinity();
    } else {
        out.found    = true;
        out.distance = range;  // conservative
    }
    return out;
}

// ===== Bug-1 wall follow (right-hand) =====
MyBugAlgorithm::WallFollowResult
MyBugAlgorithm::Bug1WallFollow(const Vec& hit_point,
                               const Vec& initial_heading_right,
                               const amp::Problem2D& problem)
{
    WallFollowResult R;
    R.closest_point_to_goal = hit_point;
    double best_d = dist(hit_point, problem.q_goal);

    Vec q  = hit_point;
    Vec hd = unit(initial_heading_right);    // tangent heading (clockwise/right-hand)
    append(R.discovery_path, q);

    const int    MIN_STEPS_BEFORE_CLOSE = 120;
    const double H_CORR_GAIN            = 0.20;

    double loop_len = 0.0;
    double angTurned = 0.0;
    double last_th  = angle(hd);

    int steps = 0;
    while (steps++ < MAX_WALK_STEPS) {
        // Boolean sensors
        Eigen::Rotation2Dd Rright(-M_PI/2.0);
        Eigen::Rotation2Dd Rleft( +M_PI/2.0);
        Vec dir_front = hd;
        Vec dir_right = Rright * hd;

        auto front = checkSensor(q, dir_front, WALL_DIST,    problem.obstacles);
        auto right = checkSensor(q, dir_right, SENSOR_RANGE, problem.obstacles);

        if (front.found) {
            hd = unit(Rleft * hd);                  // inner corner
        } else if (!right.found || right.distance > WALL_DIST * 1.0) {
            Eigen::Rotation2Dd bias(-M_PI / 18.0);  // bias toward wall
            hd = unit(bias * hd);
        } else {
            double error = WALL_DIST - right.distance;        // +ve if too close
            Eigen::Rotation2Dd tiny(-H_CORR_GAIN * error);    // steer toward/away
            hd = unit(tiny * hd);
        }

        // Step; handle residual collisions
        Vec next = q + STEP_SIZE * hd;
        if (!freeSeg(q, next, problem.obstacles)) {
            Vec mid = q + 0.5*STEP_SIZE * hd;
            if (freeSeg(q, mid, problem.obstacles)) next = mid;
            else { Eigen::Rotation2Dd Rleft(+M_PI/2.0); hd = unit(Rleft * hd); next = q + STEP_SIZE * hd; }
        }

        // Accumulate
        double step_taken = (next - q).norm();
        double th_now = angle(hd);
        double dth    = angDiff(th_now, last_th);
        angTurned    += dth;
        last_th       = th_now;

        q = next;
        append(R.discovery_path, q);
        loop_len += step_taken;

        // Track L*
        double d = dist(q, problem.q_goal);
        if (d + 1e-12 < best_d - 1e-12) { best_d = d; R.closest_point_to_goal = q; }

        // Full loop detection
        if (steps > MIN_STEPS_BEFORE_CLOSE &&
            loop_len > MIN_LOOP_ARC &&
            std::fabs(angTurned) > 2.0*M_PI*0.8 &&
            dist(q, hit_point) < HIT_RADIUS) {
            break;
        }
    }
    return R;
}

// ===== Dispatcher =====
amp::Path2D MyBugAlgorithm::plan(const amp::Problem2D& problem) {
    return (mode_ == BugMode::Bug1) ? planBug1(problem)
                                    : planBug2(problem);
}

// ===== Bug-1 plan (your existing body, unchanged except name) =====
amp::Path2D MyBugAlgorithm::planBug1(const amp::Problem2D& problem) {
    using Vec = Eigen::Vector2d;

    amp::Path2D path;
    std::vector<Vec>& W = path.waypoints;

    Vec cur = problem.q_init;
    append(W, cur);

    int outer_guard = 0;
    while (dist(cur, problem.q_goal) > GOAL_RADIUS && outer_guard++ < MAX_OUTER_ITERS) {
        // 1) Go-to-goal until contact (boolean range)
        Vec toGoal = problem.q_goal - cur;
        SensorHit front = checkSensor(cur, toGoal, toGoal.norm(), problem.obstacles);

        if (!front.found) {
            append(W, problem.q_goal);
            cur = problem.q_goal;
            break;
        }

        // move to hit
        Vec hit = cur + unit(toGoal) * front.distance;
        append(W, hit);
        cur = hit;

        // 2) Full survey (RIGHT-hand follow)
        Eigen::Rotation2Dd turn_right(-M_PI / 2.0);
        Vec initial_heading_right = turn_right * unit(toGoal);
        auto survey = Bug1WallFollow(hit, initial_heading_right, problem);

        // If L* is not a genuine improvement, Bug-1 declares blocked
        if (dist(hit, problem.q_goal) - dist(survey.closest_point_to_goal, problem.q_goal) < IMPROVE_MIN) {
            break;
        }

        // 3) VIS: append the FULL survey loop so the boundary is fully green
        const auto& loop  = survey.discovery_path;
        const Vec   Lstar = survey.closest_point_to_goal;

        for (size_t i = 1; i < loop.size(); ++i) append(W, loop[i]);

        // Find indices near current and L*
        auto nearestIdx = [&](const Vec& q)->size_t {
            size_t k = 0; double best = std::numeric_limits<double>::infinity();
            for (size_t i = 0; i < loop.size(); ++i) {
                double s = (loop[i] - q).squaredNorm();
                if (s < best) { best = s; k = i; }
            }
            return k;
        };

        const size_t startIdx = nearestIdx(W.back());
        const size_t lIdx     = nearestIdx(Lstar);

        // Arc length helpers
        auto arcLenFwd = [&](size_t i, size_t j)->double {
            if (i == j) return 0.0;
            double L = 0.0; size_t n = loop.size(); size_t k = i;
            while (k != j) { size_t kn = (k + 1) % n; L += (loop[kn] - loop[k]).norm(); k = kn; }
            return L;
        };
        auto arcLenBack = [&](size_t i, size_t j)->double {
            if (i == j) return 0.0;
            double L = 0.0; size_t n = loop.size(); size_t k = i;
            while (k != j) { size_t kp = (k + n - 1) % n; L += (loop[k] - loop[kp]).norm(); k = kp; }
            return L;
        };

        // Shorter arc to L*
        double Lf = arcLenFwd(startIdx, lIdx);
        double Lb = arcLenBack(startIdx, lIdx);

        if (Lf <= Lb) {
            size_t n = loop.size(); size_t k = startIdx;
            while (k != lIdx) { k = (k + 1) % n; append(W, loop[k]); }
        } else {
            size_t n = loop.size(); size_t k = startIdx;
            while (k != lIdx) { size_t kp = (k + n - 1) % n; append(W, loop[kp]); k = kp; }
        }

        // At L*
        append(W, Lstar);
        cur = Lstar;

        // 4) Leave boundary (same as before)
        Vec outward(0,0);
        if (W.size() >= 2) {
            Vec lastDir = unit(W.back() - *(W.end()-2));
            double th_last = angle(lastDir);
            for (int k=0; k<8; ++k) {
                double th = th_last - M_PI/2.0 + k*(M_PI/16.0);
                Vec probe = cur + 0.5*WALL_DIST*ang2vec(th);
                if (!freePoint(probe, problem.obstacles)) { outward = unit(cur - probe); break; }
            }
            if (outward.norm() == 0) outward = ang2vec(th_last + M_PI/2.0);
        } else {
            outward = unit(problem.q_goal - cur);
        }

        Vec leave = cur + LEAVE_EPS*outward + GOAL_BIAS_EPS*unit(problem.q_goal - cur);
        append(W, leave);
        cur = leave;
    }

    path.valid = (dist(cur, problem.q_goal) <= GOAL_RADIUS);
    return path;
}

// ===== Bug-2 helpers =====
bool MyBugAlgorithm::onMLineCloser(const Eigen::Vector2d& p,
                                   const amp::Problem2D& P,
                                   const HitInfo& hit) const
{
    const Eigen::Vector2d AB = hit.B - hit.A;
    const double AB2 = AB.squaredNorm();
    if (AB2 == 0.0) return false;

    const Eigen::Vector2d AP = p - hit.A;
    const double dist_perp = std::abs(AP.x()*AB.y() - AP.y()*AB.x()) / std::sqrt(AB2);
    const double t = AP.dot(AB) / AB2;

    const double MLINE_EPS = 2e-2;   // tune if your units differ
    if (dist_perp > MLINE_EPS) return false;
    if (t < -1e-6 || t > 1.0 + 1e-6) return false;

    const double rho = (p - hit.B).norm();   // distance to goal
    return (rho + 1e-12 < hit.rho_H - IMPROVE_MIN);
}

void MyBugAlgorithm::followRightUntilMLineLeave(Eigen::Vector2d& q,
                                                const HitInfo& hit,
                                                const amp::Problem2D& P,
                                                std::vector<Eigen::Vector2d>& out_trace)
{
    using Vec = Eigen::Vector2d;
    Eigen::Rotation2Dd turn_right(-M_PI/2.0), turn_left(+M_PI/2.0);

    // Start heading = right-hand tangent w.r.t. goal direction at hit
    Vec hd = unit(turn_right * unit(P.q_goal - q));
    append(out_trace, q);

    const int    MIN_STEPS_BEFORE_CLOSE = 120;
    const double H_CORR_GAIN            = 0.20;

    double loop_len = 0.0, cum_turn = 0.0, last_th = angle(hd);

    int steps = 0;
    while (steps++ < MAX_WALK_STEPS) {
        // Sensors
        Vec dir_front = hd;
        Vec dir_right = turn_right * hd;

        auto front = checkSensor(q, dir_front, WALL_DIST,    P.obstacles);
        auto right = checkSensor(q, dir_right, SENSOR_RANGE, P.obstacles);

        if (front.found) {
            hd = unit(turn_left * hd);
        } else if (!right.found || right.distance > WALL_DIST * 1.0) {
            Eigen::Rotation2Dd bias(-M_PI / 18.0);
            hd = unit(bias * hd);
        } else {
            double error = WALL_DIST - right.distance;
            Eigen::Rotation2Dd tiny(-H_CORR_GAIN * error);
            hd = unit(tiny * hd);
        }

        // Step and fix residual collisions
        Vec next = q + STEP_SIZE * hd;
        if (!freeSeg(q, next, P.obstacles)) {
            Vec mid = q + 0.5*STEP_SIZE * hd;
            if (freeSeg(q, mid, P.obstacles)) next = mid;
            else { hd = unit(turn_left * hd); next = q + STEP_SIZE * hd; }
        }

        // accumulate
        const double step_taken = (next - q).norm();
        const double th_now = angle(hd);
        const double dth    = angDiff(th_now, last_th);
        cum_turn += dth; last_th = th_now;

        q = next;
        append(out_trace, q);
        loop_len += step_taken;

        // Bug-2 leave: on M-line and closer than at hit
        if (onMLineCloser(q, P, hit)) return;

        // Loop closed near H -> give up (blocked)
        if (steps > MIN_STEPS_BEFORE_CLOSE &&
            loop_len > MIN_LOOP_ARC &&
            std::fabs(cum_turn) > 2.0*M_PI*0.8 &&
            dist(q, hit.H) < HIT_RADIUS) {
            return;
        }
    }
}

// ===== Bug-2 plan =====
amp::Path2D MyBugAlgorithm::planBug2(const amp::Problem2D& P) {
    using Vec = Eigen::Vector2d;

    amp::Path2D path;
    std::vector<Vec>& W = path.waypoints;

    Vec cur = P.q_init;
    append(W, cur);

    int outer_guard = 0;
    while (dist(cur, P.q_goal) > GOAL_RADIUS && outer_guard++ < MAX_OUTER_ITERS) {
        // 1) Go-to-goal until contact
        Vec toGoal = P.q_goal - cur;
        SensorHit front = checkSensor(cur, toGoal, toGoal.norm(), P.obstacles);

        if (!front.found) {
            append(W, P.q_goal);
            cur = P.q_goal;
            break;
        }

        // Move to hitpoint H
        Vec H = cur + unit(toGoal) * front.distance;
        append(W, H);
        cur = H;

        // 2) Follow boundary until M-line (closer) or loop close
        HitInfo hit;
        hit.H     = H;
        hit.rho_H = dist(H, P.q_goal);
        hit.A     = P.q_init;
        hit.B     = P.q_goal;
        hit.m_hat = unit(hit.B - hit.A);

        std::vector<Vec> trace;
        followRightUntilMLineLeave(cur, hit, P, trace);

        // Append traced boundary (for visualization)
        for (size_t i = 1; i < trace.size(); ++i) append(W, trace[i]);

        // 3) If not on a valid leave point, we are blocked
        if (!onMLineCloser(cur, P, hit)) {
            break;
        }

        // 4) Leave boundary immediately toward goal (reuse Bug-1 nudge)
        Vec outward(0,0);
        if (W.size() >= 2) {
            Vec lastDir = unit(W.back() - *(W.end()-2));
            double th_last = angle(lastDir);
            for (int k=0; k<8; ++k) {
                double th = th_last - M_PI/2.0 + k*(M_PI/16.0);
                Vec probe = cur + 0.5*WALL_DIST*ang2vec(th);
                if (!freePoint(probe, P.obstacles)) { outward = unit(cur - probe); break; }
            }
            if (outward.norm() == 0) outward = ang2vec(th_last + M_PI/2.0);
        } else {
            outward = unit(P.q_goal - cur);
        }

        Vec leave = cur + LEAVE_EPS*outward + GOAL_BIAS_EPS*unit(P.q_goal - cur);
        append(W, leave);
        cur = leave;
    }

    path.valid = (dist(cur, P.q_goal) <= GOAL_RADIUS);
    return path;
}
