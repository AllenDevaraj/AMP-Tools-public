#include "MyBugAlgorithm.h"
#include <algorithm>
#include <cmath>

#define _USE_MATH_DEFINES
#include <math.h>

// ===== Link-time constants =====
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

// ===== Small helpers =====
namespace {
    using Vec2 = Eigen::Vector2d;
    constexpr double EPS = 1e-9;

    inline double dist(const Vec2& a, const Vec2& b) { return (a - b).norm(); }
    inline Vec2   unit(const Vec2& v) { double n = v.norm(); return (n > 0.0) ? (v / n) : Vec2(1,0); }
    inline double angOf(const Vec2& v) { return std::atan2(v.y(), v.x()); }
    inline Vec2   fromAng(double th) { return Vec2(std::cos(th), std::sin(th)); }
    inline double wrapDiff(double a, double b) {
        double d = std::fmod(a - b + M_PI, 2.0*M_PI);
        if (d < 0) d += 2.0*M_PI;
        return d - M_PI;
    }
    inline void pushUnique(std::vector<Vec2>& poly, const Vec2& q, double tol = 1e-9) {
        if (poly.empty() || (poly.back() - q).norm() > tol) poly.push_back(q);
    }

    // ===== Boolean collision layer (planner never iterates vertices to decide) =====

    // seg/seg (free helper)
    bool segSegIntersect_free(const Vec2& p, const Vec2& q,
                              const Vec2& a, const Vec2& b) {
        const Vec2 r = q - p, s = b - a;
        const double denom = r.x()*s.y() - r.y()*s.x();
        const double num_t = (a.x()-p.x())*s.y() - (a.y()-p.y())*s.x();
        const double num_u = (a.x()-p.x())*r.y() - (a.y()-p.y())*r.x();
        if (std::abs(denom) < EPS) return false;
        const double t = num_t / denom, u = num_u / denom;
        return (t >= -EPS && t <= 1.0+EPS && u >= -EPS && u <= 1.0+EPS);
    }

    // CW vertices (normalize sign for robustness)
    std::vector<Vec2> vertsCW(const amp::Obstacle2D& obs) {
        std::vector<Vec2> V = obs.verticesCW();
        if (V.size() < 3) return V;
        double A = 0.0;
        for (size_t i=0;i<V.size();++i){
            const auto& p = V[i]; const auto& q = V[(i+1)%V.size()];
            A += p.x()*q.y() - q.x()*p.y();
        }
        if (A > 0) std::reverse(V.begin(), V.end()); // ensure CW
        return V;
    }

    // point-in-polygon (ray crossing)
    bool pointInPoly(const Vec2& P, const std::vector<Vec2>& V) {
        bool inside = false; const size_t n = V.size();
        for (size_t i=0, j=n-1; i<n; j=i++){
            const Vec2& A = V[j]; const Vec2& B = V[i];
            const bool cond = ((A.y() > P.y()) != (B.y() > P.y())) &&
                              (P.x() < (B.x()-A.x())*(P.y()-A.y())/((B.y()-A.y())+0.0) + A.x());
            if (cond) inside = !inside;
        }
        return inside;
    }

    bool freePoint(const Vec2& p, const std::vector<amp::Obstacle2D>& obs) {
        for (const auto& O: obs) {
            auto V = vertsCW(O);
            if (V.size() < 3) continue;
            if (pointInPoly(p, V)) return false;
        }
        return true;
    }

    bool freeSeg(const Vec2& a, const Vec2& b, const std::vector<amp::Obstacle2D>& obs) {
        if (!freePoint(a, obs) || !freePoint(b, obs)) return false;
        for (const auto& O: obs) {
            auto V = vertsCW(O); const size_t n = V.size(); if (n < 3) continue;
            for (size_t i=0;i<n;++i) {
                if (segSegIntersect_free(a, b, V[i], V[(i+1)%n])) return false;
            }
        }
        return true;
    }
} // anonymous namespace

// ===== Boolean “range” sensor (march + bisection) =====
MyBugAlgorithm::SensorHit
MyBugAlgorithm::checkSensor(const Vec2& pos, const Vec2& dir, double range,
                            const std::vector<amp::Obstacle2D>& obstacles)
{
    SensorHit out;
    const Vec2 u = unit(dir);
    const int N  = std::max(1, int(std::ceil(range / STEP_SIZE)));

    Vec2 a = pos;
    for (int i=0; i<N; ++i) {
        Vec2 b = pos + u * (STEP_SIZE * (i+1));
        if (!freeSeg(a, b, obstacles)) {
            // localize first contact on [a,b]
            double lo = 0.0, hi = STEP_SIZE;
            for (int it=0; it<40; ++it) {
                double mid = 0.5*(lo+hi);
                Vec2 m = a + u*mid;
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

// ===== Right-hand wall follow (clockwise) — full Bug-1 survey =====
MyBugAlgorithm::WallFollowResult
MyBugAlgorithm::performBug1WallFollow(const Vec2& hit_point,
                                      const Vec2& initial_heading_right,
                                      const amp::Problem2D& problem)
{
    WallFollowResult R;
    R.closest_point_to_goal = hit_point;
    double best_d = dist(hit_point, problem.q_goal);

    Vec2 q  = hit_point;
    Vec2 hd = unit(initial_heading_right);    // tangent heading (clockwise/right-hand)
    pushUnique(R.discovery_path, q);

    const int    MIN_STEPS_BEFORE_CLOSE = 120; // stricter than before
    const double H_CORR_GAIN            = 0.20;

    // *** new: robust loop-completion guards ***
    double loop_len = 0.0;                 // arc length walked
    double cum_turn = 0.0;                 // cumulative signed heading change (≈ 2π for a loop)
    double last_th  = angOf(hd);

    int steps = 0;
    while (steps++ < MAX_WALK_STEPS) {
        // Boolean sensors
        Eigen::Rotation2Dd Rright(-M_PI/2.0);
        Eigen::Rotation2Dd Rleft( +M_PI/2.0);
        Vec2 dir_front = hd;
        Vec2 dir_right = Rright * hd;

        auto front = checkSensor(q, dir_front, WALL_DIST,    problem.obstacles);
        auto right = checkSensor(q, dir_right, SENSOR_RANGE, problem.obstacles);

        if (front.found) {
            // inner corner
            hd = unit(Rleft * hd);
        } else if (!right.found || right.distance > WALL_DIST * 2.0) {
            // lost wall / too far
            hd = unit(Rright * hd);
        } else {
            // hold stand-off
            double error = WALL_DIST - right.distance;        // +ve if too close
            Eigen::Rotation2Dd tiny(-H_CORR_GAIN * error);    // steer toward/away
            hd = unit(tiny * hd);
        }

        // Step; if still colliding, shrink / rotate out slightly
        Vec2 next = q + STEP_SIZE * hd;
        if (!freeSeg(q, next, problem.obstacles)) {
            Vec2 mid = q + 0.5*STEP_SIZE * hd;
            if (freeSeg(q, mid, problem.obstacles)) next = mid;
            else { hd = unit(Rleft * hd); next = q + STEP_SIZE * hd; }
        }

        // accumulate arc and heading change
        double step_taken = (next - q).norm();
        double th_now = angOf(hd);
        double dth    = wrapDiff(th_now, last_th);
        cum_turn     += dth;
        last_th       = th_now;

        q = next;
        pushUnique(R.discovery_path, q);
        loop_len += step_taken;

        // Track L*
        double d = dist(q, problem.q_goal);
        if (d + 1e-12 < best_d - 1e-12) { best_d = d; R.closest_point_to_goal = q; }

        // Full loop detection: must satisfy all guards
        if (steps > MIN_STEPS_BEFORE_CLOSE &&
            loop_len > MIN_LOOP_ARC &&
            std::fabs(cum_turn) > 2.0*M_PI*0.9 &&   // ~360° of turning
            dist(q, hit_point) < HIT_RADIUS)
        {
            double dth_close = std::fabs(wrapDiff(angOf(hd), angOf(initial_heading_right)));
            if (dth_close < 15.0 * M_PI / 180.0) break;
        }
    }

    return R;
}

// ===== Bug-1 plan =====
amp::Path2D MyBugAlgorithm::plan(const amp::Problem2D& problem) {
    amp::Path2D path;
    std::vector<Vec2>& W = path.waypoints;

    Vec2 cur = problem.q_init;
    pushUnique(W, cur);

    int outer_guard = 0;
    while (dist(cur, problem.q_goal) > GOAL_RADIUS && outer_guard++ < MAX_OUTER_ITERS) {
        // 1) Go-to-goal until contact (boolean range)
        Vec2 toGoal = problem.q_goal - cur;
        SensorHit front = checkSensor(cur, toGoal, toGoal.norm(), problem.obstacles);

        if (!front.found) {
            pushUnique(W, problem.q_goal);
            cur = problem.q_goal;
            break;
        }

        // move to hit
        Vec2 hit = cur + unit(toGoal) * front.distance;
        pushUnique(W, hit);
        cur = hit;

        // 2) Full survey (RIGHT-hand follow)
        Eigen::Rotation2Dd turn_right(-M_PI / 2.0);
        Vec2 initial_heading_right = turn_right * unit(toGoal);
        auto survey = performBug1WallFollow(hit, initial_heading_right, problem);

        // If L* is not a genuine improvement, Bug-1 declares blocked
        if (dist(hit, problem.q_goal) - dist(survey.closest_point_to_goal, problem.q_goal) < IMPROVE_MIN) {
            break;
        }

        // 3) VISUALIZATION HOOK: append the FULL survey loop so the boundary is fully green
        const auto& loop  = survey.discovery_path;
        const Vec2  Lstar = survey.closest_point_to_goal;

        // append entire circumnavigation (loop[0] is the hit point you already pushed)
        for (size_t i = 1; i < loop.size(); ++i) pushUnique(W, loop[i]);

        // Find indices (nearest) of: where we are now (end of loop) and L*
        auto nearestIdx = [&](const Vec2& q)->size_t {
            size_t k = 0; double best = std::numeric_limits<double>::infinity();
            for (size_t i = 0; i < loop.size(); ++i) {
                double s = (loop[i] - q).squaredNorm();
                if (s < best) { best = s; k = i; }
            }
            return k;
        };

        const size_t startIdx = nearestIdx(W.back());  // current position (near hit after full loop)
        const size_t lIdx     = nearestIdx(Lstar);     // index of L*

        // Cyclic forward arc length from i -> j (wrap allowed)
        auto arcLenFwd = [&](size_t i, size_t j)->double {
            if (i == j) return 0.0;
            double L = 0.0;
            size_t n = loop.size();
            size_t k = i;
            while (k != j) {
                size_t kn = (k + 1) % n;
                L += (loop[kn] - loop[k]).norm();
                k = kn;
            }
            return L;
        };
        // Cyclic backward arc length from i -> j
        auto arcLenBack = [&](size_t i, size_t j)->double {
            if (i == j) return 0.0;
            double L = 0.0;
            size_t n = loop.size();
            size_t k = i;
            while (k != j) {
                size_t kp = (k + n - 1) % n;
                L += (loop[k] - loop[kp]).norm();
                k = kp;
            }
            return L;
        };

        // Choose shorter way along the boundary from where we are now to L*
        double Lf = arcLenFwd(startIdx, lIdx);
        double Lb = arcLenBack(startIdx, lIdx);

        if (Lf <= Lb) {
            // forward along loop
            size_t n = loop.size();
            size_t k = startIdx;
            while (k != lIdx) {
                k = (k + 1) % n;
                pushUnique(W, loop[k]);
            }
        } else {
            // backward along loop
            size_t n = loop.size();
            size_t k = startIdx;
            while (k != lIdx) {
                size_t kp = (k + n - 1) % n;
                pushUnique(W, loop[kp]);
                k = kp;
            }
        }

        // We are now exactly at L*
        pushUnique(W, Lstar);
        cur = Lstar;

        // 4) Leave boundary: outward (right-hand outward) + slight goal bias (unchanged)
        Vec2 outward(0,0);
        if (W.size() >= 2) {
            Vec2 lastDir = unit(W.back() - *(W.end()-2));
            double th_last = angOf(lastDir);
            for (int k=0; k<8; ++k) {
                double th = th_last - M_PI/2.0 + k*(M_PI/16.0);
                Vec2 probe = cur + 0.5*WALL_DIST*fromAng(th);
                if (!freePoint(probe, problem.obstacles)) { outward = unit(cur - probe); break; }
            }
            if (outward.norm() == 0) outward = fromAng(th_last + M_PI/2.0); // fallback
        } else {
            outward = unit(problem.q_goal - cur);
        }

        Vec2 leave = cur + LEAVE_EPS*outward + GOAL_BIAS_EPS*unit(problem.q_goal - cur);
        pushUnique(W, leave);
        cur = leave;

    }

    path.valid = (dist(cur, problem.q_goal) <= GOAL_RADIUS);
    return path;
}
// ====== BUG-2: plan (separate class) =================