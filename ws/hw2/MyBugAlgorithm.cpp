#include "MyBugAlgorithm.h"
#include <algorithm>
#include <cmath>

#define _USE_MATH_DEFINES
#include <math.h>

// ========== Link-time constants ==========
constexpr double MyBugAlgorithm::STEP_SIZE;
constexpr double MyBugAlgorithm::CONTACT_TOL;
constexpr double MyBugAlgorithm::GOAL_RADIUS;
constexpr double MyBugAlgorithm::HIT_RADIUS;
constexpr double MyBugAlgorithm::WALL_DIST;
constexpr double MyBugAlgorithm::SENSOR_RANGE;
constexpr double MyBugAlgorithm::IMPROVE_MIN;
constexpr double MyBugAlgorithm::LEAVE_EPS;
constexpr double MyBugAlgorithm::GOAL_BIAS_EPS;
constexpr int    MyBugAlgorithm::MAX_OUTER_ITERS;
constexpr int    MyBugAlgorithm::MAX_WALK_STEPS;

// ========== Small helpers ==========
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
    void pushUnique(std::vector<Vec2>& poly, const Vec2& q, double tol = 1e-9) {
        if (poly.empty() || (poly.back() - q).norm() > tol) poly.push_back(q);
    }

    // ===== Boolean collision layer (free helpers) =====

    // Segment/segment (free function, not a class member)
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

    // Ensure CW (negative signed area) for robustness in pointInPoly
    std::vector<Vec2> vertsCW(const amp::Obstacle2D& obs) {
        std::vector<Vec2> V = obs.verticesCW();   // AMP usually returns CW already
        if (V.size() < 3) return V;
        double A = 0.0;
        for (size_t i=0;i<V.size();++i){
            const auto& p = V[i]; const auto& q = V[(i+1)%V.size()];
            A += p.x()*q.y() - q.x()*p.y();
        }
        if (A > 0) std::reverse(V.begin(), V.end());
        return V;
    }

    // Point in polygon (ray crossing)
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

// ==================== Boolean “range” sensor ====================
// March forward in STEP_SIZE; when the next step would collide, bisect that
// step to localize contact to CONTACT_TOL and return the distance.
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
            // localize on [a,b]
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
        out.distance = range; // conservative fallback
    }
    return out;
}

// ========================= PLAN (Bug-1) =========================
amp::Path2D MyBugAlgorithm::plan(const amp::Problem2D& problem) {
    amp::Path2D path; 
    std::vector<Vec2>& W = path.waypoints;

    Vec2 cur = problem.q_init;
    pushUnique(W, cur);

    int outer_guard = 0;

    while (dist(cur, problem.q_goal) > GOAL_RADIUS && outer_guard++ < MAX_OUTER_ITERS) {
        // --- 1) Move toward goal until contact (boolean range sensor) ---
        Vec2 toGoal = problem.q_goal - cur;
        SensorHit front = checkSensor(cur, toGoal, toGoal.norm(), problem.obstacles);

        if (!front.found) {
            // clear shot to goal
            pushUnique(W, problem.q_goal);
            cur = problem.q_goal;
            break;
        }

        // move to hit point
        Vec2 hit = cur + unit(toGoal) * front.distance;
        pushUnique(W, hit);
        cur = hit;

        // --- 2) BUG-1: full circumnavigation (RIGHT-hand wall follow) ---
        // initial heading: rotate goal direction by -90° (clockwise) to keep wall on RIGHT
        Eigen::Rotation2Dd turn_right(-M_PI / 2.0);
        Vec2 initial_heading_right = turn_right * unit(toGoal);

        // survey loop and find L*
        auto survey = performBug1WallFollow(hit, initial_heading_right, problem);

        // If L* not meaningfully better than hit, declare blocked (Bug-1 behavior)
        if (dist(hit, problem.q_goal) - dist(survey.closest_point_to_goal, problem.q_goal) < IMPROVE_MIN) {
            break;
        }

        // --- 3) Walk the recorded loop from hit → L* ---
        // Find index of L* in the loop (nearest point)
        size_t j = 0; double best = std::numeric_limits<double>::infinity();
        for (size_t i=0;i<survey.discovery_path.size();++i) {
            double s = (survey.discovery_path[i] - survey.closest_point_to_goal).squaredNorm();
            if (s < best) { best = s; j = i; }
        }
        // Append hit→...→L*
        for (size_t i=0; i<=j && i<survey.discovery_path.size(); ++i) {
            pushUnique(W, survey.discovery_path[i]);
        }
        pushUnique(W, survey.closest_point_to_goal);
        cur = survey.closest_point_to_goal;

        // --- 4) Leave the boundary: outward (right-hand outward = to the RIGHT of motion) + slight goal bias ---
        Vec2 outward(0,0);
        if (W.size() >= 2) {
            Vec2 lastDir = unit(W.back() - *(W.end()-2));
            double th_last = angOf(lastDir);

            // Probe some right-leaning directions; choose the first that collides, and push opposite
            for (int k=0; k<8; ++k) {
                double th = th_last - M_PI/2.0 + k*(M_PI/16.0);       // directions on RIGHT side
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

// ================== Right-hand wall follow (clockwise) ==================
MyBugAlgorithm::WallFollowResult
MyBugAlgorithm::performBug1WallFollow(const Vec2& hit_point,
                                      const Vec2& initial_heading_right,
                                      const amp::Problem2D& problem)
{
    WallFollowResult R;
    R.closest_point_to_goal = hit_point;
    double best_d = dist(hit_point, problem.q_goal);

    Vec2 q  = hit_point;
    Vec2 hd = unit(initial_heading_right);    // tangent heading (clockwise / right-hand)
    pushUnique(R.discovery_path, q);

    const int    MIN_STEPS_BEFORE_CLOSE = 80;
    const double H_CORR_GAIN            = 0.30;   // small gain to hold WALL_DIST on right

    int steps = 0;
    while (steps++ < MAX_WALK_STEPS) {
        // Sensors (boolean range):
        Eigen::Rotation2Dd Rright(-M_PI/2.0);
        Eigen::Rotation2Dd Rleft( +M_PI/2.0);
        Vec2 dir_front = hd;
        Vec2 dir_right = Rright * hd;

        auto front = checkSensor(q, dir_front, WALL_DIST,    problem.obstacles);
        auto right = checkSensor(q, dir_right, SENSOR_RANGE, problem.obstacles);

        if (front.found) {
            // Inner corner approaching ⇒ turn LEFT in place (rotate heading)
            hd = unit(Rleft * hd);
        } else if (!right.found || right.distance > WALL_DIST * 2.0) {
            // Too far from wall on right ⇒ steer RIGHT to reacquire boundary
            hd = unit(Rright * hd);
        } else {
            // Hold desired stand-off to right with a small proportional rotation
            double error = WALL_DIST - right.distance;        // positive if too close
            Eigen::Rotation2Dd tiny(-H_CORR_GAIN * error);    // negative = rotate slightly toward right if too far
            hd = unit(tiny * hd);
        }

        // Advance one step along adjusted heading
        Vec2 next = q + STEP_SIZE * hd;
        // If the step would collide (rare after the logic above), reduce step once
        if (!freeSeg(q, next, problem.obstacles)) {
            Vec2 mid = q + 0.5*STEP_SIZE * hd;
            if (freeSeg(q, mid, problem.obstacles)) next = mid;
            else {
                // rotate left slightly to escape tight corner
                hd = unit(Rleft * hd);
                next = q + STEP_SIZE * hd;
            }
        }

        q = next;
        pushUnique(R.discovery_path, q);

        // Track L*
        double d = dist(q, problem.q_goal);
        if (d + 1e-12 < best_d - 1e-12) { best_d = d; R.closest_point_to_goal = q; }

        // Full loop detection: back near hit point with similar heading
        if (steps > MIN_STEPS_BEFORE_CLOSE && dist(q, hit_point) < HIT_RADIUS) {
            double dth = std::fabs(wrapDiff(angOf(hd), angOf(initial_heading_right)));
            if (dth < 15.0 * M_PI / 180.0) break;
        }
    }

    return R;
}
