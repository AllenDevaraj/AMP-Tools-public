#include "MyBugAlgorithm.h"
#include <algorithm>
#include <cmath>

#define _USE_MATH_DEFINES
#include <math.h>

// Constants
constexpr double MyBugAlgorithm::STEP_SIZE;
constexpr double MyBugAlgorithm::CONTACT_TOL;
constexpr double MyBugAlgorithm::GOAL_RADIUS;
constexpr double MyBugAlgorithm::HIT_RADIUS;
constexpr double MyBugAlgorithm::WALL_DIST;
constexpr double MyBugAlgorithm::SENSOR_RANGE;
constexpr double MyBugAlgorithm::IMPROVE_MIN;
constexpr double MyBugAlgorithm::LEAVE_EPS;
constexpr double MyBugAlgorithm::GOAL_EPS;
constexpr double MyBugAlgorithm::MIN_LOOP_ARC;
constexpr int    MyBugAlgorithm::MAX_OUTER_ITERS;
constexpr int    MyBugAlgorithm::MAX_STEPS;

// Helpers
namespace {

    using Vec = Eigen::Vector2d;
    constexpr double EPS = 1e-9;

    inline double dist(const Vec& a, const Vec& b) {
        return (a - b).norm();
    }

    inline Vec unit(const Vec& v) {
        double n = v.norm();
        if (n > 0.0) return v / n;
        return Vec(1.0, 0.0);
    }

    inline double angle(const Vec& v) {
        return std::atan2(v.y(), v.x());
    }

    inline Vec ang2vec(double th) {
        return Vec(std::cos(th), std::sin(th));
    }

    inline double angDiff(double a, double b) {
        double d = std::fmod(a - b + M_PI, 2.0 * M_PI);
        if (d < 0) d += 2.0 * M_PI;
        return d - M_PI;
    }

    inline void append(std::vector<Vec>& path, const Vec& p, double tol = 1e-9) {
        if (path.empty()) {
            path.push_back(p);
            return;
        }
        if ((path.back() - p).norm() > tol) {
            path.push_back(p);
        }
    }

    // Segment Intersection checker
    bool segIntersect(const Vec& p, const Vec& q_cur, const Vec& a, const Vec& b) {
        const Vec r = q_cur - p, s = b - a;
        const double denom = r.x()*s.y() - r.y()*s.x();
        const double num_t = (a.x()-p.x())*s.y() - (a.y()-p.y())*s.x();
        const double num_u = (a.x()-p.x())*r.y() - (a.y()-p.y())*r.x();
        if (std::abs(denom) < EPS) return false;
        const double t = num_t / denom, u = num_u / denom;
        return (t >= -EPS && t <= 1.0+EPS && u >= -EPS && u <= 1.0+EPS);
    }

    // CW ordering of the vertices
    std::vector<Vec> vertsCW(const amp::Obstacle2D& obs) {
        std::vector<Vec> V = obs.verticesCW();
        if (V.size() < 3) return V;
        double A = 0.0;
        for (size_t i=0;i<V.size();++i){
            const auto& p = V[i]; const auto& q_cur = V[(i+1)%V.size()];
            A += p.x()*q_cur.y() - q_cur.x()*p.y();
        }
        if (A > 0) std::reverse(V.begin(), V.end());
        return V;
    }

    // Point inside obstacle checker
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

    // Is point free
    bool freePoint(const Vec& p, const std::vector<amp::Obstacle2D>& obs) {
        for (const auto& O: obs) {
            auto V = vertsCW(O);
            if (V.size() < 3) continue;
            if (pointInObs(p, V)) return false;
        }
        return true;
    }

    // Is segment free
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
}

// Range sensor
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
        out.distance = range; 
    }
    return out;
}

// Bug 1 wall follow
MyBugAlgorithm::WallFollowResult
MyBugAlgorithm::Bug1WallFollow(const Vec& hit_point,
                               const Vec& initial_heading_right,
                               const amp::Problem2D& problem)
{
    const bool followRight = (side_ == FollowSide::Right);

    double th_side, th_away;
    double CORNER_SIGN, CORR_SIGN;

    if (followRight) {
        th_side   = -M_PI / 2.0;   // look toward wall (right-hand)
        th_away   = +M_PI / 2.0;   // rotate away from wall
        CORNER_SIGN = -1.0;          // -10° bias toward wall
        CORR_SIGN = -1.0;          // correction sign
    } else {
        th_side   = +M_PI / 2.0;   // look toward wall (left-hand)
        th_away   = -M_PI / 2.0;   // rotate away from wall
        CORNER_SIGN = +1.0;          // +10° bias toward wall
        CORR_SIGN = +1.0;          // correction sign
    }

    const Eigen::Rotation2Dd Rside(th_side);  // toward wall
    const Eigen::Rotation2Dd Raway(th_away);  // away from wall

    WallFollowResult R;
    R.Lstar = hit_point;
    double best_d = dist(hit_point, problem.q_goal);

    Vec q_cur  = hit_point;
    Vec hd;
    if (side_ == FollowSide::Right) {
        hd = unit(initial_heading_right);
    } else {
        hd = unit(-initial_heading_right);   
    }
    append(R.path_gone, q_cur);

    const int    MIN_STEPS_BEFORE_CLOSE = 120;
    const double P_GAIN            = 0.20;

    double loop_len = 0.0;
    double angTurned = 0.0;
    double last_th  = angle(hd);

    int steps = 0;
    while (steps++ < MAX_STEPS) {
        Vec dir_front = hd;
        Vec dir_side  = Rside * hd;

        auto front = checkSensor(q_cur, dir_front, WALL_DIST,    problem.obstacles);
        auto side  = checkSensor(q_cur, dir_side,  SENSOR_RANGE, problem.obstacles);

        if (front.found) {
            hd = unit(Raway * hd);                               // inner corner
        } else if (!side.found || side.distance > WALL_DIST) {
            Eigen::Rotation2Dd corner(CORNER_SIGN * (M_PI / 18.0));  // gone off the edge
            hd = unit(corner * hd);
        } else {
            double error = WALL_DIST - side.distance;            // maintain wall distance parallely 
            Eigen::Rotation2Dd tiny(CORR_SIGN * (-P_GAIN * error));
            hd = unit(tiny * hd);
        }

        // step; resolve residual collisions
        Vec next = q_cur + STEP_SIZE * hd;
        if (!freeSeg(q_cur, next, problem.obstacles)) {
            Vec mid = q_cur + 0.5*STEP_SIZE * hd;
            if (freeSeg(q_cur, mid, problem.obstacles)) next = mid;
            else { hd = unit(Raway * hd); next = q_cur + STEP_SIZE * hd; }
        }

        // summation of steps taken, angle turned
        double step_taken = (next - q_cur).norm();
        double th_now = angle(hd);
        double dth    = angDiff(th_now, last_th);
        angTurned    += dth;
        last_th       = th_now;

        q_cur = next;
        append(R.path_gone, q_cur);
        loop_len += step_taken;

        // update L*
        double d = dist(q_cur, problem.q_goal);
        if (d + 1e-12 < best_d - 1e-12) { best_d = d; R.Lstar = q_cur; }

        // loop close
        if (steps > MIN_STEPS_BEFORE_CLOSE &&
            loop_len > MIN_LOOP_ARC &&
            std::fabs(angTurned) > 2.0*M_PI*0.8 &&
            dist(q_cur, hit_point) < HIT_RADIUS) {
            break;
        }
    }

    R.arc_length        = loop_len;
    last_boundary_arc_  = loop_len;
    return R;
}

// Helper
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

// Bug 2 wall follow
void MyBugAlgorithm::Bug2WallFollow(Eigen::Vector2d& q_cur,
                                                const HitInfo& hit,
                                                const amp::Problem2D& P,
                                                std::vector<Eigen::Vector2d>& out_trace)
{
    const bool followRight = (side_ == FollowSide::Right);
    double th_side, th_away;   // to rotate thetas
    double CORNER_SIGN, CORR_SIGN;

    if (followRight) {
        th_side   = -M_PI / 2.0;   // toward wall
        th_away   = +M_PI / 2.0;   // away from wall
        CORNER_SIGN = -1.0;          // corner sign
        CORR_SIGN = -1.0;          // correction sign
    } else {
        th_side   = +M_PI / 2.0;
        th_away   = -M_PI / 2.0;
        CORNER_SIGN = +1.0;
        CORR_SIGN = +1.0;
    }

    const Eigen::Rotation2Dd Rside(th_side);
    const Eigen::Rotation2Dd Raway(th_away);

    // heading direction initializing
    Vec hd = unit(Rside * unit(P.q_goal - q_cur));
    append(out_trace, q_cur);

    const int    MIN_STEPS_BEFORE_CLOSE = 120;
    const double P_GAIN            = 0.20;

    double loop_len = 0.0, total_turn = 0.0, last_th = angle(hd);

    // Conditions for wall following
    int steps = 0;
    while (steps++ < MAX_STEPS) {
        Vec dir_front = hd;
        Vec dir_side  = Rside * hd;

        auto front = checkSensor(q_cur, dir_front, WALL_DIST,    P.obstacles);
        auto side  = checkSensor(q_cur, dir_side,  SENSOR_RANGE, P.obstacles);

        if (front.found) {
            hd = unit(Raway * hd); // inner corner
        } else if (!side.found || side.distance > WALL_DIST) {
            Eigen::Rotation2Dd corner(CORNER_SIGN * (M_PI / 18.0));
            hd = unit(corner * hd); // gone off the edge
        } else {
            double error = WALL_DIST - side.distance;
            Eigen::Rotation2Dd tiny(CORR_SIGN * (-P_GAIN * error));
            hd = unit(tiny * hd); // maintain wall distance parallely
        }

        Vec next = q_cur + STEP_SIZE * hd;
        if (!freeSeg(q_cur, next, P.obstacles)) {
            Vec mid = q_cur + 0.5*STEP_SIZE * hd;
            if (freeSeg(q_cur, mid, P.obstacles)) next = mid;
            else { hd = unit(Raway * hd); next = q_cur + STEP_SIZE * hd; }
        }

        const double step_taken = (next - q_cur).norm();
        const double th_now = angle(hd);
        const double dth    = angDiff(th_now, last_th);
        total_turn += dth; last_th = th_now;

        q_cur = next;
        append(out_trace, q_cur);
        loop_len += step_taken;

        if (onMLineCloser(q_cur, P, hit)) { last_boundary_arc_ = loop_len; return; }

        if (steps > MIN_STEPS_BEFORE_CLOSE &&
            loop_len > MIN_LOOP_ARC &&
            std::fabs(total_turn) > 2.0*M_PI*0.8 &&
            dist(q_cur, hit.H) < HIT_RADIUS) {
            last_boundary_arc_ = loop_len;
            return;
        }
    }
}

amp::Path2D MyBugAlgorithm::plan(const amp::Problem2D& problem) {
    if (mode_ == BugMode::Bug1) {
        return planBug1(problem);
    } else {
        return planBug2(problem);
    }
}


// Bug 1 plan
amp::Path2D MyBugAlgorithm::planBug1(const amp::Problem2D& problem) {
    using Vec = Eigen::Vector2d;

    amp::Path2D path;
    std::vector<Vec>& W = path.waypoints;

    Vec cur = problem.q_init;
    append(W, cur);

    int outer_guard = 0;
    while (dist(cur, problem.q_goal) > GOAL_RADIUS && outer_guard++ < MAX_OUTER_ITERS) {
        // Go to goal until hit
        Vec toGoal = problem.q_goal - cur;
        SensorHit front = checkSensor(cur, toGoal, toGoal.norm(), problem.obstacles);

        if (!front.found) {
            append(W, problem.q_goal);
            cur = problem.q_goal;
            break;
        }

        Vec hit = cur + unit(toGoal) * front.distance;
        append(W, hit);
        cur = hit;

        // Fix initial heading for wall follow
        Vec u = unit(toGoal);
        double th = std::atan2(u.y(), u.x());
        if (side_ == FollowSide::Right) {
            th += -M_PI / 2.0;
        }
        else {                            
            th += +M_PI / 2.0;
        }

        Vec initial_heading_side = Vec(std::cos(th), std::sin(th));
        auto survey = Bug1WallFollow(hit, initial_heading_side, problem);

        // Log waypoints of wall follow
        const auto& loop  = survey.path_gone;
        const Vec   Lstar = survey.Lstar;
        for (size_t i = 1; i < loop.size(); ++i) append(W, loop[i]);

        // Find shortest route to backtrack
        auto nearestIdx = [&](const Vec& q)->size_t {
            size_t k = 0; double best = std::numeric_limits<double>::infinity();
            for (size_t i = 0; i < loop.size(); ++i) {
                double s = (loop[i] - q).squaredNorm();
                if (s < best) { 
                    best = s; k = i; 
                }
            }
            return k;
        };

        size_t n   = loop.size();
        size_t k   = nearestIdx(W.back());  
        size_t kL  = nearestIdx(Lstar);     

        while (k != kL) {
            k = (k + 1) % n;
            append(W, loop[k]);
        }

        append(W, Lstar);
        cur = Lstar;


        // Leave towards goal
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

        Vec leave = cur + LEAVE_EPS*outward + GOAL_EPS*unit(problem.q_goal - cur);
        append(W, leave);
        cur = leave;
    }

    path.valid = (dist(cur, problem.q_goal) <= GOAL_RADIUS);
    return path;
}

// Bug 2 plan
amp::Path2D MyBugAlgorithm::planBug2(const amp::Problem2D& P) {
    using Vec = Eigen::Vector2d;

    amp::Path2D path;
    std::vector<Vec>& W = path.waypoints;

    Vec cur = P.q_init;
    append(W, cur);

    int outer_guard = 0;
    while (dist(cur, P.q_goal) > GOAL_RADIUS && outer_guard++ < MAX_OUTER_ITERS) {
        // Go to goal until hit
        Vec toGoal = P.q_goal - cur;
        SensorHit front = checkSensor(cur, toGoal, toGoal.norm(), P.obstacles);

        if (!front.found) {
            append(W, P.q_goal);
            cur = P.q_goal;
            break;
        }

        Vec H = cur + unit(toGoal) * front.distance;
        append(W, H);
        cur = H;

        // Follow boundary until M-line gets closer or loop ends
        HitInfo hit;
        hit.H     = H;
        hit.rho_H = dist(H, P.q_goal);
        hit.A     = P.q_init;
        hit.B     = P.q_goal;
        hit.m_hat = unit(hit.B - hit.A);

        std::vector<Vec> trace;
        Bug2WallFollow(cur, hit, P, trace);

        // Traced boundary for visualization
        for (size_t i = 1; i < trace.size(); ++i) append(W, trace[i]);
 
        // Leave boundary towards goal
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

        Vec leave = cur + LEAVE_EPS*outward + GOAL_EPS*unit(P.q_goal - cur);
        append(W, leave);
        cur = leave;
    }

    path.valid = (dist(cur, P.q_goal) <= GOAL_RADIUS);
    return path;
}
