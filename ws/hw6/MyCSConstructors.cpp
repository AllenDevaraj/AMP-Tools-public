#include "MyCSConstructors.h"
#include "HelpfulClass.h"

#include <queue>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <Eigen/Dense>

#define _USE_MATH_DEFINES
#include <math.h>

using MotionPlanningHelpers::CollisionChecker;
using MotionPlanningHelpers::Vec;

// =======================
// MyGridCSpace2D
// =======================
std::pair<std::size_t, std::size_t>
MyGridCSpace2D::getCellFromPoint(double x0, double x1) const {
    const double dx = (m_xmax - m_xmin) / static_cast<double>(m_nx);
    const double dy = (m_ymax - m_ymin) / static_cast<double>(m_ny);

    auto clamp_to = [](double v, double lo, double hi) {
        return std::max(lo, std::min(v, hi));
    };

    const double fx = std::floor((x0 - m_xmin) / dx);
    const double fy = std::floor((x1 - m_ymin) / dy);

    const std::size_t ix = static_cast<std::size_t>(clamp_to(fx, 0.0, static_cast<double>(m_nx - 1)));
    const std::size_t iy = static_cast<std::size_t>(clamp_to(fy, 0.0, static_cast<double>(m_ny - 1)));
    return {ix, iy};
}

// =======================
// MyPointAgentCSConstructor  (USED for Ex1)
// =======================
std::unique_ptr<amp::GridCSpace2D>
MyPointAgentCSConstructor::construct(const amp::Environment2D& env) {
    // Enforce ~0.25m cell size
    const double cell = 0.25;
    const std::size_t nx = static_cast<std::size_t>(std::ceil((env.x_max - env.x_min) / cell));
    const std::size_t ny = static_cast<std::size_t>(std::ceil((env.y_max - env.y_min) / cell));

    // Boolean occupancy grid
    auto cspace_ptr = std::make_unique<MyGridCSpace2D>(nx, ny, env.x_min, env.x_max, env.y_min, env.y_max);
    MyGridCSpace2D& cspace = *cspace_ptr;

    std::cout << "Constructing C-space for point agent (cell=0.25, 4-connected)\n";

    // Mark a cell occupied if its center collides in workspace
    for (std::size_t ix = 0; ix < nx; ++ix) {
        for (std::size_t iy = 0; iy < ny; ++iy) {
            const Eigen::Vector2d center = cspace.cellCenter(ix, iy);
            const Vec p(center.x(), center.y());
            const bool occupied = CollisionChecker::isPointInCollision(p, env.obstacles);
            cspace(ix, iy) = occupied; // true => blocked
        }
    }

    return cspace_ptr;
}

// =======================
// MyManipulatorCSConstructor  (NOW IMPLEMENTED for Ex2)
// =======================
std::unique_ptr<amp::GridCSpace2D>
MyManipulatorCSConstructor::construct(const amp::LinkManipulator2D& manipulator,
                                      const amp::Environment2D& env) {
    // Define the bounds for the C-space (joint angles from -PI to PI)
    const double theta_min = -M_PI;
    const double theta_max = M_PI;

    // Create the C-space grid for the joint angles
    auto cspace_ptr = std::make_unique<MyGridCSpace2D>(
        m_cells_per_dim, m_cells_per_dim,
        theta_min, theta_max,
        theta_min, theta_max);
    MyGridCSpace2D& cspace = *cspace_ptr;

    std::cout << "Constructing C-space for 2-link manipulator...\n";

    // Iterate through each cell of the C-space grid
    for (std::size_t ix = 0; ix < m_cells_per_dim; ++ix) {
        for (std::size_t iy = 0; iy < m_cells_per_dim; ++iy) {
            // Get the configuration (q1, q2) for the center of the current cell
            Eigen::Vector2d q_center_vec = cspace.cellCenter(ix, iy);
            amp::ManipulatorState q_center(manipulator.nLinks());
            q_center << q_center_vec.x(), q_center_vec.y();

            // Use forward kinematics to find the positions of the joints in the workspace
            Eigen::Vector2d base_pos  = manipulator.getJointLocation(q_center, 0); // Base of link 1
            Eigen::Vector2d elbow_pos = manipulator.getJointLocation(q_center, 1); // Elbow (end of link 1, base of link 2)
            Eigen::Vector2d ee_pos    = manipulator.getJointLocation(q_center, 2); // End-effector (end of link 2)

            // Check if either link is in collision with any obstacle
            bool link1_collides = CollisionChecker::isSegmentInCollision(base_pos, elbow_pos, env.obstacles);
            bool link2_collides = CollisionChecker::isSegmentInCollision(elbow_pos, ee_pos, env.obstacles);

            // If any part of the manipulator collides, mark the cell as occupied (true)
            cspace(ix, iy) = link1_collides || link2_collides;
        }
    }
    return cspace_ptr;
}


// =======================
// Concrete WaveFront implementation (no base call)
//  - 4-connected BFS from goal
//  - reconstruct to start
//  - interior points = cell centers, endpoints = exact q_init/q_goal
// =======================
amp::Path2D MyWaveFrontAlgorithm::planInCSpace(const Eigen::Vector2d& q_init,
                                               const Eigen::Vector2d& q_goal,
                                               const amp::GridCSpace2D& grid_cspace,
                                               bool /*isManipulator*/) {
    const auto* mygrid = dynamic_cast<const MyGridCSpace2D*>(&grid_cspace);
    if (!mygrid) return amp::Path2D{};

    const std::size_t NX = mygrid->width();
    const std::size_t NY = mygrid->height();

    // Map start/goal to cells
    auto [sx, sy] = mygrid->getCellFromPoint(q_init.x(), q_init.y());
    auto [gx, gy] = mygrid->getCellFromPoint(q_goal.x(), q_goal.y());

    auto is_blocked = [&](std::size_t ix, std::size_t iy) -> bool {
        return (*mygrid)(ix, iy); // true means occupied
    };
    if (is_blocked(sx, sy) || is_blocked(gx, gy)) {
        return amp::Path2D{}; // start or goal in collision
    }

    // 4-neighbor BFS from goal to compute distances
    const int INF = std::numeric_limits<int>::max();
    std::vector<int> dist(NX * NY, INF);
    std::vector<std::pair<int,int>> parent(NX * NY, {-1, -1});

    auto idx = [&](std::size_t ix, std::size_t iy) -> std::size_t {
        return iy * NX + ix;
    };

    std::queue<std::pair<std::size_t,std::size_t>> q;
    dist[idx(gx, gy)] = 0;
    q.push({gx, gy});

    const int DX[4] = {+1, -1,  0,  0};
    const int DY[4] = { 0,  0, +1, -1};

    auto in_bounds = [&](int ix, int iy) -> bool {
        return (ix >= 0 && iy >= 0 && ix < static_cast<int>(NX) && iy < static_cast<int>(NY));
    };

    while (!q.empty()) {
        auto [cx, cy] = q.front(); q.pop();
        const int cd = dist[idx(cx, cy)];
        for (int k = 0; k < 4; ++k) {
            int nx_ = static_cast<int>(cx) + DX[k];
            int ny_ = static_cast<int>(cy) + DY[k];
            if (!in_bounds(nx_, ny_)) continue;
            const std::size_t nix = static_cast<std::size_t>(nx_);
            const std::size_t niy = static_cast<std::size_t>(ny_);
            if (is_blocked(nix, niy)) continue;
            std::size_t nidx = idx(nix, niy);
            if (dist[nidx] == INF) {
                dist[nidx] = cd + 1;
                parent[nidx] = {static_cast<int>(cx), static_cast<int>(cy)};
                q.push({nix, niy});
            }
        }
    }

    // If start unreachable, return empty path
    if (dist[idx(sx, sy)] == INF) {
        return amp::Path2D{};
    }

    // Reconstruct from start -> goal
    std::vector<std::pair<std::size_t,std::size_t>> cells;
    {
        std::size_t cx = sx, cy = sy;
        cells.push_back({cx, cy});
        while (!(cx == gx && cy == gy)) {
            auto pr = parent[idx(cx, cy)];
            if (pr.first < 0) break; // safety
            cx = static_cast<std::size_t>(pr.first);
            cy = static_cast<std::size_t>(pr.second);
            cells.push_back({cx, cy});
        }
    }

    // Convert to world; interior = centers; endpoints = exact q_init/q_goal
    std::vector<Eigen::Vector2d> pts;
    pts.reserve(cells.size());
    for (auto [ix, iy] : cells) {
        pts.emplace_back(mygrid->cellCenter(ix, iy));
    }
    if (!pts.empty()) {
        pts.front() = q_init;
        pts.back()  = q_goal;
    }

    amp::Path2D path;
    path.waypoints = std::move(pts);
    return path;
}