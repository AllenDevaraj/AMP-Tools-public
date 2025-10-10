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

std::unique_ptr<amp::GridCSpace2D>
MyPointAgentCSConstructor::construct(const amp::Environment2D& env) {
    const double cell = 0.25;
    const std::size_t nx = static_cast<std::size_t>(std::ceil((env.x_max - env.x_min) / cell));
    const std::size_t ny = static_cast<std::size_t>(std::ceil((env.y_max - env.y_min) / cell));

    auto cspace_ptr = std::make_unique<MyGridCSpace2D>(nx, ny, env.x_min, env.x_max, env.y_min, env.y_max);
    MyGridCSpace2D& cspace = *cspace_ptr;

    std::cout << "Constructing C-space for point agent (cell=0.25, 4-connected)\n";

    // Center-point collision
    for (std::size_t ix = 0; ix < nx; ++ix) {
        for (std::size_t iy = 0; iy < ny; ++iy) {
            const Eigen::Vector2d center = cspace.cellCenter(ix, iy);
            const Vec p(center.x(), center.y());
            const bool occupied = CollisionChecker::isPointInCollision(p, env.obstacles);
            cspace(ix, iy) = occupied; 
        }
    }
    return cspace_ptr;
}

std::unique_ptr<amp::GridCSpace2D>
MyManipulatorCSConstructor::construct(const amp::LinkManipulator2D& manipulator,
                                      const amp::Environment2D& env) {
    const double theta_min = -M_PI;
    const double theta_max =  M_PI;

    auto cspace_ptr = std::make_unique<MyGridCSpace2D>(
        m_cells_per_dim, m_cells_per_dim,
        theta_min, theta_max,
        theta_min, theta_max);
    MyGridCSpace2D& cspace = *cspace_ptr;

    std::cout << "Constructing C-space for 2-link manipulator (wrap-around angles)\n";

    for (std::size_t ix = 0; ix < m_cells_per_dim; ++ix) {
        for (std::size_t iy = 0; iy < m_cells_per_dim; ++iy) {
            // Joint-angle center (q1, q2)
            const Eigen::Vector2d q_center = cspace.cellCenter(ix, iy);

            amp::ManipulatorState q(manipulator.nLinks());
            q << q_center.x(), q_center.y();

            // FK
            const Eigen::Vector2d base     = manipulator.getJointLocation(q, 0);
            const Eigen::Vector2d joint    = manipulator.getJointLocation(q, 1);
            const Eigen::Vector2d endeff   = manipulator.getJointLocation(q, 2);

            const bool link1_collides = CollisionChecker::isSegmentInCollision(base,  joint,  env.obstacles);
            const bool link2_collides = CollisionChecker::isSegmentInCollision(joint, endeff, env.obstacles);

            cspace(ix, iy) = link1_collides || link2_collides;
        }
    }
    return cspace_ptr;
}

amp::Path2D MyWaveFrontAlgorithm::planInCSpace(const Eigen::Vector2d& q_init,
                                               const Eigen::Vector2d& q_goal,
                                               const amp::GridCSpace2D& grid_cspace,
                                               bool isManipulator) {
    const auto* mygrid = dynamic_cast<const MyGridCSpace2D*>(&grid_cspace);
    if (!mygrid) return amp::Path2D{};

    const std::size_t NX = mygrid->width();
    const std::size_t NY = mygrid->height();

    auto [sx, sy] = mygrid->getCellFromPoint(q_init.x(), q_init.y());
    auto [gx, gy] = mygrid->getCellFromPoint(q_goal.x(), q_goal.y());

    auto is_blocked = [&](std::size_t ix, std::size_t iy) -> bool {
        return (*mygrid)(ix, iy); 
    };
    if (is_blocked(sx, sy) || is_blocked(gx, gy)) {
        return amp::Path2D{}; 
    }

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

            if (isManipulator) {
                if (nx_ < 0) nx_ += static_cast<int>(NX);
                if (ny_ < 0) ny_ += static_cast<int>(NY);
                if (nx_ >= static_cast<int>(NX)) nx_ -= static_cast<int>(NX);
                if (ny_ >= static_cast<int>(NY)) ny_ -= static_cast<int>(NY);
            } else {
                if (!in_bounds(nx_, ny_)) continue;
            }

            const std::size_t nix = static_cast<std::size_t>(nx_);
            const std::size_t niy = static_cast<std::size_t>(ny_);
            if (!isManipulator && !in_bounds(nx_, ny_)) continue;
            if (is_blocked(nix, niy)) continue;

            std::size_t nidx = idx(nix, niy);
            if (dist[nidx] == INF) {
                dist[nidx] = cd + 1;
                parent[nidx] = {static_cast<int>(cx), static_cast<int>(cy)};
                q.push({nix, niy});
            }
        }
    }


    std::vector<std::pair<std::size_t,std::size_t>> cells;
    {
        std::size_t cx = sx, cy = sy;
        cells.push_back({cx, cy});
        while (!(cx == gx && cy == gy)) {
            auto pr = parent[idx(cx, cy)];
            if (pr.first < 0) break;
            cx = static_cast<std::size_t>(pr.first);
            cy = static_cast<std::size_t>(pr.second);
            cells.push_back({cx, cy});
        }
    }

    std::vector<Eigen::Vector2d> pts;
    pts.reserve(cells.size());
    for (auto [ix, iy] : cells) {
        pts.emplace_back(mygrid->cellCenter(ix, iy));
    }

    // Endpoints:
    if (!pts.empty()) {
        pts.front() = q_init;
        pts.back()  = q_goal;
    }

    if (isManipulator) {
        const double TWO_PI = 2.0 * M_PI;
        for (std::size_t i = 1; i < pts.size(); ++i) {
            for (int j = 0; j < 2; ++j) {
                double d = pts[i][j] - pts[i-1][j];
                if (d >  M_PI) pts[i][j] -= TWO_PI;
                if (d < -M_PI) pts[i][j] += TWO_PI;
            }
        }
    }

    amp::Path2D path;
    path.waypoints = std::move(pts);
    return path;
}
