#include "CSpaceSkeleton.h"
#include "HelpfulClass.h" 

// Override this method for returning whether or not a point is in collision
std::pair<std::size_t, std::size_t> MyGridCSpace2D::getCellFromPoint(double x0, double x1) const {
    // Implement your discretization procedure here, such that the point (x0, x1) lies within the returned cell
    
    // FINAL FIX: Use our class's own member variables that we stored in the constructor.
    const double x0_range = m_x0_max - m_x0_min;
    const double x1_range = m_x1_max - m_x1_min;

    const double cell_width_x0 = x0_range / m_x0_cells;
    const double cell_width_x1 = x1_range / m_x1_cells;

    // Calculate the cell index for each dimension
    std::size_t cell_x = std::max(0, std::min((int)m_x0_cells - 1, (int)floor((x0 - m_x0_min) / cell_width_x0)));
    std::size_t cell_y = std::max(0, std::min((int)m_x1_cells - 1, (int)floor((x1 - m_x1_min) / cell_width_x1)));

    return {cell_x, cell_y};
}

// Override this method for computing all of the boolean collision values for each cell in the cspace
std::unique_ptr<amp::GridCSpace2D> MyManipulatorCSConstructor::construct(const amp::LinkManipulator2D& manipulator, const amp::Environment2D& env) {
    // Create an object of my custom cspace type and store it in a unique pointer.
    std::unique_ptr<MyGridCSpace2D> cspace_ptr = std::make_unique<MyGridCSpace2D>(m_cells_per_dim, m_cells_per_dim, 0.0, 2.0 * M_PI, 0.0, 2.0 * M_PI);
    MyGridCSpace2D& cspace = *cspace_ptr;

    // Get grid properties to help calculate angles for cell centers
    // FINAL FIX: Use the new public getter methods we defined in our MyGridCSpace2D class.
    const double cell_width_theta1 = (cspace.getX0_max() - cspace.getX0_min()) / cspace.getX0_cells();
    const double cell_width_theta2 = (cspace.getX1_max() - cspace.getX1_min()) / cspace.getX1_cells();

    // Iterate through every cell in the C-space grid
    for (std::size_t i = 0; i < cspace.getX0_cells(); ++i) {
        for (std::size_t j = 0; j < cspace.getX1_cells(); ++j) {
            // Find the angles corresponding to the center of the current cell
            double theta1 = cspace.getX0_min() + (i + 0.5) * cell_width_theta1;
            double theta2 = cspace.getX1_min() + (j + 0.5) * cell_width_theta2;
            
            amp::ManipulatorState state(2);
            state << theta1, theta2;

            bool in_collision = false;
            // Now, check each link of the manipulator for collision at this configuration
            for (uint32_t k = 0; k < manipulator.nLinks(); ++k) {
                // Get the start and end points of the current link
                Eigen::Vector2d p1 = manipulator.getJointLocation(state, k);
                Eigen::Vector2d p2 = manipulator.getJointLocation(state, k + 1);

                if (MotionPlanningHelpers::CollisionChecker::isSegmentInCollision(p1, p2, env.obstacles)) {
                    in_collision = true;
                    break; 
                }
            }
            cspace(i, j) = in_collision;
        }
    }
    return cspace_ptr;
}