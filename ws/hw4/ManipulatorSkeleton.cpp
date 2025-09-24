#include "ManipulatorSkeleton.h"
#include <cmath>

// Default constructor
MyManipulator2D::MyManipulator2D() 
    : LinkManipulator2D({1.0, 1.0}) 
{}

// Constructor for a vector of link lengths
MyManipulator2D::MyManipulator2D(const std::vector<double>& link_lengths) 
    : LinkManipulator2D(link_lengths) 
{}

// Forward Kinematics for an N-Link Manipulator
Eigen::Vector2d MyManipulator2D::getJointLocation(const amp::ManipulatorState& state, uint32_t joint_index) const {
    const amp::ManipulatorState* state_ptr = &state;
    amp::ManipulatorState temp_state;
    if (state.size() != nLinks()) {
        temp_state.resize(nLinks());
        temp_state.setZero();
        state_ptr = &temp_state;
    }
    
    Eigen::Vector2d joint_location(0.0, 0.0);
    double cumulative_angle = 0.0;
    
    // Sum the vectors of the links up to the desired joint
    for (uint32_t i = 0; i < joint_index; ++i) {
        cumulative_angle += (*state_ptr)(i);
        joint_location.x() += m_link_lengths[i] * cos(cumulative_angle);
        joint_location.y() += m_link_lengths[i] * sin(cumulative_angle);
    }

    return joint_location;
}

// Inverse Kinematics for a 3-Link Manipulator
amp::ManipulatorState MyManipulator2D::getConfigurationFromIK(const Eigen::Vector2d& end_effector_location) const {
    // If you have different implementations for 2/3/n link manipulators, you can separate them here
    if (nLinks() != 3) {
        amp::ManipulatorState state(nLinks());
        state.setZero(); 
        return state;
    }

    const double x_e = end_effector_location.x();
    const double y_e = end_effector_location.y();
    const double a1 = m_link_lengths[0];
    const double a2 = m_link_lengths[1];
    const double a3 = m_link_lengths[2];

    amp::ManipulatorState state(3);

    // Assume absolute orientation of the final link is 0 radians
    const double phi_e = 0.0;

    // Link between 2 and 3 position
    const double x_w = x_e - a3 * cos(phi_e);
    const double y_w = y_e - a3 * sin(phi_e);

    // IK to move wrist
    const double cos_theta2 = (x_w*x_w + y_w*y_w - a1*a1 - a2*a2) / (2.0 * a1 * a2);

    // Check if target is reachable
    if (cos_theta2 < -1.0 || cos_theta2 > 1.0) {
        state.setZero(); // If unreachable
        return state;
    }
    
    // Alternate configuration
    const double theta2 = -acos(cos_theta2);
    const double theta1 = atan2(y_w, x_w) - atan2(a2 * sin(theta2), a1 + a2 * cos(theta2));

    // Calculate last angle theta3
    const double theta3 = phi_e - theta1 - theta2;

    state << theta1, theta2, theta3;
    return state;
}