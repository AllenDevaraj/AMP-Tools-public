#include "ManipulatorSkeleton.h"
#include <cmath>

MyManipulator2D::MyManipulator2D() 
    : LinkManipulator2D({1.0, 1.0}) 
{}

MyManipulator2D::MyManipulator2D(const std::vector<double>& link_lengths) 
    : LinkManipulator2D(link_lengths) 
{}

// Forward Kinematics
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
    
    for (uint32_t i = 0; i < joint_index; ++i) {
        cumulative_angle += (*state_ptr)(i);
        joint_location.x() += m_link_lengths[i] * cos(cumulative_angle);
        joint_location.y() += m_link_lengths[i] * sin(cumulative_angle);
    }

    return joint_location;
}

// Inverse Kinematics
amp::ManipulatorState MyManipulator2D::getConfigurationFromIK(const Eigen::Vector2d& end_effector_location) const {
    
    // 2 Link
    if (nLinks() == 2) {
        amp::ManipulatorState state(2);
        const double x = end_effector_location.x();
        const double y = end_effector_location.y();
        const double l1 = m_link_lengths[0];
        const double l2 = m_link_lengths[1];

        const double distance_squared = x*x + y*y;
        if (distance_squared > (l1 + l2) * (l1 + l2) || distance_squared < (l1 - l2) * (l1 - l2)) {
            state.setZero();
            return state;
        }

        const double cos_theta2 = (x*x + y*y - l1*l1 - l2*l2) / (2.0 * l1 * l2);
        const double theta2 = -acos(cos_theta2); // Elbow up
        const double theta1 = atan2(y, x) - atan2(l2 * sin(theta2), l1 + l2 * cos(theta2));
        state << theta1, theta2;
        return state;
    } 
    // 3 Link
    else if (nLinks() == 3) {
        amp::ManipulatorState state(3);
        const double x = end_effector_location.x();
        const double y = end_effector_location.y();
        const double a1 = m_link_lengths[0];
        const double a2 = m_link_lengths[1];
        const double a3 = m_link_lengths[2];

        // Assumption: theta3 = 0
        const double theta3 = 0.0;
        const double l1_virtual = a1;
        const double l2_virtual = a2 + a3;

        const double distance_squared = x*x + y*y;
        if (distance_squared > (l1_virtual + l2_virtual) * (l1_virtual + l2_virtual)) {
            state.setZero(); // Target is unreachable
            return state;
        }

        const double cos_theta2 = (x*x + y*y - l1_virtual*l1_virtual - l2_virtual*l2_virtual) / (2.0 * l1_virtual * l2_virtual);
        if (cos_theta2 < -1.0 || cos_theta2 > 1.0) {
             state.setZero();
             return state;
        }

        const double theta2 = -acos(cos_theta2);
        const double theta1 = atan2(y, x) - atan2(l2_virtual * sin(theta2), l1_virtual + l2_virtual * cos(theta2));
        
        state << theta1, theta2, theta3;
        return state;
    }
    // n Link
    else {
        amp::ManipulatorState current_state(nLinks());
        current_state.setZero();
        const int max_iterations = 100;
        const double tolerance = 1e-3;
        for (int iter = 0; iter < max_iterations; ++iter) {
            if ((getJointLocation(current_state, nLinks()) - end_effector_location).norm() < tolerance) break;
            for (int i = nLinks() - 1; i >= 0; --i) {
                Eigen::Vector2d joint_pos = getJointLocation(current_state, i);
                Eigen::Vector2d end_effector_pos = getJointLocation(current_state, nLinks());
                Eigen::Vector2d vec_to_target = end_effector_location - joint_pos;
                Eigen::Vector2d vec_to_effector = end_effector_pos - joint_pos;
                double angle_to_target = atan2(vec_to_target.y(), vec_to_target.x());
                double angle_to_effector = atan2(vec_to_effector.y(), vec_to_effector.x());
                current_state(i) += angle_to_target - angle_to_effector;
            }
        }
        return current_state;
    }
}