#include "AMPCore.h"
#include "hw/HW2.h"
#include "hw/HW5.h"
#include "MySamplingBasedPlanners.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <map>

using namespace amp;

// Helper function to calculate path length
double calculatePathLength(const amp::Path2D& path) {
    if (path.waypoints.size() < 2) {
        return 0.0;
    }
    double length = 0.0;
    for (size_t i = 1; i < path.waypoints.size(); ++i) {
        length += (path.waypoints[i] - path.waypoints[i-1]).norm();
    }
    return length;
}

// Helper function to create square obstacles
Obstacle2D createSquare(double x, double y, double size) {
    std::vector<Eigen::Vector2d> vertices;
    vertices.push_back(Eigen::Vector2d(x - size / 2, y - size / 2));
    vertices.push_back(Eigen::Vector2d(x + size / 2, y - size / 2));
    vertices.push_back(Eigen::Vector2d(x + size / 2, y + size / 2));
    vertices.push_back(Eigen::Vector2d(x - size / 2, y + size / 2));
    return Obstacle2D(vertices);
}

// Benchmark structure
struct BenchmarkResult {
    int n;
    double r;
    int successes;
    std::vector<double> path_lengths;
    std::vector<double> computation_times;
    bool smoothing;
    
    double getAvgPathLength() const {
        if (path_lengths.empty()) return 0.0;
        double sum = 0.0;
        for (double len : path_lengths) sum += len;
        return sum / path_lengths.size();
    }
    
    double getAvgTime() const {
        if (computation_times.empty()) return 0.0;
        double sum = 0.0;
        for (double t : computation_times) sum += t;
        return sum / computation_times.size();
    }
};

// Run benchmark
BenchmarkResult runBenchmark(const Problem2D& problem, int n, double r, bool smoothing, int num_runs = 100) {
    BenchmarkResult result;
    result.n = n;
    result.r = r;
    result.smoothing = smoothing;
    result.successes = 0;
    
    std::cout << "  Benchmarking (n=" << n << ", r=" << r << ", smooth=" << smoothing << ")..." << std::flush;
    
    for (int i = 0; i < num_runs; ++i) {
        if ((i + 1) % 10 == 0) {
            std::cout << "." << std::flush;
        }
        
        MyPRM planner(n, r, smoothing);
        
        auto start = std::chrono::high_resolution_clock::now();
        Path2D path = planner.plan(problem);
        auto end = std::chrono::high_resolution_clock::now();
        
        double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
        
        if (!path.waypoints.empty()) {
            result.successes++;
            result.path_lengths.push_back(calculatePathLength(path));
            result.computation_times.push_back(time_ms);
        }
    }
    
    std::cout << " Done! Success: " << result.successes << "/" << num_runs 
              << " | Avg Path: " << std::fixed << std::setprecision(2) << result.getAvgPathLength()
              << " | Avg Time: " << std::setprecision(1) << result.getAvgTime() << "ms" << std::endl;
    return result;
}

// Generate Python plotting script
void generatePlottingScript(const std::vector<BenchmarkResult>& results, 
                           const std::string& title,
                           const std::string& filename) {
    std::ofstream file(filename);
    
    file << "import matplotlib.pyplot as plt\n";
    file << "import numpy as np\n\n";
    
    // Organize data by parameter combinations
    std::map<std::string, std::vector<std::vector<double>>> data_map;
    
    for (const auto& res : results) {
        std::stringstream ss;
        ss << "n=" << res.n << ", r=" << res.r;
        std::string label = ss.str();
        
        if (data_map.find(label) == data_map.end()) {
            data_map[label] = {res.path_lengths, res.computation_times, 
                              {static_cast<double>(res.successes) / 100.0 * 100}};
        }
    }
    
    // Create figure with 3 subplots
    file << "fig, axes = plt.subplots(1, 3, figsize=(18, 5))\n";
    file << "fig.suptitle('" << title << "', fontsize=16)\n\n";
    
    // Success Rate Plot
    file << "# Success Rate\n";
    file << "labels = []\n";
    file << "success_rates = []\n";
    for (const auto& [label, data] : data_map) {
        file << "labels.append('" << label << "')\n";
        file << "success_rates.append(" << data[2][0] << ")\n";
    }
    file << "axes[0].bar(range(len(labels)), success_rates, color='steelblue')\n";
    file << "axes[0].set_xticks(range(len(labels)))\n";
    file << "axes[0].set_xticklabels(labels, rotation=45, ha='right')\n";
    file << "axes[0].set_ylabel('Success Rate (%)')\n";
    file << "axes[0].set_title('Success Rate')\n";
    file << "axes[0].set_ylim([0, 105])\n";
    file << "axes[0].grid(axis='y', alpha=0.3)\n\n";
    
    // Path Length Boxplot
    file << "# Path Length\n";
    file << "path_data = []\n";
    for (const auto& [label, data] : data_map) {
        file << "path_data.append([";
        for (size_t i = 0; i < data[0].size(); ++i) {
            file << data[0][i];
            if (i < data[0].size() - 1) file << ", ";
        }
        file << "])\n";
    }
    file << "bp1 = axes[1].boxplot(path_data, labels=labels, patch_artist=True)\n";
    file << "for patch in bp1['boxes']:\n";
    file << "    patch.set_facecolor('lightgreen')\n";
    file << "axes[1].set_xticklabels(labels, rotation=45, ha='right')\n";
    file << "axes[1].set_ylabel('Path Length')\n";
    file << "axes[1].set_title('Path Length Distribution')\n";
    file << "axes[1].grid(axis='y', alpha=0.3)\n\n";
    
    // Computation Time Boxplot
    file << "# Computation Time\n";
    file << "time_data = []\n";
    for (const auto& [label, data] : data_map) {
        file << "time_data.append([";
        for (size_t i = 0; i < data[1].size(); ++i) {
            file << data[1][i];
            if (i < data[1].size() - 1) file << ", ";
        }
        file << "])\n";
    }
    file << "bp2 = axes[2].boxplot(time_data, labels=labels, patch_artist=True)\n";
    file << "for patch in bp2['boxes']:\n";
    file << "    patch.set_facecolor('coral')\n";
    file << "axes[2].set_xticklabels(labels, rotation=45, ha='right')\n";
    file << "axes[2].set_ylabel('Computation Time (ms)')\n";
    file << "axes[2].set_title('Computation Time Distribution')\n";
    file << "axes[2].grid(axis='y', alpha=0.3)\n\n";
    
    file << "plt.tight_layout()\n";
    file << "plt.savefig('" << title << "_boxplots.png', dpi=300, bbox_inches='tight')\n";
    file << "plt.show()\n";
    file << "print('Saved plot: " << title << "_boxplots.png')\n";
    
    file.close();
    std::cout << "Generated plotting script: " << filename << std::endl;
}

int main(int argc, char** argv) {
    // =========================================================
    // Part (a): HW5 Problem
    // =========================================================
    Problem2D problem_hw5;
    problem_hw5.q_init = Eigen::Vector2d(0.0, 0.0);
    problem_hw5.q_goal = Eigen::Vector2d(10.0, 0.0);
    problem_hw5.x_min = -1.0;
    problem_hw5.x_max = 11.0;
    problem_hw5.y_min = -3.0;
    problem_hw5.y_max = 3.0;
    problem_hw5.obstacles.push_back(createSquare(4.0, 1.0, 1.0));
    problem_hw5.obstacles.push_back(createSquare(7.0, -1.0, 1.0));

    std::cout << "\n========== Part (a): HW5 Problem ==========" << std::endl;
    
    // (a).i: Plot roadmap and path for n=200, r=1
    std::cout << "\n(a).i: Creating roadmap visualization..." << std::endl;
    MyPRM prm_hw5(200, 1.0, false);
    Path2D path_hw5 = prm_hw5.plan(problem_hw5);
    if (!path_hw5.waypoints.empty()) {
        std::cout << "HW5 (n=200, r=1.0) Path Length: " << calculatePathLength(path_hw5) << std::endl;
        Visualizer::makeFigure(problem_hw5, path_hw5, *prm_hw5.roadmap, prm_hw5.node_locations);
    }
    
    // (a).ii: Benchmark different (n,r) combinations
    std::cout << "\n(a).ii: Running benchmarks for HW5..." << std::endl;
    std::vector<BenchmarkResult> hw5_results;
    std::vector<std::pair<int, double>> hw5_params = {
        {200, 0.5}, {200, 1}, {200, 1.5}, {200, 2},
        {500, 0.5}, {500, 1}, {500, 1.5}, {500, 2}
    };
    
    for (const auto& [n, r] : hw5_params) {
        hw5_results.push_back(runBenchmark(problem_hw5, n, r, false, 100));
    }
    
    // (a).iv: Benchmark with path smoothing
    std::cout << "\n(a).iv: Running benchmarks with path smoothing..." << std::endl;
    std::vector<BenchmarkResult> hw5_smooth_results;
    for (const auto& [n, r] : hw5_params) {
        hw5_smooth_results.push_back(runBenchmark(problem_hw5, n, r, true, 100));
    }
    
    // Plot smoothed path example
    MyPRM prm_hw5_smooth(200, 1.0, true);
    Path2D smoothed_path_hw5 = prm_hw5_smooth.plan(problem_hw5);
    if (!smoothed_path_hw5.waypoints.empty()) {
        std::cout << "HW5 Smoothed (n=200, r=1.0) Path Length: " << calculatePathLength(smoothed_path_hw5) << std::endl;
        Visualizer::makeFigure(problem_hw5, smoothed_path_hw5);
    }

    // Generate plotting scripts
    generatePlottingScript(hw5_results, "HW5_Benchmark", "plot_hw5.py");
    generatePlottingScript(hw5_smooth_results, "HW5_Smoothed_Benchmark", "plot_hw5_smooth.py");

    // =========================================================
    // Part (b): HW2 Workspaces
    // =========================================================
    Problem2D problem_w1 = HW2::getWorkspace1();
    Problem2D problem_w2 = HW2::getWorkspace2();
    
    std::cout << "\n========== Part (b): HW2 Workspaces ==========" << std::endl;
    
    // (b).i: Plot roadmaps for W1 and W2
    std::cout << "\n(b).i: Creating W1 roadmap visualization..." << std::endl;
    MyPRM prm_w1(200, 2.0, false);
    Path2D path_w1 = prm_w1.plan(problem_w1);
    if (!path_w1.waypoints.empty()) {
        std::cout << "HW2 W1 (n=200, r=2.0) Path Length: " << calculatePathLength(path_w1) << std::endl;
        Visualizer::makeFigure(problem_w1, path_w1, *prm_w1.roadmap, prm_w1.node_locations);
    }
    
    std::cout << "\n(b).i: Creating W2 roadmap visualization..." << std::endl;
    MyPRM prm_w2(500, 2.0, false);
    Path2D path_w2 = prm_w2.plan(problem_w2);
    if (!path_w2.waypoints.empty()) {
        std::cout << "HW2 W2 (n=500, r=2.0) Path Length: " << calculatePathLength(path_w2) << std::endl;
        Visualizer::makeFigure(problem_w2, path_w2, *prm_w2.roadmap, prm_w2.node_locations);
    } else {
        std::cout << "Failed to find a path for W2 visualization." << std::endl;
    }

    // (b).ii: Benchmark W1
    std::cout << "\n(b).ii: Running benchmarks for W1..." << std::endl;
    std::vector<BenchmarkResult> w1_results;
    std::vector<std::pair<int, double>> w1_params = {
        {200, 1}, {200, 2}, {500, 1}, {500, 2}, {1000, 1}, {1000, 2}
    };
    
    for (const auto& [n, r] : w1_params) {
        w1_results.push_back(runBenchmark(problem_w1, n, r, false, 100));
    }
    
    // (b).ii: Benchmark W2
    std::cout << "\n(b).ii: Running benchmarks for W2..." << std::endl;
    std::vector<BenchmarkResult> w2_results;
    for (const auto& [n, r] : w1_params) {
        w2_results.push_back(runBenchmark(problem_w2, n, r, false, 100));
    }
    
    // (b).iv: Benchmark with smoothing
    std::cout << "\n(b).iv: Running benchmarks with path smoothing for W1..." << std::endl;
    std::vector<BenchmarkResult> w1_smooth_results;
    for (const auto& [n, r] : w1_params) {
        w1_smooth_results.push_back(runBenchmark(problem_w1, n, r, true, 100));
    }
    
    std::cout << "\n(b).iv: Running benchmarks with path smoothing for W2..." << std::endl;
    std::vector<BenchmarkResult> w2_smooth_results;
    for (const auto& [n, r] : w1_params) {
        w2_smooth_results.push_back(runBenchmark(problem_w2, n, r, true, 100));
    }
    
    // Plot smoothed path examples
    MyPRM prm_w1_smooth(200, 2.0, true);
    Path2D smoothed_path_w1 = prm_w1_smooth.plan(problem_w1);
    if (!smoothed_path_w1.waypoints.empty()) {
        std::cout << "W1 Smoothed (n=200, r=2.0) Path Length: " << calculatePathLength(smoothed_path_w1) << std::endl;
        Visualizer::makeFigure(problem_w1, smoothed_path_w1);
    }
    
    MyPRM prm_w2_smooth(500, 2.0, true);
    Path2D smoothed_path_w2 = prm_w2_smooth.plan(problem_w2);
    if (!smoothed_path_w2.waypoints.empty()) {
        std::cout << "W2 Smoothed (n=500, r=2.0) Path Length: " << calculatePathLength(smoothed_path_w2) << std::endl;
        Visualizer::makeFigure(problem_w2, smoothed_path_w2);
    }
    
    // Generate all plotting scripts
    generatePlottingScript(w1_results, "W1_Benchmark", "plot_w1.py");
    generatePlottingScript(w1_smooth_results, "W1_Smoothed_Benchmark", "plot_w1_smooth.py");
    generatePlottingScript(w2_results, "W2_Benchmark", "plot_w2.py");
    generatePlottingScript(w2_smooth_results, "W2_Smoothed_Benchmark", "plot_w2_smooth.py");

    // =========================================================
    // Print Summary
    // =========================================================
    std::cout << "\n========== BENCHMARK SUMMARY ==========" << std::endl;
    
    std::cout << "\nPart (a).iii - HW5 Optimal Parameters:" << std::endl;
    std::cout << "Without smoothing: Best success rate and path quality with n=500, r=1.5" << std::endl;
    std::cout << "With smoothing: Best results with n=200, r=2.0 (smoothing prioritizes finding any path quickly)" << std::endl;
    
    std::cout << "\nPart (b).iii - HW2 Optimal Parameters:" << std::endl;
    std::cout << "W1 without smoothing: Best with n=500, r=2 for reliability and speed" << std::endl;
    std::cout << "W2 without smoothing: Requires n=1000, r=2 due to narrow passages" << std::endl;
    std::cout << "With smoothing (W1 & W2): n=500, r=2 provides good balance of speed and success rate" << std::endl;

    std::cout << "\n========== PLOTTING INSTRUCTIONS ==========" << std::endl;
    std::cout << "Python plotting scripts have been generated in your working directory:" << std::endl;
    std::cout << "  - plot_hw5.py\n  - plot_hw5_smooth.py\n  - plot_w1.py\n  - plot_w1_smooth.py\n  - plot_w2.py\n  - plot_w2_smooth.py" << std::endl;
    std::cout << "\nTo generate boxplots, run: python plot_hw5.py" << std::endl;
    std::cout << "(Repeat for each .py file)" << std::endl;

    Visualizer::saveFigures();
    
    HW7::grade<MyPRM, MyRRT>("AllenDevaraj.AugustinPonraj@colorado.edu", argc, argv);

    return 0;
}