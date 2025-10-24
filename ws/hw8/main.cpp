#include "AMPCore.h"
#include "hw/HW8.h"
#include "MyMultiAgentPlanners.h"

#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <map>

using namespace amp;

// Benchmark structure
struct BenchmarkResult {
    int num_agents;
    int successes;
    std::vector<double> tree_sizes;
    std::vector<double> computation_times;
    
    double getAvgTreeSize() const {
        if (tree_sizes.empty()) return 0.0;
        double sum = 0.0;
        for (double size : tree_sizes) sum += size;
        return sum / tree_sizes.size();
    }
    
    double getAvgTime() const {
        if (computation_times.empty()) return 0.0;
        double sum = 0.0;
        for (double t : computation_times) sum += t;
        return sum / computation_times.size();
    }
};

BenchmarkResult runBenchmark(int num_agents, int num_runs = 100) {
    BenchmarkResult result;
    result.num_agents = num_agents;
    result.successes = 0;
    
    std::cout << "  Benchmarking (m=" << num_agents << ")..." << std::flush;
    
    for (int i = 0; i < num_runs; ++i) {
        if ((i + 1) % 10 == 0) {
            std::cout << "." << std::flush;
        }
        
        MyCentralPlanner planner;
        MultiAgentProblem2D problem = HW8::getWorkspace1(num_agents);
        
        MultiAgentPath2D path = planner.plan(problem);
        
        // Check for success (at least one agent has a non-empty path)
        if (!path.agent_paths.empty() && !path.agent_paths[0].waypoints.empty()) {
            result.successes++;
            result.tree_sizes.push_back(static_cast<double>(planner.last_tree_size)); 
            result.computation_times.push_back(planner.last_comp_time_ms);
        }
    }
    
    std::cout << " Done! Success: " << result.successes << "/" << num_runs 
              << " | Avg Tree: " << std::fixed << std::setprecision(0) << result.getAvgTreeSize()
              << " | Avg Time: " << std::setprecision(1) << result.getAvgTime() << "ms" << std::endl;
    return result;
}

// Python plotting script generator
void generatePlottingScript(const std::vector<BenchmarkResult>& results) {
    std::string filename = "plot_hw8_benchmark.py";
    std::ofstream file(filename);
    
    file << "import matplotlib.pyplot as plt\n";
    file << "import numpy as np\n\n";
    
    // Data for plots
    std::vector<std::string> labels;
    std::vector<double> avg_times;
    std::vector<double> avg_tree_sizes;
    std::vector<double> success_rates;
    std::vector<int> agent_counts;

    for (const auto& res : results) {
        labels.push_back(std::to_string(res.num_agents));
        agent_counts.push_back(res.num_agents);
        avg_times.push_back(res.getAvgTime());
        avg_tree_sizes.push_back(res.getAvgTreeSize());
        success_rates.push_back(static_cast<double>(res.successes));
    }

    // --- Figure 1: Boxplots ---
    file << "fig1, axes1 = plt.subplots(1, 3, figsize=(18, 6))\n";
    file << "fig1.suptitle('HW8 Centralized RRT Benchmark (100 Runs)', fontsize=16)\n\n";
    
    // Success Rate Bar Plot
    file << "labels = " << "['m=" + labels[0];
    for(size_t i = 1; i < labels.size(); ++i) file << "', 'm=" + labels[i];
    file << "']\n";
    file << "success_rates = [";
    for(size_t i = 0; i < success_rates.size(); ++i) file << success_rates[i] << (i == success_rates.size() - 1 ? "]\n" : ", ");
    file << "axes1[0].bar(labels, success_rates, color='steelblue')\n";
    file << "axes1[0].set_ylabel('Successes (out of 100)')\n";
    file << "axes1[0].set_title('Success Rate')\n";
    file << "axes1[0].set_ylim([0, 105])\n";
    file << "axes1[0].grid(axis='y', alpha=0.3)\n\n";

    // Computation Time Boxplot
    file << "time_data = []\n";
    for (const auto& res : results) {
        file << "time_data.append([";
        for (size_t i = 0; i < res.computation_times.size(); ++i) {
            file << res.computation_times[i] << (i == res.computation_times.size() - 1 ? "" : ", ");
        }
        file << "])\n";
    }
    file << "bp1 = axes1[1].boxplot(time_data, labels=labels, patch_artist=True)\n";
    file << "for patch in bp1['boxes']:\n";
    file << "    patch.set_facecolor('coral')\n";
    file << "axes1[1].set_ylabel('Computation Time (ms)')\n";
    file << "axes1[1].set_title('Computation Time Distribution')\n";
    file << "axes1[1].grid(axis='y', alpha=0.3)\n\n";
    
    // Tree Size Boxplot
    file << "tree_data = []\n";
    for (const auto& res : results) {
        file << "tree_data.append([";
        for (size_t i = 0; i < res.tree_sizes.size(); ++i) {
            file << res.tree_sizes[i] << (i == res.tree_sizes.size() - 1 ? "" : ", ");
        }
        file << "])\n";
    }
    file << "bp2 = axes1[2].boxplot(tree_data, labels=labels, patch_artist=True)\n";
    file << "for patch in bp2['boxes']:\n";
    file << "    patch.set_facecolor('lightgreen')\n";
    file << "axes1[2].set_ylabel('Tree Size (nodes)')\n";
    file << "axes1[2].set_title('Tree Size Distribution')\n";
    file << "axes1[2].grid(axis='y', alpha=0.3)\n\n";
    
    file << "plt.tight_layout()\n";
    file << "plt.savefig('hw8_boxplots.png', dpi=300, bbox_inches='tight')\n";
    
    // --- Figure 2: Average Plots (for Exercise 1.e) ---
    file << "\nfig2, axes2 = plt.subplots(1, 2, figsize=(14, 6))\n";
    file << "fig2.suptitle('HW8 Benchmark Averages vs. Number of Agents', fontsize=16)\n\n";
    
    file << "agent_counts = [";
    for(size_t i = 0; i < agent_counts.size(); ++i) file << agent_counts[i] << (i == agent_counts.size() - 1 ? "]\n" : ", ");
    file << "avg_times = [";
    for(size_t i = 0; i < avg_times.size(); ++i) file << avg_times[i] << (i == avg_times.size() - 1 ? "]\n" : ", ");
    file << "avg_tree_sizes = [";
    for(size_t i = 0; i < avg_tree_sizes.size(); ++i) file << avg_tree_sizes[i] << (i == avg_tree_sizes.size() - 1 ? "]\n" : ", ");

    // Avg. Computation Time vs. m
    file << "axes2[0].plot(agent_counts, avg_times, 'bo-', label='Avg. Time')\n";
    file << "axes2[0].set_xlabel('Number of Agents (m)')\n";
    file << "axes2[0].set_ylabel('Avg. Computation Time (ms)')\n";
    file << "axes2[0].set_title('Avg. Computation Time vs. Agents')\n";
    file << "axes2[0].legend()\n";
    file << "axes2[0].grid(True, alpha=0.3)\n\n";
    
    // Avg. Tree Size vs. m
    file << "axes2[1].plot(agent_counts, avg_tree_sizes, 'go-', label='Avg. Tree Size')\n";
    file << "axes2[1].set_xlabel('Number of Agents (m)')\n";
    file << "axes2[1].set_ylabel('Avg. Tree Size (nodes)')\n";
    file << "axes2[1].set_title('Avg. Tree Size vs. Agents')\n";
    file << "axes2[1].legend()\n";
    file << "axes2[1].grid(True, alpha=0.3)\n\n";

    file << "plt.tight_layout()\n";
    file << "plt.savefig('hw8_average_plots.png', dpi=300, bbox_inches='tight')\n";
    file << "plt.show()\n";
    file << "print('Saved plot: hw8_boxplots.png')\n";
    file << "print('Saved plot: hw8_average_plots.png')\n";
    
    file.close();
    std::cout << "Generated plotting script: " << filename << std::endl;
}

int main(int argc, char** argv) {
    amp::RNG::seed(amp::RNG::randiUnbounded());
    
    std::cout << "\n========== HW8 Central Planner Benchmark ==========" << std::endl;
    std::cout << "Running 100 trials for m = 2, 3, 4, 5, 6..." << std::endl;

    std::vector<BenchmarkResult> all_results;
    std::vector<int> agent_counts = {2, 3, 4, 5, 6};
    
    for (int m : agent_counts) {
        all_results.push_back(runBenchmark(m, 100));
    }
    
    // Generate plotting script
    generatePlottingScript(all_results);

    // Print summary table
    std::cout << "\n\n========== BENCHMARK SUMMARY ==========" << std::endl;
    std::cout << "(n=50000, r=0.5, p_goal=0.05, epsilon=0.25)" << std::endl;
    std::cout << std::left << std::setw(15) << "# Agents (m)" 
              << std::setw(15) << "Successes" 
              << std::setw(20) << "Avg. Tree Size" 
              << "Avg. Comp. Time (ms)" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    
    for (const auto& res : all_results) {
        std::cout << std::left << std::setw(15) << res.num_agents
                  << std::setw(15) << std::to_string(res.successes) + "/100"
                  << std::fixed << std::setprecision(0) << std::setw(20) << res.getAvgTreeSize()
                  << std::fixed << std::setprecision(2) << res.getAvgTime() << std::endl;
    }
    std::cout << "\n==============================================" << std::endl;
    std::cout << "\nTo generate boxplots, run: python plot_hw8_benchmark.py" << std::endl;

    return 0;
}