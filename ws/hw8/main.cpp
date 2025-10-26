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
// Ex 1 Centralized
struct BenchmarkResultCentral {
    int num_agents;
    int successes;
    std::vector<double> tree_sizes;
    std::vector<double> computation_times;
    
    double getAvgTreeSize() const {
        if (tree_sizes.empty()) return 0.0;
        double sum = 0.0; for (double size : tree_sizes) sum += size; return sum / tree_sizes.size();
    }
    double getAvgTime() const {
        if (computation_times.empty()) return 0.0;
        double sum = 0.0; for (double t : computation_times) sum += t; return sum / computation_times.size();
    }
};

// Ex 2 Decentralized
struct BenchmarkResultDecentral {
    int num_agents;
    int successes;
    std::vector<double> computation_times;
    
    double getAvgTime() const {
        if (computation_times.empty()) return 0.0;
        double sum = 0.0; for (double t : computation_times) sum += t; return sum / computation_times.size();
    }
};

BenchmarkResultCentral runBenchmarkCentral(int num_agents, int num_runs = 100) {
    BenchmarkResultCentral result;
    result.num_agents = num_agents;
    result.successes = 0;
    std::cout << "  Benchmarking Central (m=" << num_agents << ")..." << std::flush;
    for (int i = 0; i < num_runs; ++i) {
        if ((i + 1) % 10 == 0) std::cout << "." << std::flush;
        MyCentralPlanner planner;
        MultiAgentProblem2D problem = HW8::getWorkspace1(num_agents);
        MultiAgentPath2D path = planner.plan(problem);
        if (!path.agent_paths.empty() && !path.agent_paths[0].waypoints.empty()) {
            result.successes++;
            result.tree_sizes.push_back(static_cast<double>(planner.last_tree_size)); 
            result.computation_times.push_back(planner.last_comp_time_ms);
        }
    }
    std::cout << " Done! Success: " << result.successes << "/" << num_runs << std::endl;
    return result;
}

BenchmarkResultDecentral runBenchmarkDecentral(int num_agents, int num_runs = 100) {
    BenchmarkResultDecentral result;
    result.num_agents = num_agents;
    result.successes = 0;
    std::cout << "  Benchmarking Decentral (m=" << num_agents << ")..." << std::flush;
    for (int i = 0; i < num_runs; ++i) {
        if ((i + 1) % 10 == 0) std::cout << "." << std::flush;
        MyDecentralPlanner planner;
        MultiAgentProblem2D problem = HW8::getWorkspace1(num_agents);
        MultiAgentPath2D path = planner.plan(problem);
        if (!path.agent_paths.empty() && !path.agent_paths.back().waypoints.empty()) {
            result.successes++;
            result.computation_times.push_back(planner.last_comp_time_ms);
        }
    }
    std::cout << " Done! Success: " << result.successes << "/" << num_runs << std::endl;
    return result;
}

// Ex 1
void generatePlottingScriptCentral(const std::vector<BenchmarkResultCentral>& results) {
    std::string filename = "plot_hw8_centralized_benchmark.py";
    std::ofstream file(filename);
    file << "import matplotlib.pyplot as plt\nimport numpy as np\n\n";
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
    file << "fig1, axes1 = plt.subplots(1, 3, figsize=(18, 6))\n";
    file << "fig1.suptitle('HW8 Centralized RRT Benchmark (100 Runs)', fontsize=16)\n\n";
    file << "labels = " << "['m=" + labels[0]; for(size_t i = 1; i < labels.size(); ++i) file << "', 'm=" + labels[i]; file << "']\n";
    file << "success_rates = ["; for(size_t i = 0; i < success_rates.size(); ++i) file << success_rates[i] << (i == success_rates.size() - 1 ? "]\n" : ", ");
    file << "axes1[0].bar(labels, success_rates, color='steelblue')\n";
    file << "axes1[0].set_ylabel('Successes (out of 100)')\naxes1[0].set_title('Success Rate')\naxes1[0].set_ylim([0, 105])\naxes1[0].grid(axis='y', alpha=0.3)\n\n";
    file << "time_data = []\n";
    for (const auto& res : results) { file << "time_data.append(["; for (size_t i = 0; i < res.computation_times.size(); ++i) file << res.computation_times[i] << (i == res.computation_times.size() - 1 ? "" : ", "); file << "])\n"; }
    file << "bp1 = axes1[1].boxplot(time_data, labels=labels, patch_artist=True)\n";
    file << "for patch in bp1['boxes']:\n    patch.set_facecolor('coral')\n";
    file << "axes1[1].set_ylabel('Computation Time (ms)')\naxes1[1].set_title('Computation Time Distribution')\naxes1[1].grid(axis='y', alpha=0.3)\n\n";
    file << "tree_data = []\n";
    for (const auto& res : results) { file << "tree_data.append(["; for (size_t i = 0; i < res.tree_sizes.size(); ++i) file << res.tree_sizes[i] << (i == res.tree_sizes.size() - 1 ? "" : ", "); file << "])\n"; }
    file << "bp2 = axes1[2].boxplot(tree_data, labels=labels, patch_artist=True)\n";
    file << "for patch in bp2['boxes']:\n    patch.set_facecolor('lightgreen')\n";
    file << "axes1[2].set_ylabel('Tree Size (nodes)')\naxes1[2].set_title('Tree Size Distribution')\naxes1[2].grid(axis='y', alpha=0.3)\n\n";
    file << "plt.tight_layout()\nplt.savefig('hw8_centralized_boxplots.png', dpi=300, bbox_inches='tight')\n";
    file << "\nfig2, axes2 = plt.subplots(1, 2, figsize=(14, 6))\n";
    file << "fig2.suptitle('HW8 Centralized Averages vs. Number of Agents', fontsize=16)\n\n";
    file << "agent_counts = ["; for(size_t i = 0; i < agent_counts.size(); ++i) file << agent_counts[i] << (i == agent_counts.size() - 1 ? "]\n" : ", ");
    file << "avg_times = ["; for(size_t i = 0; i < avg_times.size(); ++i) file << avg_times[i] << (i == avg_times.size() - 1 ? "]\n" : ", ");
    file << "avg_tree_sizes = ["; for(size_t i = 0; i < avg_tree_sizes.size(); ++i) file << avg_tree_sizes[i] << (i == avg_tree_sizes.size() - 1 ? "]\n" : ", ");
    file << "axes2[0].plot(agent_counts, avg_times, 'bo-', label='Avg. Time')\n";
    file << "axes2[0].set_xlabel('Number of Agents (m)')\naxes2[0].set_ylabel('Avg. Computation Time (ms)')\naxes2[0].set_title('Avg. Computation Time vs. Agents')\naxes2[0].legend()\naxes2[0].grid(True, alpha=0.3)\n\n";
    file << "axes2[1].plot(agent_counts, avg_tree_sizes, 'go-', label='Avg. Tree Size')\n";
    file << "axes2[1].set_xlabel('Number of Agents (m)')\naxes2[1].set_ylabel('Avg. Tree Size (nodes)')\naxes2[1].set_title('Avg. Tree Size vs. Agents')\naxes2[1].legend()\naxes2[1].grid(True, alpha=0.3)\n\n";
    file << "plt.tight_layout()\nplt.savefig('hw8_centralized_average_plots.png', dpi=300, bbox_inches='tight')\n";
    file << "print('Generated " << filename << "')\n"; file.close();
}

// Ex 2
void generatePlottingScriptDecentral(const std::vector<BenchmarkResultDecentral>& results) {
    std::string filename = "plot_hw8_decentralized_benchmark.py";
    std::ofstream file(filename);
    file << "import matplotlib.pyplot as plt\nimport numpy as np\n\n";
    std::vector<std::string> labels;
    std::vector<double> avg_times;
    std::vector<double> success_rates;
    std::vector<int> agent_counts;
    for (const auto& res : results) {
        labels.push_back(std::to_string(res.num_agents));
        agent_counts.push_back(res.num_agents);
        avg_times.push_back(res.getAvgTime());
        success_rates.push_back(static_cast<double>(res.successes));
    }
    file << "fig1, axes1 = plt.subplots(1, 2, figsize=(14, 6))\n";
    file << "fig1.suptitle('HW8 Decentralized (Prioritized) RRT Benchmark (100 Runs)', fontsize=16)\n\n";
    file << "labels = " << "['m=" + labels[0]; for(size_t i = 1; i < labels.size(); ++i) file << "', 'm=" + labels[i]; file << "']\n";
    file << "success_rates = ["; for(size_t i = 0; i < success_rates.size(); ++i) file << success_rates[i] << (i == success_rates.size() - 1 ? "]\n" : ", ");
    file << "axes1[0].bar(labels, success_rates, color='steelblue')\n";
    file << "axes1[0].set_ylabel('Successes (out of 100)')\naxes1[0].set_title('Success Rate')\naxes1[0].set_ylim([0, 105])\naxes1[0].grid(axis='y', alpha=0.3)\n\n";
    file << "time_data = []\n";
    for (const auto& res : results) { file << "time_data.append(["; for (size_t i = 0; i < res.computation_times.size(); ++i) file << res.computation_times[i] << (i == res.computation_times.size() - 1 ? "" : ", "); file << "])\n"; }
    file << "bp1 = axes1[1].boxplot(time_data, labels=labels, patch_artist=True)\n";
    file << "for patch in bp1['boxes']:\n    patch.set_facecolor('coral')\n";
    file << "axes1[1].set_ylabel('Computation Time (ms)')\naxes1[1].set_title('Computation Time Distribution')\naxes1[1].grid(axis='y', alpha=0.3)\n\n";
    file << "plt.tight_layout()\nplt.savefig('hw8_decentralized_boxplots.png', dpi=300, bbox_inches='tight')\n";
    file << "\nfig2, ax2 = plt.subplots(1, 1, figsize=(7, 6))\n";
    file << "fig2.suptitle('HW8 Decentralized Averages vs. Number of Agents', fontsize=16)\n\n";
    file << "agent_counts = ["; for(size_t i = 0; i < agent_counts.size(); ++i) file << agent_counts[i] << (i == agent_counts.size() - 1 ? "]\n" : ", ");
    file << "avg_times = ["; for(size_t i = 0; i < avg_times.size(); ++i) file << avg_times[i] << (i == avg_times.size() - 1 ? "]\n" : ", ");
    file << "ax2.plot(agent_counts, avg_times, 'bo-', label='Avg. Time')\n";
    file << "ax2.set_xlabel('Number of Agents (m)')\nax2.set_ylabel('Avg. Computation Time (ms)')\nax2.set_title('Avg. Computation Time vs. Agents')\nax2.legend()\nax2.grid(True, alpha=0.3)\n\n";
    file << "plt.tight_layout()\nplt.savefig('hw8_decentralized_average_plot.png', dpi=300, bbox_inches='tight')\n";
    file << "print('Generated " << filename << "')\n";
    file << "plt.show()\n";
    file.close();
}

int main(int argc, char** argv) {
    amp::RNG::seed(amp::RNG::randiUnbounded());
    std::vector<int> agent_counts = {2, 3, 4, 5, 6};
    for (int m : agent_counts) {
        std::cout << "\nVisualizing for m = " << m << " " << std::endl;
        // Centralized (Ex 1b)
        std::cout << "Running Centralized Planner (m=" << m << ")" << std::flush;
        MyCentralPlanner central_planner;
        MultiAgentProblem2D problem_cen = HW8::getWorkspace1(m);
        MultiAgentPath2D path_cen;
        int central_attempts = 0;
        do {
            path_cen = central_planner.plan(problem_cen);
            central_attempts++;
            if (central_attempts > 1 && central_attempts % 10 == 0) std::cout << "." << std::flush; // Show progress
        } while ((path_cen.agent_paths.empty() || path_cen.agent_paths[0].waypoints.empty()) && central_attempts < 100); // 100 attempt limit

        if (central_attempts >= 100) {
            std::cout << " FAILED TO VISUALIZE (100 attempts)" << std::endl;
        } else {
            std::cout << " Found in " << central_attempts << " attempt(s)." << std::endl;
        }
        // Visualizer::makeFigure(problem_cen, path_cen);

        // Decentralized (Ex 2b)
        std::cout << "Running Decentralized Planner (m=" << m << ")..." << std::flush;
        MyDecentralPlanner decentral_planner;
        MultiAgentProblem2D problem_dec = HW8::getWorkspace1(m);
        MultiAgentPath2D path_dec;
        int decentral_attempts = 0;
        do {
            path_dec = decentral_planner.plan(problem_dec);
            decentral_attempts++;
        } while ((path_dec.agent_paths.empty() || path_dec.agent_paths.back().waypoints.empty()) && decentral_attempts < 100); // 100 attempt limit

        if (decentral_attempts >= 100) {
            std::cout << " FAILED TO VISUALIZE (100 attempts)" << std::endl;
        } else {
            std::cout << " Found in " << decentral_attempts << " attempt(s)." << std::endl;
        }
        // Visualizer::makeFigure(problem_dec, path_dec);
    }
    // Benchmark for Ex 1 (Centralized)
    std::cout << "\n            HW8 Central Planner Benchmark (Ex 1)           " << std::endl;
    std::vector<BenchmarkResultCentral> central_results;
    for (int m : agent_counts) {
        central_results.push_back(runBenchmarkCentral(m, 100));
    }
    generatePlottingScriptCentral(central_results);
    std::cout << "\n--- Centralized Benchmark Summary ---" << std::endl;
    std::cout << std::left << std::setw(15) << "# Agents (m)" << std::setw(15) << "Successes" << std::setw(20) << "Avg. Tree Size" << "Avg. Comp. Time (ms)" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    for (const auto& res : central_results) {
        std::cout << std::left << std::setw(15) << res.num_agents
                  << std::setw(15) << std::to_string(res.successes) + "/100"
                  << std::fixed << std::setprecision(0) << std::setw(20) << res.getAvgTreeSize()
                  << std::fixed << std::setprecision(2) << res.getAvgTime() << std::endl;
    }

    // Benchmark for Ex 2 (Decentralized)
    std::cout << "\n           HW8 Decentralized Planner Benchmark (Ex 2)           " << std::endl;
    std::vector<BenchmarkResultDecentral> decentral_results;
    for (int m : agent_counts) {
        decentral_results.push_back(runBenchmarkDecentral(m, 100));
    }
    generatePlottingScriptDecentral(decentral_results);
    std::cout << "\n--- Decentralized Benchmark Summary ---" << std::endl;
    std::cout << std::left << std::setw(15) << "# Agents (m)" << std::setw(15) << "Successes" << "Avg. Comp. Time (ms)" << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;
    for (const auto& res : decentral_results) {
        std::cout << std::left << std::setw(15) << res.num_agents
                  << std::setw(15) << std::to_string(res.successes) + "/100"
                  << std::fixed << std::setprecision(2) << res.getAvgTime() << std::endl;
    }

    std::cout << "\nAll visualizations and benchmarks complete." << std::endl;
    
    // Visualizer::saveFigures();
    HW8::grade<MyCentralPlanner, MyDecentralPlanner>("AllenDevaraj.AugustinPonraj@colorado.edu", argc, argv, std::make_tuple(), std::make_tuple());
    return 0;
}