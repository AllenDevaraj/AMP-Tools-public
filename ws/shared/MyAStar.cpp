#include "MyAStar.h"
#include <queue>
#include <vector>
#include <map>
#include <limits>
#include <algorithm>
#include <iostream>
#include <list>

struct NodeInfo {
    amp::Node node;
    double f_score;

    bool operator>(const NodeInfo& other) const {
        return f_score > other.f_score;
    }
};


// A* algorithm
amp::AStar::GraphSearchResult MyAStarAlgo::search(const amp::ShortestPathProblem& problem, const amp::SearchHeuristic& heuristic) {
    const std::shared_ptr<amp::Graph<double>>& graph = problem.graph;
    const amp::Node start_node = problem.init_node;
    const amp::Node goal_node = problem.goal_node;
    const std::vector<amp::Node>& all_nodes = graph->nodes();

    std::priority_queue<NodeInfo, std::vector<NodeInfo>, std::greater<NodeInfo>> open_set;
    std::map<amp::Node, amp::Node> came_from; 
    std::map<amp::Node, double> g_score;     // Cost from start to current node

    for (const auto& node : all_nodes) {
        g_score[node] = std::numeric_limits<double>::infinity();
    }

    // Start node
    g_score[start_node] = 0.0;
    double start_f_score = heuristic(start_node);
    open_set.push({start_node, start_f_score});
    
    int iterations = 0;
    
    // A* loop
    while (!open_set.empty()) {
        iterations++;
        amp::Node current = open_set.top().node;
        open_set.pop();

        if (current == goal_node) {
            std::list<amp::Node> path_list;
            amp::Node temp = current;
            path_list.push_front(temp);
            while (came_from.count(temp)) {
                temp = came_from[temp];
                path_list.push_front(temp);
            }

            GraphSearchResult result;
            result.success = true;
            result.path_cost = g_score[goal_node];
            result.node_path = path_list;
            
            
            // std::cout << "A* search completed." << std::endl;
            // std::cout << "  - Iterations: " << iterations << std::endl;
            // std::cout << "  - Path Cost: " << result.path_cost << std::endl;
            // std::cout << "  - Path: ";
            bool first = true;
            for (const auto& node : result.node_path) {
                if (!first) {
                    std::cout << " -> ";
                }
                std::cout << node;
                first = false;
            }
            std::cout << std::endl;

            return result;
        }

        const std::vector<amp::Node>& neighbors = graph->children(current);
        const std::vector<double>& edge_weights = graph->outgoingEdges(current); 

        for (size_t i = 0; i < neighbors.size(); ++i) {
            amp::Node neighbor = neighbors[i];
            double weight = edge_weights[i];
            
            double tentative_g_score = g_score[current] + weight;
            if (tentative_g_score < g_score[neighbor]) {
                came_from[neighbor] = current;
                g_score[neighbor] = tentative_g_score;
                double f_score = tentative_g_score + heuristic(neighbor);
                open_set.push({neighbor, f_score});
            }
        }
    }

    std::cout << "A* search failed to find a path from " << start_node << " to " << goal_node << "." << std::endl;
    return {false, {}, 0.0};
}