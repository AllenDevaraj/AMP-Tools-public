#include "MyAStar.h"
#include <queue>
#include <vector>
#include <map>
#include <limits>
#include <algorithm>
#include <list>

struct NodeInfo {
    amp::Node node;
    double f_score;

    bool operator>(const NodeInfo& other) const {
        return f_score > other.f_score;
    }
};

// A* algorithm (Silent version for benchmarking)
amp::AStar::GraphSearchResult MyAStarAlgo::search(const amp::ShortestPathProblem& problem, const amp::SearchHeuristic& heuristic) {
    const std::shared_ptr<amp::Graph<double>>& graph = problem.graph;
    const amp::Node start_node = problem.init_node;
    const amp::Node goal_node = problem.goal_node;
    const std::vector<amp::Node>& all_nodes = graph->nodes();

    std::priority_queue<NodeInfo, std::vector<NodeInfo>, std::greater<NodeInfo>> open_set;
    std::map<amp::Node, amp::Node> came_from; 
    std::map<amp::Node, double> g_score;

    for (const auto& node : all_nodes) {
        g_score[node] = std::numeric_limits<double>::infinity();
    }

    // Start node
    g_score[start_node] = 0.0;
    double start_f_score = heuristic(start_node);
    open_set.push({start_node, start_f_score});

    // A* loop
    while (!open_set.empty()) {
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

    // Failed to find path - return silently
    return {false, {}, 0.0};
}