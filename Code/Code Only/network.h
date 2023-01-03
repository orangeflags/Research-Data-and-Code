#pragma once
using namespace std;

class Network {
private:
	vector<vector<int>> keep_going(vector<int> path, vector<vector<int>> paths) {
		if (path.size() > 1 && is_terminal[path.back()]) {
			paths.push_back(path);
		} else {
			vector<int> neighbours = can_travel_to_nodes(path.back());
			for (int i = 0; i < neighbours.size(); i++) {
				bool append = true;
				for (int j = 0; j < path.size(); j++) if (path[j] == neighbours[i]) append = false;

				if (append) {
					path.push_back(neighbours[i]);
					paths = keep_going(path, paths);
					path.pop_back();
				}
			}
		}
		return paths;
	}

public:
	struct Node {
	public:
		int id;
		double x;
		double y;
		bool terminal;

		Node(int i_id, double i_x, double i_y, int i_terminal) {
			id = i_id;
			x = i_x;
			y = i_y;
			terminal = i_terminal == 1;
		}
	};
	struct Edge {
	public:
		int n1;
		int n2;
		double l;

		Edge(int i_n1, int i_n2, double i_l) {
			n1 = i_n1;
			n2 = i_n2;
			l = i_l;
		}
	};
	struct Demand {
	public:
		int n1;
		int n2;
		int d;

		Demand(int i_n1, int i_n2, int i_d) {
			n1 = i_n1;
			n2 = i_n2;
			d = i_d;
		}
	};

	vector<Node> nodes;
	vector<Edge> edges;
	vector<Demand> demand;

	vector<vector<double>> shortest_travel_time;
	vector<vector<double>> edges_m;

	vector<int> terminals;
	vector<bool> is_terminal;
	int num_nodes;
	
	Network() {
		num_nodes = -1;
	}

	void set_nodes(vector<vector<int>> i_nodes) {
		for (int i = 0; i < i_nodes.size(); i++) {
			Node new_node = Node(i_nodes[i][0] - 1, i_nodes[i][1], i_nodes[i][2], i_nodes[i][3]);
			nodes.push_back(new_node);

			is_terminal.push_back(false);

			if (i_nodes[i][3] == 1) {
				terminals.push_back(i);
				is_terminal[i] = true;
			}

			//terminals.push_back(i);
		}
		num_nodes = (int) i_nodes.size();
	}
	void set_edges(vector<vector<double>> i_edges) {
		for (int i = 0; i < num_nodes; i++) {
			vector<double> row;
			for (int i = 0; i < num_nodes; i++) {
				row.push_back(1024 * 1024);
			}
			edges_m.push_back(row);
		}

		for (int i = 0; i < i_edges.size(); i++) {
			Edge new_edge = Edge((int) i_edges[i][0] - 1, (int) i_edges[i][1] - 1, i_edges[i][2]);
			edges.push_back(new_edge);
			edges_m[edges[i].n1][edges[i].n2] = edges[i].l;
		}
	}
	void set_demand(vector<vector<int>> i_demand) {
		for (int i = 0; i < i_demand.size(); i++) {
			Demand new_demand = Demand(i_demand[i][0] - 1, i_demand[i][1] - 1, i_demand[i][2]);
			demand.push_back(new_demand);
		}
	}
	void calc_shortest_paths() {
		for (int i = 0; i < num_nodes; i++) {
			vector<double> row;
			for (int i = 0; i < num_nodes; i++) {
				row.push_back(1024 * 1024);
			}
			shortest_travel_time.push_back(row);
		}

		for (int i = 0; i < edges.size(); i++) {
			shortest_travel_time[edges[i].n1][edges[i].n2] = edges[i].l;
		}

		for (int i = 0; i < num_nodes; i++) {
			shortest_travel_time[i][i] = 0;
		}

		for (int k = 0; k < num_nodes; k++) {
			for (int i = 0; i < num_nodes; i++) {
				for (int j = 0; j < num_nodes; j++) {
					if (shortest_travel_time[i][j] > shortest_travel_time[i][k] + shortest_travel_time[k][j]) {
						shortest_travel_time[i][j] = shortest_travel_time[i][k] + shortest_travel_time[k][j];
					}
				}
			}
		}
	}

	vector<int> get_terminals() { return terminals; };
	vector<int> can_travel_to_nodes(int node_id) {
		vector<int> output;

		for (int i = 0; i < edges.size(); i++) {
			if (edges[i].n1 == node_id) output.push_back(edges[i].n2);
		}

		std::sort(output.begin(), output.end());

		return output;
	}
	double get_travel_time(int n1, int n2) {
		for (int i = 0; i < edges.size(); i++) if (edges[i].n1 == n1 && edges[i].n2 == n2) return edges[i].l;
		return -1;
	}
	double get_shortest_travel_time(int o_node_index, int d_node_index) {
		return shortest_travel_time[o_node_index][d_node_index];
	}

	double get_travel_time(vector<int> route) {
		double time = 0;

		for (int i = 0; i < route.size() - 1; i++) {
			time += edges_m[route[i]][route[i + 1]];
		}

		return time;
	}
	double get_travel_time_r(vector<int> route) {
		reverse(route.begin(), route.end());
		double time = 0;

		for (int i = 0; i < route.size() - 1; i++) {
			time += edges_m[route[i]][route[i + 1]];
		}

		return time;
	}
	vector<int> generate_shortest_route(int origin_node, int destination_node) {
		int i;

		int current_node = origin_node;
		vector<int> previous_node;
		for (i = 0; i < num_nodes; i++) previous_node.push_back(-1);

		// Dijkstra
		// 01 - Mark all nodes unvisited
		vector<bool> visited;
		for (i = 0; i < num_nodes; i++) visited.push_back(false);

		// 02 - Set tentative distance
		vector<double> distance;
		for (i = 0; i < num_nodes; i++) distance.push_back(1024 * 1024);
		distance[origin_node] = 0;

		while (current_node != destination_node) {
			// 03 - Consider all unvisited neighbours
			vector<int> neighbours = can_travel_to_nodes(current_node);
			for (i = 0; i < neighbours.size(); i++) {
				if (distance[neighbours[i]] > distance[current_node] + get_travel_time(current_node, neighbours[i])) {
					distance[neighbours[i]] = distance[current_node] + get_travel_time(current_node, neighbours[i]);
					previous_node[neighbours[i]] = current_node;
				}
			}

			// 04 - Mark as visited
			visited[current_node] = true;
			int closest_node = -1;
			double closest_distance = 1024 * 1024;
			for (i = 0; i < num_nodes; i++) {
				if (!visited[i] && distance[i] < closest_distance) {
					closest_distance = distance[i];
					closest_node = i;
				}
			}

			// 05 - Continue until destination reached
			current_node = closest_node;
		}

		vector<int> path;
		path.push_back(destination_node);
		while (path[path.size() - 1] != origin_node) {
			path.push_back(previous_node[path[path.size() - 1]]);
		}
		std::reverse(path.begin(), path.end());
		return path;
	}

	vector<int> get_unconnected(vector<direct_route> total_direct_routes) {
		vector<bool> connected;
		for (int i = 0; i < num_nodes; i++) connected.push_back(false);

		for (int i = 0; i < total_direct_routes.size(); i++) for (int j = 0; j < total_direct_routes[i].route.size(); j++) connected[total_direct_routes[i].route[j]] = true;

		vector<int> unconnected;
		for (int i = 0; i < num_nodes; i++) if (connected[i] == false) unconnected.push_back(i);

		return unconnected;
	}
	vector<vector<int>> get_paths_to_terminals(int node) {
		vector<vector<int>> paths;

		vector<int> path;
		path.push_back(node);
		paths = keep_going(path, paths);

		return paths;
	}
};