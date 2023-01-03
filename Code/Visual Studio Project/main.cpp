// initial_implementation.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#pragma warning(disable : 4996)

// Include List
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <random>

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include "routes.h"
#include "network.h"
#include "model.h"
#include "CPLEX_Model.h"

using namespace std;

// Global Definitions
auto main_start = chrono::high_resolution_clock::now();

string networkName;
int BUSES;
int CAP;
int MAX_ROUTES;
int MIN_NODES_PER_ROUTE;
int MAX_NODES_PER_ROUTE;
int MAX_ROUTE_LENGTH;

int data_output_type;

#define MAX_THREADS_CPLEX 8
#define MAX_TIME_CPLEX 3600

int attempts = 0;
int solves = 0;

string path_INPUT;
string path_OUTPUT;

// Function Declaration
void smaller_data_set();
void computation_time();
void other_solutions();
void proposed_method();

void alt_main_shortest(string);
void alt_main_all(string);
void alt_main_subset(int);
void improve_solution(string, int);
void add_routes(vector<direct_route>, vector<direct_route>&, Network&, float, vector<direct_route>&, float&);
vector<vector<int>> read_csv_int(string);
vector<vector<double>> read_csv_double(string);
void import_network(Network&);
void generate_routes(Network&, vector<direct_route>&);
void generate_routes_shortest(Network&, vector<direct_route>&);
void recursive_construction(vector<direct_route>&, vector<int>, double, Network&);
void route_prune(Network&, vector<direct_route>&, vector<direct_route>&);
void process_potential_transfers(vector<direct_route>&);
void route_joining(vector<transfer_route>&, vector<direct_route>&);
void intersect_routes(vector<transfer_route>&, vector<int>, vector<int>, int, int);
void remove_part(vector<transfer_route>&, vector<int>, vector<int>);
void connect_disconnected_nodes(Network&, vector<direct_route>&);
float optimise_problem(vector<direct_route>&, vector<transfer_route>&, Network&, bool);
void export_output(vector<direct_route>&, vector<transfer_route>&, Model&, IloCplex&, CPLEX_Model&, string);
void check_other_solutions(string);
void import_routeset(Network&, vector<direct_route>&, string);

int main() {
	//// Smaller Data Set
	//smaller_data_set();

	//// Computation Time (Input Size Variation)
	//computation_time();

	//// Other Mandl Solutions
	//other_solutions();

	// Proposed Method
	attempts = 0;
	proposed_method();
}

void smaller_data_set() {
	cout << "Smaller Data Set:\n\n";

	BUSES = 20;
	CAP = 50;
	MAX_ROUTES = 6;
	MIN_NODES_PER_ROUTE = 2;
	MAX_NODES_PER_ROUTE = 5;
	MAX_ROUTE_LENGTH = 50;

	ofstream file;
	auto sub_start = chrono::high_resolution_clock::now();
	auto end = chrono::high_resolution_clock::now();
	auto timetaken = chrono::duration_cast<chrono::milliseconds>(end - sub_start);

	int optimal_solution_output;

	networkName = "toy_network";

	string sets[] = { "10-a", "10-b", "10-c", "10-d", "10-e", "10-f", "10-g",
					  "11-a", "11-b", "11-c", "11-d", "11-e", "11-f", "11-g" };
	for (string set : sets) {
		path_INPUT = "../_input/" + networkName + "/";
		path_INPUT += set + "/";

		//Step 0: Optimal
		cout << networkName << "-" << set << ": Optimal\n";
		path_OUTPUT = "../_output/smaller_data_set/" + set + "/Optimal/";
		auto sub_start = chrono::high_resolution_clock::now();
		data_output_type = 1;
		alt_main_all(path_OUTPUT);
		end = chrono::high_resolution_clock::now();
		timetaken = chrono::duration_cast<chrono::milliseconds>(end - sub_start);
		file.open("../_output/smaller_data_set/Timings/Optimal.csv", std::ios_base::app);
		file << networkName << " - " << set << "; " << "Time (ms): " << timetaken.count() << ";\n";
		file.close();
		
		//Step 1:
		cout << networkName << "-" << set << ": Proposed Method Step 1\n";
		path_OUTPUT = "../_output/smaller_data_set/" + set + "/Initial/";
		attempts = 1;
		sub_start = chrono::high_resolution_clock::now();
		data_output_type = 2;
		alt_main_shortest(path_OUTPUT);
		end = chrono::high_resolution_clock::now();
		timetaken = chrono::duration_cast<chrono::milliseconds>(end - sub_start);
		file.open("../_output/smaller_data_set/Timings/Initial.csv", std::ios_base::app);
		file << networkName << " - " << set << "; " << "Time (ms): " << timetaken.count() << ";\n";
		file.close();

		//Step 2:
		cout << networkName << "-" << set << ": Proposed Method Step 2\n";
		path_OUTPUT = "../_output/smaller_data_set/" + set + "/Improved/";
		attempts = 100;
		sub_start = chrono::high_resolution_clock::now();
		data_output_type = 3;
		improve_solution("../_output/smaller_data_set/" + set + "/Initial/Solution.csv", 1);
		end = chrono::high_resolution_clock::now();
		timetaken = chrono::duration_cast<chrono::milliseconds>(end - sub_start);
		file.open("../_output/smaller_data_set/Timings/Improved.csv", std::ios_base::app);
		file << networkName << " - " << set << "; " << "Time (ms): " << timetaken.count() << ";\n";
		file.close();

		cout << "\n\n";
	}
}
void computation_time() {
	cout << "Computation Time:\n\n";
	data_output_type = 100;

	BUSES = 20;
	CAP = 50;
	MAX_ROUTES = 8;
	MIN_NODES_PER_ROUTE = 3;
	MAX_NODES_PER_ROUTE = 8;
	MAX_ROUTE_LENGTH = 50;

	networkName = "mandls_network";
	path_INPUT = "../_input/" + networkName + "/";
	path_OUTPUT = "../_output/computation_time/";

	// Incrementally Increase Number of Input Routes
	for (int i = 10; i <= 150; i += 10) {
		alt_main_subset(i);
	}
}
void other_solutions() {
	cout << "Others Solutions" << "\n\n";
	networkName = "mandls_network";
	data_output_type = 200;

	BUSES = 20;
	CAP = 50;
	MIN_NODES_PER_ROUTE = 3;
	MAX_NODES_PER_ROUTE = 8;
	MAX_ROUTE_LENGTH = 50;

	path_INPUT = "../_input/" + networkName + "/";
	path_OUTPUT = "../_output/others_solutions/";

	ofstream file;
	file.open(path_OUTPUT + "_Summary.csv");
	file << "dataset,d_0,d_1,d_2,d_u,ATT\n";
	file.close();

	string dataset_names[] = { "1979-Mandl", "1991-Baaj", "1994-Shih", "2002-Chakroborty", "2009-Fan(a)", "2009-Fan(b)", "2010-Fan", "2010-Zhang(a)", "2010-Zhang(b)", "2010-Zhang(c)", "2011-Bagloee", "2012-Chew", "2013-Kilic", "2013-Nikolic", "2013-Yan", "2014-Kechagiopoulos", "2014-Kilic(a)", "2014-Kilic(b)", "2014-Nayeem", "2014-Nikolic(a)", "2014-Nikolic(b)", "2015-Amiripour", "2015-Arbex", "2016-Buba", "2018-Buba", "2019-Jha", "2019-Moghaddam", "2020-Capali", "2020-Katsaragakis", "2021-Ahern(a)", "2021-Ahern(b)", "2021-Ahern(d)" };
	string dataset_sizes[] = { "4", "6", "7", "8", "10", "12" };

	for (string size : dataset_sizes) {
		for (string name : dataset_names) {
			MAX_ROUTES = stoi(size);

			string temp = name + "-" + size;

			cout << "\t" << temp << "\n";

			std::ifstream input_csv(path_INPUT + "others_solutions/" + temp + ".csv");
			if (input_csv) {
				check_other_solutions(temp);
			} else {
				file.open(path_OUTPUT + "_Summary.csv", std::ios_base::app);
				file << temp << ",-,-,-,-,-\n";
				file.close();
			}
		}
	}
}
void proposed_method() {
	cout << "Proposed Method:\n\n";

	BUSES = 20;
	CAP = 50;
	MAX_ROUTES = 12;
	MIN_NODES_PER_ROUTE = 3;
	MAX_NODES_PER_ROUTE = 8;
	MAX_ROUTE_LENGTH = 50;

	networkName = "mandls_network";
	path_INPUT = "../_input/" + networkName + "/";
	path_OUTPUT = "../_output/proposed_method/";

	ofstream file;
	auto sub_start = chrono::high_resolution_clock::now();
	auto end = chrono::high_resolution_clock::now();
	auto timetaken = chrono::duration_cast<chrono::milliseconds>(end - sub_start);

	//Step 1:
	cout << MAX_ROUTES << " Routes" << ": Step 1\n";
	path_OUTPUT = "../_output/proposed_method/" + to_string(MAX_ROUTES) + "-" + to_string(BUSES) + "/Initial/";
	sub_start = chrono::high_resolution_clock::now();
	data_output_type = 2;
	alt_main_shortest(path_OUTPUT);
	end = chrono::high_resolution_clock::now();
	timetaken = chrono::duration_cast<chrono::milliseconds>(end - sub_start);
	file.open("../_output/proposed_method/Timings/Initial.csv", std::ios_base::app);
	file << MAX_ROUTES << " - " << BUSES << "; " << "Time (ms): " << timetaken.count() << ";\n";
	file.close();

	//Step 2:
	cout << MAX_ROUTES << " Routes" << ": Step 2\n";
	path_OUTPUT = "../_output/proposed_method/" + to_string(MAX_ROUTES) + "-" + to_string(BUSES) + "/Improved/";
	attempts = 100;
	sub_start = chrono::high_resolution_clock::now();
	data_output_type = 3;
	improve_solution("../_output/proposed_method/" + to_string(MAX_ROUTES) + "-" + to_string(BUSES) + "/Initial/Solution.csv", 1);
	end = chrono::high_resolution_clock::now();
	timetaken = chrono::duration_cast<chrono::milliseconds>(end - sub_start);
	file.open("../_output/proposed_method/Timings/Improved.csv", std::ios_base::app);
	file << MAX_ROUTES << " - " << BUSES << "; " << "Time (ms): " << timetaken.count() << ";\n";
	file.close();
}

void alt_main_all(string INPUT_path_OUTPUT) {
	auto start = chrono::high_resolution_clock::now();
	int i;

	// 00 - Import network
	Network network;
	import_network(network);

	// 01 - Generate set of routes (and transfers)
	vector<direct_route> new_routes;
	vector<direct_route> total_direct_routes;
	vector<transfer_route> total_transfer_routes;
	generate_routes(network, new_routes);
	route_prune(network, total_direct_routes, new_routes);
	for (i = 0; i < new_routes.size(); i++) {
		if (new_routes[i].to_keep) total_direct_routes.push_back(new_routes[i]);
	}
	new_routes.clear();

	int obj_value = optimise_problem(total_direct_routes, total_transfer_routes, network, true);

	ofstream file;

	file.open(INPUT_path_OUTPUT + "Solution.csv");
	int outputted = 0;
	for (int r = 0; r < total_direct_routes.size(); r++) {
		if (total_direct_routes[r].used) {
			for (int s = 0; s < total_direct_routes[r].route.size(); s++) {
				file << total_direct_routes[r].route[s];
				file << (s == total_direct_routes[r].route.size() - 1 ? "" : ",");
			}
			outputted++;

			if (outputted < MAX_ROUTES) file << "\n";
		}
	}
	file.close();
}
void alt_main_subset(int num_input_routes) {
	auto start = chrono::high_resolution_clock::now();
	int i;

	attempts = num_input_routes;

	// 00 - Import network
	Network network;
	import_network(network);

	// 01 - Generate set of routes (and transfers)
	vector<direct_route> new_routes;
	vector<direct_route> holding_routes;
	vector<direct_route> total_direct_routes;
	vector<transfer_route> total_transfer_routes;
	generate_routes(network, new_routes);
	route_prune(network, total_direct_routes, new_routes);
	for (i = 0; i < new_routes.size(); i++) {
		if (new_routes[i].to_keep) holding_routes.push_back(new_routes[i]);
	}
	new_routes.clear();

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle(holding_routes.begin(), holding_routes.end(), std::default_random_engine(seed));

	for (i = 0; i < num_input_routes; i++) {
		total_direct_routes.push_back(holding_routes[i]);
	}
	holding_routes.clear();

	ofstream file;
	file.open(path_OUTPUT + "Timings.csv", std::ios_base::app);
	file << "Routes: " << num_input_routes << "; Choose: " << MAX_ROUTES << "; ";
	file.close();

	start = chrono::high_resolution_clock::now();
	int obj_value = optimise_problem(total_direct_routes, total_transfer_routes, network, true);

	auto end = chrono::high_resolution_clock::now();
	auto timetaken = chrono::duration_cast<chrono::milliseconds>(end - start);

	cout << "\t" << MAX_ROUTES << " from \t" << num_input_routes << " routes \t (" << timetaken.count() << " ms)\n";

	file.open(path_OUTPUT + "Timings.csv", std::ios_base::app);
	file << "Time (ms): " << timetaken.count() << "\n";
	file.close();
}
void alt_main_shortest(string INPUT_path_OUTPUT) {
	int i;

	float most_recent_obj_value;
	float second_most_recent_obj_value;
	vector<direct_route> total_direct_routes;
	vector<transfer_route> total_transfer_routes;
	vector<direct_route> new_routes;

	// 00 - Import network
	Network network;
	import_network(network);

	// 01 - Generate set of shortest routes (and transfers)
	generate_routes_shortest(network, new_routes);
	route_prune(network, total_direct_routes, new_routes);
	for (i = 0; i < new_routes.size(); i++) {
		if (new_routes[i].to_keep) total_direct_routes.push_back(new_routes[i]);
	}
	new_routes.clear();

	// Inside optimise_problem()
	// 02 - Optimise initial set of routes
	// 03 - Setup for iterations
	attempts++;

	auto start = chrono::high_resolution_clock::now();
	most_recent_obj_value = optimise_problem(total_direct_routes, total_transfer_routes, network, true);

	auto end = chrono::high_resolution_clock::now();
	auto timetaken = chrono::duration_cast<chrono::milliseconds>(end - start);


	string path = INPUT_path_OUTPUT + "Solution";
	ofstream file;
	file.open(path + ".csv");
	int outputted = 0;
	for (int r = 0; r < total_direct_routes.size(); r++) {
		if (total_direct_routes[r].used) {
			for (int s = 0; s < total_direct_routes[r].route.size(); s++) {
				file << total_direct_routes[r].route[s];
				file << (s == total_direct_routes[r].route.size() - 1 ? "" : ",");
			}
			outputted++;

			if (outputted < MAX_ROUTES) file << "\n";
		}
	}
	file.close();

	cout << "\t" << attempts << "\t" << most_recent_obj_value << "\n";

	second_most_recent_obj_value = -1;

	// 04 - Loop until new solution values are the same as previous
	while (most_recent_obj_value != second_most_recent_obj_value) {
		second_most_recent_obj_value = most_recent_obj_value;

		// 03a - Create set of routes from base set and transfers used
		for (i = 0; i < total_transfer_routes.size(); i++) {
			if (total_transfer_routes[i].used) {
				direct_route temp(total_transfer_routes[i].route, network.get_travel_time(total_transfer_routes[i].route), network.get_travel_time_r(total_transfer_routes[i].route));
				new_routes.push_back(temp);
			}
		}
		route_prune(network, total_direct_routes, new_routes);

		// Add to 'Master Route List'
		for (i = 0; i < new_routes.size(); i++) {
			if (new_routes[i].to_keep) total_direct_routes.push_back(new_routes[i]);
		}
		new_routes.clear();

		// 03c - Save solution values
		attempts++;
		start = chrono::high_resolution_clock::now();
		most_recent_obj_value = optimise_problem(total_direct_routes, total_transfer_routes, network, true);

		end = chrono::high_resolution_clock::now();
		timetaken = chrono::duration_cast<chrono::milliseconds>(end - start);

		cout << "\t" << attempts << "\t" << most_recent_obj_value << "\n";

		if (most_recent_obj_value > second_most_recent_obj_value) break;
		// 04 (as above) - Loop until new solution values are the same as previous

		path = INPUT_path_OUTPUT + "Solution";
		file.open(path + ".csv");

		int outputted = 0;
		for (int r = 0; r < total_direct_routes.size(); r++) {
			if (total_direct_routes[r].used) {
				for (int s = 0; s < total_direct_routes[r].route.size(); s++) {
					file << total_direct_routes[r].route[s];
					file << (s == total_direct_routes[r].route.size() - 1 ? "" : ",");
				}
				outputted++;

				if (outputted < MAX_ROUTES) file << "\n";
			}
		}
		file.close();
	}
}

void improve_solution(string routeset, int option) {
	// 00 - Import network
	Network network;
	import_network(network);

	// 01 - Import Solution
	vector<direct_route> solution_to_improve;
	import_routeset(network, solution_to_improve, routeset);

	// Generate additional routes to test
	vector<direct_route> potential_additions;
	generate_routes(network, potential_additions);

	vector<transfer_route> temp;

	vector<int> most_recent_adds;

	float best_obj = optimise_problem(solution_to_improve, temp, network, true);
	float initial_obj = best_obj;
	vector<direct_route> best_routes;
	for (int i = 0; i < solution_to_improve.size(); i++) best_routes.push_back(solution_to_improve[i]);

	// Test additions recursively
	add_routes(solution_to_improve, potential_additions, network, initial_obj, best_routes, best_obj);
}
void add_routes(vector<direct_route> initial_final, vector<direct_route>& potential_additions, Network& network, float most_recent_obj, vector<direct_route>& best_routes, float& best_obj) {
	for (int i = 0; i < potential_additions.size(); i++) {
		initial_final.push_back(potential_additions[i]);
		vector<transfer_route> initial_final_transfers;

		float objective_after = optimise_problem(initial_final, initial_final_transfers, network, true);

		if (objective_after < most_recent_obj) {
			vector<direct_route> new_network;
			for (int j = 0; j < initial_final.size(); j++) {
				if (initial_final[j].used) new_network.push_back(initial_final[j]);
			}

			if (objective_after < best_obj) {
				best_obj = objective_after;
				best_routes.clear();
				for (int j = 0; j < new_network.size(); j++) best_routes.push_back(new_network[j]);

				string path = path_OUTPUT + (networkName == "mandls_network" ? (to_string(MAX_ROUTES) + "-" + to_string(BUSES) + "/Improved/") : "") + to_string(attempts) + " - final routes";

				ofstream file;
				file.open(path + ".csv");
				for (int r = 0; r < best_routes.size(); r++) {
					if (best_routes[r].used) {
						for (int s = 0; s < best_routes[r].route.size(); s++) {
							file << best_routes[r].route[s];
							file << (s == best_routes[r].route.size() - 1 ? "" : ",");
						}
						file << (r == best_routes.size() - 1 ? "" : "\n");
					}
				}
				file.close();

				attempts++;
			}

			std::cout << "\t" << i << "\t" << objective_after << "\n";

			//return;
			add_routes(new_network, potential_additions, network, objective_after, best_routes, best_obj);
			return;
		}

		initial_final.pop_back();
	}
}

vector<vector<int>> read_csv_int(string file_path) {
	std::ifstream input_csv(file_path);

	vector<vector<int>> full_input;
	string line;

	//getline(input_csv, line);

	if (input_csv.is_open()) {
		while (input_csv) {
			vector<int> row;
			getline(input_csv, line);

			size_t start = 0;
			size_t end = line.find(",");
			while (end != -1) {
				try {
					row.push_back(stoi(line.substr(start, end - start)));
				} catch (exception) {
					row.push_back(stoi(line.substr(start + 3, end - (start + 3))));
				}
				start = end + 1;
				end = line.find(",", start);
			}
			row.push_back(stoi(line.substr(start, end - start)));

			full_input.push_back(row);

			if (input_csv.eof()) {
				input_csv.close();
				break;
			}
		}
	}

	return full_input;
}
vector<vector<double>> read_csv_double(string file_path) {
	std::ifstream input_csv(file_path);

	vector<vector<double>> full_input;
	string line;

	//getline(input_csv, line);

	if (input_csv.is_open()) {
		while (input_csv) {
			vector<double> row;
			getline(input_csv, line);

			size_t start = 0;
			size_t end = line.find(",");
			while (end != -1) {
				try {
					row.push_back(stod(line.substr(start, end - start)));
				} catch (exception) {
					row.push_back(stod(line.substr(start + 3, end - (start + 3))));
				}
				start = end + 1;
				end = line.find(",", start);
			}
			row.push_back(stoi(line.substr(start, end - start)));

			full_input.push_back(row);

			if (input_csv.eof()) {
				input_csv.close();
				break;
			}
		}
	}

	return full_input;
}
void import_network(Network& network) {
	vector<vector<int>> i_nodes = read_csv_int(path_INPUT + "nodes.csv");
	vector<vector<double>> i_edges = read_csv_double(path_INPUT + "links.csv");
	vector<vector<int>> i_demand = read_csv_int(path_INPUT + "demand.csv");

	network.set_nodes(i_nodes);
	network.set_edges(i_edges);
	network.set_demand(i_demand);
	network.calc_shortest_paths();
}

void generate_routes(Network& network, vector<direct_route>& total_direct_routes) {
	vector<direct_route> temp_direct_routes;

	for (int n = 0; n < network.terminals.size(); n++) {
		vector<int> route_start;
		route_start.push_back(network.terminals[n]);

		// Generate
		recursive_construction(temp_direct_routes, route_start, 0, network);

		// Prune
		route_prune(network, total_direct_routes, temp_direct_routes);

		// Add to 'Master Route List'
		for (int i = 0; i < temp_direct_routes.size(); i++) {
			if (temp_direct_routes[i].to_keep) total_direct_routes.push_back(temp_direct_routes[i]);
		}
		temp_direct_routes.clear();
	}
}
void generate_routes_shortest(Network& network, vector<direct_route>& total_direct_routes) {
	for (int o = 0; o < network.terminals.size(); o++) {
		for (int d = o + 1; d < network.terminals.size(); d++) {
			int origin_node = network.terminals[o];
			int destination_node = network.terminals[d];

			vector<int> path = network.generate_shortest_route(origin_node, destination_node);
			direct_route route(path, network.get_travel_time(path), network.get_travel_time_r(path));

			total_direct_routes.push_back(route);
		}
	}

	connect_disconnected_nodes(network, total_direct_routes);
}

void recursive_construction(vector<direct_route>& route_set, vector<int> route, double travel_time, Network& network) {
	vector<int> input_route = route;
	vector<int> can_travel_to_nodes = network.can_travel_to_nodes(route.back());
	
	for (int i = 0; i < can_travel_to_nodes.size(); i++) {
		double edge_travel_time = network.get_travel_time(route.back(), can_travel_to_nodes[i]);

		route.push_back(can_travel_to_nodes[i]);
		travel_time += edge_travel_time;

		if (route.size() > 2) {
			int s = (int) route.size() - 1;
			
			bool rep = false;
			for (int j = 1; j < route.size() - 1; j++) {
				if (route.back() == route[j]) rep = true;
			}

			if (route[s] == route[s - 2]) {
				travel_time -= edge_travel_time;
				route.pop_back();

				direct_route new_route = direct_route(route, travel_time, network.get_travel_time_r(route));
				route_set.push_back(new_route);
			} else if (rep) {
				travel_time -= edge_travel_time;
				route.pop_back();

				direct_route new_route = direct_route(route, travel_time, network.get_travel_time_r(route));
				route_set.push_back(new_route);
			} else if (route[0] == route.back()) {
				direct_route new_route = direct_route(route, travel_time, network.get_travel_time_r(route));
				route_set.push_back(new_route);

				travel_time -= edge_travel_time;
				route.pop_back();
			}
		}

		if (input_route.size() != route.size()) {
			recursive_construction(route_set, route, travel_time, network);

			travel_time -= edge_travel_time;
			route.pop_back();
		}
	}
}

void route_prune(Network& network, vector<direct_route>& total_direct_routes, vector<direct_route>& temp_direct_routes) {
	for (int i = 0; i < temp_direct_routes.size(); i++) {
		if (temp_direct_routes[i].to_keep == true) {
			// Start and End at Terminal
			if (!temp_direct_routes[i].start_end_terminal(network.get_terminals())) {
				temp_direct_routes[i].to_keep = false;
			}

			// Remove Duplicate Cycles
			if (temp_direct_routes[i].circular) {
				for (int j = 0; j < total_direct_routes.size(); j++) {
					if (total_direct_routes[j].circular) {
						if (temp_direct_routes[i].route.size() == total_direct_routes[j].route.size()) {
							if (temp_direct_routes[i].compare_are_offset_circular(total_direct_routes[j].route)) temp_direct_routes[i].to_keep = false;
							if (temp_direct_routes[i].compare_are_offset_circular_reverse(total_direct_routes[j].route)) temp_direct_routes[i].to_keep = false;
						}
					}
				}
			}

			// Remove Duplicates and Reverse in New
			for (int j = i + 1; j < temp_direct_routes.size(); j++) {
				if (temp_direct_routes[i].compare_are_duplicate(temp_direct_routes[j].route)) temp_direct_routes[j].to_keep = false;
				if (temp_direct_routes[i].compare_are_reverse(temp_direct_routes[j].route)) temp_direct_routes[i].to_keep = false;
			}

			// Remove Duplicates and Reverse from New based on Master
			for (int j = 0; j < total_direct_routes.size(); j++) {
				if (temp_direct_routes[i].compare_are_reverse(total_direct_routes[j].route)) {
					temp_direct_routes[i].to_keep = false;
				}
			}

			// Remove Shorter than 3
			if (temp_direct_routes[i].route.size() < MIN_NODES_PER_ROUTE) {
				temp_direct_routes[i].to_keep = false;
			}

			// Remove Longer than 8
			if (temp_direct_routes[i].route.size() > MAX_NODES_PER_ROUTE) {
				temp_direct_routes[i].to_keep = false;
			}

			// Remove if longer than 50
			if (MAX_ROUTE_LENGTH > 0) {
				if (temp_direct_routes[i].length > MAX_ROUTE_LENGTH) {
					temp_direct_routes[i].to_keep = false;
				}
			}
		}
	}
}

void process_potential_transfers(vector<direct_route>& total_direct_routes) {
	for (int i = 0; i < total_direct_routes.size(); i++) {
		total_direct_routes[i].transfers.clear();
		total_direct_routes[i].rev_transfers.clear();

		total_direct_routes[i].overlap.clear();
		total_direct_routes[i].subset.clear();
		total_direct_routes[i].potential_transfers.clear();

		vector<int> overlap;
		vector<int> subset;
		vector<int> potential_transfers;

		for (int j = 0; j < total_direct_routes.size(); j++) {
			if (i != j) {
				if (total_direct_routes[i].overlaps(total_direct_routes[j].route)) {
					overlap.push_back(j);
					potential_transfers.push_back(j);
				}

				//if (total_direct_routes[i].subset_of(total_direct_routes[j].route) || total_direct_routes[i].reverse_subset_of(total_direct_routes[j].route)) {
				if (total_direct_routes[i].subset_of(total_direct_routes[j].route)) {
					subset.push_back(j);
					potential_transfers.pop_back();
				}
			}
		}

		total_direct_routes[i].overlap = overlap;
		total_direct_routes[i].subset = subset;
		total_direct_routes[i].potential_transfers = potential_transfers;
	}
}
void route_joining(vector<transfer_route>& total_transfer_routes, vector<direct_route>& total_direct_routes) {
	for (int i = 0; i < total_direct_routes.size(); i++) {
		vector<int> to_check = total_direct_routes[i].potential_transfers;

		for (int j = 0; j < to_check.size(); j++) {
			if (to_check[j] > i) {
				//if (to_check[j] == 5) {
				//	int temp = 100;
				//}

				vector<transfer_route> transfer_part;

				vector<int> r1 = total_direct_routes[i].route;
				vector<int> r2 = total_direct_routes[to_check[j]].route;

				vector<int> r1_;
				vector<int> r2_;

				// r1 Forward - r2 Forward
				r1_ = r1;
				r2_ = r2;
				intersect_routes(transfer_part, r1_, r2_, (i + 1), (to_check[j] + 1));

				// r2 Forward - r1 Forward
				r1_ = r2;
				r2_ = r1;
				intersect_routes(transfer_part, r1_, r2_, (to_check[j] + 1), (i + 1));

				// r1 Forward - r2 Reverse
				r1_ = r1;
				r2_ = r2;
				std::reverse(r2_.begin(), r2_.end());
				intersect_routes(transfer_part, r1_, r2_, (i + 1), -(to_check[j] + 1));

				// r1 Reverse - r2 Forward
				r1_ = r1;
				std::reverse(r1_.begin(), r1_.end());
				r2_ = r2;
				intersect_routes(transfer_part, r1_, r2_, -(i + 1), (to_check[j] + 1));

				remove_part(transfer_part, r1, r2);

				for (int k = 0; k < transfer_part.size(); k++) {
					if (transfer_part[k].to_keep) total_transfer_routes.push_back(transfer_part[k]);
				}

				transfer_part.clear();
			}
		}
	}
}
void intersect_routes(vector<transfer_route>& transfer_part, vector<int> r1, vector<int> r2, int r1_i, int r2_i) {
	vector<int> intersection_points;
	for (int i = 0; i < r1.size(); i++) {
		for (int j = 0; j < r2.size(); j++) {
			if (r1[i] == r2[j]) {
				bool insert = true;
				for (int k = 0; k < intersection_points.size(); k++) {
					if (intersection_points[k] == r1[i]) insert = false;
				}
				if (insert) intersection_points.push_back(r1[i]);
			}
		}
	}

	for (int i = 0; i < intersection_points.size(); i++) {
		vector<int> r1_index;
		for (int j = 0; j < r1.size(); j++) {
			if (r1[j] == intersection_points[i]) r1_index.push_back(j);
		}

		vector<int> r2_index;
		for (int j = 0; j < r2.size(); j++) {
			if (r2[j] == intersection_points[i]) r2_index.push_back(j);
		}

		for (int j = 0; j < r1_index.size(); j++) {
			for (int k = 0; k < r2_index.size(); k++) {
				vector<int> route;

				for (int l = 0; l < r1_index[j]; l++) {
					route.push_back(r1[l]);
				}
				int ji = route.size();
				for (int l = r2_index[k]; l < r2.size(); l++) {
					route.push_back(r2[l]);
				}

				//transfer_route new_route = transfer_route(route, r1_i, r2_i, r1_index[j] - i);
				transfer_route new_route = transfer_route(route, r1_i, r2_i, ji);
				transfer_part.push_back(new_route);
			}
		}
	}
}
void remove_part(vector<transfer_route>& transfer_part, vector<int> r1, vector<int> r2) {
	for (int i = 0; i < transfer_part.size(); i++) {
		//// Circular Routes
		//if (transfer_part[i].route[0] == transfer_part[i].route[transfer_part[i].route.size() - 1]) {
		//	transfer_part[i].to_keep = false;
		//}

		// Mirror Routes
		int half_index = (int) ceil(transfer_part[i].route.size());
		vector<int> first_half;
		vector<int> back_half;
		for (int j = 0; j < half_index; j++) {
			first_half.push_back(transfer_part[i].route[j]);
			back_half.push_back(transfer_part[i].route[(transfer_part[i].route.size() - 1) - j]);
		}
		bool t0 = true;
		for (int j = 0; j < first_half.size(); j++) {
			if (first_half[j] != back_half[j]) t0 = false;
		}
		if (t0) transfer_part[i].to_keep = false;

		// Route Shorting
		vector<int> unique_route_total;
		for (int j = 0; j < transfer_part[i].route.size(); j++) {
			bool insert = true;
			for (int k = 0; k < unique_route_total.size(); k++) {
				if (unique_route_total[k] == transfer_part[i].route[j]) insert = false;
			}
			if (insert) unique_route_total.push_back(transfer_part[i].route[j]);
		}
		vector<int> unique_route_start;
		for (int j = 0; j < transfer_part[i].route.size() - 1; j++) {
			bool insert = true;
			for (int k = 0; k < unique_route_start.size(); k++) {
				if (unique_route_start[k] == transfer_part[i].route[j]) insert = false;
			}
			if (insert) unique_route_start.push_back(transfer_part[i].route[j]);
		}
		vector<int> unique_route_end;
		for (int j = 1; j < transfer_part[i].route.size(); j++) {
			bool insert = true;
			for (int k = 0; k < unique_route_end.size(); k++) {
				if (unique_route_end[k] == transfer_part[i].route[j]) insert = false;
			}
			if (insert) unique_route_end.push_back(transfer_part[i].route[j]);
		}

		//bool t1 = unique_route_total.size() == transfer_part[i].route.size();
		bool t1 = unique_route_total.size() == transfer_part[i].route.size();
		bool t2 = unique_route_start.size() == (transfer_part[i].route.size() - 1);
		bool t3 = unique_route_end.size() == (transfer_part[i].route.size() - 1);

		//if (!t1 && (!t2 || !t3)) transfer_part[i].to_keep = false;
		//if (!t1 && (t2 || t3)) transfer_part[i].to_keep = false;
		//if (!t1) transfer_part[i].to_keep = false;
		if (!t2 && !t3) transfer_part[i].to_keep = false;

		// Copy/Reverse of Original
		if (transfer_part[i].compare_are_duplicate(r1)) transfer_part[i].to_keep = false;
		if (transfer_part[i].compare_are_reverse(r1)) transfer_part[i].to_keep = false;

		if (transfer_part[i].compare_are_duplicate(r2)) transfer_part[i].to_keep = false;
		if (transfer_part[i].compare_are_reverse(r2)) transfer_part[i].to_keep = false;

		// Subset of Originals
		if (transfer_part[i].subset_of(r1)) transfer_part[i].to_keep = false;
		if (transfer_part[i].subset_of(r2)) transfer_part[i].to_keep = false;
		if (transfer_part[i].reverse_subset_of(r1)) transfer_part[i].to_keep = false;
		if (transfer_part[i].reverse_subset_of(r2)) transfer_part[i].to_keep = false;

		//// Copy/Reverse of New Route
		//for (int j = i + 1; j < transfer_part.size(); j++) {
		//	if (transfer_part[i].compare_are_duplicate(transfer_part[j].route)) transfer_part[i].to_keep = false;
		//	if (transfer_part[i].compare_are_reverse(transfer_part[j].route)) transfer_part[i].to_keep = false;
		//}

		//// Subset of New Route
		//for (int j = i + 1; j < transfer_part.size(); j++) {
		//	if (transfer_part[i].subset_of(transfer_part[j].route)) transfer_part[i].to_keep = false;
		//	if (transfer_part[i].reverse_subset_of(transfer_part[j].route)) transfer_part[i].to_keep = false;
		//}
	}
}

void connect_disconnected_nodes(Network& network, vector<direct_route>& total_direct_routes) {
	vector<int> to_connect = network.get_unconnected(total_direct_routes);
	if (to_connect.size() > 0) {
		for (int i = 0; i < to_connect.size(); i++) {
			vector<vector<int>> paths_to_terminals = network.get_paths_to_terminals(to_connect[i]);

			for (int j = 0; j < paths_to_terminals.size(); j++) {
				for (int k = j + 1; k < paths_to_terminals.size(); k++) {
					vector<int> new_route;
					for (int l = 0; l < paths_to_terminals[j].size(); l++) new_route.push_back(paths_to_terminals[j][l]);
					std::reverse(new_route.begin(), new_route.end());
					for (int l = 1; l < paths_to_terminals[k].size(); l++) new_route.push_back(paths_to_terminals[k][l]);
					direct_route temp = direct_route(new_route, network.get_travel_time(new_route), network.get_travel_time_r(new_route));
					total_direct_routes.push_back(temp);
				}
			}
		}
	}
}
float optimise_problem(vector<direct_route>& total_direct_routes, vector<transfer_route>& total_transfer_routes, Network& network, bool output_stream_null) {
	int i, j;
	total_transfer_routes.clear();
	process_potential_transfers(total_direct_routes);
	route_joining(total_transfer_routes, total_direct_routes);

	// 02 - Solve base set
	Model model = Model(network, total_direct_routes, total_transfer_routes, BUSES, CAP);

	if (data_output_type == 100) {
		ofstream file;
		file.open(path_OUTPUT + "Timings.csv", std::ios_base::app);
		file << "Transfers: " << total_transfer_routes.size() << "; Edges: " << model.E_input << "; ";
		file.close();
	}	

	IloEnv env;

	auto start = chrono::high_resolution_clock::now();

	CPLEX_Model cplex_model = CPLEX_Model(env, model, MAX_ROUTES);

	auto end = chrono::high_resolution_clock::now();
	auto timetaken = chrono::duration_cast<chrono::milliseconds>(end - start);

	IloCplex cplex(cplex_model.cplex_model);
	if (output_stream_null) {
		cplex.setOut(env.getNullStream());
	}
	cplex.setParam(cplex.Threads, MAX_THREADS_CPLEX);
	cplex.setParam(cplex.NodeSel, 2);
	if (MAX_TIME_CPLEX > 0) cplex.setParam(cplex.TiLim, MAX_TIME_CPLEX);

	if (!cplex.solve()) {
		env.error() << "Failed to optimize LP" << endl;
		throw(-1);
	} else {
		//if (!output_stream_null) {
		//	output_cplex_solution(model, cplex, cplex_model);
		//}
	}

	// 03a - Setup for iterations
	for (i = 0; i < total_transfer_routes.size(); i++) {
		if (cplex.getValue(cplex_model.lambda_1_r[i]) == 1) {
			double passengers_assigned = 0;
			for (j = total_transfer_routes[i].edges[0]; j < total_transfer_routes[i].edges[1]; j++) {
				passengers_assigned += cplex.getValue(cplex_model.Gamma_e[j]);
				passengers_assigned += cplex.getValue(cplex_model.Gamma_e_dash[j]);
			}

			if (passengers_assigned > 0) {
				total_transfer_routes[i].used = true;
			}
		}
	}

	for (i = 0; i < total_direct_routes.size(); i++) {
		if (cplex.getValue(cplex_model.lambda_0_r[i]) > 0.9) {
			total_direct_routes[i].used = true;
		} else {
			total_direct_routes[i].used = false;
		}
	}

	std::replace(networkName.begin(), networkName.end(), '/', '=');
	export_output(total_direct_routes, total_transfer_routes, model, cplex, cplex_model, "");
	std::replace(networkName.begin(), networkName.end(), '=', '/');

	float output = (float)cplex.getObjValue();

	cplex.clearModel();
	cplex_model.cplex_model.end();
	cplex.end();

	return output;
}

void export_output(vector<direct_route>& total_direct_routes, vector<transfer_route>& total_transfer_routes, Model& model, IloCplex& cplex, CPLEX_Model& cplex_model, string input_path_part) {
	string path = path_OUTPUT + input_path_part;

	if (data_output_type == 200) {
		path = path_OUTPUT + input_path_part + " - ";
	} else {
		if (attempts < 9) path += "0" + to_string(attempts) + " - ";
		else path += "" + to_string(attempts) + " - ";
	}

	ofstream file;
	file.open(path + "solution_summary.csv");
	auto stop = chrono::high_resolution_clock().now();
	auto duration_milliseconds = chrono::duration_cast<chrono::milliseconds>(stop - main_start);
	auto duration_seconds = chrono::duration_cast<chrono::seconds>(stop - main_start);

	file << "Time Taken (seconds): ," << duration_seconds.count() << "." << std::setfill('0') << std::setw(3) << (duration_milliseconds.count() % 1000) << "\n\n";

	file << "Solve Number: ," << solves << "\n\n";

	file << "Nodes: ," << model.S_input << "\n\n";
	file << "Direct Routes: ," << model.R_0_input << "\n";
	file << "Transfer Routes: ," << model.R_1_input << "\n\n";

	file << "Buses: ," << model.buses << "\n";
	file << "Capacity: ," << model.cap << "\n";
	file << "Total Demand: ," << model.J << "\n";
	file << "Ratio: ," << ((double)model.J / ((double)model.buses * (double)cplex_model.time_units_per_day * (double)model.cap)) << "\n\n";

	file << "pi_0: ," << cplex.getValue(cplex_model.pi_0) << "\n";
	file << "pi_1: ," << cplex.getValue(cplex_model.pi_1) << "\n";
	file << "pi_2: ," << cplex.getValue(cplex_model.pi_2) << "\n" << "\n";

	double buses_required = 0;
	for (int r = 0; r < model.R_0_input; r++) {
		for (int f = 0; f <= model.F_input; f++) {
			buses_required += cplex.getValue(cplex_model.Lambda_r_f[r][f]) * (model.f_r_f[r][f] + model.f_r_f_dash[r][f]);
		}
	}

	file << "Buses: ," << buses_required << "\n\n";

	double obj_a = 0;

	for (int e = 0; e < model.E_input; e++) {
		obj_a += cplex.getValue(cplex_model.Gamma_e[e]) * model.t_e[e];
		obj_a += cplex.getValue(cplex_model.Gamma_e_dash[e]) * model.t_e_dash[e];
	}

	double obj_b = 0;
	obj_b += cplex.getValue(cplex_model.pi_1) * model.W * model.J;

	double obj_c = 0;
	for (int r = 0; r < model.R_0_input; r++) {
		for (int f = 0; f <= model.F_input; f++) {
			obj_c += cplex.getValue(cplex_model.Lambda_r_f[r][f]) * model.u_f[f] * (model.l_r[r] + model.l_r_dash[r]) * (double)cplex_model.time_units_per_day;
		}
	}

	double obj_d = 0;
	obj_d += cplex.getValue(cplex_model.pi_2) * model.J * 256 * 256;

	file << fixed << "Objective Function Value: ," << cplex.getObjValue() << "\n\n";

	file << fixed << "Passenger Cost (Travel Time): ," << obj_a << "\n";
	file << fixed << "Passenger Cost (Waiting Time): ," << obj_b << "\n";
	file << fixed << "Operator Cost: ," << obj_c << "\n";
	file << fixed << "Unsatisfied Demand Penalty: ," << obj_d << "\n";
	file.close();

	file.open(path + "routes_direct.csv");
	for (int r = 0; r < model.R_0_input; r++) {
		for (int s = 0; s < total_direct_routes[r].route.size(); s++) {
			file << total_direct_routes[r].route[s];
			file << (s == total_direct_routes[r].route.size() - 1 ? "" : ",");
		}
		file << "\n";
	}
	file.close();

	file.open(path + "routes_transfer.csv");
	for (int r = 0; r < model.R_1_input; r++) {
		for (int s = 0; s < total_transfer_routes[r].route.size(); s++) {
			file << total_transfer_routes[r].route[s];
			file << (s == total_transfer_routes[r].route.size() - 1 ? "" : ",");
		}
		file << "\n";
	}
	file.close();

	file.open(path + "unsatisfied_demand.csv");
	for (int s1 = 0; s1 < model.S_input; s1++) {
		for (int s2 = 0; s2 < model.S_input; s2++) {
			file << cplex.getValue(cplex_model.delta_s1_s2[s1][s2]);
			file << (s2 == model.S_input - 1 ? "" : ",");
		}
		file << "\n";
	}
	file.close();

	file.open(path + "usage_frequencies_direct.csv");
	for (int r = 0; r < model.R_0_input; r++) {
		file << cplex.getValue(cplex_model.lambda_0_r[r]);
		file << ",,";

		for (int f = 0; f <= model.F_input; f++) {
			file << cplex.getValue(cplex_model.Lambda_r_f[r][f]);
			file << (f == model.F_input ? "" : ",");
		}

		file << "\n";
	}
	file.close();

	file.open(path + "usage_transfer.csv");
	for (int r = 0; r < model.R_1_input; r++) {
		file << cplex.getValue(cplex_model.lambda_1_r[r]);
		file << "\n";
	}
	file.close();

	file.open(path + "Gamma_e.csv");
	for (int e = 0; e < model.E_input; e++) {
		file << model.edges_list[e].n << ",";
		file << model.edges_list[e].r << ",";
		file << model.edges_list[e].e << ",";
		file << model.edges_list[e].s1 << ",";
		file << model.edges_list[e].s2 << ",";
		file << model.edges_list[e].r1 << ",";
		file << model.edges_list[e].r2 << ",";
		file << cplex.getValue(cplex_model.Gamma_e[e]) << ",";
		file << cplex.getValue(cplex_model.Gamma_e_dash[e]);
		file << "\n";
	}
	file.close();

	file.open(path + "edge_travel_time.csv");
	for (int e = 0; e < model.E_input; e++) {
		file << model.t_e[e] << ",";
		file << model.t_e_dash[e];
		file << "\n";
	}
	file.close();

	file.open(path + "edges_used.csv");
	for (int p = 0; p < model.P_input; p++) {
		for (int e = 0; e < model.capacity_constraints[p].edge_list.size(); e++) {
			file << model.capacity_constraints[p].edge_list[e] << ",";
		}
		file << "\n";
	}
	file << "\n";
	for (int p = 0; p < model.P_input; p++) {
		for (int e = 0; e < model.rev_capacity_constraints[p].edge_list.size(); e++) {
			file << model.rev_capacity_constraints[p].edge_list[e] << ",";
		}
		file << "\n";
	}
	file.close();
}

void check_other_solutions(string routeset) {
	Network network;
	import_network(network);

	vector<direct_route> solution_to_test;
	import_routeset(network, solution_to_test, routeset);

	vector<transfer_route> transfer_routes;

	process_potential_transfers(solution_to_test);
	route_joining(transfer_routes, solution_to_test);

	Model model = Model(network, solution_to_test, transfer_routes, BUSES, CAP);
	IloEnv env;
	CPLEX_Model cplex_model = CPLEX_Model(env, model, MAX_ROUTES);

	IloCplex cplex(cplex_model.cplex_model);
	cplex.setOut(env.getNullStream());
	cplex.setParam(cplex.Threads, MAX_THREADS_CPLEX);
	cplex.setParam(cplex.NodeSel, 2);
	if (MAX_TIME_CPLEX > 0) cplex.setParam(cplex.TiLim, MAX_TIME_CPLEX);

	if (!cplex.solve()) {
		env.error() << "Failed to optimize LP" << endl;
		throw(-1);
	} else {
		export_output(solution_to_test, transfer_routes, model, cplex, cplex_model, routeset);
	}

	double obj_a = 0;
	for (int e = 0; e < model.E_input; e++) {
		obj_a += cplex.getValue(cplex_model.Gamma_e[e]) * model.t_e[e];
		obj_a += cplex.getValue(cplex_model.Gamma_e_dash[e]) * model.t_e_dash[e];
	}

	ofstream file;
	file.open(path_OUTPUT + "_Summary.csv", std::ios_base::app);
	file << routeset << ",";
	file << cplex.getValue(cplex_model.pi_0) << ",";
	file << cplex.getValue(cplex_model.pi_1) << ",";
	file << "-" << ",";
	file << cplex.getValue(cplex_model.pi_2) << ",";
	file << (obj_a + cplex.getValue(cplex_model.pi_1) * model.W * model.J) / (model.J * (cplex.getValue(cplex_model.pi_0) + cplex.getValue(cplex_model.pi_1))) << "\n"; // << "\n";
	file.close();

	float output = (float)cplex.getObjValue();
	cplex_model.cplex_model.end();
	cplex.end();
}

void import_routeset(Network& network, vector<direct_route>& solution_to_test, string routeset) {
	string full_path = path_INPUT + routeset + ".csv";

	if (data_output_type == 3) {
		full_path = routeset;
	}

	if (data_output_type == 200) {
		full_path = path_INPUT + "others_solutions/" + routeset + ".csv";
	}

	std::ifstream input_csv(full_path);
	string line;

	if (input_csv.is_open()) {
		while (input_csv) {
			vector<int> route;
			getline(input_csv, line);

			size_t start = 0;
			size_t end = line.find(",");
			while (end != -1) {
				try {
					route.push_back(stoi(line.substr(start, end - start)));
				} catch (exception) {
					route.push_back(stoi(line.substr(start + 3, end - (start + 3))));
				}
				start = end + 1;
				end = line.find(",", start);
			}
			route.push_back(stoi(line.substr(start, end - start)));

			solution_to_test.push_back(direct_route(route, network.get_travel_time(route), network.get_travel_time_r(route)));

			if (input_csv.eof()) {
				input_csv.close();
				break;
			}
		}
	}
}