#pragma once
using namespace std;

class Model {
public:
	// More care can be put into these
	int cap = 50;
	int buses = 20;

	// Pre-processing Parameters
	int N_input = 1;

	// Scalar Parameters
	int B;
	int C;
	int J;
	int M = 1024 * 1024 * 2;
	int W = 5;

	int E_input;
	int F_input = 6;
	int P_input;
	int S_input;
	int R_0_input;
	int R_1_input;

	// Vector Parameters
	vector<double> l_r;
	vector<double> l_r_dash;
	vector<double> t_e;
	vector<double> t_e_dash;
	vector<int> u_f;
	vector<int> w_f;
	vector<int> y_p;

	// Matrix Parameters
	vector<vector<double>> f_r_f;
	vector<vector<double>> f_r_f_dash;
	vector<vector<int>> d_s1_s2;
	vector<vector<int>> z_s1_s2;
	vector<vector<int>> q_n_rn_r;
	vector<vector<int>> v_s_r;

	vector<vector<vector<int>>> e_n_r;

	// Tuple Parameters
	vector<vector<vector<vector<int>>>> o_n_s1_s2;

	// Working Variables
	struct Edge {
	public:
		int n;
		int r;
		int e;
		int s1;
		int s2;
		int r1;
		int r2;

		Edge(int _n, int _r, int _e, int _s1, int _s2, int _r1, int _r2) {
			n = _n;
			r = _r;
			e = _e;
			s1 = _s1;
			s2 = _s2;
			r1 = _r1;
			r2 = _r2;
		}
	};
	struct CapacityConstraint {
	public:
		int p;
		int r;
		vector<int> edge_list;

		CapacityConstraint(int _p, int _r, vector<int> _edge_list) {
			p = _p;
			r = _r;
			edge_list = _edge_list;
		}
	};

	vector<Edge> edges_list; // E_n_r
	vector<CapacityConstraint> capacity_constraints; // x_p
	vector<CapacityConstraint> rev_capacity_constraints; // x_p_dash

	// Functions
	void generate_edge_list_part(vector<int> route, int& insert_index, int r1, int r2, int n, int r) {
		if (route[0] == route.back()) {
			route.pop_back();

			for (int i = 0; i < route.size(); i++) {
				for (int j = 0; j < route.size(); j++) {
					if (i != j) {
						Edge edge = Edge(n, r, insert_index, route[i], route[j], r1, r2);
						edges_list.push_back(edge);
						insert_index++;
					}
				}
			}
		} else {
			for (int i = 0; i < route.size(); i++) {
				for (int j = i; j < route.size(); j++) {
					if (i != j) {
						Edge edge = Edge(n, r, insert_index, route[i], route[j], r1, r2);
						edges_list.push_back(edge);
						insert_index++;
					}
				}
			}
		}
	}
	void create_edges_list(vector<direct_route>& total_direct_routes, vector<transfer_route>& total_transfer_routes) {
		int insert_index = 0;

		for (int r = 0; r < total_direct_routes.size(); r++) {
			vector<int> route = total_direct_routes[r].route;
			int start_index = insert_index;

			generate_edge_list_part(route, insert_index, r + 1, -(r + 1), 0, r);

			total_direct_routes[r].edges[0] = start_index;
			total_direct_routes[r].edges[1] = insert_index - 1;
		}

		for (int r = 0; r < total_transfer_routes.size(); r++) {
			vector<int> route = total_transfer_routes[r].route;
			vector<int> edge;

			int start_index = insert_index;

			int r1 = abs(total_transfer_routes[r].r1) - 1;
			int r2 = abs(total_transfer_routes[r].r2) - 1;

			if (total_transfer_routes[r].r1 > 0) total_direct_routes[r1].transfers.push_back(r);
			else total_direct_routes[r1].rev_transfers.push_back(r);

			if (total_transfer_routes[r].r2 > 0) total_direct_routes[r2].transfers.push_back(r);
			else total_direct_routes[r2].rev_transfers.push_back(r);

			generate_edge_list_part(route, insert_index, total_transfer_routes[r].r1, total_transfer_routes[r].r2, 1, r);

			total_transfer_routes[r].edges[0] = start_index;
			total_transfer_routes[r].edges[1] = insert_index - 1;
		}
	}
	void get_links(vector<int>& final_needed, vector<int> route, int n, vector<Edge> route_edges) {
		if (route[0] == route[route.size() - 1]) {
			route.pop_back();

			int route_split = -1;
			for (int i = 0; i < route.size(); i++) if (route[i] == n) route_split = i;
			for (int i = 0; i < route_split; i++) {
				route.push_back(route[0]);
				route.erase(route.begin());
			}

			if (route_split == -1) return;

			for (int i = 1; i < route.size(); i++) {
				for (int j = 0; j < route_edges.size(); j++) if (route_edges[j].s1 == route[0] && route_edges[j].s2 == route[i]) final_needed.push_back(route_edges[j].e);
			}

			for (int i = (int) route.size() - 1; i > 0; i--) {
				for (int j = 1; j < i; j++) {
					for (int k = 0; k < route_edges.size(); k++) if (route_edges[k].s1 == route[j] && route_edges[k].s2 == route[i]) final_needed.push_back(route_edges[k].e);
				}
			}
		} else {
			int route_split = -1;
			for (int i = 0; i < route.size(); i++) if (route[i] == n) route_split = i;
			
			vector<int> route_start;
			for (int i = 0; i <= route_split; i++) route_start.push_back(route[i]);

			vector<int> route_end;
			for (int i = route_split + 1; i < route.size(); i++)  route_end.push_back(route[i]);

			for (int i = 0; i < route_start.size(); i++) {
				for (int j = 0; j < route_end.size(); j++) {
					for (int k = 0; k < route_edges.size(); k++) if (route_edges[k].s1 == route_start[i] && route_edges[k].s2 == route_end[j]) final_needed.push_back(route_edges[k].e);
				}
			}
		}
	}
	void get_links_transfer(vector<int>& final_needed, vector<int> route, int n, vector<Edge> route_edges, int transfer_point, int r1, int re) {
		//if (route[0] == route[route.size() - 1]) {
		//	//throw(-1);
		//} else {
			int eval_point = -1;
			for (int i = 0; i < route.size(); i++) if (route[i] == n) eval_point = i;

			vector<int> route_start;
			for (int i = 0; i <= eval_point; i++) route_start.push_back(route[i]);

			vector<int> route_end;
			for (int i = eval_point + 1; i < route.size(); i++)  route_end.push_back(route[i]);

			if (abs(r1) == abs(re)) {
				if (eval_point >= transfer_point) route_end.clear();
			} else {
				if (eval_point < transfer_point) route_end.clear();
			}

			for (int i = 0; i < route_start.size(); i++) {
				for (int j = 0; j < route_end.size(); j++) {
					for (int k = 0; k < route_edges.size(); k++) if (route_edges[k].s1 == route_start[i] && route_edges[k].s2 == route_end[j]) final_needed.push_back(route_edges[k].e);
				}
			}
		//}
	}
	void swap_edges(vector<Edge>& route_edges) {
		for (int i = 0; i < route_edges.size(); i++) {
			int temp = route_edges[i].s1;
			route_edges[i].s1 = route_edges[i].s2;
			route_edges[i].s2 = temp;
		}
	}
	void create_capacity_constraint(vector<direct_route>& total_direct_routes, vector<transfer_route>& total_transfer_routes) {
		int capacity_constratints_identified = -1;

		for (int r = 0; r < total_direct_routes.size(); r++) {
			vector<int> route = total_direct_routes[r].route;
			vector<Edge> route_edges;
			for (int i = total_direct_routes[r].edges[0]; i <= total_direct_routes[r].edges[1]; i++) {
				route_edges.push_back(edges_list[i]);
			}

			for (int n = 0; n < route.size(); n++) {
				capacity_constratints_identified += 1;
				vector<int> final_needed;

				bool process = true;
				for (int i = 0; i < n; i++) if (route[i] == route[n]) process = false;

				if (process) {
					get_links(final_needed, route, route[n], route_edges);

					for (int rt = 0; rt < total_direct_routes[r].transfers.size(); rt++) {
						int rt_index = total_direct_routes[r].transfers[rt];

						vector<int> route_transfer = total_transfer_routes[rt_index].route;
						vector<Edge> route_transfer_edges;
						for (int i = total_transfer_routes[rt_index].edges[0]; i <= total_transfer_routes[rt_index].edges[1]; i++) {
							route_transfer_edges.push_back(edges_list[i]);
						}

						get_links_transfer(final_needed, route_transfer, route[n], route_transfer_edges, total_transfer_routes[rt_index].ji, total_transfer_routes[rt_index].r1, r + 1);
					}

					for (int rt = 0; rt < total_direct_routes[r].rev_transfers.size(); rt++) {
						int rt_index = total_direct_routes[r].rev_transfers[rt];

						vector<int> route_transfer = total_transfer_routes[rt_index].route;
						std::reverse(route_transfer.begin(), route_transfer.end());

						vector<Edge> route_transfer_edges;
						for (int i = total_transfer_routes[rt_index].edges[0]; i <= total_transfer_routes[rt_index].edges[1]; i++) {
							route_transfer_edges.push_back(edges_list[i]);
						}
						swap_edges(route_transfer_edges);

						get_links_transfer(final_needed, route_transfer, route[n], route_transfer_edges, route_transfer.size() - (1 + total_transfer_routes[rt_index].ji), total_transfer_routes[rt_index].r2, r + 1);
					}
				}

				CapacityConstraint capacityConstraint = CapacityConstraint((int) capacity_constratints_identified, (int) r, final_needed);
				capacity_constraints.push_back(capacityConstraint);
			}
		}
	}
	void create_rev_capacity_constraint(vector<direct_route>& total_direct_routes, vector<transfer_route>& total_transfer_routes) {
		int rev_capacity_constratints_identified = -1;

		for (int r = 0; r < total_direct_routes.size(); r++) {
			vector<int> route = total_direct_routes[r].route;
			std::reverse(route.begin(), route.end());

			vector<Edge> route_edges;
			for (int i = total_direct_routes[r].edges[0]; i <= total_direct_routes[r].edges[1]; i++) {
				route_edges.push_back(edges_list[i]);
			}
			swap_edges(route_edges);

			for (int n = 0; n < route.size(); n++) {
				rev_capacity_constratints_identified += 1;
				vector<int> final_needed;

				bool process = true;
				for (int i = 0; i < n; i++) if (route[i] == route[n]) process = false;

				if (process) {
					get_links(final_needed, route, route[n], route_edges);

					for (int rt = 0; rt < total_direct_routes[r].rev_transfers.size(); rt++) {
						int rt_index = total_direct_routes[r].rev_transfers[rt];

						vector<int> route_transfer = total_transfer_routes[rt_index].route;
						vector<Edge> route_transfer_edges;
						for (int i = total_transfer_routes[rt_index].edges[0]; i <= total_transfer_routes[rt_index].edges[1]; i++) {
							route_transfer_edges.push_back(edges_list[i]);
						}

						get_links_transfer(final_needed, route_transfer, route[n], route_transfer_edges, total_transfer_routes[rt_index].ji, total_transfer_routes[rt_index].r1, r + 1);
					}

					for (int rt = 0; rt < total_direct_routes[r].transfers.size(); rt++) {
						int rt_index = total_direct_routes[r].transfers[rt];

						vector<int> route_transfer = total_transfer_routes[rt_index].route;
						std::reverse(route_transfer.begin(), route_transfer.end());

						vector<Edge> route_transfer_edges;
						for (int i = total_transfer_routes[rt_index].edges[0]; i <= total_transfer_routes[rt_index].edges[1]; i++) {
							route_transfer_edges.push_back(edges_list[i]);
						}
						swap_edges(route_transfer_edges);

						get_links_transfer(final_needed, route_transfer, route[n], route_transfer_edges, route_transfer.size() - (1 + total_transfer_routes[rt_index].ji), total_transfer_routes[rt_index].r2, r + 1);
					}
				}

				CapacityConstraint capacityConstraint = CapacityConstraint((int) rev_capacity_constratints_identified, (int) r, final_needed);
				rev_capacity_constraints.push_back(capacityConstraint);
			}
		}
	}

	Model(Network& network, vector<direct_route>& total_direct_routes, vector<transfer_route>& total_transfer_routes, int in_buses, int in_cap) {
		if (in_buses != 0) buses = in_buses;
		if (in_cap != 0) cap = in_cap;

		B = buses;

		// __ E_input __ P_input __
		E_input = 0;
		P_input = 0;
		for (int i = 0; i < total_direct_routes.size(); i++) {
			int len = (int) total_direct_routes[i].route.size();
			P_input += len;
			E_input += total_direct_routes[i].route[0] == total_direct_routes[i].route.back() ? (len - 2) * (len - 1) : (len * (len - 1)) / 2;
		}
		for (int i = 0; i < total_transfer_routes.size(); i++) {
			int len = (int) total_transfer_routes[i].route.size();
			E_input += total_transfer_routes[i].route[0] == total_transfer_routes[i].route.back() ? (len - 2) * (len - 1) : (len * (len - 1)) / 2;
		}

		// Edge List
		create_edges_list(total_direct_routes, total_transfer_routes);

		// Capacity Constraints
		create_capacity_constraint(total_direct_routes, total_transfer_routes);
		create_rev_capacity_constraint(total_direct_routes, total_transfer_routes);

		// __ J __ d_s1_s2 __
		for (int i = 0; i < network.num_nodes; i++) {
			vector<int> row;
			for (int j = 0; j < network.num_nodes; j++) {
				row.push_back(0);
			}
			d_s1_s2.push_back(row);
			z_s1_s2.push_back(row);
		}
		for (int i = 0; i < network.demand.size(); i++) {
			d_s1_s2[network.demand[i].n1][network.demand[i].n2] = network.demand[i].d;
			z_s1_s2[network.demand[i].n1][network.demand[i].n2] = network.demand[i].d > 0 ? 1 : 0;
			J += network.demand[i].d;
		}

		// __ y_p __
		for (int p = 0; p < P_input; p++) {
			//y_p.push_back(edges_list[p].r);
			y_p.push_back(capacity_constraints[p].r);
		}

		// __ S_input __
		S_input = network.num_nodes;

		// __ f_f __ w_f __
		w_f.push_back(M);
		w_f.push_back(60);
		w_f.push_back(30);
		w_f.push_back(20);
		w_f.push_back(15);
		w_f.push_back(10);
		w_f.push_back(5);

		u_f.push_back(0);
		u_f.push_back(1);
		u_f.push_back(2);
		u_f.push_back(3);
		u_f.push_back(4);
		u_f.push_back(6);
		u_f.push_back(12);

		// __ R_0_input __ c_r __ l_r __ u_r __ v_r __
		R_0_input = (int) total_direct_routes.size();
		for (int i = 0; i < R_0_input; i++) {
			l_r.push_back(total_direct_routes[i].travel_time_f);
			l_r_dash.push_back(total_direct_routes[i].travel_time_r);

			vector<double> bus_req;
			for (int j = 0; j < u_f.size(); j++) {
				bus_req.push_back(u_f[j] * l_r[i] / 60.0);
			}
			f_r_f.push_back(bus_req);

			vector<double> bus_req_dash;
			for (int j = 0; j < u_f.size(); j++) {
				bus_req_dash.push_back(u_f[j] * l_r_dash[i] / 60.0);
			}
			f_r_f_dash.push_back(bus_req_dash);
		}
		C = cap;

		// __ R_1_input __
		R_1_input = (int) total_transfer_routes.size();

		// __ q_n_rn_r __
		for (int i = 0; i < total_transfer_routes.size(); i++) {
			vector<int> temp;
			temp.push_back(abs(total_transfer_routes[i].r1) - 1);
			temp.push_back(abs(total_transfer_routes[i].r2) - 1);
			q_n_rn_r.push_back(temp);
		}

		// __ i_n_r __
		vector<vector<int>> n;
		for (int r = 0; r < R_0_input; r++) {
			vector<int> row;
			row.push_back(total_direct_routes[r].edges[0]);
			row.push_back(total_direct_routes[r].edges[1]);
			n.push_back(row);
		}
		e_n_r.push_back(n);
		n.clear();
		for (int r = 0; r < R_1_input; r++) {
			vector<int> row;
			row.push_back(total_transfer_routes[r].edges[0]);
			row.push_back(total_transfer_routes[r].edges[1]);
			n.push_back(row);
		}
		e_n_r.push_back(n);

		// __ t_e __ t_e_dash __
		for (int e = 0; e < E_input; e++) {
			vector<int> route;
			if (edges_list[e].n == 0) route = total_direct_routes[edges_list[e].r].route;
			else route = total_transfer_routes[edges_list[e].r].route;

			if (route[0] == route.back()) {
				vector<int> temp_route;
				for (int i = 0; i < 2; i++) {
					for (int j = 0; j < route.size() - 1; j++) {
						temp_route.push_back(route[j]);
					}
				}
				route = temp_route;
			}

			int origin = -1;
			for (int i = 0; i < route.size(); i++) if (origin == -1 && route[i] == edges_list[e].s1) origin = i;

			int destination = -1;
			for (int i = origin; i < route.size(); i++) if (destination == -1 && route[i] == edges_list[e].s2) destination = i;

			double travel_time = 0;
			for (int i = origin; i < destination; i++) {
				travel_time += network.get_travel_time(route[i], route[i + 1]);
			}

			t_e.push_back(travel_time);
		}

		for (int e = 0; e < E_input; e++) {
			vector<int> route;
			if (edges_list[e].n == 0) route = total_direct_routes[edges_list[e].r].route;
			else route = total_transfer_routes[edges_list[e].r].route;
			std::reverse(route.begin(), route.end());


			if (route[0] == route.back()) {
				vector<int> temp_route;
				for (int i = 0; i < 2; i++) {
					for (int j = 0; j < route.size() - 1; j++) {
						temp_route.push_back(route[j]);
					}
				}
				route = temp_route;
			}

			int origin = -1;
			for (int i = 0; i < route.size(); i++) if (origin == -1 && route[i] == edges_list[e].s2) origin = i;

			int destination = -1;
			for (int i = origin; i < route.size(); i++) if (destination == -1 && route[i] == edges_list[e].s1) destination = i;

			double travel_time = 0;
			for (int i = origin; i < destination; i++) {
				travel_time += network.get_travel_time(route[i], route[i + 1]);
			}

			t_e_dash.push_back(travel_time);
		}
	
		// __ o_n_s1_s2 __
		for (int n = 0; n < 2; n++) {
			vector<vector<vector<int>>> _n;

			for (int s1 = 0; s1 < S_input; s1++) {
				vector<vector<int>> _s1;

				for (int s2 = 0; s2 < S_input; s2++) {
					vector<int> _s2;

					for (int e = 0; e < E_input; e++) {
						if (edges_list[e].s1 == s1 && edges_list[e].s2 == s2 && edges_list[e].n == n) _s2.push_back(e);
					}

					_s1.push_back(_s2);
				}

				_n.push_back(_s1);
			}

			o_n_s1_s2.push_back(_n);
		}

		//// __ v_s_r __
		//for (int s = 0; s < S_input; s++) {
		//	vector<int> v_s;

		//	for (int r = 0; r < R_0_input; r++) {
		//		v_s.push_back(0);
		//	}
		//	v_s_r.push_back(v_s);
		//}
		//for (int r = 0; r < R_0_input; r++) {
		//	for (int n = 0; n < total_direct_routes[r].route.size(); n++) {
		//		v_s_r[total_direct_routes[r].route[n]][r] = 1;
		//	}
		//}

	}

	IloNumArray get_cplex_t_e(IloEnv& env) {
		IloNumArray cplex_t_e = IloNumArray(env, E_input);

		for (int e = 0; e < E_input; e++) {
			cplex_t_e[e] = t_e[e];
		}

		return cplex_t_e;
	}
	IloNumArray get_cplex_t_e_dash(IloEnv& env) {
		IloNumArray cplex_t_e_dash = IloNumArray(env, E_input);

		for (int e = 0; e < E_input; e++) {
			cplex_t_e_dash[e] = t_e_dash[e];
		}

		return cplex_t_e_dash;
	}
	IloNumArray get_cplex_o_s1_s2(IloEnv& env, int s1, int s2) {
		IloNumArray output = IloNumArray(env, E_input);

		for (int n = 0; n <= 1; n++) {
			for (int e = 0; e < o_n_s1_s2[n][s1][s2].size(); e++) {
				output[o_n_s1_s2[n][s1][s2][e]] = 1;
			}
		}

		return output;
	}
	IloNumArray get_cplex_o_e(IloEnv& env) {
		IloNumArray cplex_o_e = IloNumArray(env, E_input);

		for (int e = 0; e < E_input; e++) {
			if (edges_list[e].n == 0) cplex_o_e[e] = 1 / (double)J;
		}

		return cplex_o_e;
	}
	IloNumArray get_cplex_x_p(IloEnv& env, int p) {
		IloNumArray p_ = IloNumArray(env, E_input);

		for (int e = 0; e < capacity_constraints[p].edge_list.size(); e++) {
			p_[capacity_constraints[p].edge_list[e]] = 1;
		}

		return p_;
	}
	IloNumArray get_cplex_x_p_dash(IloEnv& env, int p) {
		IloNumArray p_ = IloNumArray(env, E_input);

		for (int e = 0; e < rev_capacity_constraints[p].edge_list.size(); e++) {
			p_[rev_capacity_constraints[p].edge_list[e]] = 1;
		}

		return p_;
	}
};