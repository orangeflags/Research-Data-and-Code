#pragma once
class direct_route {

public:
	vector<int> route;
	double travel_time_f;
	double travel_time_r;
	int length;

	int edges[2];

	vector<int> subset;
	vector<int> overlap;

	vector<int> potential_transfers;
	vector<int> transfers;
	vector<int> rev_transfers;

	vector<vector<double>> demand_satisfied;

	bool to_keep = true;
	bool circular = false;

	bool used = false;

	direct_route(vector<int> input_route, double input_travel_time_f, double input_travel_time_r) {
		route = input_route;
		travel_time_f = input_travel_time_f;
		travel_time_r = input_travel_time_r;
		length = (int) route.size();

		if (route[0] == route.back()) circular = true;
	}
	
	bool start_end_terminal(vector<int> terminals) {
		bool start = false;
		bool end = false;

		for (int i = 0; i < terminals.size(); i++) {
			if (route[0] == terminals[i]) start = true;
			if (route.back() == terminals[i]) end = true;
		}

		return start && end;
	}
	bool compare_are_duplicate(vector<int> test_route) {
		if (route.size() != test_route.size()) return false;
		for (int i = 0; i < route.size(); i++) {
			if (route[i] != test_route[i]) return false;
		}
		return true;
	}
	bool compare_are_reverse(vector<int> test_route) {
		if (route.size() != test_route.size()) return false;

		std::reverse(test_route.begin(), test_route.end());
		for (int i = 0; i < route.size(); i++) {
			if (route[i] != test_route[i]) return false;
		}
		return true;
	}
	bool compare_are_offset_circular(vector<int> test_route) {
		test_route.pop_back();
		int offset = -1;
		for (int i = 0; i < test_route.size(); i++) {
			if (route[0] == test_route[i]) offset = i;
		}
		for (int i = 0; i < route.size() - 1; i++) {
			if (route[i] != test_route[(i + offset) % test_route.size()]) return false;
		}
		return true;
	}
	bool compare_are_offset_circular_reverse(vector<int> test_route) {
		std::reverse(test_route.begin(), test_route.end());
		test_route.pop_back();
		int offset = -1;
		for (int i = 0; i < test_route.size(); i++) {
			if (route[0] == test_route[i]) offset = i;
		}
		for (int i = 0; i < route.size() - 1; i++) {
			if (route[i] != test_route[(i + offset) % test_route.size()]) return false;
		}
		return true;
	}

	bool overlaps(vector<int> test_route) {
		for (int i = 0; i < route.size(); i++) {
			for (int j = 0; j < test_route.size(); j++) {
				if (route[i] == test_route[j]) return true;
			}
		}
		return false;
	}
	bool subset_of(vector<int> test_route) {
		if (test_route.size() < route.size()) return false;

		int to_check = (int) test_route.size() - (int) route.size() + 1;

		for (int i = 0; i < to_check; i++) {
			vector<int> sub_route;
			for (int j = 0; j < route.size(); j++) {
				sub_route.push_back(test_route[j + i]);
			}
			if (compare_are_duplicate(sub_route)) return true;
			sub_route.clear();
		}

		return false;
	}
	bool reverse_subset_of(vector<int> test_route) {
		std::reverse(test_route.begin(), test_route.end());

		if (test_route.size() < route.size()) return false;

		int to_check = (int) test_route.size() - (int) route.size() + 1;

		for (int i = 0; i < to_check; i++) {
			vector<int> sub_route;
			for (int j = 0; j < route.size(); j++) {
				sub_route.push_back(test_route[j + i]);
			}
			if (compare_are_duplicate(sub_route)) return true;
			sub_route.clear();
		}

		return false;
	}
};

class transfer_route {
public:
	vector<int> route;
	int length;
	int r1;
	int r2;
	int ji;

	int edges[2];

	bool to_keep = true;
	bool circular = false;
	double travel_time;

	bool used = false;

	transfer_route(vector<int> i_route, int i_r1, int i_r2, int i_ji) {
		route = i_route;
		r1 = i_r1;
		r2 = i_r2;
		ji = i_ji;
		length = (int)route.size();

		if (route[0] == route.back()) circular = true;
	}

	bool compare_are_duplicate(vector<int> test_route) {
		if (route.size() != test_route.size()) return false;
		for (int i = 0; i < route.size(); i++) {
			if (route[i] != test_route[i]) return false;
		}
		return true;
	}
	bool compare_are_reverse(vector<int> test_route) {
		if (route.size() != test_route.size()) return false;

		std::reverse(test_route.begin(), test_route.end());
		for (int i = 0; i < route.size(); i++) {
			if (route[i] != test_route[i]) return false;
		}
		return true;
	}
	bool compare_are_offset_circular(vector<int> test_route) {
		test_route.pop_back();
		int offset = -1;
		for (int i = 0; i < test_route.size(); i++) {
			if (route[0] == test_route[i]) offset = i;
		}
		for (int i = 0; i < route.size() - 1; i++) {
			if (route[i] != test_route[(i + offset) % test_route.size()]) return false;
		}
		return true;
	}
	bool compare_are_offset_circular_reverse(vector<int> test_route) {
		std::reverse(test_route.begin(), test_route.end());
		test_route.pop_back();
		int offset = -1;
		for (int i = 0; i < test_route.size(); i++) {
			if (route[0] == test_route[i]) offset = i;
		}
		for (int i = 0; i < route.size() - 1; i++) {
			if (route[i] != test_route[(i + offset) % test_route.size()]) return false;
		}
		return true;
	}
	bool subset_of(vector<int> test_route) {
		if (test_route.size() < route.size()) return false;

		int to_check = (int) test_route.size() - (int) route.size() + 1;

		for (int i = 0; i < to_check; i++) {
			vector<int> sub_route;
			for (int j = 0; j < route.size(); j++) {
				sub_route.push_back(test_route[j + i]);
			}
			if (compare_are_duplicate(sub_route)) return true;
			sub_route.clear();
		}

		return false;
	}
	bool reverse_subset_of(vector<int> test_route) {
		std::reverse(test_route.begin(), test_route.end());

		if (test_route.size() < route.size()) return false;

		int to_check = (int) test_route.size() - (int) route.size() + 1;

		for (int i = 0; i < to_check; i++) {
			vector<int> sub_route;
			for (int j = 0; j < route.size(); j++) {
				sub_route.push_back(test_route[j + i]);
			}
			if (compare_are_duplicate(sub_route)) return true;
			sub_route.clear();
		}

		return false;
	}
};