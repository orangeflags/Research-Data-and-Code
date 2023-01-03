#pragma once
class CPLEX_Model {
private:
	int A = 1024;

public:
	IloModel cplex_model;

	IloBoolVarArray lambda_0_r;
	IloBoolVarArray lambda_1_r;
	vector<IloBoolVarArray> Lambda_r_f;

	IloNumVar pi_0;
	IloNumVar pi_1;
	IloNumVar pi_2;

	IloNumVarArray Gamma_e;
	IloNumVarArray Gamma_e_dash;
	vector<IloNumVarArray> delta_s1_s2;

	IloBoolVarArray Beta_s;

	int time_units_per_day = 10;

	CPLEX_Model(IloEnv& env, Model& model, int in_max_routes) {
		if (in_max_routes != 0) A = in_max_routes;

		int r, f, e, p, s1, s2, n, s;
		string name;

		cplex_model = IloModel(env);

		// Decision Variables
		lambda_0_r = IloBoolVarArray(env);
		for (r = 0; r < model.R_0_input; r++) {
			name = "lambda_0_" + to_string(r);
			lambda_0_r.add(IloBoolVar(env, name.c_str()));
		}

		lambda_1_r = IloBoolVarArray(env);
		for (r = 0; r < model.R_1_input; r++) {
			name = "lambda_1_" + to_string(r);
			lambda_1_r.add(IloBoolVar(env, name.c_str()));
		}

		for (r = 0; r < model.R_0_input; r++) {
			Lambda_r_f.push_back(IloBoolVarArray(env));

			for (f = 0; f <= model.F_input; f++) {
				name = "Lambda_" + to_string(r) + "_" + to_string(f);
				Lambda_r_f[r].add(IloBoolVar(env, name.c_str()));
			}
		}

		pi_0 = IloNumVar(env, 0, 1, "pi_0");
		pi_1 = IloNumVar(env, 0, 1, "pi_1");
		pi_2 = IloNumVar(env, 0, 1, "pi_2");

		Gamma_e = IloNumVarArray(env);
		Gamma_e_dash = IloNumVarArray(env);
		for (e = 0; e < model.E_input; e++) {
			name = "Gamma_" + to_string(e);
			Gamma_e.add(IloNumVar(env, 0, model.M, name.c_str()));

			name = "Gamma_" + to_string(e) + "_dash";
			Gamma_e_dash.add(IloNumVar(env, 0, model.M, name.c_str()));
		}

		for (s1 = 0; s1 < model.S_input; s1++) {
			delta_s1_s2.push_back(IloNumVarArray(env));

			for (s2 = 0; s2 < model.S_input; s2++) {
				name = "delta_" + to_string(s1) + "_" + to_string(s2);
				delta_s1_s2[s1].add(IloNumVar(env, 0, model.M, name.c_str()));
			}
		}

		// Objective
		IloObjective objective(env, 0, IloObjective::Minimize, "obj");
		
		// Obj 1a
		objective.setLinearCoefs(Gamma_e, model.get_cplex_t_e(env));
		objective.setLinearCoefs(Gamma_e_dash, model.get_cplex_t_e_dash(env));

		// Obj 1b
		objective.setLinearCoef(pi_1, model.W * model.J);

		// Pen 1
		objective.setLinearCoef(pi_2, model.J * 256 * 256);

		cplex_model.add(objective);

		// Constraints
		IloRange ct_LimitedRoutes(env, 0, A, "ct_LimitedRoutes"); // Change the second input to this function between 0 and A to change between an exact number of solutions (A) and a maximum (0).
		for (r = 0; r < model.R_0_input; r++) {
			ct_LimitedRoutes.setLinearCoef(lambda_0_r[r], 1);
		}
		cplex_model.add(ct_LimitedRoutes);

		IloRangeArray ct_Frequency1(env);
		for (r = 0; r < model.R_0_input; r++) {
			name = "ct_Frequency1_" + to_string(r);
			IloRange temp(env, 1, 1, name.c_str());
			for (f = 0; f <= model.F_input; f++) {
				temp.setLinearCoef(Lambda_r_f[r][f], 1);
			}
			ct_Frequency1.add(temp);
		}
		cplex_model.add(ct_Frequency1);

		IloRangeArray ct_Frequency2(env);
		for (r = 0; r < model.R_0_input; r++) {
			name = "ct_Frequency2_" + to_string(r);
			IloRange temp(env, 0, 0, name.c_str());
			temp.setLinearCoef(lambda_0_r[r], -1);

			for (f = 1; f <= model.F_input; f++) {
				temp.setLinearCoef(Lambda_r_f[r][f], 1);
			}
			ct_Frequency2.add(temp);
		}
		cplex_model.add(ct_Frequency2);

		IloConstraintArray ct_RouteUsage1(env);
		for (r = 0; r < model.R_1_input; r++) {
			ct_RouteUsage1.add(lambda_1_r[r] >= lambda_0_r[model.q_n_rn_r[r][0]] + lambda_0_r[model.q_n_rn_r[r][1]] - 1);
			ct_RouteUsage1[ct_RouteUsage1.getSize() - 1].setName(("ct_RouteUsage1_" + to_string(r)).c_str());
		}
		cplex_model.add(ct_RouteUsage1);

		IloConstraintArray ct_RouteUsage2(env);
		for (r = 0; r < model.R_1_input; r++) {
			ct_RouteUsage2.add(lambda_1_r[r] <= (lambda_0_r[model.q_n_rn_r[r][0]] + lambda_0_r[model.q_n_rn_r[r][1]]) / 2);
			ct_RouteUsage2[ct_RouteUsage2.getSize() - 1].setName(("ct_RouteUsage2_" + to_string(r)).c_str());
		}
		cplex_model.add(ct_RouteUsage2);

		IloRangeArray ct_DemandAssignment(env);
		for (s1 = 0; s1 < model.S_input; s1++) {
			for (s2 = 0; s2 < model.S_input; s2++) {
				name = "ct_DemandAssignment_" + to_string(s1) + "_" + to_string(s2);
				IloRange temp(env, model.d_s1_s2[s1][s2], model.d_s1_s2[s1][s2], name.c_str());
				temp.setLinearCoefs(Gamma_e, model.get_cplex_o_s1_s2(env, s1, s2));
				temp.setLinearCoefs(Gamma_e_dash, model.get_cplex_o_s1_s2(env, s2, s1));
				temp.setLinearCoef(delta_s1_s2[s1][s2], 1);
				ct_DemandAssignment.add(temp);
			}
		}
		cplex_model.add(ct_DemandAssignment);

		IloRange ct_DirectTravelProportion(env, 0, 0, "ct_DirectTravelProportion");
		ct_DirectTravelProportion.setLinearCoef(pi_0, -1);
		ct_DirectTravelProportion.setLinearCoefs(Gamma_e, model.get_cplex_o_e(env));
		ct_DirectTravelProportion.setLinearCoefs(Gamma_e_dash, model.get_cplex_o_e(env));
		cplex_model.add(ct_DirectTravelProportion);

		IloRange ct_TransferTravelProportion(env, 1, 1, "ct_TransferTravelProportion");
		ct_TransferTravelProportion.setLinearCoef(pi_0, 1);
		ct_TransferTravelProportion.setLinearCoef(pi_1, 1);
		ct_TransferTravelProportion.setLinearCoef(pi_2, 1);
		cplex_model.add(ct_TransferTravelProportion);

		IloRange ct_UnsatisfiedTravelProportion(env, 0, 0, "ct_UnsatisfiedTravelProportion");
		ct_UnsatisfiedTravelProportion.setLinearCoef(pi_2, -1);
		for (s1 = 0; s1 < model.S_input; s1++) {
			for (s2 = 0; s2 < model.S_input; s2++) {
				ct_UnsatisfiedTravelProportion.setLinearCoef(delta_s1_s2[s1][s2], 1 / (double)model.J);
			}
		}
		cplex_model.add(ct_UnsatisfiedTravelProportion);

		IloRange ct_FleetSize(env, 0, model.B, "ct_FleetSize");
		for (r = 0; r < model.R_0_input; r++) {
			for (f = 0; f <= model.F_input; f++) {
				ct_FleetSize.setLinearCoef(Lambda_r_f[r][f], model.f_r_f[r][f] + model.f_r_f_dash[r][f]);
			}
		}
		cplex_model.add(ct_FleetSize);

		IloRangeArray ct_BusCapacityForward(env);
		for (p = 0; p < model.P_input; p++) {
			name = "ct_BusCapacityForward_" + to_string(p);
			IloRange temp(env, -model.M, 0, name.c_str());
			temp.setLinearCoefs(Gamma_e, model.get_cplex_x_p(env, p));
			for (f = 0; f <= model.F_input; f++) {
				temp.setLinearCoef(Lambda_r_f[model.y_p[p]][f], -1 * model.u_f[f] * model.C * time_units_per_day);
			}
			ct_BusCapacityForward.add(temp);
		}
		cplex_model.add(ct_BusCapacityForward);

		IloRangeArray ct_BusCapacityReverse(env);
		for (p = 0; p < model.P_input; p++) {
			name = "ct_BusCapacityReverse_" + to_string(p);
			IloRange temp(env, -model.M, 0, name.c_str());
			temp.setLinearCoefs(Gamma_e_dash, model.get_cplex_x_p_dash(env, p));
			for (f = 0; f <= model.F_input; f++) {
				temp.setLinearCoef(Lambda_r_f[model.y_p[p]][f], -1 * model.u_f[f] * model.C * time_units_per_day);
			}
			ct_BusCapacityReverse.add(temp);
		}
		cplex_model.add(ct_BusCapacityReverse);

		IloConstraintArray ct_RouteLinkUsage(env);
		for (n = 0; n <= model.N_input; n++) {
			int r_max = n == 0 ? model.R_0_input : model.R_1_input;

			for (r = 0; r < r_max; r++) {
				for (e = model.e_n_r[n][r][0]; e <= model.e_n_r[n][r][1]; e++) {
					if (n == 0) {
						ct_RouteLinkUsage.add(Gamma_e[e] <= model.M * lambda_0_r[r]);
						ct_RouteLinkUsage[ct_RouteLinkUsage.getSize() - 1].setName(("ct_RouteLinkUsage_Forward_" + to_string(e)).c_str());

						ct_RouteLinkUsage.add(Gamma_e_dash[e] <= model.M * lambda_0_r[r]);
						ct_RouteLinkUsage[ct_RouteLinkUsage.getSize() - 1].setName(("ct_RouteLinkUsage_Reverse_" + to_string(e)).c_str());
					} else {
						ct_RouteLinkUsage.add(Gamma_e[e] <= model.M * lambda_1_r[r]);
						ct_RouteLinkUsage[ct_RouteLinkUsage.getSize() - 1].setName(("ct_RouteLinkUsage_Forward_" + to_string(e)).c_str());

						ct_RouteLinkUsage.add(Gamma_e_dash[e] <= model.M * lambda_1_r[r]);
						ct_RouteLinkUsage[ct_RouteLinkUsage.getSize() - 1].setName(("ct_RouteLinkUsage_Reverse_" + to_string(e)).c_str());
					}
				}
			}
		}
		cplex_model.add(ct_RouteLinkUsage);
	}
};