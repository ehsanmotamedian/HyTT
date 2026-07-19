function [results, best_sol] = solveHyTT(model, glucose_rate, params)
% solveHyTT - Optimizes the HyTT model using binary search for growth rate and calculates proteome reallocation.
%
% USAGE:
%   [results, best_sol] = solveHyTT(model, glucose_rate, params)
%
% INPUTS:
%   model               The integrated HyTT MILP model structure.
%   glucose_rate        Specific glucose uptake rate limit (mmol/gDW/h).
%
% OPTIONAL INPUTS:
%   params              Struct containing simulation and kinetic parameters:
%                       * .sigma          - Global saturation factor (Default: 0.70)
%                       * .kappa_t        - Effective translational capacity (Default: 1.3 h^-1)
%                       * .k_trans        - Translation elongation rate (Default: 28800 aa/h)
%                       * .phi_0          - Baseline ribosomal fraction (Default: 0.05)
%                       * .f_pool         - Active proteome fraction (Default: 0.71)
%                       * .total_protein  - Total cell protein content (Default: 460 mg/gDW)
%                       * .k_deg          - Protein degradation rate (Default: 0.05 1/h)
%                       * .k_deg_RNA      - mRNA degradation rate (Default: 2.1 1/h)
%                       * .f_active       - Active ribosomal fraction (Default: 0.85)
%                       * .MW_ribosome    - Ribosome molecular weight (Default: 3.0e6 Da)
%                       * .RNAP_total     - total RNA polymerase capacity (Default: 0.01 g/gDCW)
%
% OUTPUTS:
%   results             Struct containing predicted macroscopic fluxes, growth rate, and proteome fractions.
%   best_sol            The complete optimization solution object from the MILP solver.
%
% .. Author: - Ehsan Motamedian, July 2026

%% 1. Set Defaults
if nargin < 3 || isempty(params), params = struct(); end
if ~isfield(params, 'sigma'),         params.sigma = 0.7; end
if ~isfield(params, 'kappa_t'),       params.kappa_t = 1.3; end
if ~isfield(params, 'k_trans'),       params.k_trans = 8 * 3600; end
if ~isfield(params, 'phi_0'),         params.phi_0 = 0.05; end
if ~isfield(params, 'f_pool'),        params.f_pool = 0.71; end
if ~isfield(params, 'total_protein'), params.total_protein = 460; end
if ~isfield(params, 'k_deg'),         params.k_deg = 0.05; end
if ~isfield(params, 'k_deg_RNA'),     params.k_deg_RNA = 2.1; end
if ~isfield(params, 'f_active'),      params.f_active = 0.85; end
if ~isfield(params, 'MW_ribosome'),   params.MW_ribosome = 3.0e6; end
if ~isfield(params, 'k_poly'),        params.k_poly = 50 * 3600; end
if ~isfield(params, 'RNAP_total'),    params.RNAP_total = 0.01; end
params.P_total = params.total_protein / 1000;

fprintf('\n--- Solving HyTT at Glucose Input = %.2f ---\n', glucose_rate);

%% 2. Setup Constraints
ML = model.ML_data;
total_capacity_mg_ml = params.total_protein * params.f_pool;

% Apply glucose uptake as lower bound (negative for uptake)
glc_rxn = 'r_1714';
if ~ismember(glc_rxn, model.rxns), glc_rxn = 'r_1714_REV'; end
model = changeRxnBounds(model, glc_rxn, -glucose_rate, 'l');

biomassIdx = find(strcmp(model.rxns, 'r_4041'));

% Prepare MILP Problem
Problem = struct();
Problem.A = model.S;
Problem.b = model.b;
Problem.csense = model.csense;
Problem.lb = model.lb;
Problem.ub = model.ub;
Problem.vartype = model.vartype;
Problem.c = model.c;
Problem.osense = 1; 

%% 3. Binary Search for Mu
changeCobraSolver('gurobi'); changeCobraSolver('gurobi', 'MILP');
mu_low = 0; mu_high = 0.6; tol_mu = 1e-3; 
max_mu_found = 0; best_sol = [];

while (mu_high - mu_low) > tol_mu
    mu_test = (mu_low + mu_high) / 2;
    Problem.lb(biomassIdx) = mu_test - 1e-6;
    Problem.ub(biomassIdx) = mu_test + 1e-6;
    
    % Update growth-dependent matrix coefficients
    Problem.A(ML.row_trans_cap, ML.ml_rxn_cols) = ((mu_test + params.k_deg) .* ML.l_prot_arr ./ (params.sigma .* ML.MW_arr))';
    Problem.A(ML.row_trsc_cap, ML.mRNA_rxn_cols) = -((mu_test + params.k_deg_RNA) * 1.66e-9 .* ML.l_mRNA_arr)';
    
    R_active = (params.P_total * (params.phi_0 + mu_test/params.kappa_t) * params.f_active) / params.MW_ribosome * 1000;
    
    Problem.b(ML.row_trans_cap) = params.k_trans * R_active;
    Problem.b(ML.row_trsc_cap) = params.k_poly * params.RNAP_total;
    Problem.b(ML.row_ribo_demand) = total_capacity_mg_ml * (params.phi_0 + mu_test/params.kappa_t);
    
    s = solveCobraMILP(Problem, 'timeLimit', 300, 'printLevel', 0);
    
    if (s.stat == 1 || s.stat == 3)
        if s.full(ML.idx_slack_trans) <= 1e-3
            mu_low = mu_test; max_mu_found = mu_test; best_sol = s; 
        else
            mu_high = mu_test; 
        end
    else
        mu_high = mu_test; 
    end
end

%% 4. Extract Results
results = struct();
results.glucose_input = glucose_rate;
results.mu = max_mu_found;

if max_mu_found > 0 && ~isempty(best_sol)
    s_f = best_sol.full;
    
    idx_eth  = find(strcmp(model.rxns, 'r_1761'));
    idx_o2   = find(strcmp(model.rxns, 'r_1992'));
    idx_co2  = find(strcmp(model.rxns, 'r_1672'));
    idx_gluc = find(strcmp(model.rxns, 'r_1714'));
    if isempty(idx_gluc), idx_gluc = find(strcmp(model.rxns, 'r_1714_REV')); end
    idx_pool = find(strcmp(model.rxns, 'prot_pool_exchange'));
    
    results.actual_glucose = abs(s_f(idx_gluc));
    results.ethanol = s_f(idx_eth);
    results.oxygen  = abs(s_f(idx_o2));
    results.co2     = s_f(idx_co2);
    results.pool_usage_mg = abs(s_f(idx_pool));
    
    a_ribo = 0; a_meta = 0; g_cnt = 0;
    
    for i = 1:length(ML.realGeneIdx)
        if ML.realGeneIdx(i)
            g_cnt = g_cnt + 1; 
            k = ML.ml_rxn_cols(g_cnt); 
            
            if k > 0
                m_i = s_f(k) / params.sigma;
                if m_i > 1e-8
                    if ML.is_ribo_arr(g_cnt) == 1
                        a_ribo = a_ribo + m_i; 
                    else
                        a_meta = a_meta + m_i; 
                    end
                end
            end
        end
    end
    
    results.ribo_pct = (a_ribo / params.total_protein) * 100;
    results.meta_pct = (a_meta / params.total_protein) * 100;
    
    fprintf('>>> mu: %.4f | Gluc_act: %.2f | EtOH: %.2f | O2: %.2f | Meta: %.1f%% | Ribo: %.1f%%\n', ...
            results.mu, results.actual_glucose, results.ethanol, results.oxygen, results.meta_pct, results.ribo_pct);
else
    fprintf('>>> Result: INFEASIBLE or Growth is Zero.\n');
end

end