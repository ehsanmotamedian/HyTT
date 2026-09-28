% =========================================================================
% HYTT FRAMEWORK: PhaA + PhaB HETEROLOGOUS PATHWAY BURDEN
% (Case study where enzymes of a heterologous pathway
% are expressed" -- extends the single-enzyme PhaB case to a real 2-step
% pathway: acetyl-CoA -[PhaA]-> acetoacetyl-CoA -[PhaB]-> 3-hydroxybutyryl-CoA)
%
% Enzymes:
%   PhaA, acetyl-CoA acetyltransferase / beta-ketothiolase (EC 2.3.1.9)
%   PhaB, acetoacetyl-CoA reductase (EC 1.1.1.36)
%   Source (both): Cupriavidus necator H16
% Reactions modeled:
%   PhaA: 2 acetyl-CoA -> acetoacetyl-CoA + CoA
%   PhaB: acetoacetyl-CoA + NADPH + H+ -> (R)-3-hydroxybutyryl-CoA + NADP+
% Sequences (both counted directly from real protein sequences, no
% GFP-derived scaling):
%   PhaA (UniProt P14611, PDB 4O9C, 393 aa, MW = 40,532.8 Da):
%   PhaB (UniProt P14697, 246 aa, MW = 26,370 Da)
% =========================================================================
clc; clear; close all;
warning('off', 'all');

disp('=========================================================================');
disp('   HYBRID RECOMBINANT BURDEN: PhaA + PhaB HETEROLOGOUS PATHWAY (15% MASS)');
disp('=========================================================================');

% --- 1. GECKO ENVIRONMENT & BASE RECONSTRUCTION ---
geckoRoot = findGECKOroot;
adapterPath = fullfile(geckoRoot, 'tutorials', 'full_ecModel', 'YeastGEMAdapter.m');
ModelAdapterManager.setDefault(adapterPath, true);
changeCobraSolver('gurobi', 'MILP');

ecModel_base = loadEcModel('ecYeastGEM.yml');
ecModel_base = addRibosomalUsageReactions(ecModel_base);
ecModel_base = addmRNAUsageReactions(ecModel_base);

if ~ismember('EX_mRNA_pool[c]', ecModel_base.rxns)
    ecModel_base = addExchangeRxn(ecModel_base, 'mRNA_pool[c]', -100000, 0);
end
geneFeaturesTable = readtable('Yeast_GeneFeatures.xlsx', 'VariableNamingRule', 'preserve');

% --- 2. DEFINE SYSTEM BIOPHYSICAL PARAMETERS ---
params_wt = struct();
params_wt.f_pool = 0.71;
params_wt.sigma = 0.7;
params_wt.kappa_t = 1.3;
params_wt.phi_0 = 0.05;
params_wt.f_active = 0.85;
params_wt.total_protein = 460;
params_wt.MW_ribosome = 3.0e6;
params_wt.ML_function = @applyMARS_Yeast;

fixed_glucose = 11.1;
burden_fraction_total = 0.15;   % total pathway allocation, matched to the GFP/single-PhaB case
burden_fraction_A = burden_fraction_total / 2;  % split evenly between the two enzymes (0.075 each)
burden_fraction_B = burden_fraction_total / 2;

% --- 3. LAYER 0: WILD-TYPE BASELINE OPTIMIZATION ---
ML_WT = buildHyTT(ecModel_base, geneFeaturesTable, params_wt);
[results_WT, sol_WT] = solveHyTT(ML_WT, fixed_glucose, params_wt);

if isempty(sol_WT) || results_WT.mu == 0
    error('FATAL: Wild-type baseline optimization failed. Check base model bounds.');
end

% --- 4. LAYER 1: PRECURSOR DEMAND INJECTION (PhaA + PhaB Synthesis) ---
ecModel_burden = ecModel_base;

% ---- 4a. PhaA (393 aa) ----
if ~ismember('prot_PhaA[c]', ecModel_burden.mets)
    ecModel_burden = addMetabolite(ecModel_burden, 'prot_PhaA[c]', 'metName', 'Acetyl-CoA Acetyltransferase (PhaA) Protein');
end

% EXACT composition, counted directly from the real PhaA sequence above.
% Residue counts (order matches aaMets_PhaA below):
%   Gly 43, Leu 27, Lys 26, Val 38, Asp 19, Glu 20, Thr 16, Asn 13,
%   Phe 8, Ile 20, Tyr 4, Ser 19, Pro 18, His 5, Ala 65, Gln 15,
%   Met 18, Arg 14, Cys 2, Trp 3   (sum = 393)
aaMets_PhaA = {'s_1003', 's_1021', 's_1025', 's_1056', 's_0973', 's_0991', 's_1045', 's_0969', 's_1032', 's_1016', 's_1051', 's_1039', 's_1035', 's_1006', 's_0955', 's_0999', 's_1029', 's_0965', 's_0981', 's_1048'};
aaCoeffs_PhaA = [43, 27, 26, 38, 19, 20, 16, 13, 8, 20, 4, 19, 18, 5, 65, 15, 18, 14, 2, 3];

pep_bonds_A = 392;         % 393 aa -> 392 peptide bonds
atpCost_A = 4 * 393;       % 1572
gtpCost_A = 2 * 393;       % 786

rxnFormula_PhaA = '';
for i = 1:length(aaMets_PhaA)
    rxnFormula_PhaA = [rxnFormula_PhaA num2str(aaCoeffs_PhaA(i)) ' ' aaMets_PhaA{i} ' + '];
end
phaa_synth_formula = [rxnFormula_PhaA num2str(atpCost_A) ' s_0434 + ' num2str(gtpCost_A) ' s_0785 => ' num2str(pep_bonds_A) ' s_0803 + ' ...
    num2str(atpCost_A) ' s_0394 + ' num2str(gtpCost_A) ' s_0739 + ' num2str(3 * 393) ' s_1322 + prot_PhaA[c]'];

ecModel_burden = addReaction(ecModel_burden, 'r_PhaA_synthesis', ...
    'reactionFormula', phaa_synth_formula, 'reversible', false, 'lowerBound', 0, 'upperBound', 1000);

% ---- 4b. PhaB (246 aa) 
if ~ismember('prot_PhaB[c]', ecModel_burden.mets)
    ecModel_burden = addMetabolite(ecModel_burden, 'prot_PhaB[c]', 'metName', 'Acetoacetyl-CoA Reductase (PhaB) Protein');
end

aaMets_PhaB = {'s_1003', 's_1021', 's_1025', 's_1056', 's_0973', 's_0991', 's_1045', 's_0969', 's_1032', 's_1016', 's_1051', 's_1039', 's_1035', 's_1006', 's_0955', 's_0999', 's_1029', 's_0965', 's_0981', 's_1048'};
aaCoeffs_PhaB = [28, 14, 14, 23, 16, 10, 19, 11, 10, 16, 3, 15, 5, 2, 23, 10, 7, 12, 3, 5];

pep_bonds_B = 245;
atpCost_B = 4 * 246;
gtpCost_B = 2 * 246;

rxnFormula_PhaB = '';
for i = 1:length(aaMets_PhaB)
    rxnFormula_PhaB = [rxnFormula_PhaB num2str(aaCoeffs_PhaB(i)) ' ' aaMets_PhaB{i} ' + '];
end
phab_synth_formula = [rxnFormula_PhaB num2str(atpCost_B) ' s_0434 + ' num2str(gtpCost_B) ' s_0785 => ' num2str(pep_bonds_B) ' s_0803 + ' ...
    num2str(atpCost_B) ' s_0394 + ' num2str(gtpCost_B) ' s_0739 + ' num2str(3 * 246) ' s_1322 + prot_PhaB[c]'];

ecModel_burden = addReaction(ecModel_burden, 'r_PhaB_synthesis', ...
    'reactionFormula', phab_synth_formula, 'reversible', false, 'lowerBound', 0, 'upperBound', 1000);

% ---- 4c. Lock both to Biomass (The Spatial Squeeze) ----
biomass_idx = find(strcmp(ecModel_burden.rxns, 'r_4041'), 1);

phaa_mw_da = 40532.8;  % exact average MW computed from the P14611 sequence above
phab_mw_da = 26370;    % exact average MW computed from the P14697 sequence (unchanged)

coupling_coeff_A = (burden_fraction_A * params_wt.total_protein) / (phaa_mw_da * (1 - burden_fraction_total));
coupling_coeff_B = (burden_fraction_B * params_wt.total_protein) / (phab_mw_da * (1 - burden_fraction_total));

ecModel_burden.S(find(strcmp(ecModel_burden.mets, 'prot_PhaA[c]')), biomass_idx) = -coupling_coeff_A;
ecModel_burden.S(find(strcmp(ecModel_burden.mets, 'prot_PhaB[c]')), biomass_idx) = -coupling_coeff_B;
fprintf('>> Precursor Matrix Lock Engaged. PhaA Coupling: %.6e | PhaB Coupling: %.6e\n', coupling_coeff_A, coupling_coeff_B);

% --- 5. CATALYTIC FLUX: 2 Acetyl-CoA -[PhaA]-> Acetoacetyl-CoA -[PhaB]-> 3-Hydroxybutyryl-CoA ---
if ~ismember('3hb_CoA[c]', ecModel_burden.mets)
    ecModel_burden = addMetabolite(ecModel_burden, '3hb_CoA[c]', 'metName', '(R)-3-hydroxybutyryl-CoA');
end

% Real ecYeastGEM metabolite IDs (verified against the Metabolite List sheet):
%   s_0373 = acetyl-CoA, s_0367 = acetoacetyl-CoA, s_0529 = coenzyme A,
%   s_1212 = NADPH, s_1207 = NADP(+), s_0793 = H+
phaa_catalytic_formula = '2 s_0373 => s_0367 + s_0529';
ecModel_burden = addReaction(ecModel_burden, 'r_PhaA_catalysis', 'reactionFormula', phaa_catalytic_formula, ...
    'lowerBound', 0, 'upperBound', 1000);

phab_catalytic_formula = 's_0367 + s_1212 + s_0793 => 3hb_CoA[c] + s_1207';
ecModel_burden = addReaction(ecModel_burden, 'r_PhaB_catalysis', 'reactionFormula', phab_catalytic_formula, ...
    'lowerBound', 0, 'upperBound', 1000);

% Drain the truncated pathway product (no PhaC in this case study)
if ~ismember('sink_3hbCoA_repair', ecModel_burden.rxns)
    ecModel_burden = addReaction(ecModel_burden, 'sink_3hbCoA_repair', 'reactionFormula', '3hb_CoA[c] => ', ...
        'lowerBound', 0, 'upperBound', 1000);
end

% --- 6. LAYER 2 & 3: SPATIAL AND TRANSLATIONAL BURDEN ---
k_t_host = params_wt.kappa_t;
k_t_enzyme = 0.80;  % same assumed elongation-capacity factor for both heterologous enzymes
k_t_effective = 1 / ( ((1 - burden_fraction_total) / k_t_host) + (burden_fraction_total / k_t_enzyme) );

params_burden = params_wt;
params_burden.f_pool = params_wt.f_pool - burden_fraction_total;
params_burden.kappa_t = k_t_effective;

% --- 7. COMPILE AND SOLVE ---
ML_burden_compiled = buildHyTT(ecModel_burden, geneFeaturesTable, params_burden);

% Reopen all custom reactions (buildHyTT appears to cap/rescale new
% reaction bounds during matrix expansion):
ML_burden_compiled = changeRxnBounds(ML_burden_compiled, 'r_PhaA_synthesis', 1000, 'u');
ML_burden_compiled = changeRxnBounds(ML_burden_compiled, 'r_PhaB_synthesis', 1000, 'u');
ML_burden_compiled = changeRxnBounds(ML_burden_compiled, 'r_PhaA_catalysis', 1000, 'u');
ML_burden_compiled = changeRxnBounds(ML_burden_compiled, 'r_PhaB_catalysis', 1000, 'u');
ML_burden_compiled = changeRxnBounds(ML_burden_compiled, 'sink_3hbCoA_repair', 1000, 'u');

% --- 7b. BISECTION
lb_low = 0;
lb_high = 5;      
tol_flux = 0.05;
max_iter = 6;
best_feasible_lb = 0;
best_results = [];
best_sol = [];

fprintf('\n>> Bisecting for max sustainable PhaA/PhaB coupled catalytic flux (0 to %.3f mmol/gDCW/h)...\n', lb_high);
iter = 0;
while (lb_high - lb_low) > tol_flux && iter < max_iter
    iter = iter + 1;
    lb_test = (lb_low + lb_high) / 2;
    fprintf('   [iter %d/%d] testing coupled lowerBound = %.4f ...\n', iter, max_iter, lb_test);
    ML_test = changeRxnBounds(ML_burden_compiled, 'r_PhaA_catalysis', lb_test, 'l');
    ML_test = changeRxnBounds(ML_test, 'r_PhaB_catalysis', lb_test, 'l');
    [res_test, sol_test] = solveHyTT(ML_test, fixed_glucose, params_burden);
    if ~isempty(sol_test) && res_test.mu > 0
        lb_low = lb_test;
        best_feasible_lb = lb_test;
        best_results = res_test;
        best_sol = sol_test;
    else
        lb_high = lb_test;
    end
end
fprintf('>> Max sustainable coupled PhaA/PhaB flux found: %.4f mmol/gDCW/h (%d bisection iterations)\n', best_feasible_lb, iter);

% --- 7c. TWO REPORTING POINTS ---
% The bisected maximum sits right at the edge of near-zero growth (a
% "stress test" of the model's absolute ceiling), which is not a fair or
% convincing basis for the main comparison against GFP/PhaB-alone -- showing 
% a barely-alive phenotype overstates the effect. Report BOTH:
%   (a) the stress-test maximum (100% of best_feasible_lb)
%   (b) a moderate, non-edge operating point (50% of best_feasible_lb)
% and let the moderate point carry the main growth/oxygen/ethanol story.
flux_max = best_feasible_lb;
flux_moderate = best_feasible_lb / 2;

ML_moderate = changeRxnBounds(ML_burden_compiled, 'r_PhaA_catalysis', flux_moderate, 'l');
ML_moderate = changeRxnBounds(ML_moderate, 'r_PhaB_catalysis', flux_moderate, 'l');
[results_Moderate, sol_Moderate] = solveHyTT(ML_moderate, fixed_glucose, params_burden);

ML_burden_compiled = changeRxnBounds(ML_burden_compiled, 'r_PhaA_catalysis', flux_max, 'l');
ML_burden_compiled = changeRxnBounds(ML_burden_compiled, 'r_PhaB_catalysis', flux_max, 'l');
results_Burden = best_results;
sol_Burden = best_sol;

% --- 8. RESULTS REPORT (USING EXACT MASS EXTRACTION) ---
if ~isempty(sol_Burden) && results_Burden.mu > 0

    vec_WT = get_flux_vector(sol_WT);
    vec_Max = get_flux_vector(sol_Burden);
    vec_Mod = get_flux_vector(sol_Moderate);

    [~, ribo_wt, meta_wt] = extract_ML_Masses_Exact(vec_WT, ML_WT, params_wt);
    pct_meta_wt = (meta_wt / params_wt.total_protein) * 100;
    pct_ribo_wt = (ribo_wt / params_wt.total_protein) * 100;

    [~, ribo_max, meta_max] = extract_ML_Masses_Exact(vec_Max, ML_burden_compiled, params_wt);
    pct_meta_max = (meta_max / params_wt.total_protein) * 100;
    pct_ribo_max = (ribo_max / params_wt.total_protein) * 100;

    have_moderate = ~isempty(sol_Moderate) && results_Moderate.mu > 0;
    if have_moderate
        [~, ribo_mod, meta_mod] = extract_ML_Masses_Exact(vec_Mod, ML_moderate, params_wt);
        pct_meta_mod = (meta_mod / params_wt.total_protein) * 100;
        pct_ribo_mod = (ribo_mod / params_wt.total_protein) * 100;
    end

    idx_A = find(strcmp(ML_burden_compiled.rxns, 'r_PhaA_catalysis'), 1);
    idx_B = find(strcmp(ML_burden_compiled.rxns, 'r_PhaB_catalysis'), 1);
    flux_A_max = 0; flux_B_max = 0; flux_A_mod = 0; flux_B_mod = 0;
    if ~isempty(idx_A), flux_A_max = vec_Max(idx_A); end
    if ~isempty(idx_B), flux_B_max = vec_Max(idx_B); end
    if have_moderate && ~isempty(idx_A), flux_A_mod = vec_Mod(idx_A); end
    if have_moderate && ~isempty(idx_B), flux_B_mod = vec_Mod(idx_B); end

    fprintf('\n=========================================================================\n');
    fprintf('   HYTT FLASH RECONCILIATION REPORT (PhaA+PhaB Pathway, 15%% Combined Burden)\n');
    fprintf('=========================================================================\n');
    fprintf('   Phenotypic Metric      |  Wild-Type  |  Moderate (50%% max)  |  Stress-Test Max \n');
    fprintf('-------------------------------------------------------------------------------\n');
    fprintf('   Growth Rate (mu h-1)   |   %.4f    |       %.4f          |      %.4f      \n', results_WT.mu, results_Moderate.mu, results_Burden.mu);
    fprintf('   Glucose Consumption    |   %.2f    |       %.2f          |      %.2f      \n', results_WT.actual_glucose, results_Moderate.actual_glucose, results_Burden.actual_glucose);
    fprintf('   Oxygen Consumption     |   %.2f    |       %.2f          |      %.2f      \n', results_WT.oxygen, results_Moderate.oxygen, results_Burden.oxygen);
    fprintf('   Ethanol Secretion      |   %.2f    |       %.2f          |      %.2f      \n', results_WT.ethanol, results_Moderate.ethanol, results_Burden.ethanol);
    fprintf('   Metabolic Proteome (%%) |   %.2f%%    |       %.2f%%          |      %.2f%%      \n', pct_meta_wt, pct_meta_mod, pct_meta_max);
    fprintf('   Ribosomal Proteome (%%) |   %.2f%%    |       %.2f%%          |      %.2f%%      \n', pct_ribo_wt, pct_ribo_mod, pct_ribo_max);
    fprintf('   PhaA Catalytic Flux    |      --      |       %.4f          |      %.4f      \n', flux_A_mod, flux_A_max);
    fprintf('   PhaB Catalytic Flux    |      --      |       %.4f          |      %.4f      \n', flux_B_mod, flux_B_max);
    fprintf('=========================================================================\n');
    fprintf('   Use the "Moderate (50%% max)" column as the\n');
    fprintf('   main reported operating point -- it carries a real, comparable growth\n');
    fprintf('   rate (not a near-zero edge case) while still exercising genuine\n');
    fprintf('   catalytic flux through both heterologous enzymes. The "Stress-Test Max"\n');
    fprintf('   column is supporting evidence for a supplementary sensitivity statement\n');
    fprintf('   (e.g., "beyond X mmol/gDCW/h combined flux, the model predicts the\n');
    fprintf('   pathway becomes growth-limiting"), not the headline comparison.\n');
    fprintf('=========================================================================\n');

    cache_file = 'HyTT_Allocation_Results_PhaA_PhaB.mat';
    save(cache_file, 'ML_WT', 'ML_burden_compiled', 'ML_moderate', 'sol_WT', 'sol_Burden', 'sol_Moderate', 'params_wt', 'burden_fraction_total');
    fprintf('>>> Saved cache to %s\n', cache_file);
else
    warning('Model is INFEASIBLE. Check acetyl-CoA/NADPH availability or capacity constraints.');
end

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================
function vec = get_flux_vector(sol)
    vec = [];
    if isempty(sol), return; end
    if isfield(sol, 'full')
        vec = sol.full;
    elseif isfield(sol, 'x')
        vec = sol.x;
    end
end

function [masses_sorted, mass_ribo, mass_meta] = extract_ML_Masses_Exact(vec, ML_model, params)
    masses = []; mass_ribo = 0; mass_meta = 0;
    if isempty(vec), masses_sorted = []; return; end

    idx_ml = find(startsWith(ML_model.rxns, 'ml_prot_'));
    for i = 1:length(idx_ml)
        col_k = idx_ml(i);
        flux = vec(col_k);
        if flux > 1e-8
            mass = flux / params.sigma;
            masses(end+1) = mass;

            uniprot = strrep(ML_model.rxns{col_k}, 'ml_prot_', '');
            usage_rxn = ['usage_prot_', uniprot];
            idx_u = find(strcmp(ML_model.rxns, usage_rxn), 1);
            is_ribo = false;
            if ~isempty(idx_u)
                subs = ML_model.subSystems{idx_u};
                if iscell(subs), subs_str = strjoin(string(subs), ' '); else, subs_str = string(subs); end
                if contains(subs_str, 'Ribosome', 'IgnoreCase', true), is_ribo = true; end
            end
            if is_ribo, mass_ribo = mass_ribo + mass; else, mass_meta = mass_meta + mass; end
        end
    end
    masses_sorted = sort(masses, 'descend');
end
