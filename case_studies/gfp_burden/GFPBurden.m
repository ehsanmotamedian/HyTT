% =========================================================================
% HYTT FRAMEWORK: MOLECULAR MODELING OF HETEROLOGOUS RECOMBINANT BURDEN
% Three-Layer Integration: Precursor Demand, Spatial Burden, Translational Burden
% Includes Auto-Repair for Stoichiometric Dead-Ends in User Formulas
% .. Author: Ehsan Motamedian, September 2026
% =========================================================================

clc; clear; close all;

disp('=========================================================================');
disp('   HYBRID RECOMBINANT BURDEN: COMPLETE DECOUPLING & THREE-LAYER Burden');
disp('=========================================================================');

cache_file = 'HyTT_Allocation_Results.mat';
disp('>> Initializing simulation and preparing fresh data cache...');

% --- 1. GECKO ENVIRONMENT & BASE RECONSTRUCTION ---
disp('Initializing GECKO Adapter...');
geckoRoot = findGECKOroot;
adapterPath = fullfile(geckoRoot, 'tutorials', 'full_ecModel', 'YeastGEMAdapter.m');
ModelAdapterManager.setDefault(adapterPath, true);

fixed_glucose = 11.1; 
burden_fraction = 0.15; % 15% Total Proteome Allocation Target for GFP

warning('off', 'all'); 
changeCobraSolver('gurobi', 'MILP');

disp('Loading Base ecYeastGEM and expanding structures...');
ecModel_base = loadEcModel('ecYeastGEM.yml');
ecModel_base = addRibosomalUsageReactions(ecModel_base);
ecModel_base = addmRNAUsageReactions(ecModel_base);

if ~ismember('EX_mRNA_pool[c]', ecModel_base.rxns)
    ecModel_base = addExchangeRxn(ecModel_base, 'mRNA_pool[c]', -100000, 0);
end

disp('Loading Genomic Biophysical Features...');
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

% --- 3. LAYER 0: WILD-TYPE BASELINE OPTIMIZATION ---
disp('>> Compiling and optimizing Wild-Type Baseline...');
ML_WT = buildHyTT(ecModel_base, geneFeaturesTable, params_wt);
[results_WT, sol_WT] = solveHyTT(ML_WT, fixed_glucose, params_wt);

if isempty(sol_WT) || results_WT.mu == 0
    error('FATAL: Wild-type baseline optimization failed. Check base model bounds.');
end

% --- 4. LAYER 1: PRECURSOR DEMAND INJECTION ---
disp('>> Injecting Layer 1: Recombinant Sequence Precursor Cost Reactions...');
ecModel_burden = ecModel_base;

if ~ismember('prot_GFP[c]', ecModel_burden.mets)
    ecModel_burden = addMetabolite(ecModel_burden, 'prot_GFP[c]', 'metName', 'GFP Protein Species');
end

aaMets = {'s_1003', 's_1021', 's_1025', 's_1056', 's_0973', 's_0991', 's_1045', 's_0969', 's_1032', 's_1016', 's_1051', 's_1039', 's_1035', 's_1006', 's_0955', 's_0999', 's_1029', 's_0965', 's_0981', 's_1048'};
aaCoeffs = [22, 21, 20, 18, 18, 16, 16, 13, 12, 12, 11, 10, 10, 9, 8, 8, 6, 6, 2, 1];
% Energy cost: 4 ATP + 2 GTP per aa
atpCost = 952; gtpCost = 476;  % 4*238, 2*238

% Build reaction formula
rxnFormula = '';
for i = 1:length(aaMets)
    rxnFormula = [rxnFormula num2str(aaCoeffs(i)) ' ' aaMets{i} ' + '];
end
gfp_formula = [rxnFormula num2str(atpCost) ' s_0434 + ' num2str(gtpCost) ' s_0785 => 237 s_0803 + ' ...
    num2str(atpCost) ' s_0394 + ' num2str(gtpCost) ' s_0739 + ' num2str(3*238) ' s_1322 + prot_GFP[c]'];

gfp_rxn_name = 'r_GFP_synthesis';
if ismember(gfp_rxn_name, ecModel_burden.rxns)
    ecModel_burden = changeRxns(ecModel_burden, gfp_rxn_name, gfp_formula);
else
    ecModel_burden = addReaction(ecModel_burden, gfp_rxn_name, ...
        'reactionFormula', gfp_formula, 'reversible', false, 'lowerBound', 0, 'upperBound', 1000);
end

% Tie GFP production directly to the structural Biomass Equation
biomass_idx = find(strcmp(ecModel_burden.rxns, 'r_4041'), 1);
if ~isempty(biomass_idx)
    gfp_mw_da = 26900; 
    coupling_coeff = (burden_fraction * params_wt.total_protein) / (gfp_mw_da * (1 - burden_fraction));

    % Minus sign because Biomass CONSUMES the demand link
    ecModel_burden.S(find(strcmp(ecModel_burden.mets, 'prot_GFP[c]')), biomass_idx) = -coupling_coeff;
    fprintf('>> Precursor Matrix Lock Engaged. Coupling Factor: %.6e\n', coupling_coeff);
else
    error('Biomass macro-reaction r_4041 not found for precursor coupling.');
end

% --- 5. LAYER 2 & 3: SPATIAL AND TRANSLATIONAL Burden CONFIGURATION ---
disp('>> Injecting Layer 2 (Spatial Burden) & Layer 3 (Translational Burden)...');

k_t_host = params_wt.kappa_t;
k_t_GFP  = 0.80; 

fraction_host = 1 - burden_fraction;
k_t_effective = 1 / ( (fraction_host / k_t_host) + (burden_fraction / k_t_GFP) );

params_burden = params_wt;
params_burden.f_pool = params_wt.f_pool - burden_fraction; 
params_burden.kappa_t = k_t_effective;                     

fprintf('>> Baseline Host kappa_t: %.4f h^-1\n', params_wt.kappa_t);
fprintf('>> Burden Effective kappa_t: %.4f h^-1\n', params_burden.kappa_t);

% --- 6. COMPILE AND SOLVE FULL METABOLIC BURDEN SYSTEM ---
ML_burden_compiled = buildHyTT(ecModel_burden, geneFeaturesTable, params_burden);
ML_Burden = changeRxnBounds(ML_burden_compiled, gfp_rxn_name, 1000, 'u');

disp('>> Executing Spatio-Translational Optimization for Recombinant Burden...');
[results_Burden, sol_Burden] = solveHyTT(ML_Burden, fixed_glucose, params_burden);

% --- 7. CACHE GENERATION & DATA INTEGRITY VERIFICATION ---
if ~isempty(sol_Burden) && results_Burden.mu > 0
    disp('>> Optimization successful. Saving structured data variables to disk...');
    
    save(cache_file, 'ML_WT', 'ML_Burden', 'sol_WT', 'sol_Burden', 'params_wt', 'burden_fraction');
    fprintf('>>> SUCCESS: Binary cache file [%s] written with complete allocation blocks.\n', cache_file);
    
    fprintf('\n=========================================================================\n');
    fprintf('   HYTT FLASH RECONCILIATION REPORT\n');
    fprintf('=========================================================================\n');
    fprintf('   Phenotypic Metric      |   Wild-Type Baseline   |   15%% GFP Burden      \n');
    fprintf('-------------------------------------------------------------------------\n');
    fprintf('   Growth Rate (mu h-1)   |        %.4f          |        %.4f          \n', results_WT.mu, results_Burden.mu);
    fprintf('   Glucose Consumption    |        %.2f          |        %.2f          \n', results_WT.actual_glucose, results_Burden.actual_glucose);
    fprintf('   Ethanol Secretion      |        %.2f          |        %.2f          \n', results_WT.ethanol, results_Burden.ethanol);
    fprintf('   Oxygen Consumption     |        %.2f          |        %.2f          \n', results_WT.oxygen, results_Burden.oxygen);
    fprintf('=========================================================================\n');
else
    warning('Critical Failure: Recombinant burden model equations are infeasible.');
end
