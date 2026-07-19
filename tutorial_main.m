% =========================================================================
% HyTT Tutorial: Ab Initio Prediction of Resource Allocation and Overflow Metabolism
% in S. cerevisiae
%
% This script demonstrates the full pipeline of the Hybrid Translation-Transcription 
% (HyTT) framework, coupling multivariate adaptive regression splines (MARS) constraints
% with an enzyme-constrained metabolic model. It also mathematically validates the 
% strict 80S stoichiometric assembly rules.
%
% .. Author: Ehsan Motamedian, July 2026
% =========================================================================

clc; clear;
% Add subdirectories to MATLAB path
addpath(genpath(pwd));
% 0. Initialize GECKO 3 Environment (The Missing Link!)
disp('Initializing GECKO Adapter...');
geckoRoot = findGECKOroot;
adapterPath = fullfile(geckoRoot, 'tutorials', 'full_ecModel', 'YeastGEMAdapter.m');
ModelAdapterManager.setDefault(adapterPath, true);

% 1. Load the Base GECKO Model
disp('Loading Base ecYeastGEM...');
ecModel = loadEcModel('ecYeastGEM.yml');

% 2. PRE-PROCESSING: Add your custom Usage Reactions
disp('Adding Ribosomal and mRNA usage reactions...');
ecModel = addRibosomalUsageReactions(ecModel);
ecModel = addmRNAUsageReactions(ecModel);

% Make sure the mRNA pool is open for the LP/MILP to calculate
if ~ismember('EX_mRNA_pool[c]', ecModel.rxns)
    ecModel = addExchangeRxn(ecModel, 'mRNA_pool[c]', -100000, 0);
end

% 3. Load Gene Features Table
disp('Loading Gene Features...');
geneFeaturesTable = readtable('Yeast_GeneFeatures.xlsx', 'VariableNamingRule', 'preserve');

% 4. Define Biological Parameters
params = struct();
params.f_pool = 0.71;
params.sigma = 0.7;
params.kappa_t = 1.3;
params.ML_function = @applyMARS_Yeast; 

% 5. Compile the HyTT Model
ML_model = buildHyTT(ecModel, geneFeaturesTable, params);

% 6. Run the Simulation (Glucose = 11.1)
glucose_rate = 11.1;
[results, solution] = solveHyTT(ML_model, glucose_rate, params);

% 7. Display Results
if results.mu > 0
    fprintf('\n--- HyTT Model Final Report ---\n');
    fprintf('Growth Rate: %.4f h^-1\n', results.mu);
    fprintf('Actual Glucose Uptake: %.2f mmol/gDW/h\n', results.actual_glucose);
    fprintf('Ethanol Production: %.2f mmol/gDW/h\n', results.ethanol);
    fprintf('Oxygen Consumption: %.2f mmol/gDW/h\n', results.oxygen);
    fprintf('CO2 Production: %.2f mmol/gDW/h\n', results.co2);
    fprintf('Active Proteome Usage: %.2f mg/gDW\n', results.pool_usage_mg);
    fprintf('Ribosomal Proteome Fraction: %.1f%%\n', results.ribo_pct);
    fprintf('Metabolic Proteome Fraction: %.1f%%\n', results.meta_pct);
end

