function [model] = buildHyTT(model, geneFeaturesTable, params)
% buildHyTT - Integrates an ab initio proteome allocation framework with 
% machine-learning derived constraints (e.g., MARS).
%
% USAGE:
%   [model] = buildHyTT(model, geneFeaturesTable, params)
%
% INPUTS:
%   model               COBRA model structure (e.g., ecYeastGEM)
%   geneFeaturesTable   Table containing ML-derived sequence features
%
% OPTIONAL INPUTS:
%   params              Struct containing physical and biological parameters:
%                       * .f_pool         - Active proteome fraction (Default: 0.71)
%                       * .sigma          - Global saturation factor (Default: 0.70)
%                       * .total_protein  - Total cell protein content (Default: 460 mg/gDW)
%                       * .ML_function    - Function handle for ML predictions (Default: @applyMARS_Yeast)
%                       * .default_mw     - Default molecular weight for unknowns (Default: 50000 Da)
%                       * .avg_aa_mw      - Average amino acid molecular weight (Default: 110 Da)
%
% OUTPUTS:
%   model               MILP/LP ready HyTT model
%
% .. Author: - Ehsan Motamedian, June 2026
%% 1. Check Inputs & Set Defaults
if nargin < 3 || isempty(params), params = struct(); end
if ~isfield(params, 'f_pool'),        params.f_pool = 0.71; end
if ~isfield(params, 'sigma'),         params.sigma = 0.7; end
if ~isfield(params, 'total_protein'), params.total_protein = 460; end
if ~isfield(params, 'ML_function')
    params.ML_function = @applyMARS_Yeast; 
end
if ~isfield(params, 'default_mw'),     params.default_mw = 50000; end
if ~isfield(params, 'avg_aa_mw'),      params.avg_aa_mw = 110; end

fprintf('\n--- Building HyTT Model ---\n');
fprintf('Parameters: f_pool = %.2f | sigma = %.2f\n', params.f_pool, params.sigma);

%% 2. Apply Global Proteome Constraints
total_capacity_mg_ml = params.total_protein * params.f_pool;
if ismember('prot_pool_exchange', model.rxns)
    model = changeRxnBounds(model, 'prot_pool_exchange', -total_capacity_mg_ml, 'l');
end
if ismember('EX_mRNA_pool[c]', model.rxns)
    model = changeRxnBounds(model, 'EX_mRNA_pool[c]', -100000, 'l');
end

%% 3. Parse Ribosomal Data for Strict 1:1 (and 2:1 P-stalk) Stoichiometry
disp('>> Parsing Ribosomal families for strict stoichiometric assembly...');
numFamilies = 0;
uniqueFamilies = {};
family_map_k_ml = containers.Map('KeyType','char','ValueType','any');
family_map_mw   = containers.Map('KeyType','char','ValueType','any');

try
    T_ribo = readtable('RibosomalGenes.xlsx', 'VariableNamingRule', 'preserve');
    ribo_uniprots = strtrim(string(T_ribo.("UniProt AC")));
    ribo_genes = strtrim(string(T_ribo.Gene));
    baseNames = regexprep(ribo_genes, '[AB]$', '');
    uniqueFamilies = unique(baseNames);
    numFamilies = length(uniqueFamilies);
    fprintf('   Found %d unique Ribosomal families.\n', numFamilies);
catch
    warning('RibosomalGenes.xlsx not found or invalid. Ribosome stoichiometry will be skipped.');
    ribo_uniprots = [];
end

%% 4. Identify Real Genes and Add Missing Reactions
T = geneFeaturesTable;
if iscell(T.Gene), T.Gene = strtrim(string(T.Gene)); end

numGenes = length(model.ec.genes); 
realGeneIdx = true(numGenes, 1); 
missing_genes = {};

pool_met = 'prot_pool'; 
if ~ismember(pool_met, model.mets), pool_met = 'prot_pool[c]'; end

rxnIndex_ml = zeros(numGenes, 1); 
rxnIndex_mRNA = zeros(numGenes, 1); 
native_rxn_cols = zeros(numGenes, 1);

for i = 1:numGenes
    g_sys = strtrim(string(model.ec.genes{i}));
    found = any(T.Gene == g_sys);
    if ~found && isfield(model, 'genes') && isfield(model, 'geneNames')
        idx_in_model = find(strcmp(model.genes, model.ec.genes{i}), 1);
        if ~isempty(idx_in_model), found = any(T.Gene == strtrim(string(model.geneNames{idx_in_model}))); end
    end
    if ~found
        missing_genes{end+1} = char(g_sys); 
        realGeneIdx(i) = false; % Exclude missing genes
    end
end

numGenesReal = sum(realGeneIdx); 
fprintf('Total genes to expand: %d. Expanding Matrix...\n', numGenesReal);

for i = 1:numGenes
    if realGeneIdx(i)
        uniprot = char(model.ec.enzymes{i});
        rxn_mRNA = ['usage_mRNA_' uniprot]; 
        rxn_usage = ['usage_prot_' uniprot]; 
        
        exact_prot_id = ''; native_idx = 0;
        
        if ismember(rxn_usage, model.rxns)
            idx = find(strcmp(model.rxns, rxn_usage)); 
            mets = model.mets(find(model.S(:, idx))); 
            match = startsWith(mets, 'prot_') & ~strcmp(mets, 'prot_pool') & ~strcmp(mets, 'prot_pool[c]'); 
            if any(match)
                found_mets = mets(match); 
                exact_prot_id = found_mets{1}; 
            end            
            model = changeRxnBounds(model, rxn_usage, 0, 'b'); 
            native_idx = idx; 
        end
        if isempty(exact_prot_id), exact_prot_id = ['prot_' uniprot]; end

        rxn_ml = ['ml_prot_' uniprot];
        if ~ismember(rxn_ml, model.rxns)
            formula = sprintf('%.6f %s => %s', 1/params.sigma, pool_met, exact_prot_id);
            model = addReaction(model, rxn_ml, 'reactionFormula', formula, 'reversible', false, 'lowerBound', 0, 'upperBound', 100000); 
        end
        if ~ismember(['sink_mRNA_' uniprot], model.rxns)
            model = addReaction(model, ['sink_mRNA_' uniprot], 'reactionFormula', [['mRNA_' uniprot '[c]'] ' => '], 'reversible', false, 'lowerBound', 0, 'upperBound', 100000); 
        end
        if ~ismember(['sink_prot_' uniprot], model.rxns)
            model = addReaction(model, ['sink_prot_' uniprot], 'reactionFormula', [exact_prot_id ' => '], 'reversible', false, 'lowerBound', 0, 'upperBound', 100000); 
        end

        rxnIndex_mRNA(i) = find(ismember(model.rxns, rxn_mRNA), 1); 
        rxnIndex_ml(i) = find(ismember(model.rxns, rxn_ml), 1); 
        native_rxn_cols(i) = native_idx;
    end
end

%% 5. Build the MILP Matrices (A, b, csense)
nRxns = length(model.rxns); 
nMets = length(model.mets); 
A_orig = model.S; 
b_orig = model.b; 
csense_orig = repmat('E', nMets, 1);

% +3 Variables: slack_trans, slack_trsc, and v_80S
numAddVars = 11 * numGenesReal + 3; 
numAddRows = 17 * numGenesReal + 2 + numFamilies + 1; 

totalRows = nMets + numAddRows; 
totalVars = nRxns + numAddVars;

A = sparse(totalRows, totalVars); 
b = zeros(totalRows, 1); 
csense = repmat('E', totalRows, 1);

A(1:nMets, 1:nRxns) = A_orig; 
b(1:nMets) = b_orig; 
csense(1:nMets) = csense_orig; 

row_prot_pool = find(strcmp(model.mets, pool_met));
row_ribo_demand = nMets + 1; 
csense(row_ribo_demand) = 'E';

max_mRNA_exp = 5000; 
max_prot_flux = 100000;

ml_rxn_cols = zeros(numGenesReal, 1); 
mRNA_rxn_cols = zeros(numGenesReal, 1); 
l_prot_arr = zeros(numGenesReal, 1); 
l_mRNA_arr = zeros(numGenesReal, 1); 
MW_arr = zeros(numGenesReal, 1); 
is_ribo_arr = zeros(numGenesReal, 1);

currentRow = nMets + 2; 
gene_counter = 0;

for i = 1:numGenes
    if realGeneIdx(i)
        uniprot = char(model.ec.enzymes{i}); 
        gene_counter = gene_counter + 1; 
        
        rowIdx = find(T.Gene == strtrim(string(model.ec.genes{i})), 1);
        if isempty(rowIdx) && isfield(model, 'genes') && isfield(model, 'geneNames')
            idx_in_model = find(strcmp(model.genes, model.ec.genes{i}), 1);
            if ~isempty(idx_in_model)
                rowIdx = find(T.Gene == strtrim(string(model.geneNames{idx_in_model})), 1);
            end
        end
        
        MW_val = model.ec.mw(i); 
        if isnan(MW_val) || MW_val == 0, MW_val = params.default_mw; end
        
        if isempty(rowIdx)
            l_prot = round(MW_val / params.avg_aa_mw);
            l_mRNA = l_prot * 3;
            constant_adj = 0; beta = [0, 0, 0, 0]; geneRow_is_ribo = NaN;
        else
            geneRow = T(rowIdx, :);
            l_prot = geneRow.Protein_Length; l_mRNA = geneRow.mRNA_Length;
            if isnan(l_prot) || l_prot <= 0, l_prot = round(MW_val / params.avg_aa_mw); end
            if isnan(l_mRNA) || l_mRNA <= 0, l_mRNA = l_prot * 3; end
            [constant_adj, beta] = params.ML_function(geneRow);
            if ismember('Is_Ribosomal', T.Properties.VariableNames), geneRow_is_ribo = geneRow.Is_Ribosomal; else, geneRow_is_ribo = NaN; end
        end
        
        l_prot_arr(gene_counter) = l_prot; 
        l_mRNA_arr(gene_counter) = l_mRNA; 
        MW_arr(gene_counter) = MW_val; 

        baseCol = nRxns + (gene_counter-1)*11; 
        col_h = baseCol + (1:4); col_b = baseCol + (5:8); 
        col_eps_pos = baseCol + 9; col_eps_neg = baseCol + 10; 
        col_y_switch = baseCol + 11;
        
        j = rxnIndex_mRNA(i); k_ml = rxnIndex_ml(i); 
        mRNA_rxn_cols(gene_counter) = j; ml_rxn_cols(gene_counter) = k_ml; 
        
        A(row_prot_pool, k_ml) = -1 / params.sigma; 
        
        is_ribo_flag = false;
        if ~isnan(geneRow_is_ribo) && geneRow_is_ribo == 1
            is_ribo_flag = true;
        else
            idx_u = find(strcmp(model.rxns, ['usage_prot_' uniprot]), 1);
            if ~isempty(idx_u)
                subs = model.subSystems{idx_u};
                if iscell(subs), subs_str = strjoin(string(subs), ' '); else, subs_str = string(subs); end
                if contains(subs_str, 'Ribosome', 'IgnoreCase', true), is_ribo_flag = true; end
            end
        end
        
        if is_ribo_flag
            is_ribo_arr(gene_counter) = 1; 
            A(row_ribo_demand, k_ml) = 1 / params.sigma; 
            
            if numFamilies > 0
                idx_in_ribo = find(ribo_uniprots == string(uniprot), 1);
                if ~isempty(idx_in_ribo)
                    fam = char(baseNames(idx_in_ribo));
                    if isKey(family_map_k_ml, fam)
                        family_map_k_ml(fam) = [family_map_k_ml(fam), k_ml];
                        family_map_mw(fam)   = [family_map_mw(fam), MW_val];
                    else
                        family_map_k_ml(fam) = k_ml;
                        family_map_mw(fam)   = MW_val;
                    end
                end
            end
        end
        
        M1 = max(8104.542, max_mRNA_exp * l_mRNA - 8104.542);
        M2 = max(9011.816, max_mRNA_exp * l_mRNA - 9011.816);
        M3 = max(6.467589, max_mRNA_exp - 6.467589);
        M4 = max(6.467589, max_mRNA_exp - 6.467589);
        M_hing = [M1; M2; M3; M4];

        A(currentRow, col_h(1)) = 1; A(currentRow, j) = l_mRNA; b(currentRow) = -8104.542; csense(currentRow) = 'G'; currentRow = currentRow + 1;
        A(currentRow, col_h(1)) = 1; A(currentRow, j) = l_mRNA; A(currentRow, col_b(1)) = M_hing(1); b(currentRow) = M_hing(1) - 8104.542; csense(currentRow) = 'L'; currentRow = currentRow + 1;
        A(currentRow, col_h(1)) = 1; A(currentRow, col_b(1)) = -M_hing(1); b(currentRow) = 0; csense(currentRow) = 'L'; currentRow = currentRow + 1;
        
        A(currentRow, col_h(2)) = 1; A(currentRow, j) = l_mRNA; b(currentRow) = -9011.816; csense(currentRow) = 'G'; currentRow = currentRow + 1;
        A(currentRow, col_h(2)) = 1; A(currentRow, j) = l_mRNA; A(currentRow, col_b(2)) = M_hing(2); b(currentRow) = M_hing(2) - 9011.816; csense(currentRow) = 'L'; currentRow = currentRow + 1;
        A(currentRow, col_h(2)) = 1; A(currentRow, col_b(2)) = -M_hing(2); b(currentRow) = 0; csense(currentRow) = 'L'; currentRow = currentRow + 1;
        
        A(currentRow, col_h(3)) = 1; A(currentRow, j) = -1; b(currentRow) = 6.467589; csense(currentRow) = 'G'; currentRow = currentRow + 1;
        A(currentRow, col_h(3)) = 1; A(currentRow, j) = -1; A(currentRow, col_b(3)) = M_hing(3); b(currentRow) = M_hing(3) + 6.467589; csense(currentRow) = 'L'; currentRow = currentRow + 1;
        A(currentRow, col_h(3)) = 1; A(currentRow, col_b(3)) = -M_hing(3); b(currentRow) = 0; csense(currentRow) = 'L'; currentRow = currentRow + 1;
        
        A(currentRow, col_h(4)) = 1; A(currentRow, j) = 1; b(currentRow) = -6.467589; csense(currentRow) = 'G'; currentRow = currentRow + 1;
        A(currentRow, col_h(4)) = 1; A(currentRow, j) = 1; A(currentRow, col_b(4)) = M_hing(4); b(currentRow) = M_hing(4) - 6.467589; csense(currentRow) = 'L'; currentRow = currentRow + 1;
        A(currentRow, col_h(4)) = 1; A(currentRow, col_b(4)) = -M_hing(4); b(currentRow) = 0; csense(currentRow) = 'L'; currentRow = currentRow + 1;

        scale_factor = (params.sigma * MW_val) / l_prot; 
        A(currentRow, k_ml) = 1;
        for hing = 1:4, A(currentRow, col_h(hing)) = -beta(hing) * scale_factor; end
        A(currentRow, col_eps_pos) = -1 * scale_factor; 
        A(currentRow, col_eps_neg) = 1 * scale_factor; 
        A(currentRow, col_y_switch) = -constant_adj * scale_factor; 
        b(currentRow) = 0; csense(currentRow) = 'E'; currentRow = currentRow + 1;
        
        % -----------------------------------------------------------------
        % BIOPHYSICAL TEMPLATE LOCK (Direct mRNA-Protein Coupling)
        % Note: v_mRNA (j) is a negative flux (consumption). 
        % Its magnitude is (-v_mRNA).
        % -----------------------------------------------------------------
        % 1. Big-M for mRNA: If y_switch=0, mRNA magnitude must be 0 
        % -v_mRNA <= M * y_switch  =>  -v_mRNA - M * y_switch <= 0
        A(currentRow, j) = -1; A(currentRow, col_y_switch) = -max_mRNA_exp; b(currentRow) = 0; csense(currentRow) = 'L'; currentRow = currentRow + 1;
        
        % 2. Big-M for Protein: If y_switch=0, Protein must be 0
        % k_ml <= M * y_switch  =>  k_ml - M * y_switch <= 0
        A(currentRow, k_ml) = 1; A(currentRow, col_y_switch) = -max_prot_flux; b(currentRow) = 0; csense(currentRow) = 'L'; currentRow = currentRow + 1;
        
        % 3. THE ULTIMATE LOCK: Protein expression is bounded by mRNA magnitude
        % Logic: If mRNA = 0, Protein MUST be 0. 
        % k_ml <= M * (-v_mRNA)  =>  k_ml + M * v_mRNA <= 0
        A(currentRow, k_ml) = 1; A(currentRow, j) = max_prot_flux; b(currentRow) = 0; csense(currentRow) = 'L'; currentRow = currentRow + 1;
        % -----------------------------------------------------------------
    end
end

idx_slack_trans = totalVars - 2; 
idx_slack_trsc  = totalVars - 1;
idx_v80S        = totalVars;

row_trans_cap = currentRow; row_trsc_cap = currentRow + 1; 
csense(row_trans_cap:row_trsc_cap) = ['L'; 'L'];
A(row_trans_cap, idx_slack_trans) = -1; 
A(row_trsc_cap, idx_slack_trsc) = -1;
currentRow = currentRow + 2;

% --- Apply 80S Assembly Constraints ---
if numFamilies > 0
    avg_mw_scale = 25000; 
    for f = 1:length(uniqueFamilies)
        fam = char(uniqueFamilies(f));
        if isKey(family_map_k_ml, fam)
            cols = family_map_k_ml(fam);
            mws = family_map_mw(fam);
            
            if startsWith(fam, 'RPP1') || startsWith(fam, 'RPP2'), S_f = 2; else, S_f = 1; end
            
            for i = 1:length(cols)
                A(currentRow, cols(i)) = avg_mw_scale / mws(i);
            end
            A(currentRow, idx_v80S) = -S_f * params.sigma * avg_mw_scale;
            b(currentRow) = 0; csense(currentRow) = 'E';
            currentRow = currentRow + 1;
        end
    end
end

% Trim matrices to exact size
actualRows = currentRow - 1;
A = A(1:actualRows, :); b = b(1:actualRows); csense = csense(1:actualRows);

%% 6. Objective Vector (MILP_c) & Variable Types
MILP_c = zeros(totalVars, 1);

idx_pool_ml = find(strcmp(model.rxns, 'prot_pool_exchange'));
if ~isempty(idx_pool_ml), MILP_c(idx_pool_ml) = -1; end 

idx_prot_sinks = find(startsWith(model.rxns, 'sink_prot_'));
MILP_c(idx_prot_sinks) = 100; 

idx_mrna_sinks = find(startsWith(model.rxns, 'sink_mRNA_'));
MILP_c(idx_mrna_sinks) = 0; 

gene_vartype = ['C';'C';'C';'C'; 'B';'B';'B';'B'; 'C';'C'; 'B']; 
add_vartype = [repmat(gene_vartype, numGenesReal, 1); 'C'; 'C'; 'C'];
full_vartype = [repmat('C', nRxns, 1); add_vartype];

for g = 1:numGenesReal
    if native_rxn_cols(g) > 0, MILP_c(native_rxn_cols(g)) = 10000; end
    baseCol = nRxns + (g-1)*11; 
    MILP_c(baseCol + 9) = 1000000; MILP_c(baseCol + 10) = 1000000;      
end

MILP_c(idx_slack_trans) = 10000; 
MILP_c(idx_slack_trsc) = 10000;
MILP_c(idx_v80S) = 0; 

add_lb = [zeros(11*numGenesReal,1); 0; 0; 0]; 
add_ub = [Inf(11*numGenesReal,1); Inf; Inf; Inf];

%% 7. Package Everything inside the COBRA model
model.S = A; model.b = b; model.csense = csense; model.c = MILP_c;
model.lb = [model.lb; add_lb]; model.ub = [model.ub; add_ub]; 
model.vartype = full_vartype; 

model.ML_data.row_trans_cap = row_trans_cap;
model.ML_data.row_trsc_cap = row_trsc_cap;
model.ML_data.row_ribo_demand = row_ribo_demand;
model.ML_data.ml_rxn_cols = ml_rxn_cols; 
model.ML_data.mRNA_rxn_cols = mRNA_rxn_cols; 
model.ML_data.l_prot_arr = l_prot_arr;
model.ML_data.l_mRNA_arr = l_mRNA_arr;
model.ML_data.MW_arr = MW_arr;
model.ML_data.idx_slack_trans = idx_slack_trans;
model.ML_data.realGeneIdx = realGeneIdx;
model.ML_data.is_ribo_arr = is_ribo_arr; 
model.ML_data.idx_v80S = idx_v80S; 

fprintf('--- HyTT Model Successfully Compiled with Biophysical Locks ---\n');
end