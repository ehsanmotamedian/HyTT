function model = addRibosomalUsageReactions(model)
% addRibosomalUsageReactions - Incorporates core ribosomal genes into the metabolic network as mass-occupying entities.
%
% USAGE:
%   model = addRibosomalUsageReactions(model)
%
% INPUTS:
%   model               COBRA enzyme-constrained genome-scale metabolic model (ecGEM).
%
% OUTPUTS:
%   model               Expanded metabolic model containing 137 core ribosomal subunits and their usage reactions.
%
% .. Author: - Ehsan Motamedian, July 2026

% Read the ribosomal genes, UniProt IDs, and molecular weights from the attached Excel file
T = readtable('RibosomalGenes.xlsx', 'Sheet', 'Sheet1', 'VariableNamingRule', 'preserve');
% Extract the lists
ribosomal_genes = T.ORF(1:end);  % Skip header row
ribosomal_uniprots = T.("UniProt AC")(1:end);  % Skip header row
ribosomal_mw = T.("Molecular Weight (Da)")(1:end);  % Da
for i = 1:length(ribosomal_genes)
    gene = ribosomal_genes{i};
    uniprot = ribosomal_uniprots{i};  % Get the corresponding UniProt ID
    mw = ribosomal_mw(i);  % Get the molecular weight in kDa
    prot_id = ['prot_' uniprot];
    rxnID = ['usage_' prot_id];
    reactionFormula = [prot_id ' <=> prot_pool'];
    % Check if the reaction already exists to avoid warning
    if ismember(rxnID, model.rxns)
        disp(['Reaction ' rxnID ' already exists, skipping.']);
        continue;
    end
    model = addReaction(model, rxnID, ...
        'reactionFormula', reactionFormula, ...
        'reversible', true, ...
        'lowerBound', -1000, ...
        'upperBound', 0, ...
        'objective', 0, ...
        'subSystem', 'Cytoplasmic Ribosome', ...
        'geneRule', gene);
    % Manually set notes and references if the fields exist
    if isfield(model, 'rxnNotes')
        model.rxnNotes{end} = '';
    end
    if isfield(model, 'rxnReferences')
        model.rxnReferences{end} = '';
    end
    % Add to ec structure if not already present
    if ~ismember(gene, model.ec.genes)
        model.ec.genes{end+1} = gene;
        model.ec.enzymes{end+1} = uniprot;
        model.ec.mw(end+1) = mw;

        % Check and add to optional fields if they exist; otherwise initialize
        if isfield(model.ec, 'kcat')
            model.ec.kcat(end+1) = NaN;
        else
            model.ec.kcat = NaN(1, length(model.ec.genes));  % Initialize if not present (adjust size if needed)
        end

        if isfield(model.ec, 'concs')
            model.ec.concs(end+1) = NaN;
        else
            model.ec.concs = NaN(1, length(model.ec.genes));  % Initialize if not present
        end

        if isfield(model.ec, 'source')
            model.ec.source{end+1} = 'added for ribosomal';
        else
            model.ec.source = cell(1, length(model.ec.genes));  % Initialize if not present
            model.ec.source{end} = 'added for ribosomal';
        end
    end
end
end