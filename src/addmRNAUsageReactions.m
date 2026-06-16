function model = addmRNAUsageReactions(model)
% addmRNAUsageReactions - Expands the ecGEM network to incorporate transcription layer kinetics and mRNA mass balances.
%
% USAGE:
%   model = addmRNAUsageReactions(model)
%
% INPUTS:
%   model               COBRA enzyme-constrained genome-scale metabolic model (ecGEM).
%
% OUTPUTS:
%   model               Expanded metabolic model with intracellular mRNA pools, specific usage reactions, and sink balances.
%
% .. Author: - Ehsan Motamedian, 2026

% Add mRNA_pool if not exists
if ~ismember('mRNA_pool[c]', model.mets)
    model = addMetabolite(model, 'mRNA_pool[c]', 'metName', 'mRNA pool');
end

for i = 1:length(model.ec.genes)
    gene = model.ec.genes{i};
    uniprot = model.ec.enzymes{i};
    mRNA_id = ['mRNA_' uniprot '[c]'];

    % Add metabolite if not exists
    if ~ismember(mRNA_id, model.mets)
        model = addMetabolite(model, mRNA_id, 'metName', ['mRNA ' uniprot]);
    end

    rxnID = ['usage_mRNA_' uniprot];
    reactionFormula = [mRNA_id ' <=> mRNA_pool[c]'];

    % Check if reaction exists
    if ismember(rxnID, model.rxns)
        disp(['Reaction ' rxnID ' already exists, skipping.']);
        continue;
    end

    % ظرفیت پایین باز شده تا متغیر مقادیر حقیقی بیولوژیک به خود بگیرد
    model = addReaction(model, rxnID, ...
        'reactionFormula', reactionFormula, ...
        'reversible', true, ...
        'lowerBound', -100000, ...
        'upperBound', 0, ...
        'objective', 0, ...
        'subSystem', 'mRNA Usage', ...
        'geneRule', gene);

    % اضافه کردن واکنش سینک خروجی برای حفظ بالانس جرم (ایده طلایی شما)
    sinkID = ['sink_mRNA_' uniprot];
    sinkFormula = [mRNA_id ' => '];
    if ~ismember(sinkID, model.rxns)
        model = addReaction(model, sinkID, ...
            'reactionFormula', sinkFormula, ...
            'reversible', false, ...
            'lowerBound', 0, ...
            'upperBound', 100000, ...
            'objective', 0, ...
            'subSystem', 'mRNA Sink');
    end
end
end