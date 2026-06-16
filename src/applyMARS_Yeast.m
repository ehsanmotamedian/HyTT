function [constant_adj, beta] = applyMARS_Yeast(geneRow)
% =========================================================================
% ARCHITECTURAL NOTE ON THE MARS EQUATION IMPLEMENTATION:
%
% The complete MARS equation consists of two components:
%   1. A "Static" component: Based on fixed sequence features (e.g., CAI, CDS_AT).
%   2. A "Dynamic" component: Based on the mRNA Expression Level.
%
% Because the mRNA Expression Level is an active decision variable solved dynamically
% by the MILP optimizer, its corresponding hinge functions CANNOT be
% evaluated statically in this script.
%
% Therefore, this function ONLY calculates the combined static penalty
% (Intercept + non-flux hinges). The dynamic hinges for the mRNA Expression Level
% (e.g., the knots at 8104.54 and 9011.82) are structurally formulated
% as Big-M constraints directly inside the 'buildHyTT.m' file.
% =========================================================================

% 1. Extract features based on your exact dataset column names
cds_at = extractFeat(geneRow, {'CDS-AT', 'CDS_AT'});
cds_a  = extractFeat(geneRow, {'CDS-A', 'CDS_A'});
cds_aa = extractFeat(geneRow, {'CDS-AA', 'CDS_AA'});
cai    = extractFeat(geneRow, {'CAI', 'CAII'});
ox_ratio = extractFeat(geneRow, {'Oxidative_Ratio'});
a_pct  = extractFeat(geneRow, {'A_%', 'A_percent'});
l_pct  = extractFeat(geneRow, {'L_%', 'L_percent'});
m_pct  = extractFeat(geneRow, {'M_%', 'M_percent'});
cys_per_len = extractFeat(geneRow, {'Cysperproteinlength'});
mrna_prot_half = extractFeat(geneRow, {'mRNAperPRoteinhalflife(min/h)'});

% 2. The FINAL MARS Equation (Static Part: Intercept + 14 non-flux Hinges)
fixed = 0.05940566 ...
    - 0.05721368  * max(0, cds_at - 4.475043) ...
    + 0.04234268  * max(0, cds_at - 4.669261) ...
    - 0.01443783  * max(0, 7.349246 - cds_at) ...
    + 0.01506918  * max(0, cds_at - 7.349246) ...
    + 0.00109021  * max(0, 29.6332 - cds_a) ...
    - 0.001580331 * max(0, 10.09306 - cds_aa) ...
    + 0.02869147  * max(0, cai - 0.8163426) ...
    + 0.2429152   * max(0, 0.02023592 - ox_ratio) ...
    + 0.001029817 * max(0, a_pct - 7.746479) ...
    - 0.002351672 * max(0, a_pct - 10.5948) ...
    - 0.0005151684* max(0, 8.333333 - l_pct) ...
    + 0.01481261  * max(0, 0.6451613 - m_pct) ...
    - 3.175094    * max(0, 0.001831502 - cys_per_len) ...
    + 0.0003695132* max(0, 7.2 - mrna_prot_half);

constant_adj = fixed;

% 3. Beta coefficients for the dynamic mRNA flux hinges (H1 to H4)
beta = [1.440092e-05, -1.482295e-05, -0.001767903, 0.001581002];
end

% --- Helper Function ---
function val = extractFeat(row, colNames)
val = 0; % Default fallback to avoid NaNs
for i = 1:length(colNames)
    if ismember(colNames{i}, row.Properties.VariableNames)
        temp = row.(colNames{i});
        if ~isempty(temp) && ~isnan(temp)
            val = temp;
            return;
        end
    end
end
end