function [csC_perAnimal, csMean_grand, h, animalIDs, csMean_byGroup] = analyzeMotifBetaCosineSimilarity_acrossAnimals(glmRezC, glmLabelC, varargin)
%ANALYZEMOTIFBETACOSINESIMILARITY_ACROSSANIMALS  Motif-pairwise beta cosine similarity, per animal and across animals.
%
% SYNOPSIS
%   [csC_perAnimal, csMean_grand, h, animalIDs] = ...
%       analyzeMotifBetaCosineSimilarity_acrossAnimals(glmRezC, glmLabelC, ...)
%
% DESCRIPTION
%   For each session, each motif's tuning profile is its column of beta
%   (glmRezC{i,j}.beta is [P predictors x K motifs]). Cosine similarity is
%   computed between all motif-pairs (columns) within that session,
%   giving a [K x K] matrix. Per animal, these session-level [K x K]
%   matrices are averaged elementwise across that animal's valid sessions
%   -> ANIMAL-LEVEL mean matrix. The 9 animal-level matrices are then
%   averaged elementwise across animals -> the single GRAND-MEAN [K x K]
%   matrix.
%
%   Four figures are generated:
%     1) figIndiv        : one heatmap per animal (session-averaged),
%                          laid out in a roughly-square subplot grid
%                          (3x3 for 9 animals), sharing one colorbar.
%                          Each animal's matrix is reordered using the
%                          SAME leaf order as figMeanReordered below, so
%                          motif blocks line up visually across animals.
%     2) figMean          : single heatmap of the grand-mean matrix,
%                          natural motif order.
%     3) figMeanReordered : the grand-mean matrix reordered by
%                          hierarchical clustering.
%     4) figDendro        : the dendrogram for that same clustering.
%
%   Clustering (and therefore the leaf order used everywhere) is computed
%   ONLY on the grand-mean matrix -- it is not meaningful (or computed)
%   on individual animal-level matrices; individual heatmaps simply
%   inherit that single ordering for visual consistency.
%
% INPUTS
%   glmRezC   : [nAnimals x nSessions] cell array; glmRezC{i,j} is either
%               [] (missing session) or a struct with field 'beta'
%               ([P x K] ridge coefficients).
%   glmLabelC : cell array of label strings containing the animal ID
%               (e.g. 'm1045' or 'm1045_122424...'), same convention as
%               used in collect_cvR2_perAnimal.m (nAnimals x nSessions or
%               nAnimals x 1).
%
% NAME-VALUE ARGS
%   'colormapName'          : colormap for all heatmaps. Default: 'parula'.
%   'linkageMethod'         : linkage method for the dendrogram (passed to
%                             MATLAB's linkage.m). Default: 'average'.
%   'visible'               : 'on'/'off' for all 4 figures. Default: 'on'.
%   'figureScaleFactorIndiv'   : figure size multiplier for figIndiv (the
%                                3x3 grid needs more room). Default: 2.
%   'figureScaleFactorSummary' : figure size multiplier for figMean,
%                                figMeanReordered, and figDendro.
%                                Default: 1.2.
%   'figSaveDir'            : folder to save all figures as PDFs. Empty
%                             (default) = don't save.
%   'figSaveKeyword'        : extra tag inserted into saved filenames.
%   'groupDefs'             : OPTIONAL scalar struct defining animal
%                             groups, e.g.:
%                               groupDefs.fast = {'m1044','m1045','m1092','m1094'};
%                               groupDefs.slow = {'m1048','m1049','m1613','m1859','m1873'};
%                             Field names become group labels; values are
%                             cellstr of animal IDs (must not overlap
%                             across groups). When supplied, for EACH
%                             group this additionally generates:
%                               - a mean-CS heatmap, reordered using the
%                                 SAME (grand-mean) leaf order as figMean
%                                 Reordered and figIndiv -- so every
%                                 heatmap in this function is directly
%                                 comparable motif-for-motif
%                               - a dendrogram from that group's OWN
%                                 clustering (independent of the
%                                 grand-mean clustering -- this is the
%                                 only place group-specific structure is
%                                 shown; it intentionally will not match
%                                 the shared heatmap ordering)
%                             as separate figures (h.figHeatmap_byGroup.(name),
%                             h.figDendro_byGroup.(name)). figIndiv subplot
%                             titles are also tagged with each animal's
%                             group. Default: struct() (no groups; only
%                             the grand-mean-level figures are produced).
%
% OUTPUTS
%   csC_perAnimal  : [nAnimals x 1] cell array; csC_perAnimal{i} is the
%                    [K x K] animal-level mean cosine similarity matrix
%                    (empty if that animal had no valid sessions).
%   csMean_grand   : [K x K] grand-mean cosine similarity matrix, averaged
%                    across animals' csC_perAnimal entries.
%   h              : struct with figure/axes handles, the linkage result Z,
%                    the clustering leaf order, per-group figure/linkage
%                    fields (if 'groupDefs' was supplied), and opt.
%   animalIDs      : [nAnimals x 1] cellstr of animal ID labels.
%   csMean_byGroup : scalar struct, one field per group in 'groupDefs',
%                    each holding that group's [K x K] mean cosine
%                    similarity matrix (empty if 'groupDefs' not supplied
%                    or a group had no valid animals). Appended as a 5th
%                    output so existing 4-output calls remain valid.
%
% NOTE
%   Cosine similarity here is computed independently within each session
%   (over that session's own predictor set P), so it is NOT sensitive to
%   P differing slightly across sessions (e.g. a session missing a
%   continuous behavioral predictor like whisker). It IS sensitive to K
%   (motif count) being consistent across all sessions/animals being
%   averaged -- this is checked and will error with a clear message if
%   violated. Also note: matrices are averaged directly (no Fisher-z-type
%   transform), which is standard for this kind of descriptive/visualization
%   summary but worth knowing if a stricter statistical average across
%   animals is ever needed elsewhere.
%
% EXAMPLE
%   groupDefs.fast = {'m1044','m1045','m1092','m1094'};
%   groupDefs.slow = {'m1048','m1049','m1613','m1859','m1873'};
%   [csC_perAnimal, csMean_grand, h, animalIDs, csMean_byGroup] = ...
%       analyzeMotifBetaCosineSimilarity_acrossAnimals(glmRezC, glmLabelC, ...
%           'groupDefs', groupDefs, ...
%           'figSaveDir', figSaveDir, 'figSaveKeyword', 'betaCosineSim');
%
% See also: collect_cvR2_perAnimal

% -------- parse args --------
p = inputParser;
p.addParameter('colormapName', 'parula', @(s) ischar(s) || isstring(s));
p.addParameter('linkageMethod', 'average', @(s) ischar(s) || isstring(s));
p.addParameter('visible', 'on', @(s) any(strcmpi(s, {'on','off'})));
p.addParameter('figureScaleFactorIndiv', 2, @(x) isnumeric(x) && x>0);
p.addParameter('figureScaleFactorSummary', 1.2, @(x) isnumeric(x) && x>0);
p.addParameter('figSaveDir', {}, @(x) isempty(x) || ischar(x) || isstring(x) || iscell(x));
p.addParameter('figSaveKeyword', '', @(s) ischar(s) || isstring(s));
p.addParameter('groupDefs', struct(), @(x) isstruct(x) && isscalar(x));
p.parse(varargin{:});
opt = p.Results;

% normalize figSaveDir
figSaveDir = opt.figSaveDir;
if iscell(figSaveDir)
    if isempty(figSaveDir), figSaveDir = ""; else, figSaveDir = string(figSaveDir{1}); end
else
    figSaveDir = string(figSaveDir);
end
figSaveKeyword = string(opt.figSaveKeyword);

% -------- 1) per-session cosine similarity, averaged per animal --------
[nAnimals, nSessions] = size(glmRezC);

csC_perAnimal = cell(nAnimals, 1);
K = [];

for i = 1:nAnimals
    stack = [];
    for j = 1:nSessions
        s = glmRezC{i,j};
        if isempty(s) || ~isstruct(s) || ~isfield(s,'beta') || isempty(s.beta)
            continue;
        end
        Cij = motif_cosine_similarity_(s.beta);

        if isempty(K)
            K = size(Cij,1);
        else
            assert(size(Cij,1)==K, ...
                'analyzeMotifBetaCosineSimilarity_acrossAnimals:motifCountMismatch', ...
                'Motif count (K) differs across sessions/animals (found %d and %d) -- cannot average elementwise.', ...
                K, size(Cij,1));
        end
        stack = cat(3, stack, Cij);
    end

    if isempty(stack)
        csC_perAnimal{i} = [];
    else
        csC_perAnimal{i} = mean(stack, 3);
    end
end

validAnimalI = find(~cellfun(@isempty, csC_perAnimal));
assert(~isempty(validAnimalI), ...
    'No animal had any valid session with a beta matrix -- nothing to compute.');

csMean_grand = mean(cat(3, csC_perAnimal{validAnimalI}), 3);

animalIDs = extract_animal_ids_(glmLabelC, nAnimals);

% -------- optional animal grouping (e.g. fast/slow learners) --------
groupNames = fieldnames(opt.groupDefs);
hasGroups = ~isempty(groupNames);

animalGroupTag = repmat({''}, nAnimals, 1);
if hasGroups
    for gi = 1:numel(groupNames)
        for gj = gi+1:numel(groupNames)
            ov = intersect(opt.groupDefs.(groupNames{gi}), opt.groupDefs.(groupNames{gj}));
            assert(isempty(ov), ...
                'analyzeMotifBetaCosineSimilarity_acrossAnimals:groupOverlap', ...
                'Animal(s) %s appear in both group "%s" and group "%s".', ...
                strjoin(ov, ', '), groupNames{gi}, groupNames{gj});
        end
    end
    for gi = 1:numel(groupNames)
        memberIDs = opt.groupDefs.(groupNames{gi});
        animalGroupTag(ismember(animalIDs, memberIDs)) = groupNames(gi);
    end
end

% -------- shared style --------
h = struct();
h.opt = opt;
climVals = [-1 1];

% -------- clustering on the GRAND-MEAN matrix only (never on individual animals) --------
D = 1 - csMean_grand;
D = (D + D') / 2;              % force exact symmetry
D(1:K+1:end) = 0;              % force exact zero diagonal (required by squareform)

Z = linkage(squareform(D), opt.linkageMethod);

h.figDendro = figure('Color','w', 'Visible', opt.visible);
set(h.figDendro, 'Units', 'normalized');
pos = get(h.figDendro, 'Position');
pos(3:4) = pos(3:4) * opt.figureScaleFactorSummary;
set(h.figDendro, 'Position', pos);

h.axDendro = axes('Parent', h.figDendro);
axes(h.axDendro); %#ok<LAXES>
[~, ~, leafOrder] = dendrogram(Z, 0, 'Orientation', 'top');
xlabel(h.axDendro, 'Motif (leaf order)');
ylabel(h.axDendro, 'Distance (1 - cosine similarity)');
title(h.axDendro, sprintf('\\beta cosine similarity -- dendrogram (%s)', opt.linkageMethod), 'Interpreter', 'tex');

h.Z = Z;
h.leafOrder = leafOrder;

%% ===================== FIGURE 1: per-animal grid (reordered by grand-mean clustering) =====================
nRows = ceil(sqrt(nAnimals));
nCols = ceil(nAnimals / nRows);

h.figIndiv = figure('Color','w', 'Visible', opt.visible);
set(h.figIndiv, 'Units', 'normalized');
pos = get(h.figIndiv, 'Position');
pos(3:4) = pos(3:4) * opt.figureScaleFactorIndiv;
set(h.figIndiv, 'Position', pos);

h.axIndiv = gobjects(nAnimals,1);
for i = 1:nAnimals
    h.axIndiv(i) = subplot(nRows, nCols, i, 'Parent', h.figIndiv);
    if isempty(csC_perAnimal{i})
        axis(h.axIndiv(i), 'off');
        title(h.axIndiv(i), sprintf('%s (no data)', animalIDs{i}), 'Interpreter', 'none');
        continue;
    end
    imagesc(h.axIndiv(i), csC_perAnimal{i}(leafOrder, leafOrder), climVals);
    axis(h.axIndiv(i), 'square');
    colormap(h.axIndiv(i), opt.colormapName);
    set(h.axIndiv(i), 'XTick', [], 'YTick', []);
    titleStr = animalIDs{i};
    if ~isempty(animalGroupTag{i})
        titleStr = sprintf('%s (%s)', titleStr, animalGroupTag{i});
    end
    title(h.axIndiv(i), titleStr, 'Interpreter', 'none');
end

sgtitle(h.figIndiv, '\beta cosine similarity by animal (session-averaged, ordered by grand-mean clustering)', 'Interpreter', 'tex');

% one shared colorbar, positioned outside the subplot grid
lastValid = find(~cellfun(@isempty, csC_perAnimal), 1, 'last');
if ~isempty(lastValid)
    cb = colorbar(h.axIndiv(lastValid));
    cb.Position = [0.93 0.11 0.02 0.8];
    h.cbIndiv = cb;
end

%% ===================== FIGURE 2: grand-mean heatmap (natural motif order) =====================
h.figMean = figure('Color','w', 'Visible', opt.visible);
set(h.figMean, 'Units', 'normalized');
pos = get(h.figMean, 'Position');
pos(3:4) = pos(3:4) * opt.figureScaleFactorSummary;
set(h.figMean, 'Position', pos);

h.axMean = axes('Parent', h.figMean);
imagesc(h.axMean, csMean_grand, climVals);
axis(h.axMean, 'square');
colormap(h.axMean, opt.colormapName);
colorbar(h.axMean);
if K <= 30
    xticks(h.axMean, 1:K); yticks(h.axMean, 1:K);
end
xlabel(h.axMean, 'Motif index');
ylabel(h.axMean, 'Motif index');
title(h.axMean, '\beta cosine similarity -- averaged across animals', 'Interpreter', 'tex');

%% ===================== FIGURE 3: grand-mean heatmap, reordered by clustering =====================
h.figMeanReordered = figure('Color','w', 'Visible', opt.visible);
set(h.figMeanReordered, 'Units', 'normalized');
pos = get(h.figMeanReordered, 'Position');
pos(3:4) = pos(3:4) * opt.figureScaleFactorSummary;
set(h.figMeanReordered, 'Position', pos);

h.axMeanReordered = axes('Parent', h.figMeanReordered);
imagesc(h.axMeanReordered, csMean_grand(leafOrder, leafOrder), climVals);
axis(h.axMeanReordered, 'square');
colormap(h.axMeanReordered, opt.colormapName);
colorbar(h.axMeanReordered);
set(h.axMeanReordered, 'XTick', 1:K, 'XTickLabel', string(leafOrder), ...
                        'YTick', 1:K, 'YTickLabel', string(leafOrder));
xtickangle(h.axMeanReordered, 90);
xlabel(h.axMeanReordered, 'Motif (reordered)');
ylabel(h.axMeanReordered, 'Motif (reordered)');
title(h.axMeanReordered, '\beta cosine similarity -- reordered by clustering', 'Interpreter', 'tex');

%% ===================== FIGURE(S) 4+: per-group mean heatmap + dendrogram =====================
% Computed only if 'groupDefs' was supplied. Each group gets its OWN
% clustering (a group's motif structure need not match the grand-mean
% clustering), but -- consistent with the rest of this function --
% heatmap and dendrogram are always kept as separate figures.
csMean_byGroup = struct();
h.figHeatmap_byGroup = struct();
h.figDendro_byGroup  = struct();
h.Z_byGroup          = struct();
h.leafOrder_byGroup  = struct();

if hasGroups
    for gi = 1:numel(groupNames)
        gName = groupNames{gi};
        memberIDs = opt.groupDefs.(gName);
        idx = find(ismember(animalIDs, memberIDs));
        idxValid = idx(~cellfun(@isempty, csC_perAnimal(idx)));

        if isempty(idxValid)
            warning('analyzeMotifBetaCosineSimilarity_acrossAnimals:emptyGroup', ...
                'Group "%s" has no animals with valid data -- skipping.', gName);
            csMean_byGroup.(gName) = [];
            continue;
        end

        csMean_g = mean(cat(3, csC_perAnimal{idxValid}), 3);
        csMean_byGroup.(gName) = csMean_g;

        % group-specific clustering (independent of the grand-mean clustering)
        Dg = 1 - csMean_g;
        Dg = (Dg + Dg') / 2;
        Dg(1:K+1:end) = 0;
        Zg = linkage(squareform(Dg), opt.linkageMethod);

        % -- dendrogram figure --
        figDg = figure('Color','w', 'Visible', opt.visible);
        set(figDg, 'Units', 'normalized');
        posg = get(figDg, 'Position');
        posg(3:4) = posg(3:4) * opt.figureScaleFactorSummary;
        set(figDg, 'Position', posg);
        axDg = axes('Parent', figDg);
        axes(axDg); %#ok<LAXES>
        [~, ~, leafOrderG] = dendrogram(Zg, 0, 'Orientation', 'top');
        xlabel(axDg, 'Motif (leaf order)');
        ylabel(axDg, 'Distance (1 - cosine similarity)');
        title(axDg, sprintf('\\beta cosine similarity -- dendrogram (%s, %s)', opt.linkageMethod, gName), 'Interpreter', 'tex');

        % -- mean heatmap figure, reordered by the SHARED (grand-mean) leaf
        %    order -- same ordering as figMean_Reordered and figIndiv, so
        %    every heatmap in this function lines up motif-for-motif --
        figHm = figure('Color','w', 'Visible', opt.visible);
        set(figHm, 'Units', 'normalized');
        posg = get(figHm, 'Position');
        posg(3:4) = posg(3:4) * opt.figureScaleFactorSummary;
        set(figHm, 'Position', posg);
        axHm = axes('Parent', figHm);
        imagesc(axHm, csMean_g(leafOrder, leafOrder), climVals);
        axis(axHm, 'square');
        colormap(axHm, opt.colormapName);
        colorbar(axHm);
        set(axHm, 'XTick', 1:K, 'XTickLabel', string(leafOrder), ...
                  'YTick', 1:K, 'YTickLabel', string(leafOrder));
        xtickangle(axHm, 90);
        xlabel(axHm, 'Motif (reordered)');
        ylabel(axHm, 'Motif (reordered)');
        title(axHm, sprintf('\\beta cosine similarity -- averaged across %s animals (reordered)', gName), 'Interpreter', 'tex');

        h.figDendro_byGroup.(gName)  = figDg;
        h.figHeatmap_byGroup.(gName) = figHm;
        h.Z_byGroup.(gName)          = Zg;
        h.leafOrder_byGroup.(gName)  = leafOrderG;
    end
end

%% -------- save (optional) --------
if strlength(figSaveDir) > 0
    if ~isfolder(figSaveDir)
        mkdir(figSaveDir);
    end
    dateStr = char(datetime("today","Format","MMddyy"));

    saveFig_ = @(fig, tag) print(fig, fullfile(figSaveDir, ...
        strjoin(strings_nonempty_({tag, figSaveKeyword, dateStr}), "_")), ...
        '-dpdf', '-painters', '-bestfit');

    saveFig_(h.figIndiv,         "betaCosineSim_perAnimal");
    saveFig_(h.figMean,          "betaCosineSim_meanAcrossAnimals");
    saveFig_(h.figMeanReordered, "betaCosineSim_meanAcrossAnimals_reordered");
    saveFig_(h.figDendro,        "betaCosineSim_dendrogram");

    if hasGroups
        for gi = 1:numel(groupNames)
            gName = groupNames{gi};
            if isfield(h.figHeatmap_byGroup, gName) && ~isempty(h.figHeatmap_byGroup.(gName))
                saveFig_(h.figHeatmap_byGroup.(gName), sprintf("betaCosineSim_mean_%s", gName));
            end
            if isfield(h.figDendro_byGroup, gName) && ~isempty(h.figDendro_byGroup.(gName))
                saveFig_(h.figDendro_byGroup.(gName), sprintf("betaCosineSim_dendrogram_%s", gName));
            end
        end
    end
end

end % function


% ===== helper: cosine similarity between columns (motifs) of beta =====
function C = motif_cosine_similarity_(beta)
% beta: [P x K]. Returns [K x K] cosine similarity between columns.
nrm = vecnorm(beta, 2, 1);      % 1 x K
nrm(nrm == 0) = 1;              % guard: zero-norm column -> stays zero vector (similarity 0)
Bn = beta ./ nrm;
C = Bn' * Bn;
C = min(max(C, -1), 1);         % clip tiny floating-point overshoot outside [-1,1]
end


% ===== helper: pull animal ID (e.g. 'm1045') from glmLabelC row =====
function animalIDs = extract_animal_ids_(glmLabelC, nAnimals)
animalIDs = cell(nAnimals, 1);

for i = 1:nAnimals
    if size(glmLabelC,1) < i
        rowLabels = {};
    elseif size(glmLabelC,2) > 1
        rowLabels = glmLabelC(i,:);
    else
        rowLabels = glmLabelC(i,1);
    end

    idFound = '';
    for j = 1:numel(rowLabels)
        lbl = rowLabels{j};
        if (ischar(lbl) || isstring(lbl)) && ~isempty(lbl)
            tok = regexp(char(lbl), 'm\d{4}', 'match', 'once');
            if ~isempty(tok)
                idFound = tok;
                break;
            end
        end
    end

    if isempty(idFound)
        idFound = sprintf('Animal%02d', i);
    end
    animalIDs{i} = idFound;
end
end


% ===== helper: drop empty strings before joining filename parts =====
function out = strings_nonempty_(parts)
parts = string(parts);
out = parts(strlength(parts) > 0);
end