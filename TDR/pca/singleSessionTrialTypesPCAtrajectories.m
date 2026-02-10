function hFig = singleSessionTrialTypesPCAtrajectories(s, varargin)
% SINGLESESSIONTRIALTYPESPCAtrajectories
%   Plot trial-type averaged PCA trajectories for a single session (one entry of pcaRezC).
%
% Usage
%   singleSessionTrialTypesPCAtrajectories(s)
%   singleSessionTrialTypesPCAtrajectories(s, 'PCs', [1 2 3], ...
%       'trialTypes', {'hit','cr','fa','miss'}, 'PCscoreType', 'global')
%
% Inputs
%   s : struct (one session) produced by your PASS1/2 pipeline. Must contain:
%       s.trI  : trial masks (fields like goI, nogoI, hitI, missI, crI, faI)
%       s.pca  : PCA results with fields depending on PCscoreType:
%           - global       -> s.pca.globalScore_traj          [Nsel x Tsel x nPC]
%           - expertGlobal -> s.pca.expertGlobalScore_traj    [Nsel x Tsel x nPC]
%           - expertWithin -> s.pca.expertWithinMouseScore_traj [Nsel x Tsel x nPC]
%
% Name-Value
%   'PCs'              : vector of PC indices (default [1 2 3])
%   'trialTypes'       : cellstr of {'go','nogo','hit','miss','cr','fa'} (case-insensitive)
%                        (default {'hit','miss','cr','fa'})
%   'PCscoreType'      : 'global'|'expertGlobal'|'expertWithin' (default 'global')
%   'trialTypeColorC'  : Nx2 cell array: col1 trialType string, col2 = 1x3 RGB (0..1)
%                        default:
%                          {'go',[0 0 1];
%                           'nogo',[1 0 0];
%                           'hit',[46 49 146]/255;
%                           'cr',[193 39 45]/255}
%                        miss uses hit color with alpha=0.4; fa uses cr color with alpha=0.4
%   'lineWidth'        : line width (default 2)
%   'alphaMain'        : alpha for main types (default 0.8)
%   'alphaAux'         : alpha for miss/fa (default 0.4)
%   'doLegend'         : true/false (default true)
%   'doGrid'           : true/false (default true)
%
% Output
%   hFig : figure handle

% -------------------- parse --------------------
p = inputParser;
p.addParameter('PCs', [1 2 3], @(x)isnumeric(x)&&isvector(x)&&numel(x)>=1);
p.addParameter('trialTypes', {'hit','miss','cr','fa'}, @(c)iscellstr(c) || isstring(c));
p.addParameter('PCscoreType', 'global', @(s)ischar(s)||isstring(s));
p.addParameter('trialTypeColorC', defaultTrialTypeColorC(), @(c)iscell(c)&&size(c,2)==2);
p.addParameter('lineWidth', 2, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('alphaMain', 0.8, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('alphaAux',  0.4, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('doLegend', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('doGrid', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
opt = p.Results;

PCs = opt.PCs(:)';

% -------------------- validate PCs --------------------
if ~isfield(s,'pca') || isempty(s.pca)
    error('singleSessionTrialTypesPCAtrajectories:MissingPCA', 's.pca is missing/empty.');
end

[scoreTraj, scoreTrajName] = getScoreTrajByType(s, opt.PCscoreType);

if isempty(scoreTraj) || ndims(scoreTraj) ~= 3
    error('singleSessionTrialTypesPCAtrajectories:BadScoreTraj', ...
        'Selected score trajectory is missing or not [N x T x nPC]. Field: %s', scoreTrajName);
end

nPCavail = size(scoreTraj, 3);
if any(PCs < 1) || any(PCs > nPCavail)
    error('singleSessionTrialTypesPCAtrajectories:PCOutOfRange', ...
        'Requested PCs must be within 1..%d (requested: %s).', nPCavail, mat2str(PCs));
end
if numel(PCs) > 3
    error('singleSessionTrialTypesPCAtrajectories:TooManyPCs', ...
        'This plotting function supports 1, 2, or 3 PCs. Requested: %d', numel(PCs));
end

% -------------------- trialTypes parse/validate --------------------
trialTypes = cellstr(string(opt.trialTypes));
trialTypes = lower(strtrim(trialTypes));
allowed = {'go','nogo','hit','miss','cr','fa'};
bad = setdiff(unique(trialTypes), allowed);
if ~isempty(bad)
    error('singleSessionTrialTypesPCAtrajectories:BadTrialTypes', ...
        'trialTypes contained invalid entries: %s. Allowed: %s', strjoin(bad,','), strjoin(allowed,','));
end

% -------------------- build colors (with miss/fa rules) --------------------
colorMap = containers.Map('KeyType','char','ValueType','any');
ttC = opt.trialTypeColorC;
for i = 1:size(ttC,1)
    k = lower(string(ttC{i,1}));
    v = ttC{i,2};
    if isstring(v) || ischar(v), v = str2num(v); end %#ok<ST2NM>
    validateattributes(v, {'numeric'},{'size',[1,3]});
    colorMap(char(k)) = v;
end

% miss defaults to hit color; fa defaults to cr color if not explicitly given
if ~isKey(colorMap,'miss') && isKey(colorMap,'hit')
    colorMap('miss') = colorMap('hit');
end
if ~isKey(colorMap,'fa') && isKey(colorMap,'cr')
    colorMap('fa') = colorMap('cr');
end

% -------------------- figure --------------------
hdr = '';
if isfield(s,'meta') && isfield(s.meta,'header'), hdr = string(s.meta.header); end
mid = '';
if isfield(s,'meta') && isfield(s.meta,'mouseId'), mid = string(s.meta.mouseId); end

hFig = figure('Color','w'); hold on;

% -------------------- compute + plot mean trajectories --------------------
legendH = gobjects(0);
legendTxt = {};

for i = 1:numel(trialTypes)
    tt = trialTypes{i};

    mask = trialTypeMaskFromTrI(s, tt);     % logical Nsel x 1 (or [])
    if isempty(mask) || ~any(mask)
        warning('singleSessionTrialTypesPCAtrajectories:NoTrials', ...
            'No trials found for trialType=%s (session=%s). Skipping.', tt, hdr);
        continue;
    end

    traj = scoreTraj(mask, :, PCs);         % [nTr x T x nPCsel]
    mu = squeeze(nanmean(traj, 1));         % [T x nPCsel]

    if size(mu,2) == 1
        x = mu(:,1);
        alphaUse = alphaForTrialType(tt, opt.alphaMain, opt.alphaAux);
        col = colorForTrialType(tt, colorMap);
        hh = plot(x, 'LineWidth', opt.lineWidth);
        hh.Color = [col, alphaUse];

    elseif size(mu,2) == 2
        x = mu(:,1); y = mu(:,2);
        alphaUse = alphaForTrialType(tt, opt.alphaMain, opt.alphaAux);
        col = colorForTrialType(tt, colorMap);
        hh = plot(x, y, 'LineWidth', opt.lineWidth);
        hh.Color = [col, alphaUse];

    else % 3 PCs
        x = mu(:,1); y = mu(:,2); z = mu(:,3);
        alphaUse = alphaForTrialType(tt, opt.alphaMain, opt.alphaAux);
        col = colorForTrialType(tt, colorMap);
        hh = plot3(x, y, z, 'LineWidth', opt.lineWidth);
        hh.Color = [col, alphaUse];
    end

    legendH(end+1) = hh; %#ok<AGROW>
    legendTxt{end+1} = tt; %#ok<AGROW>
end

% -------------------- labels --------------------
if numel(PCs) >= 1, xlabel(sprintf('PC%d', PCs(1))); end
if numel(PCs) >= 2, ylabel(sprintf('PC%d', PCs(2))); end
if numel(PCs) == 3, zlabel(sprintf('PC%d', PCs(3))); end

titleStr = sprintf('%s | %s | %s', string(mid), string(hdr), string(opt.PCscoreType));
title(titleStr, 'Interpreter','none');

axis tight;
if opt.doGrid, grid on; end
if opt.doLegend && ~isempty(legendH)
    legend(legendH, legendTxt, 'Location','best', 'Interpreter','none');
end
view(3);

end

% ========================= helpers =========================

function ttC = defaultTrialTypeColorC()
ttC = {
    'go',   [0 0 255]./255
    'nogo', [255 0 0]./255
    'hit',  [46 49 146]./255
    'cr',   [193 39 45]./255
    };
end

function [scoreTraj, fieldName] = getScoreTrajByType(s, PCscoreType)
t = lower(string(PCscoreType));
switch t
    case "global"
        fieldName = "pca.globalScore_traj";
        if isfield(s.pca,'globalScore_traj'), scoreTraj = s.pca.globalScore_traj; else, scoreTraj = []; end
    case "expertglobal"
        fieldName = "pca.expertGlobalScore_traj";
        if isfield(s.pca,'expertGlobalScore_traj'), scoreTraj = s.pca.expertGlobalScore_traj; else, scoreTraj = []; end
    case "expertwithin"
        fieldName = "pca.expertWithinMouseScore_traj";
        if isfield(s.pca,'expertWithinMouseScore_traj'), scoreTraj = s.pca.expertWithinMouseScore_traj; else, scoreTraj = []; end
    otherwise
        error('singleSessionTrialTypesPCAtrajectories:BadPCscoreType', ...
            'PCscoreType must be one of {global, expertGlobal, expertWithin}. Got: %s', PCscoreType);
end
end

function mask = trialTypeMaskFromTrI(s, trialType)
% Return Nselx1 logical mask aligned to the rows in scoreTraj (i.e., the selected trials).
% scoreTraj is built from S.select.trialSel, so we need to map from original trial indices
% (s.trI.*I) into the local selected-trial order.
if ~isfield(s,'trI') || isempty(s.trI)
    error('singleSessionTrialTypesPCAtrajectories:MissingTrI','s.trI missing/empty.');
end
if ~isfield(s,'select') || ~isfield(s.select,'trialSel')
    error('singleSessionTrialTypesPCAtrajectories:MissingSelect','s.select.trialSel missing.');
end

tt = lower(string(trialType));

% original (full) trial mask (length = N full trials)
origMask = [];
switch tt
    case "go"
        origMask = s.trI.goI(:);
    case "nogo"
        origMask = s.trI.nogoI(:);
    case "hit"
        origMask = s.trI.hitI(:);
    case "miss"
        origMask = s.trI.missI(:);
    case "cr"
        origMask = s.trI.crI(:);
    case "fa"
        origMask = s.trI.faI(:);
end

% map to selected trials
trSel = s.select.trialSel(:);            % indices into full trial list
mask = false(numel(trSel),1);
if isempty(origMask), return; end

mask = origMask(trSel);
mask = logical(mask(:));
end

function col = colorForTrialType(tt, colorMap)
k = char(lower(string(tt)));
if isKey(colorMap, k)
    col = colorMap(k);
else
    % fallback: go->blue, nogo->red, else gray
    switch k
        case "go"
            col = [0 0 1];
        case "nogo"
            col = [1 0 0];
        otherwise
            col = [0.3 0.3 0.3];
    end
end
end

function a = alphaForTrialType(tt, alphaMain, alphaAux)
k = char(lower(string(tt)));
if any(strcmp(k, {'miss','fa'}))
    a = alphaAux;
else
    a = alphaMain;
end
end
