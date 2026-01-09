function hFig = plotGammaPETH_byTrialType(gammaC, time, trI, varargin)
% plotGammaPETH_byTrialType
% -------------------------------------------------------------------------
% Plot trial-averaged HMM state occupancy (gamma) PETHs by trial type.
%
% Figures produced:
%   go, nogo, hit, miss, cr, fa
%
% Color encodes STATE identity (consistent across all figures):
%   State 1 → Hit color
%   State 2 → Miss color
%   State 3 → CR color
%   State 4 → FA color
%
% Inputs
%   gammaC : 1×N cell, each S×T
%   time   : 1×T
%   trI    : struct with fields goI, nogoI, hitI, missI, crI, faI
%
% Name–Value options
%   'TitlePrefix' : string prepended to titles (default '')
%   'LineWidth'   : line width (default 2)
%   'YLim'        : y-limits (default [0 1])
%   'LegendLoc'   : legend location (default 'best')
%
% Output
%   hFig : struct with fields .go .nogo .hit .miss .cr .fa (figure handles)
% -------------------------------------------------------------------------

% ---------------- parse ----------------
p = inputParser;
p.addParameter('TitlePrefix', '', @(x) ischar(x) || isstring(x));
p.addParameter('LineWidth', 2, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('YLim', [0 1], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('LegendLoc', 'best', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
opt = p.Results;

% ---------------- validate ----------------
assert(iscell(gammaC) && ~isempty(gammaC), 'gammaC must be non-empty.');
S = size(gammaC{1},1);
T = size(gammaC{1},2);

time = time(:)';
assert(numel(time) == T, 'time length must match gamma T.');

N = numel(gammaC);
G = cat(3, gammaC{:}); % S×T×N

assert(S == 4, 'This plotting function assumes S=4 states (found S=%d).', S);

% ---------------- STATE COLORS (repurposed) ----------------
stateColor = [
    251 180 174;   % State 1 (Hit color)
    179 205 227;   % State 2 (Miss color)
    254 217 166;   % State 3 (CR color)
    222 203 228    % State 4 (FA color)
    ] ./ 255;

% ---------------- output struct (safe field names) ----------------
hFig = struct('go',[],'nogo',[],'hit',[],'miss',[],'cr',[],'fa',[]);

% ---------------- plot all requested groups ----------------
plotOne('goI',   'Go',    'go');
plotOne('nogoI', 'No-Go', 'nogo');

plotOne('hitI',  'Hit',  'hit');
plotOne('missI', 'Miss', 'miss');
plotOne('crI',   'CR',   'cr');
plotOne('faI',   'FA',   'fa');

% ================= nested helper =================
    function plotOne(fieldName, labelName, keyName)
        if ~isfield(trI, fieldName), return; end

        idx = trI.(fieldName);
        if isempty(idx), return; end
        idx = idx(:); % Nx1 logical

        if numel(idx) ~= N
            error('trI.%s must have length N=%d, got %d.', fieldName, N, numel(idx));
        end
        if ~any(idx), return; end

        Gsub = G(:,:,idx);                      % S×T×Ng
        mu   = mean(Gsub, 3, 'omitnan');         % S×T

        h = figure('Color','w'); hold on;
        for s = 1:S
            plot(time, mu(s,:), ...
                'LineWidth', opt.LineWidth, ...
                'Color', stateColor(s,:));
        end

        xlabel('Time (s)', 'FontSize', get(gca,'FontSize')+2, 'FontWeight','bold');
        ylabel('\gamma (s,t)', 'FontSize', get(gca,'FontSize')+2, 'FontWeight','bold');

        if strlength(string(opt.TitlePrefix)) > 0
            title(sprintf('%s | %s trials (N=%d)', opt.TitlePrefix, labelName, sum(idx)), ...
                'Interpreter','none');
        else
            title(sprintf('%s trials (N=%d)', labelName, sum(idx)));
        end

        legend({'State 1','State 2','State 3','State 4'}, 'Location', char(opt.LegendLoc));
        ylim(opt.YLim);
        grid on; box on;

        hFig.(keyName) = h;
    end
end
