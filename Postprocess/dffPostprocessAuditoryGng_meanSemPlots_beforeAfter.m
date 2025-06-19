function dffPostprocessAuditoryGng_meanSemPlots_beforeAfter( ...
    statsFileBefore, statsFileAfter, ...
    labelBefore, labelAfter)
%--------------------------------------------------------------------------
% Plot mean ± SEM peri-event traces from two sessions (before vs after).
% If their time axes differ, both traces are re-interpolated onto the
% finest common grid over the interval where they overlap.
%
% INPUTS
%   statsFileBefore : ..._evtAligned_stats.mat   (baseline)
%   statsFileAfter  : ..._evtAligned_stats.mat   (post-drug)
%   labelBefore     : legend label for baseline  (default 'beforeInj')
%   labelAfter      : legend label for post-drug (default 'afterInj')
%--------------------------------------------------------------------------

if nargin < 3, labelBefore = 'beforeInj'; end
if nargin < 4, labelAfter  = 'afterInj';  end

%% 0. Load both rezOut structures ------------------------------------------
rezBef = load(statsFileBefore,'rezOut').rezOut;
rezAft = load(statsFileAfter ,'rezOut').rezOut;

% Session ID + channel (assumes …m1234_YYMMDD_<channel>_evtAligned_stats.mat)
[~,stemBef] = fileparts(statsFileBefore);
tok = regexp(stemBef,'(m\d{1,4}_\d{6})_(\w+)_evtAligned_stats','tokens','once');
hdr     = tok{1};
channel = tok{2};

%% 1. Plot-style settings ---------------------------------------------------
colorBef   = [  0 102 204]/255;    % blue
colorAft   = [204  51  51]/255;    % red
plotColors = {colorBef,colorAft};

figPos = [300 200 680 480];   lw = 1.2;

%% 2. Region / analysis labels ---------------------------------------------
regions = fieldnames(rezBef);
if isempty(regions)
    warning('No regions found in %s',statsFileBefore); return; end
analysisLabels = fieldnames(rezBef.(regions{1}));

%% 3. Figure-save directory -------------------------------------------------
matfilesDir = fileparts(statsFileBefore);     % …/Matfiles
sessionDir  = fileparts(matfilesDir);         % …/task  (or session root)
figSaveDir  = fullfile(sessionDir,'Figures_compare_beforeAfter');
if ~isfolder(figSaveDir), mkdir(figSaveDir); end

%% 4. Plot loop -------------------------------------------------------------
for r = 1:numel(regions)
    reg = regions{r};

    for a = 1:numel(analysisLabels)
        lbl = analysisLabels{a};

        % -- fetch traces --------------------------------------------------
        [mu1,se1,t1,ok1] = grabTrace(rezBef,reg,lbl);
        [mu2,se2,t2,ok2] = grabTrace(rezAft,reg,lbl);
        if ~(ok1 && ok2), continue; end

        % -- align if time axes differ ------------------------------------
        if ~isequal(t1,t2)
            [muAligned,seAligned,tCommon,success] = ...
                alignTwoTraces(mu1,se1,t1,mu2,se2,t2);
            if ~success
                warning('%s-%s: no overlap between time axes – skipped.',reg,lbl);
                continue
            end
            muPlot  = muAligned;
            sePlot  = seAligned;
            tAxis   = tCommon;
        else
            muPlot = [mu1; mu2];
            sePlot = [se1; se2];
            tAxis  = t1;
        end

        % -- plot ----------------------------------------------------------
        hFig = plotMeanSemColorC(muPlot,sePlot,tAxis,plotColors, ...
            {labelBefore,labelAfter});
        set(hFig,'Visible','off','Position',figPos);

        ax = gca;
        ax.LineWidth = lw; ax.FontSize = 11; ax.TickDir = 'out';
        ax.XTick = ceil(ax.XLim(1)/0.5)*0.5 : 0.5 : floor(ax.XLim(2)/0.5)*0.5;

        xlabel('Time (s)'); ylabel('\DeltaF/F');
        title(sprintf('%s – %s',upper(reg),lbl),'Interpreter','none');

        % -- save ----------------------------------------------------------
        figName = sprintf('%s_%s_%s_%s_%svs%s', ...
            hdr,channel,upper(reg),lbl,labelBefore,labelAfter);
        print(hFig,fullfile(figSaveDir,figName),'-dpdf','-vector','-bestfit');
        fprintf('Saved %s.pdf  (%s | %s)\n',figName,upper(reg),lbl);
        close(hFig);
    end
end
end
%% =======================================================================
function [mu,se,ts,ok] = grabTrace(rez,reg,lbl)
% Return mean, sem, time axis + flag
mu = []; se = []; ts = []; ok = false;
if isfield(rez.(reg).(lbl),'mean')
    mu = rez.(reg).(lbl).mean;
    se = rez.(reg).(lbl).sem;
    ts = rez.(reg).(lbl).ts;
    ok = ~isempty(mu);
end
end
% -----------------------------------------------------------------------
function [muOut,seOut,tCommon,success] = ...
    alignTwoTraces(mu1,se1,t1,mu2,se2,t2)

% --- sanitise & sort --------------------------------------------------
[t1,idx1] = sort(t1(:));   mu1 = mu1(idx1);   se1 = se1(idx1);
[t2,idx2] = sort(t2(:));   mu2 = mu2(idx2);   se2 = se2(idx2);

good1 = ~isnan(t1) & ~isnan(mu1);   t1 = t1(good1); mu1 = mu1(good1); se1 = se1(good1);
good2 = ~isnan(t2) & ~isnan(mu2);   t2 = t2(good2); mu2 = mu2(good2); se2 = se2(good2);

if isempty(t1) || isempty(t2)
    muOut=[]; seOut=[]; tCommon=[]; success=false; return
end

% --- find true overlap -----------------------------------------------
ovStart = max( min(t1), min(t2) );
ovEnd   = min( max(t1), max(t2) );

% enlarge by ±½ sample to avoid rounding cut-offs
dt1 = median(diff(t1));   dt2 = median(diff(t2));
dt  = min(dt1,dt2);
ovStart = ovStart + 0.5*dt;     % trim inside by half-sample
ovEnd   = ovEnd   - 0.5*dt;

if ovStart >= ovEnd          % still nothing
    muOut=[]; seOut=[]; tCommon=[]; success=false; return
end

% --- build common grid -----------------------------------------------
tCommon = ovStart : dt : ovEnd;

muOut = [interp1(t1,mu1,tCommon,'linear','extrap') ; ...
    interp1(t2,mu2,tCommon,'linear','extrap')];
seOut = [interp1(t1,se1,tCommon,'linear','extrap') ; ...
    interp1(t2,se2,tCommon,'linear','extrap')];

success = true;
end
% -----------------------------------------------------------------------
function out = ternary(cond,a,b)
if cond, out = a; else, out = b; end
end
