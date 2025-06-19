function dffPostprocessAuditoryGng_meanSemPlots_6OHDA(filePath, channel)
% Produce mean±SEM PETH plots for every region / analysis pair in rezOut.
%
% INPUTS
%   filePath : session folder that was processed by the previous pipeline
%   channel  : 'green' | 'red' | … (whatever you used upstream)
%
% J. Park – 2025-06-04
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 0. Load stats produced by *dffPostprocessAuditoryGng_align_meanSem_6OHDA*
% -------------------------------------------------------------------------
hdr  = regexp(filePath,'m\d{1,4}_\d{6}','match','once');
statsFile = fullfile(filePath,"Matfiles",sprintf('%s_%s_evtAligned_stats.mat',hdr,channel));
load(statsFile,'rezOut');

figSaveDir = fullfile(filePath, 'Figure');
if ~isfolder(figSaveDir), mkdir(figSaveDir); end

% -------------------------------------------------------------------------
% 1. User-configurable options
% -------------------------------------------------------------------------
leftColor  = [204, 51,  51]/255;  % red(ish)
rightColor = [  0,102,204]/255;   % blue(ish)
plotColors = {leftColor, rightColor};

figPos     = [300 200 680 480];   % default figure size
lw         = 1.2;                 % line-width override

% -------------------------------------------------------------------------
% 2. Auto-discover region & analysis labels
% -------------------------------------------------------------------------
regions   = fieldnames(rezOut);
hemis     = {'L','R'};

% Use the first region to discover which analyses were computed
if isempty(regions)
    warning('rezOut is empty – nothing to plot.'); return; end
analysisLabels = fieldnames(rezOut.(regions{1}));

% -------------------------------------------------------------------------
% 3. Loop through regions and analyses
% -------------------------------------------------------------------------
% ---------------------------------------------------------------
% PLOT LOOP  (replace your current for-loops with this block)
% ---------------------------------------------------------------
for r = 1:numel(regions)
    reg = regions{r};

    for a = 1:numel(analysisLabels)
        lbl = analysisLabels{a};

        % -- fetch data for both hemispheres ---------------------
        valid = true(1,2);
        for h = 1:2
            hemi = hemis{h};
            mf = sprintf('mean%s',hemi);
            sf = sprintf('sem%s',hemi);
            tf = sprintf('ts%s',hemi);

            if isfield(rezOut.(reg).(lbl),mf)
                meanDat{h,1} = rezOut.(reg).(lbl).(mf);
                semDat{h,1}  = rezOut.(reg).(lbl).(sf);
                tsDat{h,1}   = rezOut.(reg).(lbl).(tf);
                if isempty(meanDat{h}); valid(h) = false; end
            else
                valid(h) = false;
            end
        end
        if ~all(valid), continue; end

        % -- same time axis check --------------------------------
        if ~isequal(tsDat{1},tsDat{2})
            warning('%s-%s time axes differ; using left hemi.',reg,lbl);
        end
        tAxis = tsDat{1};

        % -- plotting -------------------------------------------
        hFig = plotMeanSemColorC(cell2mat(meanDat), ...
            cell2mat(semDat), ...
            tAxis, plotColors, ...
            {'Left','Right'});
        % make invisible & size it
        set(hFig,'Visible','off','Position',figPos);

        % axis cosmetics
        ax = gca;
        ax.LineWidth = lw;
        ax.FontSize  = 11;
        ax.TickDir   = 'out';
        ax.XTick     = ceil(ax.XLim(1)/0.5)*0.5 : 0.5 : floor(ax.XLim(2)/0.5)*0.5;

        xlabel('Time (s)');
        ylabel('\DeltaF/F');
        title(sprintf('%s – %s',upper(reg),lbl),'Interpreter','none');

        % -- save to PDF (vector) -------------------------------
        figName = sprintf('%s_%s_%s_%s_%s',hdr,channel,upper(reg),lbl,'LRhemi6OHDA');
        print(hFig,fullfile(figSaveDir,figName),'-dpdf','-painters','-bestfit');
        fprintf('Saved %s.pdf  (%s | %s)\n',figName,upper(reg),lbl);   % <-- progress update
        close(hFig);   % free memory
    end
end

end
