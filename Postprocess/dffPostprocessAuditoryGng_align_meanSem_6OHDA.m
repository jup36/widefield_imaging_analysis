function dffPostprocessAuditoryGng_align_meanSem_6OHDA(filePath, channel)
% Compute peri‑event averages (mean & SEM) for every region/hemisphere.
% Saves one struct per region with sub‑fields per analysis.

%% ------------------------------------------------------------------
% 0. Load input produced by previous pipeline stage
% -------------------------------------------------------------------
hdr  = regexp(filePath,'m\d{1,4}_\d{6}','match','once');
load(fullfile(filePath,"Matfiles",sprintf('%s_%s_dff_evtAligned_regionMask.mat',hdr,channel)), ...
     'rez','voluntaryHitLickI');

%% ------------------------------------------------------------------
% 1. User settings
% -------------------------------------------------------------------
regions   = {'m1','m2','ss','v1','rs'};   % regions available in rez
hemis     = {'L','R'};                    % left/right
stepDT    = 0.01;                         % interpolation step

% Each row: {rezField  outLabel   trialIndex(optional)}
analyses = { ...
    {'stimOnDffC'     ,'stimOn'           }, ...
    {'stimOnDffC'     ,'stimOnVol' ,voluntaryHitLickI}, ...
    {'hitLickFst'     ,'hitFst'           }, ...
    {'hitLickFst'     ,'hitFstVol',voluntaryHitLickI}, ...
    {'hitLickLast'    ,'hitLast'          }, ...
    {'hitLickToneAlign','hitTone'         }, ...
    {'firstWater'     ,'water'            }, ...
    {'firstWater'     ,'waterVol' ,voluntaryHitLickI}};

outFile = fullfile(filePath,'Matfiles',sprintf('%s_%s_evtAligned_stats.mat',hdr,channel));

%% ------------------------------------------------------------------
% 2. Loop through regions / analyses
% -------------------------------------------------------------------
rezOut = struct();
for r = 1:numel(regions)
    reg = regions{r};

    for a = 1:numel(analyses)
        srcField = analyses{a}{1};   % field in rez
        lbl      = analyses{a}{2};   % label in output
        if numel(analyses{a})==3
            index   = analyses{a}{3};
        else
            index   = true(size(rez.(srcField).([reg hemis{1}]),1), 1);
        end

        % Pre‑allocate output sub‑struct
        if ~isfield(rezOut,reg), rezOut.(reg) = struct(); end
        rezOut.(reg).(lbl) = struct();

        for h = 1:2
            hemi = hemis{h};
            sigCell  = rez.(srcField).([reg hemi])(index,1);
            timeCell = rez.(srcField).([reg hemi])(index,2);
            if isempty(sigCell)   % no trials
                rezOut.(reg).(lbl).(['mean' hemi]) = [];
                rezOut.(reg).(lbl).(['sem'  hemi]) = [];
                rezOut.(reg).(lbl).(['ts'   hemi]) = [];
                continue
            end
            % interpolate & stack
            [aligned,itpTs] = temporalAlignInterp1(sigCell,timeCell,stepDT);
            [mu,~,se]      = meanstdsem(cell2mat(aligned));

            rezOut.(reg).(lbl).(['mean' hemi]) = mu;
            rezOut.(reg).(lbl).(['sem'  hemi]) = se;
            rezOut.(reg).(lbl).(['ts'   hemi]) = itpTs;  % store the exact axis used
        end
    end
end

%% ------------------------------------------------------------------
% 3. Save
% -------------------------------------------------------------------
save(outFile,'rezOut');
fprintf('Saved peri‑event stats → %s\n',outFile);
end
