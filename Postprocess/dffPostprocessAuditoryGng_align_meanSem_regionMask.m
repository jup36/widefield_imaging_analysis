function dffPostprocessAuditoryGng_align_meanSem_regionMask(filePath, channel)
%--------------------------------------------------------------------------
% Compute peri-event averages (mean & SEM) for every cortical region.
%   • Uses evt-aligned traces produced by dffPostprocess_auditory_gng_…  
%   • Interpolates each trial to a common time base (stepDT)  
%   • Saves rezOut.(region).(label).mean / .sem / .ts
%--------------------------------------------------------------------------

%% 0. Load input produced by previous stage
hdr = regexp(filePath,'m\d{1,4}_\d{6}','match','once');
load(fullfile(filePath,"Matfiles", ...
     sprintf('%s_%s_dff_evtAligned_regionMask.mat',hdr,channel)), ...
     'rez','voluntaryHitLickI');

%% 1. Settings -----------------------------------------------------------
regions = {'m1','m2','ss','v1','rs'};     % region names in rez
stepDT  = 0.01;                            % time-axis resolution (s)

% analyses: {rezField  outputLabel  logicalIndex(optional)}
analyses = { ...
    {'stimOnDffC'      ,'stimOn'           }, ...
    {'stimOnDffC'      ,'stimOnVol' ,voluntaryHitLickI}, ...
    {'hitLickFst'      ,'hitFst'           }, ...
    {'hitLickFst'      ,'hitFstVol',voluntaryHitLickI}, ...
    {'hitLickLast'     ,'hitLast'          }, ...
    {'hitLickToneAlign','hitTone'          }, ...
    {'firstWater'      ,'water'            }, ...
    {'firstWater'      ,'waterVol' ,voluntaryHitLickI}};

outFile = fullfile(filePath,"Matfiles", ...
          sprintf('%s_%s_evtAligned_stats.mat',hdr,channel));

%% 2. Loop through regions / analyses -----------------------------------
rezOut = struct();

for r = 1:numel(regions)
    reg = regions{r};

    for a = 1:numel(analyses)
        srcField = analyses{a}{1};
        lbl      = analyses{a}{2};

        % trial-selection mask
        if numel(analyses{a}) == 3
            useTrials = analyses{a}{3};
        else
            useTrials = true(size(rez.(srcField).(reg),1),1);
        end

        sigCell  = rez.(srcField).(reg)(useTrials,1);
        timeCell = rez.(srcField).(reg)(useTrials,2);

        % ensure field exists
        if ~isfield(rezOut,reg), rezOut.(reg) = struct(); end
        rezOut.(reg).(lbl) = struct();

        if isempty(sigCell) || all(cellfun(@isempty,sigCell))
            rezOut.(reg).(lbl).mean = [];
            rezOut.(reg).(lbl).sem  = [];
            rezOut.(reg).(lbl).ts   = [];
            continue
        end

        % interpolate every trial to common axis
        [aligned,itpTs] = temporalAlignInterp1(sigCell,timeCell,stepDT);

        % aggregate
        [mu,~,se] = meanstdsem(cell2mat(aligned));

        % store
        rezOut.(reg).(lbl).mean = mu;
        rezOut.(reg).(lbl).sem  = se;
        rezOut.(reg).(lbl).ts   = itpTs;
    end
end

%% 3. Save ---------------------------------------------------------------
save(outFile,'rezOut');
fprintf('Saved peri-event stats → %s\n',outFile);
end

