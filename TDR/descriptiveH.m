function HstatC = descriptiveH(HsY3d, Htime, trI)
%DESCRIPTIVEH  Compute descriptive statistics (mean & SEM) of raw motif activity.
%
% SYNOPSIS
%   HstatC = descriptiveH(HsY3d, Htime, trI)
%
% DESCRIPTION
%   This function summarizes raw motif activity (H) extracted from the
%   CNMF/seqNMF pipeline. It computes **trial-averaged motif activity**
%   separately for:
%       • Go trials (trI.goI)
%       • No-Go trials (trI.nogoI)
%
%   For each motif k:
%       • Extracts H across all trials → [N x K x T]
%       • Selects trials belonging to Go or NoGo
%       • Computes mean and SEM across those trials
%       • Packs results into a struct with fields:
%           - motifIdx
%           - mean.periCueGoH
%           - mean.periCueNoGoH
%           - sem.periCueGoH
%           - sem.periCueNoGoH
%
% INPUTS
%   HsY3d  : [N x K x T] 3D array of motif activity
%       N = # trials, K = # motifs, T = # time bins
%
%   Htime  : [1 x T] vector of time points (e.g., Hs.winCtrs)
%       (Stored for clarity but not used inside the function.)
%
%   trI    : struct containing logical trial-type masks:
%               trI.goI   — [N x 1] true for Go trials
%               trI.nogoI — [N x 1] true for NoGo trials
%
% OUTPUT
%   HstatC : 1 x K cell array
%       Each cell contains a structure:
%
%       HstatC{k}.motifIdx        — motif index
%       HstatC{k}.mean.periCueGoH — [1 x T]
%       HstatC{k}.mean.periCueNoGoH
%       HstatC{k}.sem.periCueGoH  — [1 x T]
%       HstatC{k}.sem.periCueNoGoH
%
% DEPENDENCIES
%   meanstdsem.m
%       Expected signature: [m, sd, se] = meanstdsem(data)
%       Note: This function is assumed to operate along dim 1.
%
% EXAMPLE
%   HsY3d = Hs.Y3;
%   Htime = Hs.winCtrs;
%   Hstat = descriptiveH(HsY3d, Htime, trI);
%
%   % Plot motif #5 Go vs NoGo
%   figure; hold on;
%   plot(Hstat{5}.mean.periCueGoH);
%   plot(Hstat{5}.mean.periCueNoGoH);
%   legend('Go','NoGo');
%
% -------------------------------------------------------------------------

% Initialize output cell
HstatC = cell(1, size(HsY3d, 2));

% Loop over motifs
for k = 1:size(HsY3d, 2)

    % Extract motif activity for Go and NoGo trials
    % squeeze → [nTrials_type x T]
    [mean_go_H, ~, sem_go_H] = meanstdsem(squeeze(HsY3d(trI.goI,   k, :))); 
    [mean_nogo_H, ~, sem_nogo_H] = meanstdsem(squeeze(HsY3d(trI.nogoI, k, :))); 
    
    % Pack into struct
    Hstat = struct();
    Hstat.motifIdx = k; 
    Hstat.winCtrs = Htime; 

    % Mean fields
    Hstat.mean = struct( ...
        'periCueGoH',   mean_go_H, ...
        'periCueNoGoH', mean_nogo_H);

    % SEM fields
    Hstat.sem = struct( ...
        'periCueGoH',   sem_go_H, ...
        'periCueNoGoH', sem_nogo_H);

    % Assign to cell array
    HstatC{1, k} = Hstat; 
end

end
