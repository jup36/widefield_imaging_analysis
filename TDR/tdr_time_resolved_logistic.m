function out = tdr_time_resolved_logistic(figSaveDir, tbytDat_hAligned, y, varargin)
% Time-resolved TDR with ridge logistic regression (Hit vs CR).
% Returns standardized weights (Wz), original-units weights (W),
% projections, d', and plotting.

% ---- params
p = inputParser;
p.addParameter('Epoch', [-0.5 4], @(v)isnumeric(v)&&numel(v)==2);
p.addParameter('Win', 0.100, @(v)isnumeric(v)&&isscalar(v)&&v>0);      % 100 ms
p.addParameter('Step', 0.050, @(v)isnumeric(v)&&isscalar(v)&&v>0);     % 50 ms
p.addParameter('Lambda', 1, @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('SmoothSigma', 0, @(v)isnumeric(v)&&isscalar(v)&&v>=0); % (unused here)
p.addParameter('DoPlots', true, @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('printPlots', true, @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('minSD', 1e-3, @(v)isnumeric(v)&&isscalar(v)&&v>0);     % variance floor for unscaling
p.addParameter('minPerClass', 5, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
p.addParameter('cdWeightCscale', [-0.15 0.15], @(v)isnumeric(v)&&numel(v)==2);
p.addParameter('dPrimeYlim', [0.2 2], @(v)isnumeric(v)&&numel(v)==2); 
p.addParameter('projCdYlim', [-1.5 2], @(v)isnumeric(v)&&numel(v)==2); 

p.parse(varargin{:});
prm = p.Results;

% ---- basics
N = size(tbytDat_hAligned, 2);
assert(numel(y)==N, 'Label vector y must match #trials.');
K = size(tbytDat_hAligned{1,1},1); % #motifs

% windows
tStart = prm.Epoch(1); tEnd = prm.Epoch(2);
win = prm.Win; step = prm.Step;
ctr = (tStart + win/2) : step : (tEnd - win/2);
nW = numel(ctr);
winBounds = [ctr(:)-win/2, ctr(:)+win/2]; % nW x 2

% trials cache
Htrials = cell(1,N);
ttrials = cell(1,N);
for n = 1:N
    Htrials{n} = tbytDat_hAligned{1,n};     % K x Tn
    ttrials{n} = tbytDat_hAligned{2,n}(:)'; % 1 x Tn (s)
end

% ---- pass 1: assemble windowed features
Xw = cell(1,nW);         % N x K (NaNs allowed)
validw = cell(1,nW);     % N x 1
for i = 1:nW
    lb = winBounds(i,1); ub = winBounds(i,2);
    X = nan(N, K);
    for n = 1:N
        tn = ttrials{n};
        idx = find(tn >= lb & tn < ub);
        if ~isempty(idx)
            Hn = Htrials{n};
            X(n,:) = mean(Hn(:,idx), 2, 'omitnan')';
        end
    end
    Xw{i} = X;
    validw{i} = ~any(isnan(X),2);
end

% ---- global stats across all windows/trials
Xstack = cat(1, Xw{:});                  % (N*nW) x K with NaNs
muG = mean(Xstack, 1, 'omitnan');
sdG = std( Xstack, 0, 1, 'omitnan');     % 0 => normalize by N-1
sdG(sdG < prm.minSD) = prm.minSD;        % variance floor

% ---- pass 2: fit per-window models with global z-scoring
Wz = nan(K, nW);    % standardized weights (for plotting/comparison)
W  = nan(K, nW);    % original-units weights (for interpretability)
b  = nan(1, nW);
proj = nan(N, nW);
dprime = nan(1, nW);
motifRank = cell(1, nW);

for i = 1:nW
    X = Xw{i};
    valid = validw{i};
    v_hit = valid & (y(:)==1);
    v_cr  = valid & (y(:)==0);
    if sum(v_hit) < prm.minPerClass || sum(v_cr) < prm.minPerClass
        continue
    end

    % global z-score
    Xz_valid = (X(valid,:) - muG) ./ sdG;

    % ridge logistic
    mdl = fitclinear(Xz_valid, y(valid), 'Learner','logistic', ...
        'Regularization','ridge', 'Lambda', prm.Lambda, ...
        'Solver','lbfgs', 'ClassNames',[0,1]);

    w_z = mdl.Beta;                   % K x 1 (standardized)
    b_i = mdl.Bias;

    % store standardized weights
    Wz(:,i) = w_z;

    % unscale (stable due to sdG floor)
    w     = w_z ./ sdG(:);
    b_unz = b_i - muG * (w_z ./ sdG(:));

    % projections (all trials; invalid remain NaN)
    s = nan(N,1);
    s(valid) = X(valid,:) * w + b_unz;

    % align sign (Hit positive)
    if mean(s(v_hit),'omitnan') < mean(s(v_cr),'omitnan')
        w   = -w;  w_z = -w_z;
        s   = -s;  b_unz = -b_unz;
    end

    W(:,i)     = w;
    proj(:,i)  = s;
    m1 = mean(s(v_hit),'omitnan'); v1 = var(s(v_hit), 'omitnan');
    m0 = mean(s(v_cr), 'omitnan'); v0 = var(s(v_cr),  'omitnan');
    dprime(i)  = (m1 - m0) / sqrt(0.5*(v1+v0) + eps);
    [~,ord]    = sort(abs(w),'descend'); motifRank{i} = ord(:); % row
end

% % ===== NEW: Global sign alignment to response-period reference (2–4 s) =====
% flipWsign = false(1, nW);                        % track which windows are flipped
% respWins = find(ctr >= 2 & ctr < 4 & any(~isnan(Wz),1));  % windows in [2,4) with some data
% if ~isempty(respWins)
%     refVec = mean(Wz(:, respWins), 2, 'omitnan'); % reference = mean Wz in response period
%     if any(isfinite(refVec)) && norm(refVec(~isnan(refVec))) > 0
%         % Normalize reference safely
%         nr = norm(refVec(~isnan(refVec)));
%         refUnit = refVec; refUnit(~isnan(refUnit)) = refVec(~isnan(refVec)) / nr;
% 
%         for i = 1:nW
%             wi = Wz(:,i);
%             if all(isnan(wi)), continue; end
%             % cosine similarity with reference (ignore NaNs)
%             mask = isfinite(wi) & isfinite(refUnit);
%             if ~any(mask), continue; end
%             cosSim = (wi(mask)' * refUnit(mask)) / (norm(wi(mask)) * norm(refUnit(mask)));
%             if isfinite(cosSim) && cosSim < 0
%                 % flip signs to match the response-period reference
%                 Wz(:,i)   = -Wz(:,i);
%                 W(:,i)    = -W(:,i);
%                 proj(:,i) = -proj(:,i);
%                 b(i)      = -b(i);
%                 flipWsign(i) = true;
%             end
%         end
%     end
% end
% ===========================================================================

% ---- pack output
out = struct('Wz', Wz, 'W', W, 'b', b, ...
    'winCtrs', ctr, 'winBounds', winBounds, ...
    'proj', proj, 'dprime', dprime, ...
    'motifRank', cell2mat(motifRank), ...
    'params', prm);

% ---- plots
if prm.DoPlots
    fb = @(x, y1, y2, a, c) patch([x(:)' fliplr(x(:)')], ...
        [y1(:)' fliplr(y2(:)')], ...
        c, 'EdgeColor','none', 'FaceAlpha',a);

    h_proj = figure('Name','Coding projection (Hit vs CR)'); hold on;
    m_hit = mean(out.proj(y==1,:),1,'omitnan'); se_hit = std(out.proj(y==1,:),0,1,'omitnan')/sqrt(sum(y==1));
    m_cr  = mean(out.proj(y==0,:),1,'omitnan'); se_cr  = std(out.proj(y==0,:),0,1,'omitnan')/sqrt(sum(y==0));

    col_hit = [0.8 0 0];   % red
    col_cr  = [0 0 0.8];   % blue

    % plot projection score 
    fb(out.winCtrs, m_hit-se_hit, m_hit+se_hit, 0.3, col_hit);
    plot(out.winCtrs, m_hit, 'Color', col_hit, 'LineWidth', 2);

    fb(out.winCtrs, m_cr -se_cr,  m_cr +se_cr,  0.3, col_cr);
    plot(out.winCtrs, m_cr, 'Color', col_cr, 'LineWidth', 2);

    xline(0,'k:'); xline(2,'k:'); xlabel('Time (s)'); ylabel('Projection'); ylim(prm.projCdYlim); legend({'Hit \pm SE','Hit','CR \pm SE','CR'});

    % plot dPrime score
    h_dprime = figure('Name','Time-resolved d'''); plot(out.winCtrs, out.dprime,'LineWidth',2);
    xline(0,'k:'); xline(2,'k:'); ylim(prm.dPrimeYlim); xlabel('Time (s)'); ylabel('d'''); box on;

    % imagesc Wz
    h_Wz = figure('Name','Motif weights heatmap (standardized)');
    imagesc(out.winCtrs, 1:K, out.Wz); axis xy;
    xlabel('Time (s)'); ylabel('Motif #'); colorbar; clim(prm.cdWeightCscale); 
    title('Coding weights (standardized units; + = Hit)'); hold on; ylims = ylim; plot([0 0; 2 2],[ylims; ylims],'k:');

    % Print (Optional)
    if prm.printPlots
        if ~exist(figSaveDir, 'dir'); mkdir(figSaveDir); end
        header = extract_date_animalID_header(figSaveDir);
        timestampStr = datestr(now, 'mmddyy_HHMMSS');  % Date and time string
        figSaveName_proj = sprintf('proj_score_CD_%s_%s', header, timestampStr);
        print(h_proj, fullfile(figSaveDir, figSaveName_proj),'-dpdf','-painters','-bestfit')
        figSaveName_dprime = sprintf('dPrime_CD_%s_%s', header, timestampStr);
        print(h_dprime, fullfile(figSaveDir, figSaveName_dprime),'-dpdf','-painters','-bestfit')
        figSaveName_Wz = sprintf('wZ_CD_%s_%s', header, timestampStr);
        print(h_Wz, fullfile(figSaveDir, figSaveName_Wz),'-dpdf','-painters','-bestfit')
    end
end
end
