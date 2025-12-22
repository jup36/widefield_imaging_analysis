function hmmInput = buildHMMseq_fromHsY3(filePath, fileKeyword, varargin)
% buildHMMseq_fromHsY3
% Minimal prep for Gaussian HMM from motif activity Hs.Y3 (N x K x T).
%
% Outputs:
%   hmmInput.seqC   : 1xN cell, each [D x T] (trial sequence)
%   hmmInput.trI    : trial masks (goI, nogoI, etc.)
%   hmmInput.time   : 1xT time vector
%   hmmInput.pca    : PCA params if used

p = inputParser;
p.addParameter('doPCA', true, @(v)islogical(v)||ismember(v,[0 1]));
p.addParameter('nPC', 8, @(v)isnumeric(v)&&isscalar(v)&&v>=1);
p.addParameter('zscoreH', true, @(v)islogical(v)||ismember(v,[0 1]));
p.parse(varargin{:});
prm = p.Results;

% --- load like your xcorr pipeline ---
header      = extract_date_animalID_header(filePath);
keyword_beh = '_alignedPupilOrofacial.mat';

filePath_matfiles = cell2mat(GrabFiles_sort_trials('Matfiles', 0, {filePath}));
filePath_H        = cell2mat(GrabFiles_sort_trials([header '*' fileKeyword], 0, {filePath_matfiles}));
filePath_B        = cell2mat(GrabFiles_sort_trials([header '*' keyword_beh], 0, {filePath_matfiles}));

load(filePath_H, 'tbytDat_hAligned');
load(filePath_B, 'tbytDat');

trI = trialTypeInfoAuditoryGngTbytDat(tbytDat);

% --- stack trials: N x K x T ---
Hs = stack_trials_H(tbytDat_hAligned, 'zscore', prm.zscoreH);
Y3 = Hs.Y3;                % N x K x T
time = Hs.winCtrs(:)';     % 1 x T

[N,K,T] = size(Y3);

% --- reshape into (N*T) x K to fit PCA if requested ---
X = reshape(permute(Y3,[1 3 2]), N*T, K); % (N*T) x K
X(~isfinite(X)) = 0;

pcaInfo = struct();
if prm.doPCA
    [coeff, score, latent, ~, explained, mu] = pca(X, 'NumComponents', prm.nPC);
    Xred = score;  % (N*T) x nPC

    pcaInfo.coeff = coeff;
    pcaInfo.mu = mu;
    pcaInfo.latent = latent;
    pcaInfo.explained = explained;

    D = prm.nPC;
else
    Xred = X;
    D = K;
end

% --- rebuild trial sequences as cell array, each D x T ---
Xred3 = permute(reshape(Xred, [N, T, D]), [1 3 2]);  % N x D x T

seqC = cell(1,N);
for n = 1:N
    seqC{n} = squeeze(Xred3(n,:,:));  % D x T
end

hmmInput = struct();
hmmInput.seqC = seqC;
hmmInput.trI = trI;
hmmInput.time = time;
hmmInput.pca = pcaInfo;
hmmInput.D = D;
hmmInput.N = N;
hmmInput.T = T;
end
