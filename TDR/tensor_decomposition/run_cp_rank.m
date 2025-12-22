function cpRez = run_cp_rank(W_all, ranks, varargin)
% RUN_CP_RANK  Fit CP/PARAFAC to W_all for multiple ranks.
%
%   cpRez = run_cp_rank(W_all, ranks, 'tol',1e-6, 'maxiters',500, 'printitn',0)
%
% INPUT
%   W_all : [S x K x T] tensor of decoder weights
%   ranks : vector of candidate CP ranks, e.g. 1:8
%
% NAME-VALUE OPTIONS
%   'tol'      : convergence tolerance for cp_als (default 1e-6)
%   'maxiters' : max iterations (default 500)
%   'printitn' : cp_als print frequency (0 = silent; default 0)
%
% OUTPUT
%   cpRez : struct with fields
%       .R        : vector of ranks
%       .models   : cell{r} ktensor models
%       .factors  : cell{r} factor cells, each {U_sess, U_motif, U_time}
%       .fit      : 1 x numel(ranks) vector of variance explained (0–1)
%       .Xnorm    : Frobenius norm of data tensor

    % ---- parse options ----
    p = inputParser;
    p.addParameter('tol', 1e-6, @(x)isnumeric(x)&&isscalar(x));
    p.addParameter('maxiters', 500, @(x)isnumeric(x)&&isscalar(x));
    p.addParameter('printitn', 0, @(x)isnumeric(x)&&isscalar(x));
    p.parse(varargin{:});
    opt = p.Results;

    % ---- wrap data as Tensor Toolbox object ----
    Xten  = tensor(W_all);
    Xnorm = norm(Xten);

    % ---- preallocate results ----
    ranks   = ranks(:)';               % row vector
    nR      = numel(ranks);
    models  = cell(1, nR);
    factors = cell(1, nR);
    fit     = nan(1, nR);

    % ---- loop over ranks ----
    for ii = 1:nR
        R = ranks(ii);
        fprintf('\n[run_cp_rank] Fitting CP-ALS with rank R = %d\n', R);

        [M_cp, U, out_cp] = cp_als(Xten, R, ...
            'tol', opt.tol, ...
            'maxiters', opt.maxiters, ...
            'printitn', opt.printitn);

        resid   = norm(Xten - full(M_cp));
        fitFrac = 1 - (resid / Xnorm)^2;

        models{ii}  = M_cp;
        factors{ii} = U;
        fit(ii)     = fitFrac;

        fprintf('  -> Fit = %.2f%%%% variance explained (iters = %d)\n', ...
            100*fitFrac, out_cp.iters);
    end

    % ---- pack output ----
    cpRez = struct();
    cpRez.R       = ranks;
    cpRez.models  = models;
    cpRez.factors = factors;
    cpRez.fit     = fit;
    cpRez.Xnorm   = Xnorm;
end
