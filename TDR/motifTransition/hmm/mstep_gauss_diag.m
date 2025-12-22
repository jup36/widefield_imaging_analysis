function model = mstep_gauss_diag(stats, opt)
S = numel(stats.sumGamma1);

piNew = stats.sumGamma1 / (sum(stats.sumGamma1) + eps);

Anew = stats.sumXi ./ (sum(stats.sumXi, 2) + eps); % row-normalize

muNew = stats.sumX ./ (stats.sumGamma' + eps);

Ex2   = stats.sumXX ./ (stats.sumGamma' + eps);
varNew = Ex2 - muNew.^2;
varNew(~isfinite(varNew) | varNew < opt.varFloor) = opt.varFloor;

model = struct('pi',piNew,'A',Anew,'mu',muNew,'var',varNew);
end
