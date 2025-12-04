function fh = visualize_convolution_alignment(basis, step, whichCols, nBins, t0_sec)
% VISUALIZE_CONVOLUTION_ALIGNMENT
%   Visualize kernels (before/after centering) and their convolved outputs
%   to verify that early vs. late kernels produce responses at correct
%   physical times relative to an event.
%
% Inputs
%   basis     : struct with fields
%               .B [nLag x nB]  (unshifted basis bank)
%               .lagsSec [nLag x 1] (lag axis, must include 0)
%   step      : bin size in seconds (e.g., 0.05)
%   whichCols : indices of kernels to visualize (e.g., [1 size(B,2)])
%   nBins     : total time bins in the demo event train (e.g., 300)
%   t0_sec    : time of the event in seconds (e.g., 5.0)
%
% Output
%   fh        : figure handle

if nargin < 3 || isempty(whichCols)
    whichCols = [1 size(basis.B,2)]; % first & last by default
end
if nargin < 4 || isempty(nBins), nBins = 300; end
if nargin < 5 || isempty(t0_sec), t0_sec = 5.0; end

B    = basis.B;          % unshifted bank
lags = basis.lagsSec(:);     % lag vector, MUST include 0
[nLag, nB] = size(B); %#ok<ASGLU>

% --- find 0-lag index and center index for conv(...,'same') ---
idx0      = find(lags==0,1,'first');
if isempty(idx0)
    error('visualize_convolution_alignment: lagsSec must contain 0 exactly.');
end
centerIdx = ceil((nLag+1)/2);
shift     = centerIdx - idx0;

% --- make a centered (shifted) bank for convolution ---
B_centered = circshift(B, shift, 1);  % shift rows = time

% --- build a simple event train with a single impulse at t0_sec ---
t         = (0:nBins-1)*step;
[~, t0bin]= min(abs(t - t0_sec));
E         = zeros(1, nBins); 
E(t0bin)  = 1; % impulse at t0

% --- Convolve with unshifted vs centered (for comparison) ---
% We'll only plot the *centered* results for correctness, but compute both.
F_unshift = zeros(nBins, numel(whichCols));
F_center  = zeros(nBins, numel(whichCols));

for ii = 1:numel(whichCols)
    c  = whichCols(ii);
    ku = B(:,c);
    kc = B_centered(:,c);
    F_unshift(:,ii) = conv(E, ku, 'same');  % misaligned if used in practice
    F_center(:,ii)  = conv(E, kc, 'same');  % correctly aligned
end

% --- Expected peak times (approx) for intuition: use basis centers ---
% The "physical lag location" of the kernel is where it has its max on lags
lag_peak = zeros(1, numel(whichCols));
for ii = 1:numel(whichCols)
    [~,pk] = max(B(:, whichCols(ii)));
    lag_peak(ii) = lags(pk);  % in seconds
end
t_expected = t0_sec + lag_peak; % expected peak time of the convolved output

% ----------- Plot -----------
fh = figure('Color','w','Name','Convolution Alignment Diagnostics');

% 1) Original kernels vs physical lag
subplot(3,1,1); hold on
for ii = 1:numel(whichCols)
    plot(lags, B(:,whichCols(ii)), 'LineWidth', 1.5);
end
xlabel('Lag (s)'); ylabel('Amplitude'); grid on
legend(compose('Kernel %d (unshifted)', whichCols), 'Location','best');
title('Original kernels (unshifted) vs. physical lag');

% 2) Centered kernels (index space; center index = 0-lag anchor)
subplot(3,1,2); hold on
L = size(B_centered,1);
ctrAxis = ( (1:L) - centerIdx ) * step; % index mapped to seconds around center
for ii = 1:numel(whichCols)
    plot(ctrAxis, B_centered(:,whichCols(ii)), 'LineWidth', 1.5);
end
xline(0,'--k','0-lag (conv anchor)');
xlabel('Centered index \times step (s)'); ylabel('Amplitude'); grid on
legend(compose('Kernel %d (centered)', whichCols), 'Location','best');
title('Centered kernels (0-lag at center index used by conv)');

% 3) Convolved outputs vs absolute time
subplot(3,1,3); hold on
cols = lines(numel(whichCols));
for ii = 1:numel(whichCols)
    plot(t, F_center(:,ii), 'Color', cols(ii,:), 'LineWidth', 1.5);
    xline(t_expected(ii), ':', sprintf('~t0+%.2fs', lag_peak(ii)), ...
        'Color', cols(ii,:), 'LabelOrientation','horizontal');
end
xline(t0_sec, '--k', 'event (t0)');
xlabel('Time (s)'); ylabel('Convolved output'); grid on
legend(compose('Conv with centered kernel %d', whichCols), 'Location','best');
title('Convolved outputs vs. time (event at t0)');

% Annotation to explain expected ordering
txt = sprintf(['Expected: earlier kernel (e.g., col %d) should peak near t0 + %.2fs,\n' ...
               'later kernel (e.g., col %d) should peak near t0 + %.2fs.'], ...
               whichCols(1), lag_peak(1), whichCols(end), lag_peak(end));
annotation('textbox',[0.12 0.02 0.8 0.08], 'String', txt, ...
           'EdgeColor','none','Interpreter','none');

end
