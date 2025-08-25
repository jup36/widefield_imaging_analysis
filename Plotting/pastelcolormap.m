function cmap = pastelcolormap(N)
% Generates N pastel colors using HSV space
    hues = linspace(0, 1, N+1); hues(end) = [];  % avoid wrapping hue
    sat = 0.35; val = 0.95;
    cmap = hsv2rgb([hues(:), repmat(sat, N, 1), repmat(val, N, 1)]);
end
