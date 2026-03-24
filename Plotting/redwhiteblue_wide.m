function cmap = redwhiteblue_wide(n, gamma)

if nargin < 1, n = 256; end
if nargin < 2, gamma = 2; end

x = linspace(-1,1,n);

% nonlinear warping
xw = sign(x).*abs(x).^gamma;

r = zeros(n,1);
g = zeros(n,1);
b = zeros(n,1);

for i = 1:n
    if xw(i) < 0
        % blue → white
        t = (xw(i)+1);
        r(i) = t;
        g(i) = t;
        b(i) = 1;
    else
        % white → red
        t = (1-xw(i));
        r(i) = 1;
        g(i) = t;
        b(i) = t;
    end
end

cmap = [r g b];
end