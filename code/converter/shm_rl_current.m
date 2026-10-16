function I = shm_rl_current(th, s, pf, delta, xs)
%SHM_RL_CURRENT Exact periodic steady-state currents of a three-phase R-L load.
%   I = SHM_RL_CURRENT(th, s, pf, delta, xs) for the quarter-wave pattern th
%   (rad) with step signs s, a star-connected R-L load with isolated neutral
%   and power factor pf, and DC-link offset delta (positive levels scaled by
%   1 + delta, negative by 1 - delta). Units: voltage Vdc/2, |Z1| = 1. Returns
%   the phase currents (columns a, b, c) at the angles xs (rad, ascending in
%   [0, 2*pi)). The current is piecewise exponential and computed exactly.

[e, V] = segments(th, s, delta);
R = pf;  tau = tan(acos(pf));
a = exp(-diff(e) / tau);
i = zeros(1, 3);
for k = 1:numel(a), i = i * a(k) + (V(k, :) / R) * (1 - a(k)); end
i = i / (1 - exp(-2*pi / tau));
xs = xs(:);  I = zeros(numel(xs), 3);
k = 1;  tk = e(1);
for j = 1:numel(xs)
    while k < numel(e) - 1 && xs(j) >= e(k + 1)
        ak = exp(-(e(k + 1) - tk) / tau);
        i = i * ak + (V(k, :) / R) * (1 - ak);  tk = e(k + 1);  k = k + 1;
    end
    ak = exp(-(xs(j) - tk) / tau);
    I(j, :) = i * ak + (V(k, :) / R) * (1 - ak);
end
end

function [e, V] = segments(th, s, delta)
th = sort(th(:)).';
base = [th, pi - th, pi + th, 2*pi - th];
e = unique(mod([base, base + 2*pi/3, base + 4*pi/3, 0], 2*pi));
e = [e(:); 2*pi];
m = 0.5 * (e(1:end-1) + e(2:end));
L = [pole_level(th, s, m), pole_level(th, s, m - 2*pi/3), pole_level(th, s, m - 4*pi/3)];
L(L > 0) = L(L > 0) * (1 + delta);  L(L < 0) = L(L < 0) * (1 - delta);
V = L - repmat(mean(L, 2), 1, 3);
end

function lv = pole_level(th, s, x)
x = mod(x, 2*pi);  sg = ones(size(x));
k = x >= pi;  x(k) = x(k) - pi;  sg(k) = -1;
k = x > pi/2;  x(k) = pi - x(k);
c = cumsum(s(:).');  lv = zeros(size(x));
for j = 1:numel(th), lv(x >= th(j)) = c(j); end
lv = sg .* lv;
end
