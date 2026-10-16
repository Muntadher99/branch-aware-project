function d = shm_wavediff(th1, th2, s)
%SHM_WAVEDIFF Share of the period (%) in which two patterns differ.
%   d = SHM_WAVEDIFF(th1, th2, s) for angle sets th1, th2 (rad, first quarter)
%   with step signs s. Exact for piecewise-constant patterns.

e = unique([0, sort(th1(:)).', sort(th2(:)).', pi/2]);
x = (e(1:end-1) + e(2:end)) / 2;
d = 100 * sum(diff(e) .* (level(th1, s, x) ~= level(th2, s, x))) / (pi/2);
end

function L = level(th, s, x)
th = sort(th(:)).';  c = cumsum(s(1:numel(th)));
L = zeros(size(x));
for j = 1:numel(th), L(x >= th(j)) = c(j); end
end
