function X = shm_x_timing(D, E2b, opts)
%SHM_X_TIMING Stage TIM: operations, storage and time of a table read.
%   X = SHM_X_TIMING(D, E2b, opts) implements the read of the plain table
%   (binary search, linear interpolation) and of the branch-aware table
%   (search over the switch points and the branch rows, pchip evaluated by
%   Horner's rule from precomputed coefficients), counts their operations,
%   checks that they reproduce the angles of Table 1, and measures the median
%   time per read over NREP reads (machine dependent).
%
%   Options: R (100), NREP (20000)

if nargin < 3, opts = struct(); end
O.R = 100;  O.NREP = 20000;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end
P = D.P;  N = P.N;  grid = D.grid(:);

sel = D.sel(:).';  own = D.owns(:);
len = arrayfun(@(k) sum(own == k), sel);  use = sel(len > 0);  len = len(len > 0);
nb = max(2, round(O.R * len / sum(len)));
B = cell(1, numel(use));  nset = 0;
for u = 1:numel(use)
    k = use(u);  og = find(own == k);
    lo = max(1, og(1) - 1);  hi = min(numel(grid), og(end) + 1);
    ex = find(isfinite(D.THD(:, k)));  lo = max(lo, ex(1));  hi = min(hi, ex(end));
    idx = unique(round(linspace(lo, hi, nb(u))));  idx = idx(isfinite(D.THD(idx, k)));
    x = grid(idx);  Y = reshape(D.TH(idx, :, k), numel(idx), N);
    B{u} = struct('x', x, 'c', pchip_coefs(x, Y));
    nset = nset + numel(idx);
end
f = find(isfinite(own));  runs = [own(f(1)), grid(f(1))];
for g = f(1):f(end) - 1
    if isfinite(own(g)) && isfinite(own(g + 1)) && own(g) ~= own(g + 1)
        a = own(g);  b = own(g + 1);
        da = D.THD(g, a) - D.THD(g, b);  db = D.THD(g + 1, a) - D.THD(g + 1, b);
        if ~isfinite(D.THD(g + 1, a)), xc = grid(g);
        elseif ~isfinite(D.THD(g, b)), xc = grid(g + 1);
        elseif all(isfinite([da db])) && da ~= db, xc = grid(g) + da / (da - db) * (grid(g + 1) - grid(g));
        else, xc = grid(g + 1); end
        runs(end+1, :) = [b, xc]; %#ok<AGROW>
    end
end
sw = runs(2:end, 2);  rb = arrayfun(@(k) find(use == k, 1), runs(:, 1));
fg = find(isfinite(D.env));  ri = unique(fg(round(linspace(1, numel(fg), O.R))));
px = grid(ri);  PY = D.envTH(ri, :);

T = E2b.test;  Mt = T.Ma(:);
Yb = zeros(numel(Mt), N);  Ya = Yb;
for t = 1:numel(Mt)
    Yb(t, :) = ba_read(Mt(t), sw, rb, B);  Ya(t, :) = plain_read(Mt(t), px, PY);
end
X.check_deg = struct('branch_aware', max(abs(Yb(:) - T.Yb(:))) * 180/pi, 'plain', max(abs(Ya(:) - T.Ya(:))) * 180/pi);

nmax = max(cellfun(@(b) numel(b.x), B));
X.ops = struct( ...
    'plain_comparisons', ceil(log2(numel(px))), 'plain_mult', 1 + N, 'plain_add', 2 + 2 * N, ...
    'ba_comparisons', ceil(log2(numel(sw) + 1)) + ceil(log2(nmax)), 'ba_mult', 3 * N, 'ba_add', 1 + 3 * N, ...
    'note', ['plain: t = (Ma - x_i)/h_i, theta = y_i + t*(y_{i+1} - y_i); branch-aware: s = Ma - x_i, ' ...
             'theta = ((c3*s + c2)*s + c1)*s + c0 per angle']);
X.storage = struct('plain_numbers', N * numel(px), 'ba_sets', nset, 'ba_numbers', N * nset + 2 * size(runs, 1) - 1, ...
    'ba_numbers_with_coefficients', 4 * N * nset + 2 * size(runs, 1) - 1, 'switch_points', numel(sw));

mq = Mt(1) + (Mt(end) - Mt(1)) * rand(O.NREP, 1);
tb = zeros(O.NREP, 1);  tp = tb;
for r = 1:O.NREP
    t0 = tic;  ba_read(mq(r), sw, rb, B);  tb(r) = toc(t0);
    t0 = tic;  plain_read(mq(r), px, PY);  tp(r) = toc(t0);
end
X.time_us = struct('branch_aware_median', 1e6 * median(tb), 'branch_aware_p99', 1e6 * prctile_(tb, 99), ...
    'plain_median', 1e6 * median(tp), 'plain_p99', 1e6 * prctile_(tp, 99), ...
    'platform', platform(), 'note', 'interpreted code: an upper bound for compiled controller code');
X.corrector_ms = NaN;
if isfield(E2b, 'corr') && ~isempty(E2b.corr)
    c = E2b.corr(E2b.corr(:, 1) == O.R & E2b.corr(:, 2) == 2, :);
    if ~isempty(c), X.corrector_ms = c(1, 7); end
end
X.name = 'TIM_online_cost';
fprintf(['  deployed read reproduces Table 1 to %.1e deg (branch-aware), %.1e deg (plain); per read: branch-aware ' ...
         '%d comparisons + %d mult + %d add, %.2f us here; plain %d + %d + %d, %.2f us\n'], ...
    X.check_deg.branch_aware, X.check_deg.plain, X.ops.ba_comparisons, X.ops.ba_mult, X.ops.ba_add, ...
    X.time_us.branch_aware_median, X.ops.plain_comparisons, X.ops.plain_mult, X.ops.plain_add, X.time_us.plain_median);
end

function th = ba_read(m, sw, rb, B)
lo = 1;  hi = numel(sw) + 1;
while hi - lo > 0
    md = floor((lo + hi) / 2);
    if md <= numel(sw) && m >= sw(md), lo = md + 1; else, hi = md; end
end
b = B{rb(lo)};  x = b.x;
i = 1;  j = numel(x);
if numel(x) == 1, i = 1; m = x(1);
elseif m <= x(1), i = 1; elseif m >= x(end), i = numel(x) - 1;
else
    while j - i > 1, k = floor((i + j) / 2); if m >= x(k), i = k; else, j = k; end, end
end
s = m - x(i);  c = b.c(:, :, i);
th = ((c(1, :) * s + c(2, :)) * s + c(3, :)) * s + c(4, :);
end

function th = plain_read(m, x, Y)
i = 1;  j = numel(x);
if m <= x(1), i = 1; elseif m >= x(end), i = numel(x) - 1;
else
    while j - i > 1, k = floor((i + j) / 2); if m >= x(k), i = k; else, j = k; end, end
end
t = (m - x(i)) / (x(i + 1) - x(i));
th = Y(i, :) + t * (Y(i + 1, :) - Y(i, :));
end

function C = pchip_coefs(x, Y)
n = numel(x);  N = size(Y, 2);  C = zeros(4, N, max(1, n - 1));
if n == 1, C(4, :, 1) = Y(1, :); return; end
h = diff(x(:));
for c = 1:N
    y = Y(:, c);  del = diff(y) ./ h;  d = zeros(n, 1);
    if n == 2
        d(:) = del(1);
    else
        for k = 2:n-1
            if sign(del(k-1)) * sign(del(k)) > 0
                w1 = 2*h(k) + h(k-1);  w2 = h(k) + 2*h(k-1);
                d(k) = (w1 + w2) / (w1 / del(k-1) + w2 / del(k));
            end
        end
        d(1) = endslope(h(1), h(2), del(1), del(2));
        d(n) = endslope(h(n-1), h(n-2), del(n-1), del(n-2));
    end
    for k = 1:n-1
        c3 = (d(k) - 2*del(k) + d(k+1)) / h(k)^2;
        c2 = (3*del(k) - 2*d(k) - d(k+1)) / h(k);
        C(:, c, k) = [c3; c2; d(k); y(k)];
    end
end
end

function d = endslope(h1, h2, del1, del2)
d = ((2*h1 + h2) * del1 - h1 * del2) / (h1 + h2);
if sign(d) ~= sign(del1), d = 0;
elseif sign(del1) ~= sign(del2) && abs(d) > abs(3*del1), d = 3*del1; end
end

function v = prctile_(x, p)
x = sort(x(:));  v = x(max(1, min(numel(x), round(p / 100 * numel(x)))));
end

function s = platform()
if exist('OCTAVE_VERSION', 'builtin') ~= 0, s = ['Octave ' OCTAVE_VERSION]; else, s = ['MATLAB ' version]; end
s = [s ', ' computer];
end
