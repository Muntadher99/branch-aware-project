function X = shm_x_independent(D, E2b, opts)
%SHM_X_INDEPENDENT Stage IND: independent search at the 200 midpoints.
%   X = SHM_X_INDEPENDENT(D, E2b, opts) runs, at every midpoint of Table 1, a
%   multistart with fmincon's interior-point algorithm from NSTART scrambled
%   Halton starts, then two sweeps along the midpoints that pass the best
%   SWEEP_K optima to the next point. Reports whether the reference optimum is
%   found, how many starts reach it, and the gap of the best optimum found to
%   the reference. Every lower optimum is traced over the grid (SHM_TRACE) and
%   added to a union branch set, against which Table 1 is recomputed.
%
%   Options: NSTART (250), ALG ('interior-point'), TOL_DEG (0.05),
%            MATCH_DEG (0.05), LOWER_PP (0.01), SAME_DEG (0.10), SWEEP_K (5),
%            POINTS ([] = all), verbose

if nargin < 3, opts = struct(); end
O.NSTART    = 250;
O.ALG       = 'interior-point';
O.TOL_DEG   = 0.05;
O.MATCH_DEG = 0.05;
O.LOWER_PP  = 0.01;
O.SAME_DEG  = 0.10;
O.SWEEP_K   = 5;
O.POINTS    = [];
O.verbose   = true;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

F = shm_model();  P = D.P;  N = P.N;
grid = D.grid(:);  G = numel(grid);  K = size(D.THD, 2);
ref = E2b.ref;  Mm = ref.Ma(:);  env = ref.env(:);  own = ref.own(:);
pts = find(isfinite(env));
if ~isempty(O.POINTS), pts = intersect(pts, O.POINTS(:)); end
so = struct('algorithm', O.ALG);
tol = O.TOL_DEG * pi/180;  mtol = O.MATCH_DEG * pi/180;
shift = rand(1, N);
t0 = tic;  F.nfe_reset();

hdr = {'Ma', 'reference_THD49_pct', 'best_found_THD49_pct', 'best_minus_reference_pp', 'converged_starts', ...
       'distinct_optima', 'found_once', 'reference_found', 'starts_reaching_reference', ...
       'pct_starts_reaching_reference', 'reference_owner', 'best_with_sweeps_THD49_pct', ...
       'best_with_sweeps_minus_reference_pp', 'reference_found_with_sweeps'};
rows = nan(numel(pts), numel(hdr));
REP = cell(numel(pts), 1);  RT = REP;
if O.verbose
    fprintf('  independent check: fmincon %s, %d Halton starts at each of %d midpoints\n', O.ALG, O.NSTART, numel(pts));
end
for ii = 1:numel(pts)
    t = pts(ii);  m = Mm(t);
    H = halton(O.NSTART, N, 1 + (ii - 1) * O.NSTART, shift);
    span = (P.ub - P.lb) - (N - 1) * P.minGap;
    X0 = P.lb + sort(H, 2) * span + repmat(P.minGap * (0:N-1), O.NSTART, 1);
    [Y, I] = F.solve(X0, m, P, so);
    ok = I.ok;  Y = Y(ok, :);  th = I.thd(ok);
    [th, o] = sort(th);  Y = Y(o, :);
    rep = zeros(0, N);  rt = zeros(0, 1);  cnt = zeros(0, 1);
    for j = 1:size(Y, 1)
        if ~isempty(rep)
            [dm, q] = min(max(abs(rep - repmat(Y(j, :), size(rep, 1), 1)), [], 2));
            if dm <= tol, cnt(q) = cnt(q) + 1; continue; end
        end
        rep(end+1, :) = Y(j, :);  rt(end+1, 1) = th(j);  cnt(end+1, 1) = 1; %#ok<AGROW>
    end
    best = NaN;  if ~isempty(rt), best = rt(1); end
    a = ANm_at(ref, t, own(t), N);
    match = false(size(rt));
    if ~isempty(rt)
        match = max(abs(rep - repmat(a, size(rep, 1), 1)), [], 2) <= mtol | abs(rt - env(t)) <= 1e-7;
    end
    hits = sum(cnt(match));
    gap = 100 * (best - env(t));
    rows(ii, 1:11) = [m, 100 * env(t), 100 * best, gap, sum(ok), numel(rt), sum(cnt == 1), double(any(match)), ...
                      hits, 100 * hits / max(sum(ok), 1), own(t)];
    REP{ii} = rep;  RT{ii} = rt;
    if O.verbose && (mod(ii, 20) == 0 || ii == numel(pts))
        fprintf('    %3d/%d midpoints, %.0f s: reference found at %.1f%%, lower optimum at %d\n', ii, numel(pts), ...
            toc(t0), 100 * mean(rows(1:ii, 8)), sum(rows(1:ii, 4) < -O.LOWER_PP));
    end
end

PR = repmat({zeros(0, N)}, numel(pts), 1);  PT = repmat({zeros(0, 1)}, numel(pts), 1);
nsweep = 0;
if O.SWEEP_K > 0
    for dirn = [1 -1]
        ord = 1:numel(pts);  if dirn < 0, ord = fliplr(ord); end
        carry = zeros(0, N);
        for ii = ord
            if ~isempty(carry)
                [Y, I] = F.solve(carry, Mm(pts(ii)), P, so);  nsweep = nsweep + size(carry, 1);
                PR{ii} = [PR{ii}; Y(I.ok, :)];  PT{ii} = [PT{ii}; I.thd(I.ok)];
            end
            cr = distinct([REP{ii}; PR{ii}], [RT{ii}; PT{ii}], tol);
            carry = cr(1:min(O.SWEEP_K, size(cr, 1)), :);
        end
    end
end
lowers = zeros(0, N + 2);
for ii = 1:numel(pts)
    t = pts(ii);
    [cr, ct] = distinct([REP{ii}; PR{ii}], [RT{ii}; PT{ii}], tol);
    b = NaN;  f = 0;
    if ~isempty(ct)
        a = ANm_at(ref, t, own(t), N);
        b = ct(1);  f = double(any(max(abs(cr - repmat(a, size(cr, 1), 1)), [], 2) <= mtol | abs(ct - env(t)) <= 1e-7));
    end
    g = 100 * (b - env(t));
    rows(ii, 12:14) = [100 * b, g, f];
    if isfinite(g) && g < -O.LOWER_PP
        lowers(end+1, :) = [t, g, cr(1, :)]; %#ok<AGROW>
    end
end
if O.verbose && O.SWEEP_K > 0
    fprintf('    sweeps (%d solves): reference found at %.1f%% of the midpoints, lower optimum at %d\n', nsweep, ...
        100 * mean(rows(:, 14)), size(lowers, 1));
end
X.name = 'IND_independent_check';
X.header = hdr;  X.data = rows;
X.nfe = F.nfe();  X.time_s = toc(t0);
fnd = rows(:, 8) == 1;
X.summary = struct('points', numel(pts), 'starts_per_point', O.NSTART, 'algorithm', O.ALG, ...
    'converged_pct', 100 * sum(rows(:, 5)) / (O.NSTART * numel(pts)), ...
    'found_pct', 100 * mean(fnd), 'not_found', sum(~fnd), ...
    'hit_pct_median', median(rows(:, 10)), 'hit_pct_min', min(rows(:, 10)), 'hit_pct_max', max(rows(:, 10)), ...
    'points_hit_below_1pct', sum(rows(:, 10) < 1), ...
    'distinct_median', median(rows(:, 6)), 'distinct_max', max(rows(:, 6)), ...
    'once_pct_median', median(100 * rows(:, 7) ./ max(rows(:, 6), 1)), ...
    'multistart_lower_points', sum(rows(:, 4) < -O.LOWER_PP), ...
    'not_found_gap_min_pp', min([inf; rows(~fnd, 4)]), 'not_found_gap_max_pp', max([-inf; rows(~fnd, 4)]), ...
    'sweep_k', O.SWEEP_K, 'sweep_solves', nsweep, 'found_with_sweeps_pct', 100 * mean(rows(:, 14)), ...
    'not_found_with_sweeps', sum(rows(:, 14) == 0), ...
    'not_found_with_sweeps_gap_max_pp', max([-inf; rows(rows(:, 14) == 0, 13)]), ...
    'lower_points', size(lowers, 1), 'lower_gap_max_pp', abs(min([0; lowers(:, 2)])));

uo = unique(own(pts(isfinite(own(pts))))).';
X.by_owner = zeros(numel(uo), 5);
for q = 1:numel(uo)
    s = rows(:, 11) == uo(q);
    X.by_owner(q, :) = [uo(q), sum(s), 100 * mean(rows(s, 8)), median(rows(s, 10)), 100 * mean(rows(s, 14))];
end
X.by_owner_header = {'branch', 'midpoints_owned', 'found_pct', 'median_pct_starts', 'found_with_sweeps_pct'};

newTH = nan(G, N, 0);  newTHD = nan(G, 0);  newsrc = zeros(0, 2);
same = O.SAME_DEG * pi/180;
for j = 1:size(lowers, 1)
    t = lowers(j, 1);  x = lowers(j, 3:end);
    [~, gi] = min(abs(grid - Mm(t)));
    if on_any(D.TH, gi, x, same) || on_any(newTH, gi, x, same), continue; end
    B = shm_trace(x, Mm(t), grid, P, D.opts.trace);
    if isnan(B.lo), continue; end
    if on_any(D.TH, B.lo, B.TH(B.lo, :), same), continue; end
    newTH(:, :, end+1) = B.TH;  newTHD(:, end+1) = B.thd;  newsrc(end+1, :) = [Mm(t), lowers(j, 2)]; %#ok<AGROW>
    if O.verbose
        fprintf('    new branch from Ma %.5f (%.4f pp below the reference): Ma %.4f..%.4f\n', Mm(t), lowers(j, 2), ...
            grid(B.lo), grid(B.hi));
    end
end
X.new = struct('TH', newTH, 'THD', newTHD, 'from', newsrc);

THu = [D.THD, newTHD];
[envu, ownu] = min(THu, [], 2);  ownu(~isfinite(envu)) = NaN;
v = isfinite(ownu(1:end-1)) & isfinite(ownu(2:end));
gain = 100 * (D.env(:) - envu);  gain(~isfinite(gain)) = 0;
newown = unique(ownu(ownu > K)).';
U.branches = size(THu, 2);  U.new_branches = size(newTHD, 2);
U.owners = numel(unique(ownu(isfinite(ownu))));  U.owner_changes = sum(ownu(1:end-1) ~= ownu(2:end) & v);
U.new_owners = numel(newown);
U.gain_max_pp = max(gain);  U.gain_mean_pp = mean(gain(isfinite(D.env)));
U.new_owned_points = sum(ismember(ownu, newown));
og = find(ismember(ownu, newown));
if isempty(og), U.new_owned_Ma = [NaN NaN]; else, U.new_owned_Ma = [grid(og(1)), grid(og(end))]; end
U.paper_owners = D.stats.n_owners;  U.paper_owner_changes = D.stats.owner_changes;
X.union = U;

T = E2b.test;  [~, it] = ismember(round(T.Ma(:) * 1e6), round(Mm * 1e6));
envT = T.env(:);  envU = envT;
for k = 1:size(newTHD, 2)
    for q = 1:numel(it)
        g = find(grid < Mm(it(q)), 1, 'last');
        if isempty(g) || g >= G || ~all(isfinite(newTHD([g g+1], k))), continue; end
        x0 = 0.5 * (newTH(g, :, k) + newTH(g + 1, :, k));
        [x, I] = F.solve(x0, Mm(it(q)), P);
        if I.ok && max(abs(x - sort(x0))) <= pi/180, envU(q) = min(envU(q), I.thd); end
    end
end
Ra = F.report(sort(T.Ya, 2), T.Ma(:), P);  Rb = F.report(sort(T.Yb, 2), T.Ma(:), P);
ea0 = 100 * (Ra.THD49 - envT);  eb0 = 100 * (Rb.THD49 - envT);
ea1 = 100 * (Ra.THD49 - envU);  eb1 = 100 * (Rb.THD49 - envU);
X.table1 = struct('R', T.R, 'points', numel(envT), 'ref_lowered_points', sum(envU < envT - 1e-12), ...
    'ref_drop_max_pp', 100 * max(envT - envU), ...
    'plain_mean_pp', [mean(ea0), mean(ea1)], 'plain_max_pp', [max(ea0), max(ea1)], ...
    'ba_mean_pp', [mean(eb0), mean(eb1)], 'ba_max_pp', [max(eb0), max(eb1)], ...
    'note', 'first value: reference of the paper, second: union reference');
X.opts = O;
if O.verbose
    s = X.summary;
    fprintf(['  multistart: reference found at %d of %d midpoints (%.1f%%); starts reaching it: median %.2f%%, min %.2f%%;\n' ...
             '  with the sweeps: found at %.1f%% of the midpoints;\n' ...
             '  lower optimum at %d midpoints (by up to %.4f pp); union: %d owners, %d owner changes (paper %d, %d);\n' ...
             '  Table 1 at R = %d with the union reference: plain mean %.3f -> %.3f pp, BA mean %.4f -> %.4f pp (%.0f s)\n'], ...
        s.points - s.not_found, s.points, s.found_pct, s.hit_pct_median, s.hit_pct_min, s.found_with_sweeps_pct, s.lower_points, ...
        s.lower_gap_max_pp, U.owners, U.owner_changes, U.paper_owners, U.paper_owner_changes, T.R, ...
        X.table1.plain_mean_pp(1), X.table1.plain_mean_pp(2), X.table1.ba_mean_pp(1), X.table1.ba_mean_pp(2), X.time_s);
end
end

function [rep, rt] = distinct(Y, th, tol)
[th, o] = sort(th(:));  Y = Y(o, :);
rep = zeros(0, size(Y, 2));  rt = zeros(0, 1);
for j = 1:size(Y, 1)
    if ~isempty(rep) && min(max(abs(rep - repmat(Y(j, :), size(rep, 1), 1)), [], 2)) <= tol, continue; end
    rep(end+1, :) = Y(j, :);  rt(end+1, 1) = th(j); %#ok<AGROW>
end
end

function a = ANm_at(ref, t, k, N)
a = nan(1, N);
if isfinite(k), a = reshape(ref.ANm(t, :, k), 1, N); end
end

function tf = on_any(TH, gi, x, same)
tf = false;
if isempty(TH) || size(TH, 3) == 0, return; end
R = reshape(TH(gi, :, :), size(TH, 2), size(TH, 3));
d = max(abs(R - repmat(x(:), 1, size(R, 2))), [], 1);
tf = any(d <= same);
end

function U = halton(K, d, skip, shift)
p = primes(200);  p = p(1:d);
idx = (skip:skip + K - 1).';
U = zeros(K, d);
for j = 1:d
    b = p(j);  n = idx;  f = 1 / b;  v = zeros(K, 1);
    while any(n > 0)
        v = v + f * mod(n, b);  n = floor(n / b);  f = f / b;
    end
    U(:, j) = v;
end
U = mod(U + repmat(shift, K, 1), 1);
end
