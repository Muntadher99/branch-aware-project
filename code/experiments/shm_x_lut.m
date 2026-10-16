function E = shm_x_lut(D, opts)
%SHM_X_LUT Stage E2: lookup tables on the modulation-index grid.
%   E = SHM_X_LUT(D, opts) compares, at R stored sets, the plain table (linear
%   interpolation of the per-Ma optimum), the nearest row, the plain table
%   with one corrector solve, the continuity-constrained path (no angle step
%   above CONT_DEG) and the branch-aware table (pchip inside the owning
%   branch), against the envelope over all branches. Also returns the split
%   into resolution and structural error, the corrector basin, and the lowest
%   mean THD49 of a path with at most k jumps (E.pareto).
%
%   Options: R ([10 25 50 100]), R_BASIN (25), CONT_DEG (5), CORRECT (true),
%            solve, verbose

if nargin < 2, opts = struct(); end
O.R        = [10 25 50 100];
O.R_BASIN  = 25;
O.CONT_DEG = 5;
O.CORRECT  = true;
O.solve    = struct();
O.verbose  = true;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

F = shm_model();  P = D.P;
grid = D.grid;  G = numel(grid);  N = P.N;
env = D.env;  envTH = D.envTH;  own = D.own;
fin = isfinite(env);

chg = [false; own(2:end) ~= own(1:end-1) & isfinite(own(2:end)) & isfinite(own(1:end-1))];
chgMa = 0.5 * (grid(find(chg) - 1) + grid(chg));
dist2chg = inf(G, 1);
for g = 1:G
    if ~isempty(chgMa), dist2chg(g) = min(abs(grid(g) - chgMa)); end
end

[cpath, cTH, cthd, nbreak] = continuity_path(D, O.CONT_DEG);
ccost = 100 * (cthd - env);
[PF, kenv] = continuity_pareto(D, O.CONT_DEG);

hdr = {'arm','rows_R','stored_sets','test_points','THD_excess_mean_pp','THD_excess_max_pp', ...
       'pct_points_over_0p1pp','angle_err_mean_deg','angle_err_max_deg','fundamental_err_mean_pct'};
rows = {};
split = zeros(0, 11);
splitpts = {};
basin = struct();

if O.verbose
    fprintf('  continuity-constrained path (<= %g deg per step): mean cost %.4f pp, worst %.4f pp, %d breaks\n', ...
        O.CONT_DEG, mean(ccost(fin & isfinite(ccost))), max(ccost(fin & isfinite(ccost))), nbreak);
    if ~isempty(PF)
        fprintf('  the per-Ma optimum jumps %d times; mean cost of a path allowed k jumps, k = 0..%d:\n    %s pp\n', ...
            kenv, PF(end, 1), mat2str(round(PF(:, 2).' * 1e3) / 1e3));
    end
    fprintf('\n  %-8s %5s %6s %12s %11s %10s %10s\n', 'arm', 'R', 'sets', 'excess mean', 'max', 'ang mean', 'ang max');
end

for R = O.R
    fg = find(fin);
    ri = unique(fg(round(linspace(1, numel(fg), R))));
    te = setdiff(fg, ri);
    Ya = interp1(grid(ri), envTH(ri, :), grid(te), 'linear');
    [ea, aa, fa] = score(F, P, Ya, grid(te), env(te), envTH(te, :));
    rows(end+1, :) = mkrow('LUT', R, numel(ri), numel(te), ea, aa, fa); %#ok<AGROW>
    prow(O, rows(end, :));

    dm = abs(repmat(grid(te), 1, numel(ri)) - repmat(grid(ri).', numel(te), 1));
    [~, jn] = min(dm, [], 2);
    Yn = envTH(ri(jn), :);
    [en, an, fn0] = score(F, P, Yn, grid(te), env(te), envTH(te, :));
    rows(end+1, :) = mkrow('NEAREST', R, numel(ri), numel(te), en, an, fn0); %#ok<AGROW>
    prow(O, rows(end, :));

    br = zeros(numel(te), 1);
    for t = 1:numel(te)
        lo = ri(find(ri < te(t), 1, 'last'));  hi = ri(find(ri > te(t), 1, 'first'));
        if isempty(lo) || isempty(hi), br(t) = NaN; continue; end
        seg = own(lo:hi);
        br(t) = any(seg ~= seg(1));
    end
    s0 = br == 0;  s1 = br == 1;
    a0 = abs(ea(s0));  a1 = abs(ea(s1));
    split(end+1, :) = [R, mean(ea(s0)), mean(ea(s1)), sum(s0), sum(s1), ...
        pctl(a0, 50), pctl(a1, 50), pctl(a0, 10), pctl(a1, 10), pctl(a0, 90), pctl(a1, 90)]; %#ok<AGROW>
    splitpts{end+1} = [grid(te), br, ea]; %#ok<AGROW>

    if O.CORRECT
        [Yc, Ic] = F.solve(Ya, grid(te), P, O.solve);
        bad = ~Ic.ok;  Yc(bad, :) = Ya(bad, :);
        [eb, ab, fb] = score(F, P, Yc, grid(te), env(te), envTH(te, :));
        rows(end+1, :) = mkrow('LUT+SQP', R, numel(ri), numel(te), eb, ab, fb); %#ok<AGROW>
        prow(O, rows(end, :));
        if R == O.R_BASIN
            landed = nan(numel(te), 1);
            for t = 1:numel(te)
                landed(t) = nearest_branch(D, te(t), Yc(t, :));
            end
            basin.Ma       = grid(te);
            basin.dist     = dist2chg(te);
            basin.excess   = eb;
            basin.excess_lut = ea;
            basin.onowner  = landed == own(te);
            basin.landed   = landed;
            basin.owner    = own(te);
            basin.R        = R;
            basin.ok       = Ic.ok;
        end
    end

    cf = find(isfinite(cthd));
    if numel(cf) < 2, continue; end
    rc = unique(cf(round(linspace(1, numel(cf), R))));
    tc = te(ismember(te, cf) & te >= min(rc) & te <= max(rc));
    Ycn = interp1(grid(rc), cTH(rc, :), grid(tc), 'linear');
    [ec, ac, fc] = score(F, P, Ycn, grid(tc), env(tc), envTH(tc, :));
    rows(end+1, :) = mkrow('CONT', R, numel(rc), numel(tc), ec, ac, fc); %#ok<AGROW>
    prow(O, rows(end, :));

    [Yd, nset] = branch_aware(D, R, te);
    ok = all(isfinite(Yd), 2);
    [ed, ad, fd] = score(F, P, Yd(ok, :), grid(te(ok)), env(te(ok)), envTH(te(ok), :));
    rows(end+1, :) = mkrow('BRANCH', R, nset, sum(ok), ed, ad, fd); %#ok<AGROW>
    prow(O, rows(end, :));
end

E.name   = 'E2_lookup_table';
E.header = hdr;
E.data   = rows;
E.split  = split;
E.split_header = {'R','same_owner_excess_mean_pp','straddling_excess_mean_pp','n_same','n_straddle', ...
    'same_owner_abs_excess_median_pp','straddling_abs_excess_median_pp','same_owner_abs_excess_p10_pp', ...
    'straddling_abs_excess_p10_pp','same_owner_abs_excess_p90_pp','straddling_abs_excess_p90_pp'};
E.split_points = splitpts;
E.basin  = basin;
E.cont   = struct('path', cpath, 'TH', cTH, 'thd', cthd, 'cost_pp', ccost, 'breaks', nbreak, ...
                  'deg', O.CONT_DEG);
E.pareto = PF;
E.pareto_header = {'max_jumps','mean_cost_pp','worst_cost_pp','jumps_used'};
E.kenv   = kenv;
E.chgMa  = chgMa;
E.note   = sprintf(['Deployed lookup table against the envelope of all traced branches. ' ...
    'LUT = per-Ma optimum stored at R rows, linear interpolation. NEAREST = the same rows held ' ...
    'without interpolation (its cost is the fundamental error). LUT+SQP = the same followed ' ...
    'by a warm-started fmincon SQP corrector. CONT = continuity-constrained table (no angle ' ...
    'moves more than %g deg between neighbouring grid points). BRANCH = branch-aware storage ' ...
    'with pchip inside each branch and the owning branch known. Excess is THD49 above the ' ...
    'envelope in percentage points; angle error is against the envelope optimum.'], O.CONT_DEG);
end

function [ex, ae, fe] = score(F, P, Y, Ma, envt, envTH)
Y  = sort(Y, 2);
R  = F.report(Y, Ma, P);
ex = 100 * (R.THD49 - envt(:));
ae = max(abs(Y - envTH), [], 2) * 180/pi;
fe = 100 * R.b1err ./ Ma(:);
end

function r = mkrow(arm, R, nset, nte, ex, ae, fe)
r = {arm, R, nset, nte, mean(ex), max(ex), 100*mean(ex > 0.1), mean(ae), max(ae), mean(fe)};
end

function prow(O, r)
if O.verbose
    fprintf('  %-8s %5d %6d %12.4f %11.4f %10.4f %10.4f\n', r{1}, r{2}, r{3}, r{5}, r{6}, r{8}, r{9});
end
end

function k = nearest_branch(D, g, x)
k = NaN;
X = squeeze(D.TH(g, :, :));
if size(X, 1) ~= numel(x), X = X.'; end
d = max(abs(X - repmat(x(:), 1, size(X, 2))), [], 1) * 180/pi;
[dm, j] = min(d);
if isfinite(dm) && dm <= 0.25, k = j; end
end

function [Y, nset] = branch_aware(D, R, te)
grid = D.grid;  N = D.P.N;
sel  = D.sel(:).';
own  = D.owns;
len  = arrayfun(@(k) sum(own == k), sel);
use  = sel(len > 0);  len = len(len > 0);
nb   = max(2, round(R * len / sum(len)));
Y    = nan(numel(te), N);
nset = 0;
tab  = cell(1, numel(use));
for u = 1:numel(use)
    k  = use(u);
    og = find(own == k);
    lo = max(1, og(1) - 1);  hi = min(numel(grid), og(end) + 1);
    ex = find(isfinite(D.THD(:, k)));
    lo = max(lo, ex(1));  hi = min(hi, ex(end));
    idx = unique(round(linspace(lo, hi, nb(u))));
    idx = idx(isfinite(D.THD(idx, k)));
    tab{u} = struct('Ma', grid(idx), 'TH', D.TH(idx, :, k));
    nset = nset + numel(idx);
end
for t = 1:numel(te)
    k = own(te(t));
    u = find(use == k, 1);
    if isempty(u), continue; end
    T = tab{u};
    if numel(T.Ma) >= 3
        Y(t, :) = interp1(T.Ma, T.TH, grid(te(t)), 'pchip');
    elseif numel(T.Ma) == 2
        Y(t, :) = interp1(T.Ma, T.TH, grid(te(t)), 'linear', 'extrap');
    end
end
end

function [path, TH, thd, nbreak] = continuity_path(D, degmax)
G = numel(D.grid);  K = size(D.THD, 2);  N = D.P.N;
lim = degmax * pi/180;  BIG = 1e3;
C = D.THD;  alive = isfinite(C);  C(~alive) = inf;
cost = inf(G, K);  prev = zeros(G, K);
cost(1, :) = C(1, :);
for g = 1:G-1
    for j = 1:K
        if ~alive(g+1, j), continue; end
        best = inf;  bi = 0;
        for i = 1:K
            if ~isfinite(cost(g, i)), continue; end
            if i == j
                c = cost(g, i);
            else
                d = max(abs(D.TH(g+1, :, j) - D.TH(g, :, i)));
                if d <= lim, c = cost(g, i); else, c = cost(g, i) + BIG; end
            end
            if c < best, best = c; bi = i; end
        end
        if bi == 0
            fin = cost(g, :);  [best, bi] = min(fin);
            if ~isfinite(best), best = 0; bi = 0; end
            best = best + BIG;
        end
        cost(g+1, j) = best + C(g+1, j);
        prev(g+1, j) = bi;
    end
end
path = nan(G, 1);
g = find(any(isfinite(cost), 2), 1, 'last');
if isempty(g), g = 0; j = 0; else, [~, j] = min(cost(g, :)); end
while g >= 1
    if j == 0 || ~isfinite(cost(g, j))
        [cm, j] = min(cost(g, :));
        if ~isfinite(cm), g = g - 1;  j = 0;  continue; end
    end
    path(g) = j;
    j = prev(g, j);
    g = g - 1;
end
TH = nan(G, N);  thd = nan(G, 1);  nbreak = 0;  last = 0;
for g = 1:G
    if isfinite(path(g))
        TH(g, :) = D.TH(g, :, path(g));
        thd(g)   = D.THD(g, path(g));
        if last > 0 && path(g) ~= path(last) && max(abs(TH(g, :) - TH(last, :))) > lim
            nbreak = nbreak + 1;
        end
        last = g;
    end
end
end

function [PF, kenv] = continuity_pareto(D, degmax)
lim = degmax * pi/180;
N = D.P.N;  K = size(D.THD, 2);
gi = find(isfinite(D.env));  T = numel(gi);
PF = zeros(0, 4);  kenv = 0;
if T < 2, return; end
C = D.THD(gi, :);  C(~isfinite(C)) = inf;
env = D.env(gi);  own = D.own(gi);
for t = 1:T-1
    if own(t) ~= own(t+1) && max(abs(D.TH(gi(t+1), :, own(t+1)) - D.TH(gi(t), :, own(t)))) > lim
        kenv = kenv + 1;
    end
end
nk = kenv + 1;
cost = inf(K, nk);  cost(:, 1) = C(1, :).';
bi = zeros(T, K, nk);  bb = zeros(T, K, nk);
for t = 1:T-1
    ii = find(isfinite(C(t, :)));  jj = find(isfinite(C(t+1, :)));
    Xa = reshape(D.TH(gi(t), :, ii), N, []);  Xb = reshape(D.TH(gi(t+1), :, jj), N, []);
    free = false(numel(ii), numel(jj));
    for q = 1:numel(jj)
        free(:, q) = (max(abs(bsxfun(@minus, Xa, Xb(:, q))), [], 1).' <= lim) | (ii(:) == jj(q));
    end
    new = inf(K, nk);
    for b = 1:nk
        S = repmat(cost(ii, b), 1, numel(jj));  S(~free) = inf;
        [m1, p1] = min(S, [], 1);
        if b > 1
            J = repmat(cost(ii, b - 1), 1, numel(jj));  J(free) = inf;
            [m2, p2] = min(J, [], 1);
        else
            m2 = inf(1, numel(jj));  p2 = ones(1, numel(jj));
        end
        use2 = m2 < m1;
        m = m1;  m(use2) = m2(use2);
        p = p1;  p(use2) = p2(use2);
        new(jj, b) = m(:) + C(t+1, jj).';
        bi(t+1, jj, b) = ii(p);
        bb(t+1, jj, b) = b - double(use2);
    end
    cost = new;
end
PF = nan(nk, 4);
for k = 0:kenv
    sub = cost(:, 1:k+1);
    [cm, lin] = min(sub(:));
    if ~isfinite(cm), continue; end
    [j, b] = ind2sub(size(sub), lin);
    path = zeros(T, 1);
    for t = T:-1:1
        path(t) = j;
        if t > 1
            jn = bi(t, j, b);  b = bb(t, j, b);  j = jn;
        end
    end
    used = 0;
    for t = 1:T-1
        if path(t) ~= path(t+1) && max(abs(D.TH(gi(t+1), :, path(t+1)) - D.TH(gi(t), :, path(t)))) > lim
            used = used + 1;
        end
    end
    thd = C(sub2ind(size(C), (1:T).', path));
    c = 100 * (thd - env);
    PF(k + 1, :) = [k, mean(c), max(c), used];
end
end

function v = pctl(x, p)
x = sort(x(isfinite(x)));
if isempty(x), v = NaN; return; end
if numel(x) == 1, v = x; return; end
v = interp1(1:numel(x), x, 1 + (numel(x) - 1) * p / 100);
end
