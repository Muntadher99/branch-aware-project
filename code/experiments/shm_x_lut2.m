function E = shm_x_lut2(R0, opts)
%SHM_X_LUT2 Stage E2b: the table comparison at the 200 midpoints (Table 1).
%   E = SHM_X_LUT2(R, opts) scores every arm at the midpoints between grid
%   points against the lowest THD49 over all branches there (each branch
%   evaluated by one local solve from the mean of its neighbouring grid
%   solutions, accepted if no angle moves more than MOVE deg). Arms: plain,
%   nearest row, plain + corrector, continuous, sweep, branch-aware pchip,
%   branch-aware linear, branch-aware + corrector. Scored by THD49 excess (pp)
%   and fundamental error (%).
%
%   Options: R ([10 25 50 100]), CONT_DEG (5), MOVE (1.0), R_BASIN (25),
%            R_TEST (100), verbose

if nargin < 2, opts = struct(); end
O.R        = [10 25 50 100];
O.CONT_DEG = 5;
O.MOVE     = 1.0;
O.R_BASIN  = 25;
O.R_TEST   = 100;
O.verbose  = true;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

D = R0.D;  F = shm_model();  P = D.P;  N = P.N;
grid = D.grid(:);  G = numel(grid);  K = size(D.THD, 2);
Mm = 0.5 * (grid(1:end-1) + grid(2:end));  nm = numel(Mm);

THm = nan(nm, K);  ANm = nan(nm, N, K);
tic;
for k = 1:K
    g = find(isfinite(D.THD(1:end-1, k)) & isfinite(D.THD(2:end, k)));
    if isempty(g), continue; end
    f = find(isfinite(D.THD(:, k)));
    T0 = 0.5 * (D.TH(g, :, k) + D.TH(g + 1, :, k));
    Tp = T0;
    for i = 1:numel(g)
        idx = f(abs(f - g(i)) <= 3);
        if numel(idx) >= 3, Tp(i, :) = interp1(grid(idx), reshape(D.TH(idx, :, k), numel(idx), N), Mm(g(i)), 'pchip'); end
    end
    [Th, I] = F.solve(Tp, Mm(g), P);
    mv = max(abs(Th - sort(Tp, 2)), [], 2) * 180/pi;
    okk = I.ok & mv <= O.MOVE;
    if any(~okk)
        r = find(~okk);
        [Th2, I2] = F.solve(T0(r, :), Mm(g(r)), P);
        mv2 = max(abs(Th2 - sort(T0(r, :), 2)), [], 2) * 180/pi;
        ok2 = I2.ok & mv2 <= O.MOVE;
        Th(r(ok2), :) = Th2(ok2, :);  I.thd(r(ok2)) = I2.thd(ok2);  okk(r(ok2)) = true;
    end
    THm(g(okk), k) = I.thd(okk);  ANm(g(okk), :, k) = Th(okk, :);
end
[envm, ownm] = min(THm, [], 2);
fin = isfinite(envm);
sel = D.sel(:).';
[~, j] = min(THm(:, sel), [], 2);  ownsm = sel(j).';
ownsm(~any(isfinite(THm(:, sel)), 2)) = NaN;
storedgap = 100 * (min(THm(:, sel), [], 2) - envm);
own = D.owns(:);  lab = nan(nm, 1);
for i = 1:nm
    a = own(i);  b = own(i + 1);
    if isfinite(a) && isfinite(b) && a ~= b
        da = D.THD(i, a) - D.THD(i, b);  db = D.THD(i + 1, a) - D.THD(i + 1, b);
        if ~isfinite(D.THD(i + 1, a))
            lab(i) = b;
        elseif ~isfinite(D.THD(i, b))
            lab(i) = a;
        elseif all(isfinite([da db])) && da ~= db
            xc = grid(i) + da / (da - db) * (grid(i + 1) - grid(i));
            if Mm(i) < xc, lab(i) = a; else, lab(i) = b; end
        else
            lab(i) = a;
        end
    elseif isfinite(a), lab(i) = a;
    else, lab(i) = b;
    end
end
ownsm = lab;
if O.verbose
    fprintf('  reference at %d midpoints in %.0f s; stored subset above the full envelope: mean %.4f, worst %.4f pp\n', ...
        sum(fin), toc, mean(storedgap(fin & isfinite(storedgap))), max(storedgap(fin)));
end

gs = find(isfinite(D.env), 1);
SW = nan(G, N);  SW(gs, :) = D.envTH(gs, :);
swfail = 0;
for g = gs + 1:G
    [x, I] = F.solve(SW(g - 1, :), grid(g), P);
    if I.ok, SW(g, :) = x; else, SW(g, :) = SW(g - 1, :);  swfail = swfail + 1; end
end
Rs = F.report(sort(SW(gs:end, :), 2), grid(gs:end), P);
swcost = 100 * (Rs.THD49(:) - D.env(gs:end));
swjump = sum(max(abs(diff(SW(gs:end, :))), [], 2) * 180/pi > O.CONT_DEG);
if O.verbose
    fprintf('  sweep from Ma = %.3f: cost over the grid mean %.3f, worst %.3f pp, %d jumps > %g deg, %d failed solves\n', ...
        grid(gs), mean(swcost), max(swcost), swjump, O.CONT_DEG, swfail);
end
cTH = R0.E2.cont.TH;

hdr = {'arm', 'R', 'stored_sets', 'test_points', 'excess_mean_pp', 'excess_max_pp', 'pct_over_0p1pp', ...
       'fund_err_mean_pct', 'fund_err_max_pct', 'worst_at_Ma'};
rows = {};  worst = struct();  split = zeros(0, 11);  corr = zeros(0, 9);  basin = struct();  test = struct();
own_full = D.own(:);  h0 = grid(2) - grid(1);
chg = find(own_full(2:end) ~= own_full(1:end-1) & isfinite(own_full(2:end)) & isfinite(own_full(1:end-1)));
chgMa = 0.5 * (grid(chg) + grid(chg + 1));
te = find(fin);  Mt = Mm(te);  Et = envm(te);
Rlist = O.R;
for R = Rlist
    fg = find(isfinite(D.env));
    ri = unique(fg(round(linspace(1, numel(fg), R))));
    Ya = interp1(grid(ri), D.envTH(ri, :), Mt, 'linear');
    [sa, ea] = sc(F, P, Ya, Mt, Et);
    rows(end+1, :) = mk('plain', R, numel(ri), sa); %#ok<AGROW>
    br = nan(numel(te), 1);
    for t = 1:numel(te)
        lo = ri(find(ri <= te(t), 1, 'last'));  hi = ri(find(ri >= te(t) + 1, 1, 'first'));
        if isempty(lo) || isempty(hi), continue; end
        sg = own_full(lo:hi);  sg = sg(isfinite(sg));
        br(t) = any(sg ~= sg(1));
    end
    s0 = br == 0;  s1 = br == 1;  a0 = abs(ea(s0));  a1 = abs(ea(s1));
    split(end+1, :) = [R, mean(ea(s0)), mean(ea(s1)), sum(s0), sum(s1), pctl(a0, 50), pctl(a1, 50), ...
        pctl(a0, 10), pctl(a1, 10), pctl(a0, 90), pctl(a1, 90)]; %#ok<AGROW>
    [~, jn] = min(abs(Mt - grid(ri).'), [], 2);
    rows(end+1, :) = mk('nearest', R, numel(ri), sc(F, P, D.envTH(ri(jn), :), Mt, Et)); %#ok<AGROW>
    tc = tic;  [Yc, Ic] = F.solve(Ya, Mt, P);  tcall = toc(tc) / numel(Mt);
    Yc(~Ic.ok, :) = Ya(~Ic.ok, :);
    [sc1, ec] = sc(F, P, Yc, Mt, Et);
    rows(end+1, :) = mk('plain+corr', R, numel(ri), sc1); %#ok<AGROW>
    corr(end+1, :) = [R, 1, numel(Mt), sum(~Ic.ok), mean(ec(Ic.ok)), max(ec(Ic.ok)), 1e3 * tcall, ...
        mean(Ic.iter), mean(Ic.iter(Ic.ok))]; %#ok<AGROW>
    if R == O.R_BASIN
        land = nan(numel(te), 1);
        for t = 1:numel(te)
            X = reshape(ANm(te(t), :, :), N, K);
            dd = max(abs(X - repmat(sort(Yc(t, :)).', 1, K)), [], 1) * 180/pi;
            [dm, jb] = min(dd);  if dm <= 0.25, land(t) = jb; end
        end
        dist = inf(numel(te), 1);
        if ~isempty(chgMa), dist = min(abs(repmat(Mt, 1, numel(chgMa)) - repmat(chgMa(:).', numel(Mt), 1)), [], 2); end
        on = land == ownm(te);  near = dist <= 2 * h0;
        basin = struct('R', R, 'Ma', Mt, 'onowner', on, 'near', near, 'excess', ec, 'ok', Ic.ok(:), ...
            'onowner_far_pct', 100 * mean(on(~near)), 'onowner_near_pct', 100 * mean(on(near)), ...
            'n_far', sum(~near), 'n_near', sum(near), 'excess_off_owner_pp', mean(ec(~on)), ...
            'fails', sum(~Ic.ok));
    end
    cf = find(all(isfinite(cTH), 2));
    rc = unique(cf(round(linspace(1, numel(cf), R))));
    ok = Mt >= grid(rc(1)) & Mt <= grid(rc(end));
    Yn = interp1(grid(rc), cTH(rc, :), Mt(ok), 'linear');
    rows(end+1, :) = mk('continuous', R, numel(rc), sc(F, P, Yn, Mt(ok), Et(ok))); %#ok<AGROW>
    sf = (gs:G).';
    rs = unique(sf(round(linspace(1, numel(sf), R))));
    ok = Mt >= grid(rs(1)) & Mt <= grid(rs(end));
    Ys = interp1(grid(rs), SW(rs, :), Mt(ok), 'linear');
    rows(end+1, :) = mk('sweep', R, numel(rs), sc(F, P, Ys, Mt(ok), Et(ok))); %#ok<AGROW>
    [Yp, nset, okb] = branch_aware(D, R, Mt, ownsm(te), 'pchip');
    [s1, w1] = sc(F, P, Yp(okb, :), Mt(okb), Et(okb));
    rows(end+1, :) = mk('BA pchip', R, nset, s1); %#ok<AGROW>
    [Yl, ~, okl] = branch_aware(D, R, Mt, ownsm(te), 'linear');
    rows(end+1, :) = mk('BA linear', R, nset, sc(F, P, Yl(okl, :), Mt(okl), Et(okl))); %#ok<AGROW>
    tc = tic;  [Ybc, Ib] = F.solve(Yp(okb, :), Mt(okb), P);  tcall = toc(tc) / numel(okb);
    Ybc(~Ib.ok, :) = Yp(okb(~Ib.ok), :);
    [sb, eb] = sc(F, P, Ybc, Mt(okb), Et(okb));
    rows(end+1, :) = mk('BA+corr', R, nset, sb); %#ok<AGROW>
    corr(end+1, :) = [R, 2, numel(okb), sum(~Ib.ok), mean(eb(Ib.ok)), max(eb(Ib.ok)), 1e3 * tcall, ...
        mean(Ib.iter), mean(Ib.iter(Ib.ok))]; %#ok<AGROW>
    if R == O.R_TEST
        ti = te(okb);  Yo = nan(numel(ti), N);
        for t = 1:numel(ti), Yo(t, :) = ANm(ti(t), :, ownm(ti(t))); end
        test = struct('R', R, 'Ma', Mt(okb), 'env', Et(okb), 'Yopt', Yo, 'Ya', Ya(okb, :), 'Yb', Yp(okb, :));
    end
    tt = te(okb);  [~, o] = sort(w1, 'descend');
    W = zeros(min(5, numel(o)), 8);
    for i = 1:size(W, 1)
        y = sort(Yp(okb(o(i)), :));  yo = sort(ANm(tt(o(i)), :, ownm(tt(o(i)))));
        W(i, :) = [Mm(tt(o(i))), w1(o(i)), ownsm(tt(o(i))), min(diff(y)) * 180/pi, y(1) * 180/pi, ...
            90 - y(end) * 180/pi, sum(diff(yo) < P.minGap + 1e-3*pi/180), min([inf; abs(Mm(tt(o(i))) - chgMa(:))]) / h0];
    end
    worst.(sprintf('R%d', R)) = W;
    if O.verbose
        r = rows(end-7:end, :);
        for i = 1:size(r, 1)
            fprintf('  R=%3d %-11s sets %3d  excess mean %8.4f  max %7.3f pp  >0.1pp %5.1f%%  fund mean %.4f max %.4f %%\n', ...
                R, r{i, 1}, r{i, 3}, r{i, 5}, r{i, 6}, r{i, 7}, r{i, 8}, r{i, 9});
        end
    end
end
ns = unique(cell2mat(rows(strcmp(rows(:, 1), 'BA pchip'), 3))).';
for R = setdiff(ns, Rlist)
    ri = unique(fg(round(linspace(1, numel(fg), R))));
    Ya = interp1(grid(ri), D.envTH(ri, :), Mt, 'linear');
    rows(end+1, :) = mk('plain', R, numel(ri), sc(F, P, Ya, Mt, Et)); %#ok<AGROW>
    if O.verbose, fprintf('  R=%3d plain (matched storage) excess mean %.4f max %.3f pp\n', R, rows{end, 5}, rows{end, 6}); end
end

E.name = 'E2b_lookup_table_midpoints';  E.header = hdr;  E.data = rows;
E.worst = worst;  E.worst_header = {'Ma', 'excess_pp', 'owner', 'min_gap_deg', 'theta1_deg', '90_minus_thetaN_deg', ...
    'active_gaps_of_the_optimum', 'distance_to_owner_change_grid_steps'};
E.split = split;
E.split_header = {'R','same_owner_excess_mean_pp','straddling_excess_mean_pp','n_same','n_straddle', ...
    'same_owner_abs_excess_median_pp','straddling_abs_excess_median_pp','same_owner_abs_excess_p10_pp', ...
    'straddling_abs_excess_p10_pp','same_owner_abs_excess_p90_pp','straddling_abs_excess_p90_pp'};
E.corr = corr;
E.corr_header = {'R', 'start (1 plain table, 2 branch-aware table)', 'calls', 'failed', 'excess_mean_converged_pp', ...
    'excess_max_converged_pp', 'ms_per_call', 'iterations_mean', 'iterations_mean_converged'};
E.basin = basin;
E.test = test;
ba = rows(strcmp(rows(:, 1), 'BA pchip') & cell2mat(rows(:, 2)) == O.R_TEST, :);
if ~isempty(ba)
    lo = min(Mt);  hi = max(Mt);  span = hi - lo;
    hw = 2 * (ba{1, 9} / 100) * lo;
    hm = 4 * (ba{1, 8} / 100) / mean(1 ./ Mt);
    E.nearest_rows = struct('for_worst', ceil(span / hw) + 1, 'for_mean', ceil(span / hm) + 1, ...
        'target_worst_pct', ba{1, 9}, 'target_mean_pct', ba{1, 8}, 'R', O.R_TEST);
else
    E.nearest_rows = struct('for_worst', NaN, 'for_mean', NaN, 'target_worst_pct', NaN, 'target_mean_pct', NaN, 'R', O.R_TEST);
end
E.ref = struct('Ma', Mm, 'env', envm, 'own', ownm, 'owns', ownsm, 'stored_gap_pp', storedgap, 'THm', THm, 'ANm', ANm);
E.sweep = struct('TH', SW, 'cost_pp', swcost, 'jumps', swjump, 'failed', swfail);
E.opts = O;
end

function [s, ex] = sc(F, P, Y, Ma, env)
Y = sort(Y, 2);
R = F.report(Y, Ma, P);
ex = 100 * (R.THD49(:) - env(:));
fe = 100 * R.b1err(:) ./ Ma(:);
[~, iw] = max(ex);
s = [numel(ex), mean(ex), max(ex), 100 * mean(ex > 0.1), mean(fe), max(fe), Ma(iw)];
end

function r = mk(arm, R, nset, s)
r = {arm, R, nset, s(1), s(2), s(3), s(4), s(5), s(6), s(7)};
end

function [Y, nset, ok] = branch_aware(D, R, Mq, ownq, meth)
grid = D.grid(:);  N = D.P.N;
sel  = D.sel(:).';  own = D.owns(:);
len  = arrayfun(@(k) sum(own == k), sel);
use  = sel(len > 0);  len = len(len > 0);
nb   = max(2, round(R * len / sum(len)));
nset = 0;  tab = cell(1, numel(use));
for u = 1:numel(use)
    k  = use(u);  og = find(own == k);
    lo = max(1, og(1) - 1);  hi = min(numel(grid), og(end) + 1);
    ex = find(isfinite(D.THD(:, k)));  lo = max(lo, ex(1));  hi = min(hi, ex(end));
    idx = unique(round(linspace(lo, hi, nb(u))));
    idx = idx(isfinite(D.THD(idx, k)));
    tab{u} = struct('Ma', grid(idx), 'TH', reshape(D.TH(idx, :, k), numel(idx), N));
    nset = nset + numel(idx);
end
Y = nan(numel(Mq), N);
for t = 1:numel(Mq)
    u = find(use == ownq(t), 1);
    if isempty(u), continue; end
    T = tab{u};
    if numel(T.Ma) >= 3 && strcmp(meth, 'pchip')
        Y(t, :) = interp1(T.Ma, T.TH, Mq(t), 'pchip', 'extrap');
    elseif numel(T.Ma) >= 2
        Y(t, :) = interp1(T.Ma, T.TH, Mq(t), 'linear', 'extrap');
    end
end
ok = find(all(isfinite(Y), 2));
end

function v = pctl(x, p)
x = sort(x(isfinite(x)));
if isempty(x), v = NaN; return; end
if numel(x) == 1, v = x; return; end
v = interp1(1:numel(x), x, 1 + (numel(x) - 1) * p / 100);
end
