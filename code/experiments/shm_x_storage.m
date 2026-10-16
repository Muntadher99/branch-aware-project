function X = shm_x_storage(D, E2, E2b, opts)
%SHM_X_STORAGE Stage STO: Table 1 at equal storage, with published baselines.
%   X = SHM_X_STORAGE(D, E2, E2b, opts). A branch-aware table of nominal size R
%   stores its angle sets plus 2*runs - 1 numbers of switch metadata; every
%   other arm gets ceil(numbers/15) rows. Arms, scored at the 200 midpoints:
%   branch-aware (pchip, linear, + corrector), plain, nearest row, plain +
%   corrector, sweep, continuous, single continuous branch, and a learned
%   angle generator (SHM_TANHLM, no more parameters than the table's numbers)
%   with and without one corrector solve.
%
%   Options: R ([10 25 50 100]), SEEDS (1:5), NN_WIDTHS ([2 4 8 16 32 64]),
%            NN_RESTARTS (2), NN_CORR (true), verbose

if nargin < 4, opts = struct(); end
O.R           = [10 25 50 100];
O.SEEDS       = 1:5;
O.NN_WIDTHS   = [2 4 8 16 32 64];
O.NN_RESTARTS = 2;
O.NN_CORR     = true;
O.verbose     = true;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

F = shm_model();  P = D.P;  N = P.N;
grid = D.grid(:);  G = numel(grid);
ref = E2b.ref;  Mm = ref.Ma(:);  envm = ref.env(:);
te = find(isfinite(envm));  Mt = Mm(te);  Et = envm(te);  labt = ref.owns(te);
fg = find(isfinite(D.env));

own = D.owns(:);  fo = own(isfinite(own));
runs = 1 + sum(fo(2:end) ~= fo(1:end-1));
meta = 2 * runs - 1;

cov = sum(isfinite(D.THD), 1);
cand = find(cov == max(cov));
mth = arrayfun(@(k) mean(D.THD(isfinite(D.THD(:, k)), k)), cand);
[~, q] = min(mth);  kb = cand(q);
cb = find(isfinite(D.THD(:, kb)));
sb_cost = 100 * (D.THD(cb, kb) - D.env(cb));
X.single = struct('branch', kb, 'Ma_lo', grid(cb(1)), 'Ma_hi', grid(cb(end)), 'coverage_pct', 100 * numel(cb) / G, ...
    'cost_mean_pp', mean(sb_cost), 'cost_max_pp', max(sb_cost), 'candidates', numel(cand));

SW = E2b.sweep.TH;  sf = find(all(isfinite(SW), 2));
CT = E2.cont.TH;    cf = find(all(isfinite(CT), 2));

hdr = {'arm', 'level_R', 'stored_sets', 'stored_numbers', 'test_points', 'coverage_pct', 'excess_mean_pp', ...
       'excess_max_pp', 'pct_over_0p1pp', 'fund_err_mean_pct', 'fund_err_max_pct', 'excess_mean_sd_over_seeds_pp', ...
       'note'};
rows = cell(0, numel(hdr));
nn = struct('R', {}, 'Req', {}, 'widths', {}, 'mean', {}, 'max', {}, 'mean_corr', {}, 'max_corr', {}, ...
            'fails_corr', {}, 'nparam', {});
for R = O.R
    [Yb, nset, okb] = branch_aware(D, R, Mt, labt, 'pchip');
    [Yl, ~, okl] = branch_aware(D, R, Mt, labt, 'linear');
    numbers = N * nset + meta;  Req = ceil(numbers / N);
    rows(end+1, :) = arm('BA pchip', R, nset, numbers, sc(F, P, Yb(okb, :), Mt(okb), Et(okb)), numel(te), 0, ''); %#ok<AGROW>
    rows(end+1, :) = arm('BA linear', R, nset, numbers, sc(F, P, Yl(okl, :), Mt(okl), Et(okl)), numel(te), 0, ''); %#ok<AGROW>
    [Yc, Ic] = F.solve(Yb(okb, :), Mt(okb), P);  Yc(~Ic.ok, :) = Yb(okb(~Ic.ok), :);
    rows(end+1, :) = arm('BA + corrector', R, nset, numbers, sc(F, P, Yc, Mt(okb), Et(okb)), numel(te), 0, ...
        sprintf('%d corrector solves did not converge (BA angles kept)', sum(~Ic.ok))); %#ok<AGROW>
    ri = unique(fg(round(linspace(1, numel(fg), Req))));
    Ya = interp1(grid(ri), D.envTH(ri, :), Mt, 'linear');
    rows(end+1, :) = arm('plain', R, numel(ri), N * numel(ri), sc(F, P, Ya, Mt, Et), numel(te), 0, ''); %#ok<AGROW>
    [~, jn] = min(abs(Mt - grid(ri).'), [], 2);
    rows(end+1, :) = arm('nearest row', R, numel(ri), N * numel(ri), sc(F, P, D.envTH(ri(jn), :), Mt, Et), ...
        numel(te), 0, 'fundamental quantised to the row spacing'); %#ok<AGROW>
    [Yac, Iac] = F.solve(Ya, Mt, P);  Yac(~Iac.ok, :) = Ya(~Iac.ok, :);
    rows(end+1, :) = arm('plain + corrector', R, numel(ri), N * numel(ri), sc(F, P, Yac, Mt, Et), numel(te), 0, ...
        sprintf('%d corrector solves did not converge', sum(~Iac.ok))); %#ok<AGROW>
    rs = unique(sf(round(linspace(1, numel(sf), Req))));
    rows(end+1, :) = covered('sweep [8]', R, rs, SW, grid, Mt, Et, F, P, numel(te), N); %#ok<AGROW>
    rc = unique(cf(round(linspace(1, numel(cf), Req))));
    rows(end+1, :) = covered('continuous (<= 5 deg)', R, rc, CT, grid, Mt, Et, F, P, numel(te), N); %#ok<AGROW>
    rb = unique(cb(round(linspace(1, numel(cb), Req))));
    rows(end+1, :) = covered(sprintf('single branch %d [22]', kb), R, rb, D.TH(:, :, kb), grid, Mt, Et, F, P, ...
        numel(te), N); %#ok<AGROW>
    w = O.NN_WIDTHS(17 * O.NN_WIDTHS + 45 <= N * numel(ri));
    if isempty(w), w = min(O.NN_WIDTHS); end
    S = zeros(numel(O.SEEDS), 7);  Sc = S;  fails = 0;  np = 0;  wch = zeros(1, numel(O.SEEDS));
    for si = 1:numel(O.SEEDS)
        M = shm_tanhlm(grid(ri), D.envTH(ri, :), struct('WIDTHS', w, 'SEED', O.SEEDS(si), 'RESTARTS', O.NN_RESTARTS));
        Yn = M.predict(Mt);
        for t = 1:size(Yn, 1), Yn(t, :) = F.feasible(Yn(t, :), P); end
        S(si, :) = sc(F, P, Yn, Mt, Et);  np = max(np, M.nparam);  wch(si) = M.width;
        if O.NN_CORR
            [Ync, Inc] = F.solve(Yn, Mt, P);  Ync(~Inc.ok, :) = Yn(~Inc.ok, :);
            Sc(si, :) = sc(F, P, Ync, Mt, Et);  fails = fails + sum(~Inc.ok);
        end
    end
    rows(end+1, :) = arm('learned generator', R, numel(ri), np, mean(S, 1), numel(te), std(S(:, 2)), ...
        sprintf('tanh MLP, widths chosen %s of %s, %d seeds; max over seeds %.3f pp', mat2str(wch), mat2str(w), ...
        numel(O.SEEDS), max(S(:, 3)))); %#ok<AGROW>
    rows{end, 8} = max(S(:, 3));
    if O.NN_CORR
        rows(end+1, :) = arm('learned generator + corrector', R, numel(ri), np, mean(Sc, 1), numel(te), std(Sc(:, 2)), ...
            sprintf('%d of %d corrector solves did not converge', fails, numel(O.SEEDS) * numel(te))); %#ok<AGROW>
        rows{end, 8} = max(Sc(:, 3));
    end
    nn(end+1) = struct('R', R, 'Req', numel(ri), 'widths', wch, 'mean', S(:, 2).', 'max', S(:, 3).', ...
        'mean_corr', Sc(:, 2).', 'max_corr', Sc(:, 3).', 'fails_corr', fails, 'nparam', np); %#ok<AGROW>
    if O.verbose
        fprintf('  storage level R = %3d: branch-aware %d sets + %d switch numbers = %d numbers; other arms %d rows\n', ...
            R, nset, meta, numbers, numel(ri));
        r = rows(end - (9 + O.NN_CORR):end, :);
        for i = 1:size(r, 1)
            fprintf('    %-30s numbers %5d  covered %5.1f%%  excess mean %8.4f  max %7.3f pp  fund max %.4f %%\n', ...
                r{i, 1}, r{i, 4}, r{i, 6}, r{i, 7}, r{i, 8}, r{i, 11});
        end
    end
end
X.name = 'STO_equal_storage';
X.header = hdr;  X.data = rows;  X.nn = nn;
X.meta = struct('owner_runs', runs, 'switch_numbers', meta);
X.opts = O;
end

function r = arm(name, R, nset, numbers, s, ntot, sd, note)
r = {name, R, nset, numbers, s(1), 100 * s(1) / ntot, s(2), s(3), s(4), s(5), s(6), sd, note};
end

function r = covered(name, R, rows_idx, TH, grid, Mt, Et, F, P, ntot, N)
ok = Mt >= grid(rows_idx(1)) & Mt <= grid(rows_idx(end));
Y = interp1(grid(rows_idx), TH(rows_idx, :), Mt(ok), 'linear');
r = arm(name, R, numel(rows_idx), N * numel(rows_idx), sc(F, P, Y, Mt(ok), Et(ok)), ntot, 0, '');
end

function s = sc(F, P, Y, Ma, env)
Y = sort(Y, 2);
R = F.report(Y, Ma, P);
ex = 100 * (R.THD49(:) - env(:));
fe = 100 * R.b1err(:) ./ Ma(:);
if isempty(ex), s = [0 NaN NaN NaN NaN NaN NaN]; return; end
[~, iw] = max(ex);
s = [numel(ex), mean(ex), max(ex), 100 * mean(ex > 0.1), mean(fe), max(fe), Ma(iw)];
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
