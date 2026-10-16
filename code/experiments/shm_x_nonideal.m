function E = shm_x_nonideal(D, opts)
%SHM_X_NONIDEAL Stage E9: the tables under non-ideal switching (Table 4).
%   E = SHM_X_NONIDEAL(D, opts) evaluates, for every condition of
%   SHM_NONIDEAL_CONDITIONS, the per-Ma optimum, the plain table and the
%   branch-aware table with the switching-function model SHM_NONIDEAL_WAVE.
%   A table's excess is taken over the optimum under the same condition.
%
%   Options: R (100), STEP (1), TEST (midpoints of stage E2b), Cnom (470e-6),
%            Rload (64), Lload (0.18), Vdc (800), f1 (50), M (2^17), Tvj (25),
%            verbose

if nargin < 2, opts = struct(); end
O.R = 100;  O.STEP = 1;  O.Cnom = 470e-6;  O.Rload = 64;  O.Lload = 0.18;
O.Vdc = 800;  O.f1 = 50;  O.M = 2^17;  O.Tvj = 25;  O.verbose = true;
O.TEST = [];
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

grid = D.grid;  N = D.P.N;
dev = shm_device(O.Tvj);
K   = shm_nonideal_conditions();
cfg = struct('Vdc', O.Vdc, 'f1', O.f1, 'R', O.Rload, 'L', O.Lload, 'M', O.M);

fg = find(isfinite(D.env));
ri = unique(fg(round(linspace(1, numel(fg), O.R))));
if isempty(O.TEST)
    te = setdiff(fg, ri);
    te = te(1:O.STEP:end);
    Mt = grid(te);  Mt = Mt(:);
    Ya = sort(interp1(grid(ri), D.envTH(ri, :), Mt, 'linear'), 2);
    Yb = sort(branch_aware(D, O.R, te), 2);
    Yo = D.envTH(te, :);
else
    T = O.TEST;  j = 1:O.STEP:numel(T.Ma);
    Mt = T.Ma(j);  Mt = Mt(:);
    Ya = sort(T.Ya(j, :), 2);  Yb = sort(T.Yb(j, :), 2);  Yo = sort(T.Yopt(j, :), 2);
    O.R = T.R;
end
okb = all(isfinite(Yb), 2);
go = fg(1:O.STEP:end);

hdr = {'condition', 'tag', 'opt_points', 'test_points', 'opt_dTHD_max_pp', 'opt_dTHD_mean_pp', ...
       'opt_b1_err_max_pct', 'lut_mean_pp', 'lut_max_pp', 'branch_mean_pp', 'branch_max_pp', 'np_pp_max_V'};
rows = cell(0, numel(hdr));
pts  = zeros(0, 6);
t0 = tic;
if O.verbose
    fprintf('  %d optimum points, %d test points of a %d-row table, %d conditions\n', ...
        numel(go), numel(Mt), O.R, numel(K));
end
for k = 1:numel(K)
    NI = struct('drops', K(k).drops, 'td', K(k).td, 'eps_dc', K(k).eps_dc, 'Ctot', 2*O.Cnom, 'tq', K(k).tq);
    to = nan(numel(go), 1);  dto = to;  bo = to;  np = to;
    for j = 1:numel(go)
        W = shm_nonideal_wave(D.envTH(go(j), :), grid(go(j)), NI, cfg, dev);
        to(j) = W.THD49;  dto(j) = W.dthd;  bo(j) = W.b1_err_pct;  np(j) = W.np_ripple_pp;
    end
    ta = nan(numel(Mt), 1);  tb = ta;  tm = ta;  dtm = ta;
    for j = 1:numel(Mt)
        W = shm_nonideal_wave(Yo(j, :), Mt(j), NI, cfg, dev);  tm(j) = W.THD49;  dtm(j) = W.dthd;
        W = shm_nonideal_wave(Ya(j, :), Mt(j), NI, cfg, dev);  ta(j) = W.THD49;
        if okb(j)
            W = shm_nonideal_wave(Yb(j, :), Mt(j), NI, cfg, dev);  tb(j) = W.THD49;
        end
    end
    ea = ta - tm;  eb = tb - tm;
    rows(end+1, :) = {K(k).name, K(k).tag, numel(go), numel(Mt), max(abs(dto)), mean(abs(dto)), ...
        max(abs(bo)), mean(ea), max(ea), mean(eb(okb)), max(eb(okb)), max(np)}; %#ok<AGROW>
    pts = [pts; [k * ones(numel(Mt), 1), Mt, tm, dtm, ta, tb]]; %#ok<AGROW>
    if O.verbose
        fprintf(['  %-26s optimum: |dTHD49| <= %.3f pp, |b1 err| <= %.2f %%;  table: %.3f pp mean, %.2f pp worst;  ' ...
                 'branch-aware: %.4f pp mean, %.3f pp worst  (%.1f min)\n'], strrep(K(k).name, '$\mu$', 'u'), ...
                 rows{end, 5}, rows{end, 7}, rows{end, 8}, rows{end, 9}, rows{end, 10}, rows{end, 11}, toc(t0)/60);
    end
end

E.name   = 'E9_nonideal';
E.header = hdr;
E.data   = rows;
E.points_header = {'condition', 'Ma', 'THD49_optimum_pct', 'dTHD49_optimum_pp', 'THD49_table_pct', ...
                   'THD49_branch_aware_pct'};
E.points = pts;
E.conditions = K;
E.device = dev.name;
E.opts = O;
E.note = sprintf(['Switching-function model (shm_nonideal_wave), device %s, load %g ohm + %g mH per phase, ' ...
    'C1 = C2 = %g uF, Vdc = %g V. Table: %d rows, linear interpolation; branch-aware: %d sets, pchip inside ' ...
    'the owning branch. Excess = THD49 of the table''s angles minus THD49 of the per-Ma optimum under the same ' ...
    'condition. opt_dTHD: the optimum''s THD49 minus its closed form.'], dev.name, O.Rload, 1e3 * O.Lload, ...
    1e6 * O.Cnom, O.Vdc, numel(ri), O.R);
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
