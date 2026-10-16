function X = shm_x_branchswitch(D, opts)
%SHM_X_BRANCHSWITCH Stage BSW: branch switching under R-L loads, exact (Table 3).
%   X = SHM_X_BRANCHSWITCH(D, opts)
%     static   at every owner change of the stored subset and of the full set,
%              the deviation of the phase current from the new steady state
%              after a switch at the start of a period, for each power factor
%              PF and offset DELTA (SHM_RL_CURRENT), with L/R and the deviation
%              left one period later
%     dynamic  Ma ramps over the range in RAMP_PERIODS periods and Ma dithering
%              around a switch point, with hysteresis H, for the branch-aware
%              and the plain table (R rows), at power factor PF_DYN
%
%   Options: PF ([0.995 0.95 0.85 0.75 0.6 0.5]), DELTA ([0 0.02 0.05]),
%            NTS (3600), R (100), RAMP_PERIODS ([250 50 10]), DITHER_AMP (0.002),
%            DITHER_HZ (1), DITHER_PERIODS (100), H ([0 0.0025 0.005]),
%            PF_DYN (0.75), verbose

if nargin < 2, opts = struct(); end
O.PF = [0.995 0.95 0.85 0.75 0.6 0.5];
O.DELTA = [0 0.02 0.05];
O.NTS = 3600;
O.R = 100;
O.RAMP_PERIODS = [250 50 10];
O.DITHER_AMP = 0.002;  O.DITHER_HZ = 1;  O.DITHER_PERIODS = 100;
O.H = [0 0.0025 0.005];
O.PF_DYN = 0.75;
O.verbose = true;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

F = shm_model();  P = D.P;  s = P.s(:).';
grid = D.grid(:);
xs = (0:O.NTS-1).' * (2*pi / O.NTS);

sets = {'stored', D.owns(:); 'full', D.own(:)};
X.static = struct('set', {}, 'Ma', {}, 'from', {}, 'to', {}, 'jump_deg', {}, 'wave_diff_pct', {}, ...
    'PF', {}, 'delta', {}, 'dev_period_start_pct', {}, 'dev_worst_pct', {}, 'dev_best_pct', {});
for si = 1:size(sets, 1)
    own = sets{si, 2};
    ch = find(own(1:end-1) ~= own(2:end) & isfinite(own(1:end-1)) & isfinite(own(2:end)));
    for g = ch(:).'
        a = own(g);  b = own(g + 1);
        [xc, thA, thB] = at_switch(F, P, D, g, a, b);
        if isempty(thA), continue; end
        jd = max(abs(thA - thB)) * 180/pi;  wd = shm_wavediff(thA, thB, P.s);
        for pf = O.PF
            for dl = O.DELTA
                iA = shm_rl_current(thA, s, pf, dl, xs);  iB = shm_rl_current(thB, s, pf, dl, xs);
                I1 = F.b1(thB, P);
                dv = 100 * max(abs(iA - iB), [], 2) / I1;
                X.static(end+1) = struct('set', sets{si, 1}, 'Ma', xc, 'from', a, 'to', b, 'jump_deg', jd, ...
                    'wave_diff_pct', wd, 'PF', pf, 'delta', dl, 'dev_period_start_pct', dv(1), ...
                    'dev_worst_pct', max(dv), 'dev_best_pct', min(dv)); %#ok<AGROW>
            end
        end
    end
end
st = X.static;
X.static_summary = struct('set', {}, 'PF', {}, 'delta', {}, 'changes', {}, 'start_max_pct', {}, ...
    'start_median_pct', {}, 'worst_max_pct', {}, 'best_max_pct', {}, 'tau_ms', {}, 'after_one_period_pct', {});
for si = 1:size(sets, 1)
    for pf = O.PF
        for dl = O.DELTA
            q = strcmp({st.set}, sets{si, 1}) & [st.PF] == pf & [st.delta] == dl;
            if ~any(q), continue; end
            X.static_summary(end+1) = struct('set', sets{si, 1}, 'PF', pf, 'delta', dl, 'changes', sum(q), ...
                'start_max_pct', max([st(q).dev_period_start_pct]), ...
                'start_median_pct', median([st(q).dev_period_start_pct]), ...
                'worst_max_pct', max([st(q).dev_worst_pct]), 'best_max_pct', max([st(q).dev_best_pct]), ...
                'tau_ms', 1e3 * tan(acos(pf)) / (2*pi*50), ...
                'after_one_period_pct', max([st(q).dev_period_start_pct]) * exp(-2*pi / tan(acos(pf)))); %#ok<AGROW>
        end
    end
end
if O.verbose
    for r = X.static_summary
        if r.delta ~= O.DELTA(1) && r.delta ~= O.DELTA(end), continue; end
        fprintf(['  %-6s PF %.3f delta %.2f: %2d owner changes; at the period start max %.2f%% (median %.2f%%) of I1, ' ...
                 'worst instant %.2f%%, best %.2f%%; L/R %.2f ms, left after one period %.3f%%\n'], ...
            r.set, r.PF, r.delta, r.changes, r.start_max_pct, r.start_median_pct, r.worst_max_pct, r.best_max_pct, ...
            r.tau_ms, r.after_one_period_pct);
    end
end

T = build_tables(D, O.R);
pf = O.PF_DYN;
X.dynamic = struct('scenario', {}, 'table', {}, 'H', {}, 'periods', {}, 'branch_switches', {}, ...
    'dev_max_pct', {}, 'thd_max_pct', {}, 'thd_mean_pct', {});
lo = min(T.runs(:, 2));  hi = max(T.runs(:, 3));
for np = O.RAMP_PERIODS
    Mtraj = linspace(lo, hi, np).';
    X.dynamic = [X.dynamic, run_traj(F, P, T, Mtraj, 'ramp', np, pf, O.H(min(2, numel(O.H))))];
end
[~, q] = max([st(strcmp({st.set}, 'stored') & [st.PF] == O.PF(1) & [st.delta] == 0).wave_diff_pct]);
sw = st(strcmp({st.set}, 'stored') & [st.PF] == O.PF(1) & [st.delta] == 0);
if ~isempty(sw)
    m0 = sw(q).Ma;
    k = (0:O.DITHER_PERIODS-1).';
    Mtraj = m0 + O.DITHER_AMP * sin(2*pi * O.DITHER_HZ * k / 50);
    for h = O.H
        X.dynamic = [X.dynamic, run_traj(F, P, T, Mtraj, sprintf('dither at Ma %.4f', m0), numel(k), pf, h)];
    end
end
if O.verbose
    for r = X.dynamic
        fprintf('  %-22s %-12s H %.4f: %3d periods, %3d branch switches, deviation max %.2f%% of I1, THD49 max %.2f%% mean %.3f%%\n', ...
            r.scenario, r.table, r.H, r.periods, r.branch_switches, r.dev_max_pct, r.thd_max_pct, r.thd_mean_pct);
    end
end
X.tables = struct('R', O.R, 'runs', T.runs, 'switch_points', T.runs(2:end, 2));
X.name = 'BSW_branch_switch';
X.opts = O;
end

function [xc, thA, thB] = at_switch(F, P, D, g, a, b)
grid = D.grid(:);  thA = [];  thB = [];
da = D.THD(g, a) - D.THD(g, b);  db = D.THD(g + 1, a) - D.THD(g + 1, b);
if ~isfinite(D.THD(g + 1, a)), xc = grid(g);
elseif ~isfinite(D.THD(g, b)), xc = grid(g + 1);
elseif all(isfinite([da db])) && da ~= db, xc = grid(g) + da / (da - db) * (grid(g + 1) - grid(g));
else, xc = grid(g + 1); end
thA = one(F, P, D, a, g, xc);  thB = one(F, P, D, b, g, xc);
if isempty(thA) || isempty(thB), thA = []; thB = []; end
end

function th = one(F, P, D, k, g, x)
grid = D.grid(:);  th = [];
f = [g g+1];  f = f(isfinite(D.THD(f, k)));
if isempty(f), return; end
if numel(f) == 2
    x0 = interp1(grid(f), reshape(D.TH(f, :, k), 2, P.N), x, 'linear');
else
    x0 = D.TH(f, :, k);
end
[y, I] = F.solve(x0, x, P);
if I.ok && max(abs(y - sort(x0))) <= pi/180, th = y; else, th = sort(x0); end
end

function T = build_tables(D, R)
grid = D.grid(:);  N = D.P.N;
sel = D.sel(:).';  own = D.owns(:);
len = arrayfun(@(k) sum(own == k), sel);
use = sel(len > 0);  len = len(len > 0);
nb = max(2, round(R * len / sum(len)));
T.tab = cell(1, numel(use));  T.use = use;
for u = 1:numel(use)
    k = use(u);  og = find(own == k);
    lo = max(1, og(1) - 1);  hi = min(numel(grid), og(end) + 1);
    ex = find(isfinite(D.THD(:, k)));  lo = max(lo, ex(1));  hi = min(hi, ex(end));
    idx = unique(round(linspace(lo, hi, nb(u))));  idx = idx(isfinite(D.THD(idx, k)));
    T.tab{u} = struct('Ma', grid(idx), 'TH', reshape(D.TH(idx, :, k), numel(idx), N));
end
f = find(isfinite(own));  runs = [own(f(1)), grid(f(1)), NaN];
for g = f(1):f(end) - 1
    if isfinite(own(g)) && isfinite(own(g + 1)) && own(g) ~= own(g + 1)
        a = own(g);  b = own(g + 1);
        da = D.THD(g, a) - D.THD(g, b);  db = D.THD(g + 1, a) - D.THD(g + 1, b);
        if ~isfinite(D.THD(g + 1, a)), xc = grid(g);
        elseif ~isfinite(D.THD(g, b)), xc = grid(g + 1);
        elseif all(isfinite([da db])) && da ~= db, xc = grid(g) + da / (da - db) * (grid(g + 1) - grid(g));
        else, xc = grid(g + 1); end
        runs(end, 3) = xc;  runs(end+1, :) = [b, xc, NaN]; %#ok<AGROW>
    end
end
runs(end, 3) = grid(f(end));
T.runs = runs;
fg = find(isfinite(D.env));
ri = unique(fg(round(linspace(1, numel(fg), R))));
T.plain = struct('Ma', grid(ri), 'TH', D.envTH(ri, :));
end

function [th, run] = ba_read(T, m, run, H)
r = T.runs;
if isempty(run) || m < r(run, 2) - H || m > r(run, 3) + H
    run = find(m >= r(:, 2) & m <= r(:, 3), 1);
    if isempty(run), [~, run] = min(min(abs(m - r(:, 2)), abs(m - r(:, 3)))); end
end
u = find(T.use == r(run, 1), 1);  tb = T.tab{u};
if numel(tb.Ma) >= 3, th = interp1(tb.Ma, tb.TH, m, 'pchip', 'extrap');
elseif numel(tb.Ma) == 2, th = interp1(tb.Ma, tb.TH, m, 'linear', 'extrap');
else, th = tb.TH(1, :); end
end

function out = run_traj(F, P, T, Mtraj, name, np, pf, H)
s = P.s(:).';  tau = tan(acos(pf));  a = exp(-2*pi / tau); %#ok<NASGU>
out = struct('scenario', {}, 'table', {}, 'H', {}, 'periods', {}, 'branch_switches', {}, ...
    'dev_max_pct', {}, 'thd_max_pct', {}, 'thd_mean_pct', {});
for tb = {'branch-aware', 'plain'}
    run = [];  prev = [];  i = [];  nsw = 0;  dev = 0;  thd = zeros(numel(Mtraj), 1);
    for p = 1:numel(Mtraj)
        m = Mtraj(p);
        if strcmp(tb{1}, 'branch-aware')
            [th, run] = ba_read(T, m, run, H);
            if ~isempty(prev) && run ~= prev && T.runs(run, 1) ~= T.runs(prev, 1), nsw = nsw + 1; end
            prev = run;
        else
            th = interp1(T.plain.Ma, T.plain.TH, min(max(m, T.plain.Ma(1)), T.plain.Ma(end)), 'linear');
        end
        th = F.feasible(th, P);
        i0 = shm_rl_current(th, s, pf, 0, 0);
        I1 = F.b1(th, P);
        if isempty(i), i = i0; end
        dev = max(dev, 100 * max(abs(i - i0)) / I1);
        i = i0 + (i - i0) * a;
        thd(p) = 100 * F.thd(th, P);
    end
    hh = H;  if strcmp(tb{1}, 'plain'), hh = NaN; end
    out(end+1) = struct('scenario', sprintf('%s (%d periods)', name, np), 'table', tb{1}, 'H', hh, 'periods', np, ...
        'branch_switches', nsw, 'dev_max_pct', dev, 'thd_max_pct', max(thd), 'thd_mean_pct', mean(thd)); %#ok<AGROW>
end
end
