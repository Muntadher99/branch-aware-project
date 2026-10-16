function X = shm_x_compliance(D, E2b, opts)
%SHM_X_COMPLIANCE Stage CMP: per-order compliance of the optimum and the tables.
%   X = SHM_X_COMPLIANCE(D, E2b, opts) evaluates q_n = |b_n| / (L_n * b_1) for
%   the limits L_n of SHM_MODEL (EN 50160 for orders 5-25, 0.2 + 32.5/n %
%   above) on the grid optimum and, at the 200 midpoints, for the optimum, the
%   plain table, the branch-aware table and both with one corrector solve.
%   Reports the share compliant, q_max, and where and on which order a pattern
%   fails; also THD40.

if nargin < 3, opts = struct(); end
O.verbose = true;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

F = shm_model();  P = D.P;  P40 = F.defaults(P.N, 40);
grid = D.grid(:);
ords = P.ORD(:).';

f = find(isfinite(D.env));
[Cg, Qg, Jg] = qmax(F, P, D.envTH(f, :), grid(f));
X.grid = summarise(Cg, Qg, Jg, grid(f), ords);
X.grid.thd40_max_pct = 100 * max(F.thd(sort(D.envTH(f, :), 2), P40));

T = E2b.test;  Ma = T.Ma(:);
Yc_b = T.Yb;  [Yb2, Ib] = F.solve(T.Yb, Ma, P);  Yc_b(Ib.ok, :) = Yb2(Ib.ok, :);
Yc_a = T.Ya;  [Ya2, Ia] = F.solve(T.Ya, Ma, P);  Yc_a(Ia.ok, :) = Ya2(Ia.ok, :);
arms = {'optimum', T.Yopt; 'plain', T.Ya; 'plain+corrector', Yc_a; 'branch-aware', T.Yb; ...
        'branch-aware+corrector', Yc_b};
[Co, ~] = qmax(F, P, T.Yopt, Ma);
X.midpoints = struct('arm', {}, 'R', {}, 'n', {}, 'compliant_pct', {}, 'qmax_median', {}, 'qmax_max', {}, ...
    'fail_orders', {}, 'fail_Ma', {}, 'lost_vs_optimum', {}, 'gained_vs_optimum', {}, 'thd40_max_pct', {});
for a = 1:size(arms, 1)
    [C, Q, J] = qmax(F, P, arms{a, 2}, Ma);
    s = summarise(C, Q, J, Ma, ords);
    X.midpoints(end+1) = struct('arm', arms{a, 1}, 'R', T.R, 'n', numel(Ma), 'compliant_pct', s.compliant_pct, ...
        'qmax_median', s.qmax_median, 'qmax_max', s.qmax_max, 'fail_orders', s.fail_orders, 'fail_Ma', s.fail_Ma, ...
        'lost_vs_optimum', sum(Co & ~C), 'gained_vs_optimum', sum(~Co & C), ...
        'thd40_max_pct', 100 * max(F.thd(sort(arms{a, 2}, 2), P40)));
end
X.corrector_failed = struct('plain', sum(~Ia.ok), 'branch_aware', sum(~Ib.ok));
X.limits = struct('orders', ords, 'L_pct', 100 * P.LIM(:).');
X.name = 'CMP_grid_code';
if O.verbose
    g = X.grid;
    fprintf('  optimum on the grid: compliant at %.1f%% (fails at Ma %s on order(s) %s), q_max median %.3f, max %.3f\n', ...
        g.compliant_pct, rangestr(g.fail_Ma), mat2str(g.fail_orders), g.qmax_median, g.qmax_max);
    for a = 1:numel(X.midpoints)
        m = X.midpoints(a);
        fprintf('  %-24s compliant %5.1f%%  q_max median %.3f max %7.2f  lost vs optimum %3d  orders %s\n', ...
            m.arm, m.compliant_pct, m.qmax_median, m.qmax_max, m.lost_vs_optimum, mat2str(m.fail_orders));
    end
end
end

function [C, Q, J] = qmax(F, P, TH, Ma)
TH = sort(TH, 2);
[Bo, ~, b1] = F.harmonics(TH, P);
q = abs(Bo) ./ (repmat(P.LIM(:).', size(Bo, 1), 1) .* repmat(abs(b1), 1, size(Bo, 2)));
[Q, J] = max(q, [], 2);
C = Q <= 1 + 1e-4;
end

function s = summarise(C, Q, J, Ma, ords)
s.n = numel(C);
s.compliant_pct = 100 * mean(C);
s.qmax_median = median(Q);  s.qmax_max = max(Q);
s.fail_Ma = Ma(~C);  s.fail_Ma = s.fail_Ma(:).';
s.fail_orders = unique(ords(J(~C)));
end

function s = rangestr(x)
if isempty(x), s = 'none'; else, s = sprintf('%.4f..%.4f (%d points)', min(x), max(x), numel(x)); end
end
