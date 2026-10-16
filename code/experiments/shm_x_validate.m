function E = shm_x_validate(D, opts)
%SHM_X_VALIDATE Stage E6: checks independent of how the branch set was built.
%   E = SHM_X_VALIDATE(D, opts)
%     1. a fresh multistart (NSTART starts) at NV modulation indices between
%        the grid points, compared with the stored branches there
%     2. the spectrum above the design band (THD to the 99th and 999th order,
%        weighted THD, largest harmonic above the 49th)
%     3. per-order limits on the grid, for 0.2 + 32.5/n % and 0.2 + 12.5/n %
%        above the 25th order
%     4. sensitivity of THD49 to a tighter upper angle bound (BOUND)
%
%   Options: NV (61), NSTART (200), SPEC_MA, BOUND (true), BOUND_MA, BOUND_UB,
%            BOUND_NSTART (150), BOUND_KEEP (10), CORR_DEG (1.0), solve, verbose

if nargin < 2, opts = struct(); end
O.NV       = 61;
O.NSTART   = 200;
O.SPEC_MA  = [0.70 0.80 0.90 1.00 1.10 1.20];
O.BOUND    = true;
O.BOUND_MA = [0.80 0.90 1.00 1.10 1.20];
O.BOUND_UB = [];
O.BOUND_NSTART = 150;
O.BOUND_KEEP   = 10;
O.CORR_DEG     = 1.0;
O.solve    = struct();
O.verbose  = true;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

F = shm_model();  P = D.P;  grid = D.grid;

lo = grid(1);  hi = grid(end);  h = median(diff(grid));
pts = linspace(lo + h/2, hi - h/2, O.NV).';
ref = nan(O.NV, 1);  dep = ref;  depi = ref;  owner = ref;
src = zeros(O.NV, 1);
brej = zeros(O.NV, 1);
cdeg = O.CORR_DEG * pi/180;
if O.verbose
    fprintf('  dense envelope: %d off-grid points, %d fresh starts each\n', O.NV, O.NSTART);
end
for i = 1:O.NV
    m = pts(i);
    [~, I] = F.solve(F.starts(O.NSTART, P), m, P, O.solve);
    if any(I.ok), ref(i) = min(I.thd(I.ok)); end
    bi = inf;  bv = inf;  kv = NaN;  sv = 0;  ri = 0;
    for k = D.sel(:).'
        th = branch_at(D, k, m, h, F, P);
        if isempty(th), continue; end
        R = F.report(th, m, P);
        [xc, Ic] = F.solve(th, m, P, O.solve);
        if ~Ic.ok
            [xc, Ic] = F.solve(F.feasible(th + 1e-4 * cos(2.4 * (1:numel(th))), P), m, P, O.solve);
        end
        on = Ic.ok && max(abs(xc - th)) <= cdeg;
        if R.b1err < 1e-3 && R.THD49 < bi
            bi = R.THD49;  ri = double(~Ic.ok) + 2 * double(Ic.ok && ~on);
        end
        if on,                  v = Ic.thd;   s = 1;
        elseif R.b1err < 1e-3,  v = R.THD49;  s = 2;
        else,                   v = inf;      s = 0;
        end
        if v < bv, bv = v;  kv = k;  sv = s; end
    end
    if isfinite(bi), depi(i) = bi;  brej(i) = ri; end
    if isfinite(bv), dep(i) = bv;  owner(i) = kv;  src(i) = sv; end
    if O.verbose && mod(i, max(1, round(O.NV / 6))) == 0
        fprintf('    Ma %.4f   multistart %.4f%%   stored branches %.4f%% (interpolated %.4f%%)\n', ...
            m, 100*ref(i), 100*dep(i), 100*depi(i));
    end
end
gap = 100 * (dep - ref);   g = gap(isfinite(gap));
gapi = 100 * (depi - ref); gi_ = gapi(isfinite(gapi));
E.dense = struct('Ma', pts, 'ref', ref, 'dep', dep, 'gap_pp', gap, 'owner', owner, ...
    'dep_interp', depi, 'gap_interp_pp', gapi, ...
    'mean_gap_pp', mean(g), 'median_gap_pp', median(g), 'worst_gap_pp', max(g), ...
    'best_gap_pp', min(g), 'n_over_0p1', sum(g > 0.1), 'n_over_1', sum(g > 1), ...
    'n_better', sum(g < -1e-6), 'n_uncovered', sum(~isfinite(dep)), ...
    'worst_gap_interp_pp', max(gi_), 'mean_gap_interp_pp', mean(gi_), ...
    'source', src, 'best_corrector', brej, 'n_interp_used', sum(src == 2), ...
    'n_best_noconv', sum(brej == 1), 'n_best_left', sum(brej == 2));
if O.verbose
    fprintf('  dense gap: mean %.4f pp, median %.4f pp, worst %.4f pp; %d of %d over 0.1 pp; stored better at %d\n', ...
        mean(g), median(g), max(g), sum(g > 0.1), numel(g), sum(g < -1e-6));
    fprintf('  without the corrector (the table''s own interpolation): worst %.4f pp\n', max(gi_));
    fprintf(['  corrector on the best interpolated branch: did not converge at %d, left the branch at %d ' ...
             'of %d points; interpolated value used at %d\n'], sum(brej == 1), sum(brej == 2), O.NV, sum(src == 2));
end

sp = zeros(0, 9);
for m = O.SPEC_MA
    [dg, gi] = min(abs(grid - m));
    if dg > 1e-9 || ~isfinite(D.env(gi)), continue; end
    S = F.spectrum(D.envTH(gi, :), P, 999);
    sp(end+1, :) = [grid(gi), 100*S.THD49, 100*S.THD99, 100*S.THDfull, 100*S.WTHD49, ...
                    100*S.WTHDfull, 100*S.OOB_max, S.OOB_order, 100*S.b3_ratio]; %#ok<AGROW>
end
E.spectrum = sp;
E.spectrum_header = {'Ma','THD49_pct','THD99_pct','THDfull_pct','WTHD49_pct','WTHDfull_pct', ...
                     'largest_above_49_pct','its_order','b3_pct'};

fin = find(isfinite(D.env));
Rc = F.report(D.envTH(fin, :), grid(fin), P);
Pa = P;  o = P.ORD;  Pa.LIM(o > 25) = (0.2 + 12.5 ./ o(o > 25)) / 100;
Ra = F.report(D.envTH(fin, :), grid(fin), Pa);
E.compliance = struct('Ma', grid(fin), 'MAXQ', Rc.MAXQ, 'compliant', Rc.compliant, ...
    'frac_compliant', mean(Rc.compliant), 'MAXQ_median', median(Rc.MAXQ), 'MAXQ_worst', max(Rc.MAXQ), ...
    'MAXQ_alt', Ra.MAXQ, 'frac_compliant_alt', mean(Ra.compliant), 'MAXQ_worst_alt', max(Ra.MAXQ), ...
    'note', 'limits above the 25th: 0.2+32.5/n % (default) and 0.2+12.5/n % (alt)');
if O.verbose
    fprintf(['  per-order limits met at %.1f%% of the envelope (worst ratio %.2f); ' ...
             'with 0.2+12.5/n above the 25th: %.1f%% (worst %.2f)\n'], ...
        100 * mean(Rc.compliant), max(Rc.MAXQ), 100 * mean(Ra.compliant), max(Ra.MAXQ));
end

E.bound = [];
if O.BOUND
    ub = O.BOUND_UB(:).';
    if isempty(ub), ub = [89 87 85 80]; end
    ub = sort(unique([P.ub * 180/pi, ub(ub < P.ub * 180/pi - 1e-9)]), 'descend');
    nb = numel(ub);  nm = numel(O.BOUND_MA);
    B = nan(nm, nb);  TOP = B;  SRC = zeros(nm, nb);
    Mb = nan(nm, 1);
    for i = 1:nm
        [~, gi] = min(abs(grid - O.BOUND_MA(i)));
        m = grid(gi);  Mb(i) = m;
        cand = zeros(0, 3);
        t = D.THD(gi, :);  ok = find(isfinite(t));
        X = zeros(numel(ok), P.N);
        for r = 1:numel(ok), X(r, :) = D.TH(gi, :, ok(r)); end
        cand = [cand; t(ok).', max(X, [], 2) * 180/pi, ones(numel(ok), 1)]; %#ok<AGROW>
        [~, o] = sort(t(ok));
        Xc = X(o(1:min(O.BOUND_KEEP, numel(o))), :);
        for j = 1:nb
            if isempty(Xc), break; end
            if ub(j) * pi/180 >= P.ub - 1e-12, continue; end
            Pb = P;  Pb.ub = ub(j) * pi/180;
            X0 = Xc;
            for r = 1:size(X0, 1), X0(r, :) = F.feasible(X0(r, :), Pb); end
            [Xs, I] = F.solve(X0, m, Pb, O.solve);
            cand = [cand; I.thd(I.ok), max(Xs(I.ok, :), [], 2) * 180/pi, 2 * ones(sum(I.ok), 1)]; %#ok<AGROW>
            Xc(I.ok, :) = Xs(I.ok, :);
        end
        for j = 1:nb
            Pb = P;  Pb.ub = ub(j) * pi/180;
            [Xm, I] = F.solve(F.starts(O.BOUND_NSTART, Pb), m, Pb, O.solve);
            cand = [cand; I.thd(I.ok), max(Xm(I.ok, :), [], 2) * 180/pi, 3 * ones(sum(I.ok), 1)]; %#ok<AGROW>
        end
        for j = 1:nb
            f = find(cand(:, 2) <= ub(j) + 1e-6);
            if isempty(f), continue; end
            [B(i, j), q] = min(cand(f, 1));
            TOP(i, j) = cand(f(q), 2);  SRC(i, j) = cand(f(q), 3);
        end
        if O.verbose
            fprintf('  bound study Ma %.2f: best THD49 %s %% for ub = %s deg\n', m, ...
                mat2str(round(100*B(i,:)*1e3)/1e3), mat2str(ub));
        end
    end
    E.bound = struct('Ma', Mb, 'ub_deg', ub, 'thd_pct', 100*B, 'top_angle_deg', TOP, ...
        'source', SRC, 'cost_pp', 100*(B - repmat(B(:, 1), 1, nb)), ...
        'note', ['Best THD49 with the upper angle bound tightened. Candidates: the traced branches as ' ...
                 'they are (source 1), the best branches carried down in the bound by warm-started ' ...
                 'fmincon (source 2), a fresh multistart at every bound (source 3). The best at bound ' ...
                 'ub is the minimum over every candidate whose top angle is <= ub. cost = rise over ' ...
                 'the loosest bound, pp.']);
end
end

function th = branch_at(D, k, m, h, F, P)
th = [];
grid = D.grid;  ex = find(isfinite(D.THD(:, k)));
if isempty(ex), return; end
a = grid(ex(1));  b = grid(ex(end));
if m < a - h/2 - 1e-12 || m > b + h/2 + 1e-12, return; end
if numel(ex) == 1
    th = D.TH(ex, :, k);
elseif m < a || m > b
    if m < a, s = ex(1:2); else, s = ex(end-1:end); end
    th = interp1(grid(s), D.TH(s, :, k), m, 'linear', 'extrap');
elseif numel(ex) == 2
    th = interp1(grid(ex), D.TH(ex, :, k), m, 'linear');
else
    th = interp1(grid(ex), D.TH(ex, :, k), m, 'pchip');
end
th = F.feasible(th, P);
end
