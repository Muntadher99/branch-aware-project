function E = shm_x_regress(D, opts)
%SHM_X_REGRESS Stage E3: regression inside a branch at a matched budget, d = 1.
%   E = SHM_X_REGRESS(D, opts) gives every method the same K training samples
%   of each of NB branches and scores it on the remaining samples by angle
%   error (deg) and THD49 excess (pp). Predicted angles that break the ordering
%   or the minimum gap are counted and repaired before the THD is evaluated.
%
%   Methods: linear, spline, pchip, RBF and GP (SHM_GPR), tanhLM (SHM_TANHLM),
%            MLP (SHM_MLP), ANFIS (SHM_ANFIS)
%   Options: NB (3), K ([6 9 15 25 41]), SEEDS (1:5), METHODS, mlp, anfis,
%            tanhlm, gpr, KEEP_SEEDS (false), verbose

if nargin < 2, opts = struct(); end
O.NB      = 3;
O.K       = [6 9 15 25 41];
O.SEEDS   = 1:5;
O.mlp     = struct();
O.anfis   = struct();
O.tanhlm  = struct();
O.gpr     = struct();
O.METHODS = {'linear', 'spline', 'pchip', 'MLP', 'ANFIS'};
O.verbose = true;
O.KEEP_SEEDS = false;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

F = shm_model();  P = D.P;  N = P.N;  grid = D.grid;

sel = D.sel(:).';
len = arrayfun(@(k) sum(D.owns == k), sel);
[len, o] = sort(len, 'descend');  sel = sel(o);
use = [];  hull = {};
for u = 1:numel(sel)
    k = sel(u);
    [i1, i2] = longest_run(isfinite(D.THD(:, k)));
    if isempty(i1), continue; end
    idx = (i1:i2).';
    if numel(idx) < min(O.K) + 5, continue; end
    use(end+1) = k; %#ok<AGROW>
    hull{end+1} = idx; %#ok<AGROW>
    if numel(use) >= O.NB, break; end
end

hdr = {'method','branch','K_train','n_test','angle_MAE_deg','angle_MAE_sd_deg','angle_max_err_deg', ...
       'THD_excess_mean_pp','THD_excess_max_pp','pct_infeasible','stored_parameters','seed_CV_pct'};
rows = {};
curve = [];
seedrows = cell(0, 10);
meth = O.METHODS;

for u = 1:numel(use)
    k = use(u);  idx = hull{u};  n = numel(idx);
    Ma = grid(idx);  Y = D.TH(idx, :, k);  tb = D.THD(idx, k);
    if O.verbose
        fprintf('  branch %d: Ma %.4f..%.4f, %d samples\n', k, Ma(1), Ma(end), n);
    end
    for K = O.K
        if K > n - 5, continue; end
        tr = unique(round(linspace(1, n, K)));  te = setdiff(1:n, tr);
        for m = 1:numel(meth)
            switch meth{m}
                case {'linear', 'spline', 'pchip'}
                    Yp = interp1(Ma(tr), Y(tr, :), Ma(te), meth{m});
                    r = evaluate(F, P, Yp, Y(te, :), Ma(te), tb(te));
                    if strcmp(meth{m}, 'linear'), np = K * N; else, np = 4 * (K - 1) * N; end
                    rows(end+1, :) = {meth{m}, k, K, numel(te), r(1), 0, r(2), r(3), r(4), r(5), np, 0}; %#ok<AGROW>
                    curve(end+1, :) = [K, m, r(1)]; %#ok<AGROW>
                case {'MLP', 'tanhLM'}
                    rs = zeros(numel(O.SEEDS), 5);  np = 0;
                    for si = 1:numel(O.SEEDS)
                        if strcmp(meth{m}, 'MLP')
                            mo = O.mlp;  mo.SEED = O.SEEDS(si);
                            M = shm_mlp(Ma(tr), Y(tr, :), mo);
                        else
                            mo = O.tanhlm;  mo.SEED = O.SEEDS(si);
                            M = shm_tanhlm(Ma(tr), Y(tr, :), mo);
                        end
                        Yp = M.predict(Ma(te));
                        rs(si, :) = evaluate(F, P, Yp, Y(te, :), Ma(te), tb(te));
                        if isfield(M, 'nparam'), np = M.nparam; end
                        if O.KEEP_SEEDS
                            wd = NaN;  if isfield(M, 'width'), wd = M.width; end
                            seedrows(end+1, :) = {meth{m}, k, K, O.SEEDS(si), rs(si, 1), rs(si, 2), rs(si, 3), ...
                                rs(si, 4), rs(si, 5), wd}; %#ok<AGROW>
                        end
                    end
                    cv = 100 * std(rs(:, 1)) / max(mean(rs(:, 1)), eps);
                    rows(end+1, :) = {meth{m}, k, K, numel(te), mean(rs(:,1)), std(rs(:,1)), max(rs(:,2)), ...
                        mean(rs(:,3)), max(rs(:,4)), mean(rs(:,5)), np, cv}; %#ok<AGROW>
                    curve(end+1, :) = [K, m, mean(rs(:, 1))]; %#ok<AGROW>
                case {'GP', 'RBF'}
                    go = O.gpr;  if strcmp(meth{m}, 'RBF'), go.KERNEL = 'phs'; end
                    M = shm_gpr(Ma(tr), Y(tr, :), go);
                    Yp = M.predict(Ma(te));
                    r = evaluate(F, P, Yp, Y(te, :), Ma(te), tb(te));
                    rows(end+1, :) = {meth{m}, k, K, numel(te), r(1), 0, r(2), r(3), r(4), r(5), M.nparam, 0}; %#ok<AGROW>
                    curve(end+1, :) = [K, m, r(1)]; %#ok<AGROW>
                case 'ANFIS'
                    M = shm_anfis(Ma(tr), Y(tr, :), O.anfis);
                    Yp = M.predict(Ma(te));
                    r = evaluate(F, P, Yp, Y(te, :), Ma(te), tb(te));
                    rows(end+1, :) = {'ANFIS', k, K, numel(te), r(1), 0, r(2), r(3), r(4), r(5), M.nparam, 0}; %#ok<AGROW>
                    curve(end+1, :) = [K, m, r(1)]; %#ok<AGROW>
            end
            if O.verbose
                r0 = rows(end, :);
                fprintf('    K=%3d %-7s MAE %8.4f deg  max %7.3f deg  THD excess %8.4f pp  params %6d\n', ...
                    K, r0{1}, r0{5}, r0{7}, r0{8}, r0{11});
            end
        end
    end
end

if isempty(curve), curve = zeros(0, 3); end
Ks = unique(curve(:, 1)).';
avg = nan(numel(Ks), numel(meth));
for i = 1:numel(Ks)
    for m = 1:numel(meth)
        s = curve(:, 1) == Ks(i) & curve(:, 2) == m;
        if any(s), avg(i, m) = mean(curve(s, 3)); end
    end
end

E.name    = 'E3a_regression_1D';
E.header  = hdr;
E.data    = rows;
E.curveK  = Ks(:);
E.curve   = avg;
E.methods = meth;
E.branches = use;
if O.KEEP_SEEDS
    E.seed_rows = seedrows;
    E.seed_header = {'method','branch','K_train','seed','angle_MAE_deg','angle_max_err_deg', ...
        'THD_excess_mean_pp','THD_excess_max_pp','pct_infeasible','hidden_width'};
end
E.note = sprintf(['Matched-budget regression inside a branch: every method sees the same K samples ' ...
    'of the same branch and is scored on the same held-out samples of its traced range. %d branches, ' ...
    'the most-owning first. MLP figures are the mean over %d seeds; seed_CV_pct is the seed-to-seed coefficient ' ...
    'of variation of the MAE. stored_parameters is what a deployed implementation would hold.'], ...
    numel(use), numel(O.SEEDS));
end

function [i1, i2] = longest_run(tf)
i1 = [];  i2 = [];  best = 0;  s = 0;
tf = tf(:).';
for k = 1:numel(tf)
    if tf(k)
        if s == 0, s = k; end
        if k - s + 1 > best, best = k - s + 1;  i1 = s;  i2 = k; end
    else
        s = 0;
    end
end
end

function r = evaluate(F, P, Yp, Yt, Ma, tbr)
ae  = abs(Yp - Yt) * 180/pi;
inf_ = any(diff(Yp, 1, 2) < P.minGap - 1e-12, 2) | Yp(:, 1) < P.lb | Yp(:, end) > P.ub;
Yf = Yp;
for i = 1:size(Yf, 1), Yf(i, :) = F.feasible(Yf(i, :), P); end
R  = F.report(Yf, Ma, P);
ex = 100 * (R.THD49 - tbr(:));
r  = [mean(ae(:)), max(ae(:)), mean(ex), max(ex), 100 * mean(inf_)];
end
