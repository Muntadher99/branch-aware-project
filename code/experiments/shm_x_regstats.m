function X = shm_x_regstats(D, S, V, cfg, Eold)
%SHM_X_REGSTATS Stage RGS: seed spread, model selection and ranking (Table 2).
%   X = SHM_X_REGSTATS(D, S, V, cfg, Eold) repeats stages E3, E3b and E3c with
%   the same settings and random streams, keeping every seed (the means are
%   checked against the cached stages Eold), repeats the small tanh MLP with a
%   wider width search (cfg.RGS.EXT), and summarises per dimension, budget and
%   branch whether a learned model beats the best interpolant.

if nargin < 5, Eold = struct(); end
X = struct();
rng(1000 + cfg.seed.E3);   o = cfg.E3;  o.KEEP_SEEDS = true;
X.d1 = shm_x_regress(D, o);
rng(1000 + cfg.seed.E3b);  o = cfg.E3b;  o.P = D.P;  o.KEEP_SEEDS = true;
X.d2 = shm_x_regress_nd(S, o);
rng(1000 + cfg.seed.E3c);  o = cfg.E3c;  o.P = D.P;  o.KEEP_SEEDS = true;
X.d3 = shm_x_regress_nd(V, o);
X.check = struct('d1', same_rows(X.d1, Eold, 'E3'), 'd2', same_rows(X.d2, Eold, 'E3b'), ...
                 'd3', same_rows(X.d3, Eold, 'E3c'));
fprintf('  identical to the cached Table 2 stages: d=1 %s, d=2 %s, d=3 %s\n', yn(X.check.d1), yn(X.check.d2), yn(X.check.d3));

ext = cfg.RGS.EXT;
cvr = def(ext, 'CV_RESTARTS', def(ext, 'RESTARTS', 2));
hd = sprintf('widths %s, %d restarts for the final fit, %d per cross-validation fit', ...
    mat2str(ext.WIDTHS), def(ext, 'RESTARTS', 2), cvr);
t2 = tic;
fprintf('  wider search for the small MLP (%s), d = 1 ...\n', hd);
o = cfg.E3;  o.METHODS = {'tanhLM'};  o.tanhlm = ext;  o.KEEP_SEEDS = true;  o.verbose = true;
X.ext1 = shm_x_regress(D, o);
fprintf('  wider search, d = 2 ... (%.0f min so far)\n', toc(t2) / 60);
o = cfg.E3b;  o.P = D.P;  o.METHODS = {'tanhLM'};  o.tanhlm = ext;  o.KEEP_SEEDS = true;  o.verbose = true;
X.ext2 = shm_x_regress_nd(S, o);
fprintf('  wider search, d = 3 ... (%.0f min so far)\n', toc(t2) / 60);
o = cfg.E3c;  o.P = D.P;  o.METHODS = {'tanhLM'};  o.tanhlm = ext;  o.KEEP_SEEDS = true;  o.verbose = true;
X.ext3 = shm_x_regress_nd(V, o);
fprintf('  wider search done in %.0f min\n', toc(t2) / 60);

X.summary = [summ(X.d1.seed_rows, X.d1.data, 1, X.ext1), summ(X.d2.seed_rows, X.d2.data, 2, X.ext2), ...
             summ(X.d3.seed_rows, X.d3.data, 3, X.ext3)];
X.cases = [cases(X.d1.data, X.d1.seed_rows, 1), cases(X.d2.data, X.d2.seed_rows, 2), cases(X.d3.data, X.d3.seed_rows, 3)];
X.selection = struct( ...
    'tanhLM', sprintf(['one hidden layer of tanh units; width from %s by cross-validation on the training ' ...
        'samples only (leave-one-out for n <= 12, else %d folds); best of %d random initialisations; LM ' ...
        'with weight decay %g; extended search: widths %s, best of %d initialisations for the final fit and %d ' ...
        'for each cross-validation fit'], mat2str(def(cfg.E3, 'WIDTHS', [2 4 8])), ...
        5, 2, 1e-4, mat2str(ext.WIDTHS), ext.RESTARTS, cvr), ...
    'MLP', 'two hidden layers of 128 sinusoidal units, Adam, early stopping on a validation split of the training samples', ...
    'GP', 'squared-exponential ARD kernel, hyperparameters by maximum marginal likelihood (ML-II), noiseless', ...
    'ANFIS', 'first-order Sugeno, Gaussian membership functions on a grid partition, hybrid training');
X.name = 'RGS_regression_statistics';
for i = 1:numel(X.cases)
    c = X.cases(i);
    fprintf(['  d=%d: learned model better than the best interpolant (mean over seeds) in %d of %d (branch, budget) ' ...
             'cases by angle error, %d by THD49 excess; best single seed: %d and %d\n'], c.d, c.learned_wins_mae, ...
             c.n, c.learned_wins_thd, c.bestseed_wins_mae, c.bestseed_wins_thd);
end
end

function v = def(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end

function s = yn(tf)
if tf, s = 'yes'; else, s = 'NO'; end
end

function tf = same_rows(E, Eold, name)
tf = false;
if ~isfield(Eold, name) || isempty(Eold.(name)), return; end
A = E.data;  B = Eold.(name).data;
if ~isequal(size(A), size(B)), return; end
if isempty(A), tf = true; return; end
num = cellfun(@isnumeric, A(1, :));
a = cell2mat(A(:, num));  b = cell2mat(B(:, num));
k = ~isnan(b);
tf = isequal(A(:, ~num), B(:, ~num)) && isequal(isnan(a), isnan(b)) && ...
     all(abs(a(k) - b(k)) <= 1e-9 * max(1, abs(b(k))));
end

function T = summ(seedrows, rows, d, ext)
if d == 1
    mcol = 1; bcol = 3; ecol = 5; tcol = 8;
    smcol = 1; sbcol = 3; secol = 5; stcol = 7; sbr = 4;
else
    mcol = 2; bcol = 4; ecol = 7; tcol = 10;
    smcol = 2; sbcol = 4; secol = 7; stcol = 9; sbr = 6;
end
T = struct('d', {}, 'method', {}, 'budget', {}, 'branches', {}, 'mae_mean', {}, 'mae_seed_sd', {}, ...
    'mae_best_seed', {}, 'mae_worst_seed', {}, 'thd_mean', {}, 'thd_seed_sd', {}, 'thd_best_seed', {}, ...
    'thd_worst_seed', {});
if isempty(rows), return; end
meth = unique(rows(:, mcol)).';  bud = unique(cell2mat(rows(:, bcol))).';
T = struct('d', {}, 'method', {}, 'budget', {}, 'branches', {}, 'mae_mean', {}, 'mae_seed_sd', {}, ...
    'mae_best_seed', {}, 'mae_worst_seed', {}, 'thd_mean', {}, 'thd_seed_sd', {}, 'thd_best_seed', {}, ...
    'thd_worst_seed', {});
for m = meth
    for b = bud
        s = strcmp(rows(:, mcol), m{1}) & cell2mat(rows(:, bcol)) == b;
        if ~any(s), continue; end
        e = cell2mat(rows(s, ecol));  t = cell2mat(rows(s, tcol));
        r = struct('d', d, 'method', m{1}, 'budget', b, 'branches', sum(s), 'mae_mean', mean(e), 'mae_seed_sd', 0, ...
            'mae_best_seed', mean(e), 'mae_worst_seed', mean(e), 'thd_mean', mean(t), 'thd_seed_sd', 0, ...
            'thd_best_seed', mean(t), 'thd_worst_seed', mean(t));
        ss = strcmp(seedrows(:, smcol), m{1}) & cell2mat(seedrows(:, sbcol)) == b;
        if any(ss)
            [r.mae_seed_sd, r.mae_best_seed, r.mae_worst_seed] = seedstats(seedrows(ss, :), sbr, secol);
            [r.thd_seed_sd, r.thd_best_seed, r.thd_worst_seed] = seedstats(seedrows(ss, :), sbr, stcol);
        end
        T(end+1) = r; %#ok<AGROW>
    end
end
if ~isempty(ext) && isfield(ext, 'seed_rows') && ~isempty(ext.seed_rows)
    er = ext.data;
    for b = bud
        s = cell2mat(er(:, bcol)) == b;
        if ~any(s), continue; end
        ss = cell2mat(ext.seed_rows(:, sbcol)) == b;
        r = struct('d', d, 'method', 'tanhLM (wide search)', 'budget', b, 'branches', sum(s), ...
            'mae_mean', mean(cell2mat(er(s, ecol))), 'mae_seed_sd', 0, 'mae_best_seed', NaN, 'mae_worst_seed', NaN, ...
            'thd_mean', mean(cell2mat(er(s, tcol))), 'thd_seed_sd', 0, 'thd_best_seed', NaN, 'thd_worst_seed', NaN);
        [r.mae_seed_sd, r.mae_best_seed, r.mae_worst_seed] = seedstats(ext.seed_rows(ss, :), sbr, secol);
        [r.thd_seed_sd, r.thd_best_seed, r.thd_worst_seed] = seedstats(ext.seed_rows(ss, :), sbr, stcol);
        T(end+1) = r; %#ok<AGROW>
    end
end
end

function [sd, best, worst] = seedstats(sr, seedcol, vcol)
seeds = cell2mat(sr(:, seedcol));  v = cell2mat(sr(:, vcol));
us = unique(seeds).';  avg = zeros(size(us));
for i = 1:numel(us), avg(i) = mean(v(seeds == us(i))); end
if isempty(avg), sd = NaN; best = NaN; worst = NaN; return; end
sd = std(avg);  best = min(avg);  worst = max(avg);
end

function C = cases(rows, seedrows, d)
interp = {'linear', 'spline', 'pchip', 'RBF'};  learned = {'GP', 'tanhLM', 'MLP', 'ANFIS'};
if d == 1, mcol = 1; brc = 2; bcol = 3; ecol = 5; tcol = 8; smcol = 1; sbr = 2; sbcol = 3; secol = 5; stcol = 7;
else,      mcol = 2; brc = 3; bcol = 4; ecol = 7; tcol = 10; smcol = 2; sbr = 3; sbcol = 4; secol = 7; stcol = 9; end
C = struct('d', d, 'n', 0, 'learned_wins_mae', 0, 'learned_wins_thd', 0, 'bestseed_wins_mae', 0, ...
           'bestseed_wins_thd', 0, 'list', zeros(0, 8));
if isempty(rows), return; end
key = [cell2mat(rows(:, brc)), cell2mat(rows(:, bcol))];
uk = unique(key, 'rows');  C.n = size(uk, 1);
for i = 1:size(uk, 1)
    s = key(:, 1) == uk(i, 1) & key(:, 2) == uk(i, 2);
    mi = s & ismember(rows(:, mcol), interp);  ml = s & ismember(rows(:, mcol), learned);
    if ~any(mi) || ~any(ml), continue; end
    bi_e = min(cell2mat(rows(mi, ecol)));  bi_t = min(cell2mat(rows(mi, tcol)));
    bl_e = min(cell2mat(rows(ml, ecol)));  bl_t = min(cell2mat(rows(ml, tcol)));
    ss = cell2mat(seedrows(:, sbr)) == uk(i, 1) & cell2mat(seedrows(:, sbcol)) == uk(i, 2) & ...
         ismember(seedrows(:, smcol), learned);
    bs_e = bl_e;  bs_t = bl_t;
    if any(ss), bs_e = min(bs_e, min(cell2mat(seedrows(ss, secol))));  bs_t = min(bs_t, min(cell2mat(seedrows(ss, stcol)))); end
    C.learned_wins_mae = C.learned_wins_mae + (bl_e < bi_e);
    C.learned_wins_thd = C.learned_wins_thd + (bl_t < bi_t);
    C.bestseed_wins_mae = C.bestseed_wins_mae + (bs_e < bi_e);
    C.bestseed_wins_thd = C.bestseed_wins_thd + (bs_t < bi_t);
    C.list(end+1, :) = [uk(i, :), bi_e, bl_e, bs_e, bi_t, bl_t, bs_t];
end
end
