function E = shm_x_select2(axes_, T, opts)
%SHM_X_SELECT2 Stage E4: selection of the owning branch over d parameters.
%   E = SHM_X_SELECT2({Ma}, T, opts)               d = 1, T is nM x K
%   E = SHM_X_SELECT2({Ma, delta}, T, opts)        d = 2
%   E = SHM_X_SELECT2({Ma, delta, gap}, T, opts)   d = 3
%   T holds the THD49 of each stored branch (NaN where absent); the label is
%   the branch with the lowest THD49. The scored points are split once into NF
%   disjoint folds; per fold the training set is a q^d lattice in the
%   normalised box. Selectors: nearest training label, k-NN, CART tree, softmax
%   classifier (SHM_CLSNET, one fit per seed), fixed branch. Score: THD49
%   penalty over the best branch (pp); 95% intervals by paired bootstrap.
%
%   Options: Q ([3 4 5 6]), NF (5), SEEDS (1:5), TEST_SEED (12345),
%            MIN_ALIVE (2), NBOOT (2000), KLIST ([1 3 5 7]), DEPTHS (1:8),
%            NAME, clsnet, verbose

if nargin < 3, opts = struct(); end
d = numel(axes_);
O.Q         = [3 4 5 6];
O.NF        = 5;
O.SEEDS     = 1:5;
O.TEST_SEED = 12345;
O.MIN_ALIVE = 2;
O.NBOOT     = 2000;
O.NAME      = sprintf('selection2_%dD', d);
O.clsnet    = struct();
O.KLIST     = [1 3 5 7];
O.DEPTHS    = 1:8;
O.verbose   = true;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

sz = zeros(1, d);
for a = 1:d, sz(a) = numel(axes_{a}); end
K  = size(T, d + 1);
Tf = reshape(T, prod(sz), K);
[bestT, own] = min(Tf, [], 2);
alive = sum(isfinite(Tf), 2);
good  = isfinite(bestT) & alive >= O.MIN_ALIVE;
sub = cell(1, d);
[sub{:}] = ind2sub([sz, 1], (1:prod(sz)).');
X = zeros(prod(sz), d);
for a = 1:d, X(:, a) = axes_{a}(sub{a}); end
X = X(good, :);  y = own(good);  Tg = Tf(good, :);  Bg = bestT(good);
n = numel(y);
mn = min(X, [], 1);  rg = max(X, [], 1) - mn;  rg(rg == 0) = 1;
Xs = (X - mn) ./ rg;

sv = rng;  rng(O.TEST_SEED);
fold = zeros(n, 1);  pm = randperm(n);
for f = 1:O.NF, fold(pm(f:O.NF:n)) = f; end
bootI = randi(n, n, O.NBOOT);
rng(sv);

names = {'lookup', 'kNN', 'tree', 'classifier', 'fixed'};
hdr = {'dimension', 'q', 'mean_training_labels', 'selector', 'penalty_mean_pp', ...
       'reduction_vs_lookup_pct', 'ci95_lo_pct', 'ci95_hi_pct', 'seed_sd_of_reduction_pct', ...
       'accuracy_pct', 'chosen_k_or_depth'};
rows = {};
if O.verbose
    fprintf('  %s: %d scored points, %d branches, %d disjoint folds, %d classifier seeds\n', ...
        O.NAME, n, K, O.NF, numel(O.SEEDS));
end
if n < 2 * O.NF
    fprintf('  %s: too few scored points for %d folds; not scored\n', O.NAME, O.NF);
    E.name = O.NAME;  E.header = hdr;  E.data = rows;  E.n_scored = n;  E.opts = O;
    return
end
for q = O.Q
    L = lattice(axes_, q);  Ls = (L - mn) ./ rg;
    pen = nan(n, 4);  penC = nan(n, numel(O.SEEDS));  hit = nan(n, 4);  hitC = nan(n, numel(O.SEEDS));
    ntr = zeros(O.NF, 1);  kch = zeros(O.NF, 1);  dch = zeros(O.NF, 1);
    for f = 1:O.NF
        te = find(fold == f);  cand = find(fold ~= f);
        pick = unique(cand(nearest_row(Ls, Xs(cand, :))));
        ntr(f) = numel(pick);
        Xtr = Xs(pick, :);  ytr = y(pick);  Xte = Xs(te, :);
        P = Tg(te, :);  B = Bg(te);
        yh = ytr(nearest_row(Xte, Xtr));
        [pen(te, 1), hit(te, 1)] = score(yh, y(te), P, B);
        [kb, yh] = knn_fit_predict(Xtr, ytr, Xte, O.KLIST);  kch(f) = kb;
        [pen(te, 2), hit(te, 2)] = score(yh, y(te), P, B);
        [db, yh] = tree_fit_predict(Xtr, ytr, Xte, O.DEPTHS);  dch(f) = db;
        [pen(te, 3), hit(te, 3)] = score(yh, y(te), P, B);
        [pen(te, 4), hit(te, 4)] = score(mode(ytr) * ones(numel(te), 1), y(te), P, B);
        for si = 1:numel(O.SEEDS)
            if numel(unique(ytr)) == 1
                yh = ytr(1) * ones(numel(te), 1);
            else
                co = O.clsnet;  co.SEED = O.SEEDS(si);
                M = shm_clsnet(X(pick, :), ytr, co);
                yh = M.predict(X(te, :));
            end
            [penC(te, si), hitC(te, si)] = score(yh, y(te), P, B);
        end
    end
    mL = mean(pen(:, 1));
    pCm = mean(penC, 2);
    sets = {pen(:, 1), pen(:, 2), pen(:, 3), pCm, pen(:, 4)};
    accs = [mean(hit(:, 1)), mean(hit(:, 2)), mean(hit(:, 3)), mean(hitC(:)), mean(hit(:, 4))];
    seedred = 100 * (mL - mean(penC, 1)) / max(mL, eps);
    for j = 1:numel(names)
        p = sets{j};
        red = 100 * (mL - mean(p)) / max(mL, eps);
        bL = mean(pen(bootI, 1), 1);  bX = mean(p(bootI), 1);
        br = 100 * (bL - bX) ./ max(bL, eps);
        ci = prctile_(br, [2.5 97.5]);
        ssd = NaN;  ch = NaN;
        if strcmp(names{j}, 'classifier'), ssd = std(seedred); end
        if strcmp(names{j}, 'kNN'), ch = mode(kch); end
        if strcmp(names{j}, 'tree'), ch = mode(dch); end
        if strcmp(names{j}, 'lookup'), red = 0;  ci = [0 0]; end
        rows(end+1, :) = {d, q, mean(ntr), names{j}, mean(p), red, ci(1), ci(2), ssd, accs(j), ch}; %#ok<AGROW>
    end
    if O.verbose
        r = rows(end-4:end, :);
        fprintf('   q=%d (%.0f labels): lookup %.4f | kNN %.4f (%+.1f%% [%+.1f, %+.1f], k=%d) | tree %.4f (%+.1f%% [%+.1f, %+.1f], depth %d) | classifier %.4f (%+.1f%% [%+.1f, %+.1f], seed sd %.1f) | fixed %.3f\n', ...
            q, mean(ntr), r{1,5}, r{2,5}, r{2,6}, r{2,7}, r{2,8}, r{2,11}, r{3,5}, r{3,6}, r{3,7}, r{3,8}, r{3,11}, ...
            r{4,5}, r{4,6}, r{4,7}, r{4,8}, r{4,9}, r{5,5});
    end
end
E.name = O.NAME;  E.header = hdr;  E.data = rows;  E.n_scored = n;  E.opts = O;
end

function [c, h] = score(yh, yte, P, B)
n = numel(yh);  c = zeros(n, 1);
for i = 1:n
    t = P(i, yh(i));
    if ~isfinite(t), f = P(i, :);  t = max(f(isfinite(f))); end
    c(i) = 100 * (t - B(i));
end
h = 100 * double(yh(:) == yte(:));
end

function [kb, yh] = knn_fit_predict(Xtr, ytr, Xte, KL)
n = size(Xtr, 1);  KL = KL(KL < n);  if isempty(KL), KL = 1; end
Dtr = sqd(Xtr, Xtr);  Dtr(1:n+1:end) = inf;
err = zeros(size(KL));
for i = 1:numel(KL), err(i) = mean(knn_vote(Dtr, ytr, KL(i)) ~= ytr); end
[~, b] = min(err);  kb = KL(b);
yh = knn_vote(sqd(Xte, Xtr), ytr, kb);
end

function yh = knn_vote(D2, ytr, k)
[ds, ix] = sort(D2, 2);
cls = unique(ytr);  yh = zeros(size(D2, 1), 1);
for i = 1:size(D2, 1)
    nb = ytr(ix(i, 1:k));  w = 1 ./ (sqrt(ds(i, 1:k)) + 1e-9);
    s = zeros(numel(cls), 1);
    for c = 1:numel(cls), s(c) = sum(w(nb == cls(c))); end
    [~, j] = max(s);  yh(i) = cls(j);
end
end

function [db, yh] = tree_fit_predict(Xtr, ytr, Xte, DEP)
n = size(Xtr, 1);  err = zeros(size(DEP));
nf = min(5, n);  fo = mod((0:n-1).', nf) + 1;
for di = 1:numel(DEP)
    e = 0;
    for f = 1:nf
        tr = fo ~= f;  te = ~tr;
        if ~any(tr), continue; end
        t = cart(Xtr(tr, :), ytr(tr), DEP(di));
        e = e + sum(cart_predict(t, Xtr(te, :)) ~= ytr(te));
    end
    err(di) = e;
end
[~, b] = min(err);  db = DEP(b);
t = cart(Xtr, ytr, db);
yh = cart_predict(t, Xte);
end

function node = cart(X, y, depth)
node.leaf = true;  node.label = mode(y);
if depth == 0 || numel(unique(y)) == 1 || numel(y) < 2, return; end
cls = unique(y);  nc = numel(cls);  n = numel(y);
Y = zeros(n, nc);  for c = 1:nc, Y(:, c) = (y == cls(c)); end
best = 1 - sum((sum(Y, 1) / n).^2);  bj = 0;  bt = 0;
for j = 1:size(X, 2)
    [xs, o] = sort(X(:, j));  Cl = cumsum(Y(o, :), 1);  tot = Cl(end, :);
    nl = (1:n-1).';  Cl = Cl(1:n-1, :);  Cr = tot - Cl;  nr = n - nl;
    g = (nl .* (1 - sum((Cl ./ nl).^2, 2)) + nr .* (1 - sum((Cr ./ nr).^2, 2))) / n;
    valid = xs(1:n-1) < xs(2:n);
    g(~valid) = inf;
    [gm, i] = min(g);
    if gm < best - 1e-12, best = gm;  bj = j;  bt = 0.5 * (xs(i) + xs(i + 1)); end
end
if bj == 0, return; end
l = X(:, bj) <= bt;
node.leaf = false;  node.j = bj;  node.t = bt;
node.L = cart(X(l, :), y(l), depth - 1);
node.R = cart(X(~l, :), y(~l), depth - 1);
end

function yh = cart_predict(t, X)
yh = zeros(size(X, 1), 1);
for i = 1:size(X, 1)
    nd = t;
    while ~nd.leaf
        if X(i, nd.j) <= nd.t, nd = nd.L; else, nd = nd.R; end
    end
    yh(i) = nd.label;
end
end

function D2 = sqd(A, B)
D2 = sum(A.^2, 2) + sum(B.^2, 2).' - 2 * (A * B.');
D2(D2 < 0) = 0;
end

function L = lattice(axes_, q)
d = numel(axes_);  v = cell(1, d);
for a = 1:d
    ax = axes_{a}(:);
    v{a} = ax(unique(round(linspace(1, numel(ax), min(q, numel(ax))))));
end
g = cell(1, d);  [g{:}] = ndgrid(v{:});
L = zeros(numel(g{1}), d);
for a = 1:d, L(:, a) = g{a}(:); end
end

function j = nearest_row(A, B)
[~, j] = min(sqd(A, B), [], 2);  j = j(:);
end

function p = prctile_(x, q)
x = sort(x(:));  n = numel(x);
p = zeros(size(q));
for i = 1:numel(q)
    r = q(i) / 100 * (n - 1) + 1;  lo = floor(r);  hi = min(ceil(r), n);
    p(i) = x(lo) + (r - lo) * (x(hi) - x(lo));
end
end
