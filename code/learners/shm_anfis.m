function M = shm_anfis(X, Y, opts)
%SHM_ANFIS First-order Sugeno ANFIS regressor (no toolbox).
%   M = SHM_ANFIS(X, Y, opts); Yq = M.predict(Xq). Grid partition with MF
%   Gaussian membership functions per input (MF^d rules), linear consequents
%   by least squares, membership parameters trained by Adam with early
%   stopping on a held-out part of the training data.
%
%   Options: MF (per input; from the data if empty), EPOCHS (150), LR (0.02),
%            RIDGE (1e-8), VAL_FRAC (0.15), verbose

if nargin < 3, opts = struct(); end
[n, d] = size(X);
O.MF = [];  O.EPOCHS = 150;  O.LR = 0.02;  O.RIDGE = 1e-8;  O.VAL_FRAC = 0.15;  O.verbose = false;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

nv = max(2, round(O.VAL_FRAC * n));
vi = unique(round(linspace(2, n - 1, nv)));
ti = setdiff(1:n, vi);
if numel(ti) < 3, ti = 1:n;  vi = []; end

if isempty(O.MF)
    mf = 2;
    while mf < 7 && (mf + 1)^d * (d + 1) <= max(numel(ti), 4), mf = mf + 1; end
    O.MF = mf;
end
mf = O.MF;

xlo = min(X, [], 1);  xsc = max(X, [], 1) - xlo;  xsc(xsc == 0) = 1;
Xs  = (X - repmat(xlo, n, 1)) ./ repmat(xsc, n, 1);
ymu = mean(Y(ti, :), 1);  ysd = std(Y(ti, :), 0, 1);  ysd(ysd == 0) = 1;
Ys  = (Y - repmat(ymu, n, 1)) ./ repmat(ysd, n, 1);

C = repmat(linspace(0, 1, mf).', 1, d);
S = log(0.5 / max(mf - 1, 1)) * ones(mf, d);
idx = rule_index(mf, d);
th  = [C(:); S(:)];
np  = numel(th);

loss = @(t) fit_loss(t, Xs(ti, :), Ys(ti, :), idx, mf, d, O.RIDGE);
m1 = zeros(np, 1);  m2 = zeros(np, 1);  b1 = 0.9;  b2 = 0.999;  h = 1e-5;
best = inf;  bestth = th;
for ep = 1:O.EPOCHS
    g = zeros(np, 1);
    for k = 1:np
        e = zeros(np, 1);  e(k) = h;
        g(k) = (loss(th + e) - loss(th - e)) / (2*h);
    end
    m1 = b1*m1 + (1-b1)*g;  m2 = b2*m2 + (1-b2)*g.^2;
    th = th - O.LR * (m1/(1-b1^ep)) ./ (sqrt(m2/(1-b2^ep)) + 1e-8);
    if isempty(vi)
        vl = loss(th);
    else
        P  = anfis_consequents(th, Xs(ti, :), Ys(ti, :), idx, mf, d, O.RIDGE);
        Yv = anfis_forward(th, P, Xs(vi, :), idx, mf, d);
        vl = mean(mean((Yv - Ys(vi, :)).^2));
    end
    if vl < best, best = vl;  bestth = th; end
end
th = bestth;
P  = anfis_consequents(th, Xs(ti, :), Ys(ti, :), idx, mf, d, O.RIDGE);

M.predict = @(Xq) predict_(Xq, th, P, idx, mf, d, xlo, xsc, ymu, ysd);
M.nparam  = np + numel(P);
M.mf      = mf;
M.rules   = size(idx, 1);
M.best_val = best;
if O.verbose
    fprintf('    ANFIS: %d inputs, %d MF each, %d rules, %d parameters\n', d, mf, size(idx, 1), M.nparam);
end
end

function Yq = predict_(Xq, th, P, idx, mf, d, xlo, xsc, ymu, ysd)
nq = size(Xq, 1);
Xs = (Xq - repmat(xlo, nq, 1)) ./ repmat(xsc, nq, 1);
Yq = anfis_forward(th, P, Xs, idx, mf, d) .* repmat(ysd, nq, 1) + repmat(ymu, nq, 1);
end

function idx = rule_index(mf, d)
g = cell(1, d);
[g{:}] = ndgrid(1:mf);
idx = zeros(mf^d, d);
for j = 1:d, idx(:, j) = g{j}(:); end
end

function W = anfis_firing(th, Xs, idx, mf, d)
C = reshape(th(1:mf*d), mf, d);  S = exp(reshape(th(mf*d+1:end), mf, d));
n = size(Xs, 1);  R = size(idx, 1);
W = ones(n, R);
for j = 1:d
    mu = exp(-0.5 * ((repmat(Xs(:, j), 1, mf) - repmat(C(:, j).', n, 1)) ./ repmat(S(:, j).', n, 1)).^2);
    W = W .* mu(:, idx(:, j));
end
W = W ./ repmat(sum(W, 2) + 1e-300, 1, R);
end

function Phi = anfis_design(th, Xs, idx, mf, d)
W = anfis_firing(th, Xs, idx, mf, d);
n = size(Xs, 1);  R = size(idx, 1);
Xa = [ones(n, 1), Xs];
Phi = zeros(n, R * (d + 1));
for r = 1:R
    Phi(:, (r-1)*(d+1) + (1:d+1)) = Xa .* repmat(W(:, r), 1, d + 1);
end
end

function P = anfis_consequents(th, Xs, Ys, idx, mf, d, ridge)
Phi = anfis_design(th, Xs, idx, mf, d);
A = Phi.' * Phi;
P = (A + ridge * max(trace(A) / size(A, 1), 1) * eye(size(A, 1))) \ (Phi.' * Ys);
end

function Y = anfis_forward(th, P, Xs, idx, mf, d)
Y = anfis_design(th, Xs, idx, mf, d) * P;
end

function L = fit_loss(th, Xs, Ys, idx, mf, d, ridge)
P = anfis_consequents(th, Xs, Ys, idx, mf, d, ridge);
E = anfis_forward(th, P, Xs, idx, mf, d) - Ys;
L = mean(E(:).^2);
end
