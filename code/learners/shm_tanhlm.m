function M = shm_tanhlm(X, Y, opts)
%SHM_TANHLM One-hidden-layer tanh network trained by Levenberg-Marquardt.
%   M = SHM_TANHLM(X, Y, opts); Yq = M.predict(Xq). Inputs scaled to [-1, 1],
%   targets detrended and standardised. The hidden width is chosen from WIDTHS
%   by cross-validation on the training samples (leave-one-out for n <= 12,
%   otherwise KFOLD folds) and refitted on all of them; each fit keeps the
%   best of RESTARTS initialisations (CV_RESTARTS for the validation fits).
%
%   Options: WIDTHS ([2 4 8]), RESTARTS (2), CV_RESTARTS ([] = RESTARTS),
%            ITERS (200), L2 (1e-4), KFOLD (5), SEED (1), verbose

if nargin < 3, opts = struct(); end
O.WIDTHS   = [2 4 8];
O.RESTARTS = 2;
O.CV_RESTARTS = [];
O.ITERS    = 200;
O.L2       = 1e-4;
O.KFOLD    = 5;
O.SEED     = 1;
O.verbose  = false;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

[n, d] = size(X);  m = size(Y, 2);
lo = min(X, [], 1);  sc = max(X, [], 1) - lo;  sc(sc == 0) = 1;
Xn = 2 * (X - lo) ./ sc - 1;
A  = [ones(n, 1), Xn];
if n > d + 1, Bt = A \ Y; else, Bt = [mean(Y, 1); zeros(d, m)]; end
Rz = Y - A * Bt;
sy = std(Rz, 0, 1);  sy(sy < 1e-12) = 1;
Z  = Rz ./ sy;

rng_saved = rng;  rng(O.SEED);
if n <= 12, folds = num2cell(1:n); else
    p = randperm(n);  folds = cell(1, O.KFOLD);
    for k = 1:O.KFOLD, folds{k} = p(k:O.KFOLD:n); end
end
cvr = O.RESTARTS;  if ~isempty(O.CV_RESTARTS), cvr = O.CV_RESTARTS; end
cvs = inf(1, numel(O.WIDTHS));
for wi = 1:numel(O.WIDTHS)
    H = O.WIDTHS(wi);  e = 0;
    for k = 1:numel(folds)
        te = folds{k};  tr = setdiff(1:n, te);
        if numel(tr) < 2, continue; end
        th = fit(Xn(tr, :), Z(tr, :), H, O, cvr);
        Zp = fwd(th, Xn(te, :), H, d, m);
        e = e + sum(sum(((Zp - Z(te, :)) .* sy).^2));
    end
    cvs(wi) = e;
end
[~, wb] = min(cvs);  H = O.WIDTHS(wb);
th = fit(Xn, Z, H, O, O.RESTARTS);
rng(rng_saved);

M.width   = H;
M.cv      = cvs;
M.predict = @(Xq) [ones(size(Xq, 1), 1), 2 * (Xq - lo) ./ sc - 1] * Bt + ...
                  fwd(th, 2 * (Xq - lo) ./ sc - 1, H, d, m) .* sy;
M.nparam  = numel(th) + numel(Bt);
end

function thb = fit(X, Z, H, O, nr)
[n, d] = size(X);  m = size(Z, 2);
P = H * d + H + m * H + m;
best = inf;  thb = [];
for r = 1:nr
    th = [randn(H * d, 1) / sqrt(d); 0.1 * randn(H, 1); randn(m * H, 1) / sqrt(H); zeros(m, 1)];
    mu = 1e-2;
    f = obj(th, X, Z, H, O.L2);
    for it = 1:O.ITERS
        [rv, J] = resjac(th, X, Z, H, d, m, O.L2);
        g = J' * rv;  JJ = J' * J;
        improved = false;
        for tries = 1:12
            dt = -(JJ + mu * eye(P)) \ g;
            fn = obj(th + dt, X, Z, H, O.L2);
            if fn < f
                th = th + dt;  mu = max(mu / 3, 1e-9);  improved = true;
                rel = (f - fn) / max(f, eps);  f = fn;
                break
            end
            mu = mu * 5;
        end
        if ~improved || rel < 1e-7, break; end
    end
    if f < best, best = f;  thb = th; end
end
end

function f = obj(th, X, Z, H, l2)
[n, d] = size(X);  m = size(Z, 2);
E = fwd(th, X, H, d, m) - Z;
f = sum(E(:).^2) + l2 * sum(th.^2);
end

function Zp = fwd(th, X, H, d, m)
[W1, b1, W2, b2] = unpack(th, H, d, m);
Zp = tanh(X * W1' + b1') * W2' + b2';
end

function [rv, J] = resjac(th, X, Z, H, d, m, l2)
[W1, b1, W2, b2] = unpack(th, H, d, m);
n = size(X, 1);
S = tanh(X * W1' + b1');
E = S * W2' + b2' - Z;
dS = 1 - S.^2;
P = numel(th);
J = zeros(n * m, P);
for o = 1:m
    rows = (o - 1) * n + (1:n);
    G = dS .* W2(o, :);
    for j = 1:d
        J(rows, (j - 1) * H + (1:H)) = G .* X(:, j);
    end
    J(rows, H * d + (1:H)) = G;
    J(rows, H * d + H + (0:H-1) * m + o) = S;
    J(rows, H * d + H + m * H + o) = 1;
end
rv = [E(:); sqrt(l2) * th];
J  = [J; sqrt(l2) * eye(P)];
end

function [W1, b1, W2, b2] = unpack(th, H, d, m)
k = 0;
W1 = reshape(th(k + (1:H*d)), H, d);  k = k + H*d;
b1 = th(k + (1:H));                   k = k + H;
W2 = reshape(th(k + (1:m*H)), m, H);  k = k + m*H;
b2 = th(k + (1:m));
end
