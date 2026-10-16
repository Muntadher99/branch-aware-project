function M = shm_gpr(X, Y, opts)
%SHM_GPR Gaussian-process and radial-basis-function regression.
%   M = SHM_GPR(X, Y, opts); Yq = M.predict(Xq). Zero-mean GP on the residual
%   of a per-output linear trend, outputs standardised, one ARD kernel shared
%   by all outputs, hyperparameters from the training data by marginal
%   likelihood (OBJ 'ml') or leave-one-out error (OBJ 'loo'). KERNEL 'se',
%   'matern52', or 'phs' (cubic polyharmonic RBF interpolant, no
%   hyperparameters).
%
%   Options: OBJ ('ml'), KERNEL ('se'), NOISE_FLOOR (1e-6), STARTS, verbose

if nargin < 3, opts = struct(); end
O.OBJ         = 'ml';
O.KERNEL      = 'se';
O.NOISE_FLOOR = 1e-6;
O.STARTS      = [0.1 0.3 1.0];
O.verbose     = false;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

[n, d] = size(X);  m = size(Y, 2);
lo = min(X, [], 1);  sc = max(X, [], 1) - lo;  sc(sc == 0) = 1;
Xn = (X - lo) ./ sc;
A  = [ones(n, 1), Xn];
if n > d + 1, Bt = A \ Y; else, Bt = [mean(Y, 1); zeros(d, m)]; end
Rz = Y - A * Bt;
sy = std(Rz, 0, 1);  sy(sy < 1e-12) = 1;
Z  = Rz ./ sy;

if strcmp(O.KERNEL, 'phs')
    M = phs_fit(Xn, Y, lo, sc);  return
end
best = inf;  hb = [];
fix = [0, -30];
for l0 = O.STARTS
    h0 = log(l0) * ones(1, d);
    if strcmp(O.OBJ, 'ml'), f = @(hl) nlml([hl, fix], Xn, Z, O.NOISE_FLOOR, O.KERNEL);
    else, f = @(hl) loo([hl, fix], Xn, Z, O.NOISE_FLOOR, O.KERNEL); end
    [hl, v] = fminsearch(f, h0, optimset('MaxFunEvals', 400 * (d + 2), 'MaxIter', 400 * (d + 2), ...
        'TolX', 1e-6, 'TolFun', 1e-10, 'Display', 'off'));
    if v < best, best = v;  hb = [hl, fix]; end
end
[Kc, ~] = kmat(hb, Xn, Xn, O.NOISE_FLOOR, true, O.KERNEL);
L  = chol(Kc, 'lower');
al = L' \ (L \ Z);

M.hyp     = hb;
M.predict = @(Xq) predict_(Xq, lo, sc, Bt, Xn, al, hb, sy, O.NOISE_FLOOR, O.KERNEL);
M.nparam  = n * d + n * m + numel(Bt) + numel(hb);
M.nll     = best;
end

function v = nlml(h, X, Z, floor_, ker)
[K, ok] = kmat(h, X, X, floor_, true, ker);
if ~ok, v = 1e12; return; end
[L, p] = chol(K, 'lower');
if p > 0, v = 1e12; return; end
al = L' \ (L \ Z);
m = size(Z, 2);  n = size(Z, 1);
v = 0.5 * sum(sum(Z .* al)) + m * sum(log(diag(L))) + 0.5 * n * m * log(2*pi);
v = v + 0.5 * sum((h(1:end-2) - log(0.3)).^2) / 4;
end

function v = loo(h, X, Z, floor_, ker)
[K, ok] = kmat(h, X, X, floor_, true, ker);
if ~ok, v = 1e12; return; end
[L, p] = chol(K, 'lower');
if p > 0, v = 1e12; return; end
Ki = L' \ (L \ eye(size(K, 1)));
r  = (Ki * Z) ./ diag(Ki);
v  = sum(r(:).^2) + 1e-3 * sum((h(1:end-2) - log(0.3)).^2);
end

function [K, ok] = kmat(h, A, B, floor_, addnoise, ker)
d = size(A, 2);  ell = exp(h(1:d));  sf2 = exp(2 * h(d + 1));  sn2 = exp(2 * h(d + 2));
ok = all(isfinite([ell, sf2, sn2])) && all(ell > 1e-4) && all(ell < 1e3);
D2 = zeros(size(A, 1), size(B, 1));
for j = 1:d
    D2 = D2 + ((A(:, j) - B(:, j).') / ell(j)).^2;
end
if strcmp(ker, 'se')
    K = sf2 * exp(-0.5 * D2);
else
    r = sqrt(5 * D2);
    K = sf2 * (1 + r + r.^2 / 3) .* exp(-r);
end
if addnoise, K = K + (sn2 + floor_ * sf2) * eye(size(A, 1)); end
end

function Yq = predict_(Xq, lo, sc, Bt, Xn, al, hb, sy, floor_, ker)
Xqn = (Xq - lo) ./ sc;
Ks  = kmat(hb, Xqn, Xn, floor_, false, ker);
Yq  = [ones(size(Xqn, 1), 1), Xqn] * Bt + (Ks * al) .* sy;
end

function M = phs_fit(Xn, Y, lo, sc)
[n, d] = size(Xn);
Phi = rdist(Xn, Xn).^3;  Pm = [ones(n, 1), Xn];
A = [Phi, Pm; Pm', zeros(d + 1)];
C = A \ [Y; zeros(d + 1, size(Y, 2))];
M.hyp = [];  M.nll = NaN;
M.nparam = n * d + numel(C);
M.predict = @(Xq) phs_eval((Xq - lo) ./ sc, Xn, C);
end

function Yq = phs_eval(Xq, Xn, C)
n = size(Xn, 1);
Yq = [rdist(Xq, Xn).^3, ones(size(Xq, 1), 1), Xq] * C;
end

function R = rdist(A, B)
R = sqrt(max(sum(A.^2, 2) + sum(B.^2, 2).' - 2 * (A * B.'), 0));
end
