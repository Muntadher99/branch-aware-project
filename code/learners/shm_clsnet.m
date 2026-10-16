function M = shm_clsnet(X, y, opts)
%SHM_CLSNET Softmax classifier of the owning branch (no toolbox).
%   M = SHM_CLSNET(X, y, opts) with operating points X (n x d) and branch
%   labels y; [yh, Pr] = M.predict(Xq). Random Fourier features, tanh hidden
%   layer, cross-entropy loss, Adam. Inputs are min-max scaled internally.
%
%   Options: HID (64), NFEAT (64), OMEGA (4), EPOCHS (4000), LR (3e-3),
%            SEED (1), verbose

if nargin < 3, opts = struct(); end
O.HID      = 64;
O.NFEAT    = 64;
O.OMEGA    = 4;
O.EPOCHS   = 4000;
O.LR       = 3e-3;
O.SEED     = 1;
O.verbose  = false;
fn = fieldnames(opts);
for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

y = y(:);
cls = unique(y);
C   = numel(cls);
yi  = zeros(size(y));
for i = 1:C, yi(y == cls(i)) = i; end

mn = min(X, [], 1);
rg = max(X, [], 1) - mn;
rg(rg == 0) = 1;
Xn = (X - repmat(mn, size(X,1), 1)) ./ repmat(rg, size(X,1), 1);

rng_saved = rng;  rng(O.SEED);
d  = size(Xn, 2);
Om = O.OMEGA * randn(d, O.NFEAT);
ph = 2*pi*rand(1, O.NFEAT);
Z  = feats(Xn, Om, ph);
m  = size(Z, 2);
H  = O.HID;
P  = {randn(m,H)/sqrt(m), zeros(1,H), randn(H,H)/sqrt(H), zeros(1,H), ...
      randn(H,C)/sqrt(H), zeros(1,C)};
rng(rng_saved);

n  = size(Z, 1);
Y  = zeros(n, C);
for i = 1:n, Y(i, yi(i)) = 1; end

Ms = cell(1,6); Vs = cell(1,6);
for k = 1:6, Ms{k} = zeros(size(P{k})); Vs{k} = Ms{k}; end
b1 = 0.9; b2 = 0.999; eps0 = 1e-8;

for t = 1:O.EPOCHS
    A1 = tanh(Z*P{1}  + repmat(P{2}, n, 1));
    A2 = tanh(A1*P{3} + repmat(P{4}, n, 1));
    L  = A2*P{5} + repmat(P{6}, n, 1);
    L  = L - repmat(max(L, [], 2), 1, C);
    E  = exp(L);
    Pr = E ./ repmat(sum(E, 2), 1, C);
    dL = (Pr - Y) / n;
    g5 = A2'*dL;  g6 = sum(dL, 1);
    dA2 = (dL*P{5}') .* (1 - A2.^2);
    g3 = A1'*dA2; g4 = sum(dA2, 1);
    dA1 = (dA2*P{3}') .* (1 - A1.^2);
    g1 = Z'*dA1;  g2 = sum(dA1, 1);
    G  = {g1, g2, g3, g4, g5, g6};
    a  = O.LR * 0.5 * (1 + cos(pi*t/O.EPOCHS));
    for k = 1:6
        Ms{k} = b1*Ms{k} + (1-b1)*G{k};
        Vs{k} = b2*Vs{k} + (1-b2)*G{k}.^2;
        P{k}  = P{k} - a*(Ms{k}/(1-b1^t)) ./ (sqrt(Vs{k}/(1-b2^t)) + eps0);
    end
    if O.verbose && mod(t, max(1,round(O.EPOCHS/5))) == 0
        acc = 100*mean(argmaxrow(Pr) == yi);
        fprintf('      epoch %5d  train acc %.1f%%\n', t, acc);
    end
end

M.predict = @(Xq) doPredict(Xq, P, Om, ph, mn, rg, cls);
M.classes = cls;
M.nparam  = sum(cellfun(@numel, P));
M.opts    = O;
end

function [yh, Pr] = doPredict(Xq, P, Om, ph, mn, rg, cls)
n  = size(Xq, 1);
Xn = (Xq - repmat(mn, n, 1)) ./ repmat(rg, n, 1);
Z  = feats(Xn, Om, ph);
A1 = tanh(Z*P{1}  + repmat(P{2}, n, 1));
A2 = tanh(A1*P{3} + repmat(P{4}, n, 1));
L  = A2*P{5} + repmat(P{6}, n, 1);
L  = L - repmat(max(L, [], 2), 1, size(L,2));
E  = exp(L);
Pr = E ./ repmat(sum(E, 2), 1, size(E,2));
yh = cls(argmaxrow(Pr));
yh = yh(:);
end

function Z = feats(Xn, Om, ph)
A = Xn*Om + repmat(ph, size(Xn,1), 1);
Z = [cos(A), sin(A), Xn];
end

function j = argmaxrow(A)
[~, j] = max(A, [], 2);
end
