function M = shm_mlp(X, Y, opts)
%SHM_MLP Multilayer perceptron regressor (no toolbox).
%   M = SHM_MLP(X, Y, opts) fits a network from inputs X (n x d) to targets Y
%   (n x m); Yq = M.predict(Xq). Full-batch Adam with cosine decay, targets
%   detrended and standardised per column, early stopping on a held-out part
%   of the training data.
%
%   Options: HID ([128 128]), ACT ('sin'), W0 (12), EPOCHS (8000), LR (1e-2),
%            L2 (1e-7), FOURIER (0), DETREND (true), VAL_FRAC (0.15),
%            PATIENCE (800), ENSEMBLE (1), SEED (1), verbose

if nargin < 3, opts = struct(); end
D.HID      = [128 128];
D.EPOCHS   = 8000;
D.ENSEMBLE = 1;
D.LR       = 1e-2;
D.L2       = 1e-7;
D.FOURIER  = 0;
D.ACT      = 'sin';
D.W0       = 12;
D.DETREND  = true;
D.VAL_FRAC = 0.15;
D.PATIENCE = 800;
D.SEED     = 1;
D.verbose  = false;
fn = fieldnames(opts);
for i = 1:numel(fn), D.(fn{i}) = opts.(fn{i}); end

if D.ENSEMBLE > 1
    k = D.ENSEMBLE;  o = D;  o.ENSEMBLE = 1;  ms = cell(1,k);
    for e = 1:k, o.SEED = D.SEED + e - 1;  ms{e} = shm_mlp(X, Y, o); end
    M.predict = @(Xq) ensavg(ms, Xq);
    M.members = ms;
    M.nparam  = sum(cellfun(@(z) z.nparam, ms));
    M.best_val = mean(cellfun(@(z) z.best_val, ms));
    M.opts = D;
    return;
end

rng_saved = rng;  rng(D.SEED);

xlo = min(X,[],1);  xhi = max(X,[],1);
xsc = max(xhi - xlo, eps);
if D.DETREND
    Vt = [ones(size(X,1),1), X];
    w = warning('off','all');  Ctr = Vt \ Y;  warning(w);
else
    Ctr = zeros(size(X,2)+1, size(Y,2));
end
trend = @(Z) [ones(size(Z,1),1), Z] * Ctr;
Yd    = Y - trend(X);

ymu = mean(Yd,1);   ysd = max(std(Yd,0,1), eps);

feat = @(Z) mkfeat(Z, xlo, xsc, D.FOURIER);
Xf   = feat(X);
Yn   = (Yd - ymu) ./ ysd;

[n, ~]  = size(Xf);
din     = size(Xf,2);
dout    = size(Yn,2);
sz      = [din, D.HID(:).', dout];
L       = numel(sz) - 1;

W = cell(1,L);  b = cell(1,L);
mW = cell(1,L); vW = cell(1,L); mB = cell(1,L); vB = cell(1,L);
siren = strcmpi(D.ACT, 'sin');
for l = 1:L
    if siren
        if l == 1, c = 1/sz(l); else, c = sqrt(6/sz(l))/D.W0; end
        W{l} = (2*rand(sz(l+1), sz(l)) - 1) * c;
    else
        W{l} = randn(sz(l+1), sz(l)) * sqrt(1/sz(l));
    end
    b{l} = zeros(sz(l+1), 1);
    mW{l} = zeros(size(W{l})); vW{l} = zeros(size(W{l}));
    mB{l} = zeros(size(b{l})); vB{l} = zeros(size(b{l}));
end

nv = max(2, round(D.VAL_FRAC*n));
vi = round(linspace(2, n-1, nv));
ti = setdiff(1:n, vi);
Xt = Xf(ti,:).';  Yt = Yn(ti,:).';
Xv = Xf(vi,:).';  Yv = Yn(vi,:).';

bestW = W; bestB = b; bestL = Inf; bestEp = 0;
hist  = nan(D.EPOCHS, 2);
beta1 = 0.9; beta2 = 0.999; epsA = 1e-8;

for ep = 1:D.EPOCHS
    lr = D.LR * 0.5 * (1 + cos(pi*(ep-1)/D.EPOCHS));

    a = cell(1,L+1);  a{1} = Xt;  zc = cell(1,L);
    for l = 1:L
        z = W{l}*a{l} + b{l};
        if siren && l == 1, z = D.W0 * z; end
        zc{l} = z;
        if l < L
            if siren, a{l+1} = sin(z); else, a{l+1} = tanh(z); end
        else
            a{l+1} = z;
        end
    end
    E  = a{L+1} - Yt;
    tl = mean(E(:).^2);

    dl = (2/numel(E)) * E;
    for l = L:-1:1
        gW = dl*a{l}.' + D.L2*W{l};
        gB = sum(dl, 2);
        if l > 1
            if siren
                dz = cos(zc{l-1});
                if l-1 == 1, dz = D.W0 * dz; end
                dl = (W{l}.'*dl) .* dz;
            else
                dl = (W{l}.'*dl) .* (1 - a{l}.^2);
            end
        end
        mW{l} = beta1*mW{l} + (1-beta1)*gW;
        vW{l} = beta2*vW{l} + (1-beta2)*gW.^2;
        mB{l} = beta1*mB{l} + (1-beta1)*gB;
        vB{l} = beta2*vB{l} + (1-beta2)*gB.^2;
        mh = mW{l}/(1-beta1^ep);  vh = vW{l}/(1-beta2^ep);
        W{l} = W{l} - lr * mh ./ (sqrt(vh) + epsA);
        mh = mB{l}/(1-beta1^ep);  vh = vB{l}/(1-beta2^ep);
        b{l} = b{l} - lr * mh ./ (sqrt(vh) + epsA);
    end

    av = fwdcore(Xv, W, b, L, siren, D.W0);
    ev = av - Yv;  vl = mean(ev(:).^2);
    hist(ep,:) = [tl, vl];
    if vl < bestL*(1 - 1e-6)
        bestL = vl; bestW = W; bestB = b; bestEp = ep;
    elseif ep - bestEp > D.PATIENCE
        hist = hist(1:ep,:);
        break;
    end
end

W = bestW; b = bestB;
rng(rng_saved);

M.predict  = @(Xq) fwdcore(feat(Xq).', W, b, L, siren, D.W0).' .* ysd + ymu + trend(Xq);
M.W = W;  M.b = b;  M.L = L;
M.hist = hist;  M.best_val = bestL;  M.best_epoch = bestEp;
M.opts = D;
M.nparam = sum(cellfun(@numel, W)) + sum(cellfun(@numel, b));
if D.verbose
    fprintf('    MLP %s: %d params, best val MSE %.3e at epoch %d\n', ...
            mat2str(sz), M.nparam, bestL, bestEp);
end
end

function Y = ensavg(ms, Xq)
Y = ms{1}.predict(Xq);
for e = 2:numel(ms), Y = Y + ms{e}.predict(Xq); end
Y = Y / numel(ms);
end

function a = fwdcore(a, W, b, L, siren, w0)
for l = 1:L
    z = W{l}*a + b{l};
    if siren && l == 1, z = w0 * z; end
    if l < L
        if siren, a = sin(z); else, a = tanh(z); end
    else
        a = z;
    end
end
end

function Z = mkfeat(X, xlo, xsc, nf)
Xs = 2*(X - xlo)./xsc - 1;
Z  = Xs;
if nf > 0
    t = Xs(:,1);
    F = zeros(size(X,1), 2*nf);
    for k = 1:nf
        F(:,2*k-1) = sin(k*pi*t);
        F(:,2*k)   = cos(k*pi*t);
    end
    Z = [Z, F];
end
end
