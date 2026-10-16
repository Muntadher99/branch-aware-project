function V = shm_extend(D, opts)
%SHM_EXTEND Extend the stored branches over more operating parameters.
%   V = SHM_EXTEND(D, opts) continues each stored branch from its trace at
%   r = 1, theta_gap = 0.25 deg, first in the DC-link ratio r at fixed Ma, then
%   in the minimum pulse width. A step that moves an angle by more than
%   JUMP_DEG or does not converge ends the branch in that direction.
%
%   Options: MA (0.70:0.01:1.20), R (linspace(1, 1.3, 11)), GAP (0.25),
%            MAXB (8), JUMP_DEG (0.75), solve, verbose
%   Output:  Ma, r, delta, gap, vol(k).TH / .thd / .ok / .branch, nsolve, time_s

if nargin < 2, opts = struct(); end
O.MA       = 0.70:0.01:1.20;
O.R        = linspace(1, 1.3, 11);
O.GAP      = 0.25;
O.MAXB     = 8;
O.JUMP_DEG = 0.75;
O.solve    = struct();
O.verbose  = true;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

F = shm_model();
P0 = D.P;  N = P0.N;
t0 = tic;
MaV = O.MA(:);  rV = O.R(:);  gV = O.GAP(:);
nM = numel(MaV);  nR = numel(rV);  nG = numel(gV);
jmp = O.JUMP_DEG * pi/180;

sel = D.sel(:).';
own = D.owns;
score = arrayfun(@(k) sum(own == k), sel);
[~, o] = sort(score, 'descend');
sel = sel(o);
sel = sel(1:min(O.MAXB, numel(sel)));

V.Ma = MaV;  V.r = rV;  V.delta = (rV - 1) ./ (rV + 1);  V.gap = gV;
V.vol = repmat(struct('TH', [], 'thd', [], 'ok', [], 'branch', 0), 1, numel(sel));
nsolve = 0;

for v = 1:numel(sel)
    b  = sel(v);
    TH = nan(nM, nR, nG, N);  T = nan(nM, nR, nG);  OK = false(nM, nR, nG);
    for i = 1:nM
        [dg, gi] = min(abs(D.grid - MaV(i)));
        if dg > 1e-9, continue; end
        x = D.TH(gi, :, b);
        if any(~isfinite(x)), continue; end
        xr_prev = [];  xr = x;
        for j = 1:nR
            Pr = F.set_r(P0, rV(j));
            if j > 1
                pred = xr;
                if ~isempty(xr_prev), pred = xr + (xr - xr_prev); end
                [xn, ok] = step(F, pred, xr, MaV(i), Pr, jmp, O);  nsolve = nsolve + 1;
                if ~ok, break; end
                xr_prev = xr;  xr = xn;
            end
            xg_prev = [];  xg = xr;
            for k = 1:nG
                Pg = F.set_gap(Pr, gV(k));
                if k > 1 || abs(gV(1) - P0.gapdeg) > 1e-12
                    pred = xg;
                    if ~isempty(xg_prev), pred = xg + (xg - xg_prev); end
                    pred = F.feasible(pred, Pg);
                    [xn, ok] = step(F, pred, F.feasible(xg, Pg), MaV(i), Pg, jmp, O);  nsolve = nsolve + 1;
                    if ~ok, break; end
                    xg_prev = xg;  xg = xn;
                end
                TH(i, j, k, :) = reshape(xg, 1, 1, 1, N);
                T(i, j, k)     = F.thd(xg, Pg);
                OK(i, j, k)    = true;
            end
        end
    end
    V.vol(v).TH = TH;  V.vol(v).thd = T;  V.vol(v).ok = OK;  V.vol(v).branch = b;
    if O.verbose
        fprintf('    branch %2d: %5.1f%% of the %dx%dx%d box, mean THD49 %.3f%%\n', ...
            b, 100*mean(OK(:)), nM, nR, nG, 100*mean(T(OK)));
    end
end
V.nsolve = nsolve;
V.time_s = toc(t0);
end

function [x, ok] = step(F, pred, cur, Ma, P, jmp, O)
[x, I] = F.solve(pred, Ma, P, O.solve);
ok = I.ok && max(abs(x - pred)) <= jmp && max(abs(x - cur)) <= 3*jmp;
if ok, return; end
[x, I] = F.solve(cur, Ma, P, O.solve);
ok = I.ok && max(abs(x - cur)) <= jmp;
end
