function E = shm_x_folds(D, opts)
%SHM_X_FOLDS Stage E1b: the interior ends of the branches that own the envelope.
%   E = SHM_X_FOLDS(D, opts) follows the smallest eigenvalue of the reduced
%   Hessian of the Lagrangian over the last NBACK grid points before each end
%   (a fold shows as a value falling to zero) and tries to continue the branch
%   past its end with steps of 1/REFINE grid step; continued solutions are
%   compared with every other branch.
%
%   Options: NBACK (4), REFINE (64), ACCEPT (0.5 deg), TOL (0.05 deg),
%            EXTRA (3), verbose

if nargin < 2, opts = struct(); end
O.NBACK  = 4;
O.REFINE = 64;
O.ACCEPT = 0.5;
O.TOL    = 0.05;
O.EXTRA  = 3;
O.verbose = true;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

F = shm_model();  P = D.P;  N = P.N;  grid = D.grid(:);  G = numel(grid);  h0 = grid(2) - grid(1);
owners = unique(D.own(isfinite(D.own))).';
hdr = {'branch', 'end', 'Ma_end', 'eig_rel_at_end', 'eig_rel_trend', 'kind', 'continued_steps', ...
       'joins_branch', 'active_constraints_at_end', 'eig_rel_at_refined_end', 'Ma_refined_end'};
rows = {};
for k = owners
    f = find(isfinite(D.THD(:, k)));
    if numel(f) < 2, continue; end
    ev = nan(numel(f), 1);
    for i = 1:numel(f), ev(i) = redeig(F, P, D.TH(f(i), :, k), grid(f(i))); end
    scale = median(ev(isfinite(ev)));
    for side = [-1 1]
        if side < 0, e = f(1); else, e = f(end); end
        if (side < 0 && e == 1) || (side > 0 && e == G), continue; end
        back = e - side * (0:O.NBACK - 1);
        back = back(back >= f(1) & back <= f(end));
        er = arrayfun(@(g) redeig(F, P, D.TH(g, :, k), grid(g)), back) / scale;
        [~, nact] = redeig(F, P, D.TH(e, :, k), grid(e));
        trend = er(1) - er(end);
        if er(1) < 0.05, kind = 'fold'; elseif er(1) < 0.2 && trend < 0, kind = 'fold?'; else, kind = 'other'; end
        th = D.TH(e, :, k);  thp = D.TH(e - side, :, k);  m = grid(e);  mp = grid(e - side);
        step = side * h0 / O.REFINE;  nst = 0;  joined = NaN;  okc = true;
        while okc && abs(m - grid(e)) < O.EXTRA * h0 - 1e-12 && m + step >= grid(1) && m + step <= grid(end)
            pred = th + (th - thp) * (step / (m - mp));
            [x, I] = F.solve(pred, m + step, P);
            if ~I.ok || max(abs(x - sort(pred))) * 180/pi > O.ACCEPT, okc = false; break; end
            thp = th;  mp = m;  th = x;  m = m + step;  nst = nst + 1;
            [dg, gi] = min(abs(grid - m));
            if dg < 1e-9
                Xg = squeeze(D.TH(gi, :, :));  if size(Xg, 1) ~= N, Xg = Xg.'; end
                dd = max(abs(Xg - th(:)), [], 1) * 180/pi;  dd(k) = inf;
                [dm, jb] = min(dd);
                if dm <= O.TOL, joined = jb; break; end
            end
        end
        erl = redeig(F, P, th, m) / scale;
        if erl < 0.05 && ~strcmp(kind, 'fold'), kind = 'fold'; end
        rows(end+1, :) = {k, ifelse(side < 0, 'low', 'high'), grid(e), er(1), trend, kind, nst, joined, nact, erl, m}; %#ok<AGROW>
        if O.verbose
            fprintf('  branch %2d %4s end Ma %.4f: eig %.3f of median (trend %+.3f) | refined: %3d steps to Ma %.5f, eig %.4f -> %-5s | joins %s | active %d\n', ...
                k, rows{end, 2}, grid(e), er(1), trend, nst, m, erl, kind, num2str(joined), nact);
        end
    end
end
E.name = 'E1c_owner_branch_ends';  E.header = hdr;  E.data = rows;  E.opts = O;
end

function [emin, nact] = redeig(F, P, th, Ma)
th = sort(th(:));  N = numel(th);  s = P.s(:);
[~, g] = F.objgrad(th, P);
n = P.ORD(:);
dB = -(4/pi) * sin(n * th.') .* (ones(numel(n), 1) * s.');
bn = (4 ./ (n*pi)) .* (cos(n * th.') * s);
d2 = -(4/pi) * (n .* bn).' * cos(n * th.');
H  = 2 * (dB.' * dB) + 2 * diag(d2 .* s.');
if P.delta ~= 0 && ~isempty(P.ORDE)
    e = P.ORDE(:);
    dA = -P.delta * (4/pi) * cos(e * th.') .* (ones(numel(e), 1) * s.');
    an = -P.delta * (4 ./ (e*pi)) .* (sin(e * th.') * s);
    d2a = P.delta * (4/pi) * (e .* an).' * sin(e * th.');
    H = H + 2 * (dA.' * dA) + 2 * diag(d2a .* s.');
end
c  = -(4/pi) * (s .* sin(th));
Hc = diag(-(4/pi) * (s .* cos(th)));
A = [];
for i = 1:N - 1
    if th(i + 1) - th(i) - P.minGap < 1e-7, r = zeros(1, N); r(i) = -1; r(i + 1) = 1; A = [A; r]; end %#ok<AGROW>
end
if th(1) - P.lb < 1e-7, r = zeros(1, N); r(1) = 1; A = [A; r]; end
if P.ub - th(N) < 1e-7, r = zeros(1, N); r(N) = 1; A = [A; r]; end
nact = size(A, 1);
Cc = [c.'; A];
mult = Cc.' \ g;
HL = H - mult(1) * Hc;
Z = null(Cc);
if isempty(Z), emin = NaN; return; end
emin = min(eig((Z.' * HL * Z + (Z.' * HL * Z).') / 2));
end

function s = ifelse(c, a, b)
if c, s = a; else, s = b; end
end
