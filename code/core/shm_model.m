function F = shm_model()
%SHM_MODEL Three-level SHM-PWM waveform model and constrained solver.
%   F = SHM_MODEL() returns a struct of function handles:
%     P = F.defaults(N, nmax, r, gapdeg)  problem data: N angles (15), orders
%                                         up to nmax (49), DC-link ratio r (1),
%                                         minimum pulse width in deg (0.25)
%     [TH, I] = F.solve(TH0, Ma, P, o)    one constrained solve per row of TH0
%     X = F.starts(n, P)                  n random feasible angle sets
%     t = F.thd(TH, P)                    THD49 of each row of TH
%     [Bo, Be, b1] = F.harmonics(TH, P)   odd and even harmonics, fundamental
%     R = F.report(TH, Ma, P)             THD49, b1 error, constraint
%                                         violation, per-order compliance
%
%   Quarter-wave three-level pattern with angles 0 < th_1 < ... < th_N < pi/2
%   and level steps s_k = (-1)^(k-1). In units of Vdc/2,
%       b_n = 4/(n*pi) * sum_k s_k cos(n*th_k),   n odd, and Ma = b_1.
%   THD49 is the root-sum-square of the odd non-triplen orders 5..49 relative
%   to b_1. A DC-link imbalance delta = (r-1)/(r+1) adds the even orders
%       a_n = -delta * 4/(n*pi) * sum_k s_k sin(n*th_k).
%
%   The solver minimises sum(b_n^2) (+ sum(a_n^2)) subject to b_1 = Ma,
%   th_(k+1) - th_k >= theta_gap and theta_gap/2 <= th <= 90 - theta_gap/2 deg,
%   with fmincon (SQP by default; o.algorithm = 'interior-point' for the
%   independent search) and analytic gradients.

F.defaults   = @defaults;
F.set_gap    = @set_gap;
F.set_r      = @set_r;
F.harmonics  = @harmonics;
F.b1         = @b1_of;
F.thd        = @thd_of;
F.report     = @report;
F.spectrum   = @spectrum;
F.solve      = @solve;
F.solve_defaults = @solve_defaults;
F.starts     = @random_starts;
F.feasible   = @make_feasible;
F.objgrad    = @objgrad;
F.eqcon      = @eqcon;
F.backend    = @backend;
F.nfe        = @nfe_get;
F.nfe_reset  = @nfe_reset;
F.delta_of_r = @(r) (r - 1)./(r + 1);
F.r_of_delta = @(d) (1 + d)./(1 - d);
end

function P = defaults(N, nmax, r, gapdeg)
if nargin < 1 || isempty(N),      N = 15;      end
if nargin < 2 || isempty(nmax),   nmax = 49;   end
if nargin < 3 || isempty(r),      r = 1;       end
if nargin < 4 || isempty(gapdeg), gapdeg = 0.25; end
P.N      = N;
P.nmax   = nmax;
P.s      = (-1).^(0:N-1);
P.ORD    = odd_nontriplen(nmax);
P.ORDE   = even_nontriplen(nmax);
P.LIM    = grid_limits(P.ORD);
P.LIME   = 0.25 * grid_limits(P.ORDE);
P = set_r(P, r);
P = set_gap(P, gapdeg);
end

function P = set_gap(P, gapdeg)
P.gapdeg = gapdeg;
P.minGap = gapdeg * pi/180;
P.lb     = P.minGap / 2;
P.ub     = pi/2 - P.minGap / 2;
end

function P = set_r(P, r)
P.r     = r;
P.delta = (r - 1) / (r + 1);
end

function o = odd_nontriplen(nMax)
a = 5:2:nMax;  o = a(mod(a,3) ~= 0);
end

function o = even_nontriplen(nMax)
a = 2:2:nMax;  o = a(mod(a,3) ~= 0);
end

function L = grid_limits(orders)
L = zeros(size(orders));
for i = 1:numel(orders)
    switch orders(i)
        case 2,  L(i) = 0.02;
        case 4,  L(i) = 0.01;
        case 5,  L(i) = 0.06;   case 7,  L(i) = 0.05;
        case 11, L(i) = 0.035;  case 13, L(i) = 0.03;
        case 17, L(i) = 0.02;   case 19, L(i) = 0.015;
        case 23, L(i) = 0.015;  case 25, L(i) = 0.015;
        otherwise, L(i) = (0.2 + 32.5/orders(i)) / 100;
    end
end
end

function b = b1_of(TH, P)
b = (4/pi) * (cos(TH) * P.s(:));
end

function [Bo, Be, b1] = harmonics(TH, P)
[K, N] = size(TH);
s  = P.s(:);
no = P.ORD(:).';
A  = reshape(no, 1, [], 1) .* reshape(TH, K, 1, N);
Bo = reshape(reshape(cos(A), [], N) * s, K, []) .* repmat(4 ./ (no*pi), K, 1);
if P.delta == 0 || isempty(P.ORDE)
    Be = zeros(K, numel(P.ORDE));
else
    ne = P.ORDE(:).';
    Ae = reshape(ne, 1, [], 1) .* reshape(TH, K, 1, N);
    Be = -P.delta * reshape(reshape(sin(Ae), [], N) * s, K, []) .* repmat(4 ./ (ne*pi), K, 1);
end
b1 = b1_of(TH, P);
end

function t = thd_of(TH, P)
[Bo, Be, b1] = harmonics(TH, P);
t = sqrt(sum(Bo.^2, 2) + sum(Be.^2, 2)) ./ abs(b1);
end

function R = report(TH, Ma, P)
TH = sort(TH, 2);
[Bo, Be, b1] = harmonics(TH, P);
K = size(TH, 1);
if numel(Ma) == 1, Ma = Ma * ones(K, 1); end
gv = max([zeros(K,1), P.minGap - diff(TH, 1, 2)], [], 2);
bv = max(max(0, P.lb - TH(:,1)), max(0, TH(:,end) - P.ub));
R.THD49 = sqrt(sum(Bo.^2, 2) + sum(Be.^2, 2)) ./ abs(b1);
R.b1    = b1;
R.b1err = abs(b1 - Ma(:));
R.viol  = max(gv, bv);
Q = [abs(Bo) ./ repmat(P.LIM, K, 1), abs(Be) ./ repmat(P.LIME, K, 1)] ./ repmat(abs(b1), 1, numel(P.LIM) + numel(P.LIME));
R.MAXQ  = max(Q, [], 2);
R.compliant = R.MAXQ <= 1 + 1e-6;
end

function S = spectrum(TH, P, nmax)
if nargin < 3 || isempty(nmax), nmax = 999; end
TH = sort(TH, 2);
K  = size(TH, 1);
Q  = P;  Q.ORD = odd_nontriplen(nmax);  Q.ORDE = even_nontriplen(nmax);
[Bo, Be, b1] = harmonics(TH, Q);
a1 = abs(b1);
oo = Q.ORD;  oe = Q.ORDE;
i49o = oo <= 49;  i99o = oo <= 99;  i49e = oe <= 49;  i99e = oe <= 99;
sq = @(M) sum(M.^2, 2);
S.THD49    = sqrt(sq(Bo(:, i49o)) + sq(Be(:, i49e))) ./ a1;
S.THD99    = sqrt(sq(Bo(:, i99o)) + sq(Be(:, i99e))) ./ a1;
S.THDfull  = sqrt(sq(Bo) + sq(Be)) ./ a1;
wo = repmat(1 ./ oo, K, 1);  we = repmat(1 ./ oe, K, 1);
S.WTHD49   = sqrt(sq(Bo(:, i49o) .* wo(:, i49o)) + sq(Be(:, i49e) .* we(:, i49e))) ./ a1;
S.WTHDfull = sqrt(sq(Bo .* wo) + sq(Be .* we)) ./ a1;
R = abs(Bo(:, ~i49o)) ./ repmat(a1, 1, sum(~i49o));
if isempty(R)
    S.OOB_max = zeros(K,1);  S.OOB_order = nan(K,1);
else
    [S.OOB_max, j] = max(R, [], 2);
    ob = oo(~i49o);  S.OOB_order = ob(j(:)).';  S.OOB_order = S.OOB_order(:);
end
S.b3_ratio = abs((4/(3*pi)) * (cos(3*TH) * P.s(:))) ./ a1;
end

function [f, g] = objgrad(x, P, fscale)
if nargin < 3, fscale = 1; end
nfe_bump(1);
x = x(:);  s = P.s(:);
n = P.ORD(:);
Cn = cos(n * x.');
bn = (4 ./ (n*pi)) .* (Cn * s);
f  = sum(bn.^2);
if nargout > 1
    dB = -(4/pi) * sin(n * x.') .* repmat(s.', numel(n), 1);
    g  = 2 * (dB.' * bn);
end
if P.delta ~= 0 && ~isempty(P.ORDE)
    e  = P.ORDE(:);
    an = -P.delta * (4 ./ (e*pi)) .* (sin(e * x.') * s);
    f  = f + sum(an.^2);
    if nargout > 1
        dA = -P.delta * (4/pi) * cos(e * x.') .* repmat(s.', numel(e), 1);
        g  = g + 2 * (dA.' * an);
    end
end
f = fscale * f;
if nargout > 1, g = fscale * g; end
end

function [c, ceq, gc, gceq] = eqcon(x, Ma, P)
x = x(:);  s = P.s(:);
c    = [];
ceq  = (4/pi) * sum(s .* cos(x)) - Ma;
gc   = [];
gceq = -(4/pi) * (s .* sin(x));
end

function TH = random_starts(K, P)
N    = P.N;
span = (P.ub - P.lb) - (N - 1) * P.minGap;
U    = sort(rand(K, N), 2) * span;
TH   = P.lb + U + repmat(P.minGap * (0:N-1), K, 1);
end

function th = make_feasible(th, P)
th = sort(th(:)).';  N = numel(th);
th = min(max(th, P.lb), P.ub);
for k = 2:N
    th(k) = max(th(k), th(k-1) + P.minGap);
end
if th(N) > P.ub
    th(N) = P.ub;
    for k = N-1:-1:1
        th(k) = min(th(k), th(k+1) - P.minGap);
    end
end
end

function [TH, I] = solve(TH0, Ma, P, o)
if nargin < 4, o = struct(); end
o = solve_defaults(o);

[K, N] = size(TH0);
if numel(Ma) == 1, Ma = Ma * ones(K, 1); end
TH = nan(K, N);
I.ok = false(K,1);  I.thd = nan(K,1);  I.b1err = nan(K,1);  I.viol = nan(K,1);
I.iter = zeros(K,1);  I.fev = zeros(K,1);  I.flag = zeros(K,1);  I.fo = nan(K,1);

A  = zeros(N-1, N);
for k = 1:N-1, A(k,k) = 1;  A(k,k+1) = -1;  end
bb = -P.minGap * ones(N-1, 1);
lb = P.lb * ones(N, 1);  ub = P.ub * ones(N, 1);
be = backend();
lasterr_msg = '';

for k = 1:K
    x0 = make_feasible(TH0(k,:), P).';
    m  = Ma(k);
    [x, flag, it, fe, msg, fo] = run_one(be, x0, m, P, A, bb, lb, ub, o);
    if flag == -99
        Pi = P;  Pi.minGap = P.minGap + 1e-7;  Pi.lb = P.lb + 1e-7;  Pi.ub = P.ub - 1e-7;
        [x, flag, it, fe, msg, fo] = run_one(be, make_feasible(x0.', Pi).', m, P, A, bb, lb, ub, o);
    end
    if flag == -99
        lasterr_msg = msg;
        if errcount(1) >= 25
            error('shm_model:solve', ['25 consecutive solves raised an error, which points ' ...
                'to a setup problem rather than hard starts. Last error: %s'], msg);
        end
    else
        errcount(-1);
    end
    x = sort(x(:)).';
    R = report(x, m, P);
    TH(k,:)    = x;
    I.thd(k)   = R.THD49;
    I.b1err(k) = R.b1err;
    I.viol(k)  = R.viol;
    I.iter(k)  = it;
    I.fev(k)   = fe;
    I.flag(k)  = flag;
    I.fo(k)    = fo;
    conv = flag > 0 || strcmp(be, 'sqp') || (flag == 0 && fo < 1e-6);
    I.ok(k)    = conv && flag ~= -99 && R.b1err < o.b1tol && R.viol < 1e-9 && isfinite(R.THD49);
end
I.lasterr = lasterr_msg;
end

function o = solve_defaults(o)
if ~isfield(o, 'maxit'),   o.maxit   = 600;   end
if ~isfield(o, 'tol'),     o.tol     = 1e-10; end
if ~isfield(o, 'ctol'),    o.ctol    = 1e-10; end
if ~isfield(o, 'steptol'), o.steptol = 1e-12; end
if ~isfield(o, 'fscale'),  o.fscale  = 1e3;   end
if ~isfield(o, 'b1tol'),   o.b1tol   = 1e-8;  end
if ~isfield(o, 'algorithm'), o.algorithm = 'sqp'; end
end

function [x, flag, it, fe, msg, fo] = run_one(be, x0, m, P, A, bb, lb, ub, o)
msg = '';  fo = NaN;
try
    switch be
        case 'fmincon'
            opt = fmincon_options(o);
            [x, ~, flag, out] = fmincon(@(x) objgrad(x, P, o.fscale), x0, A, bb, [], [], ...
                                        lb, ub, @(x) eqcon(x, m, P), opt);
            it = out.iterations;  fe = out.funcCount;
            if isfield(out, 'firstorderopt') && ~isempty(out.firstorderopt), fo = out.firstorderopt; end
        otherwise
            phi = {@(x) objgrad(x, P, o.fscale), @(x) grad_only(x, P, o.fscale)};
            geq = {@(x) eq_only(x, m, P), @(x) eq_jac(x, P)};
            hin = {@(x) bb - A*x, @(x) -A};
            [x, ~, info, it, fe] = sqp(x0, phi, geq, hin, lb, ub, o.maxit, o.tol);
            flag = double(info == 101 || info == 104);
            if info == 103, flag = 0; end
    end
catch err
    x = x0;  flag = -99;  it = 0;  fe = 0;  msg = err.message;  fo = NaN;
end
end

function n = errcount(k)
persistent c
if isempty(c), c = 0; end
if k > 0, c = c + 1; elseif k < 0, c = 0; end
n = c;
end

function g = grad_only(x, P, fscale)
[~, g] = objgrad(x, P, fscale);
end

function v = eq_only(x, Ma, P)
v = (4/pi) * sum(P.s(:) .* cos(x(:))) - Ma;
end

function J = eq_jac(x, P)
J = (-(4/pi) * (P.s(:) .* sin(x(:)))).';
end

function opt = fmincon_options(o)
persistent cache key
k = sprintf('%g_%g_%g_%g_%s', o.maxit, o.tol, o.ctol, o.steptol, o.algorithm);
if isempty(cache) || ~strcmp(key, k)
    cache = optimoptions('fmincon', ...
        'Algorithm',                 o.algorithm, ...
        'SpecifyObjectiveGradient',  true, ...
        'SpecifyConstraintGradient', true, ...
        'Display',                   'off', ...
        'MaxIterations',             o.maxit, ...
        'MaxFunctionEvaluations',    40 * o.maxit, ...
        'OptimalityTolerance',       o.tol, ...
        'ConstraintTolerance',       o.ctol, ...
        'StepTolerance',             o.steptol);
    key = k;
end
opt = cache;
end

function be = backend()
persistent b
if isempty(b)
    isoct = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    hasf  = exist('fmincon', 'file') == 2 || exist('fmincon', 'file') == 6;
    if hasf && ~isoct
        try
            hasf = license('test', 'Optimization_Toolbox') == 1;
        catch
        end
    end
    if hasf
        b = 'fmincon';
    elseif isoct
        b = 'sqp';
    else
        error('shm_model:toolbox', ['fmincon was not found. This package uses MATLAB''s ' ...
            'Optimization Toolbox (fmincon, SQP) to generate every solution.']);
    end
end
be = b;
end

function n = nfe_get()
n = nfe_bump(0);
end

function nfe_reset()
nfe_bump(-1);
end

function n = nfe_bump(k)
persistent c
if isempty(c), c = 0; end
if k < 0, c = 0; else, c = c + k; end
n = c;
end
