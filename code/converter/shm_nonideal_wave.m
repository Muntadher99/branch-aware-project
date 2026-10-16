function W = shm_nonideal_wave(th, Ma, NI, cfg, dev)
%SHM_NONIDEAL_WAVE Switching-function model of the three-phase NPC converter.
%   W = SHM_NONIDEAL_WAVE(th, Ma, NI, cfg, dev) evaluates one quarter-wave
%   angle set (rad) under the non-idealities NI: conduction drops (drops),
%   dead time td (s), static neutral-point offset eps_dc, neutral-point ripple
%   from the capacitance Ctot (F), timer resolution tq (s). The phase currents
%   are the exact periodic steady state of the R-L load, iterated with the
%   voltage (iters, 4).
%
%   cfg: Vdc (800), f1 (50), R (64), L (0.18), M (2^17 samples per period)
%   dev: device data (default SHM_DEVICE())
%   Output: THD49 (%), THD49_odd, b1, b1_err_pct, even_max, np_ripple_pp,
%           dthd (pp against the closed form), I1_pk, ITHD49

if nargin < 3 || isempty(NI),  NI = struct(); end
if nargin < 4 || isempty(cfg), cfg = struct(); end
if nargin < 5 || isempty(dev), dev = shm_device(); end
cfg = fillc(cfg, 'Vdc', 800);  cfg = fillc(cfg, 'f1', 50);
cfg = fillc(cfg, 'R', 64);     cfg = fillc(cfg, 'L', 0.18);
cfg = fillc(cfg, 'M', 2^17);
NI  = fillc(NI, 'drops', false);  NI = fillc(NI, 'td', 0);
NI  = fillc(NI, 'eps_dc', 0);     NI = fillc(NI, 'Ctot', Inf);
NI  = fillc(NI, 'tq', 0);         NI = fillc(NI, 'iters', 4);

M = cfg.M;  f1 = cfg.f1;  dt = 1/(f1*M);  w1 = 2*pi*f1;  E2 = cfg.Vdc/2;
ang = sort(th(:)).';
if NI.tq > 0
    d = w1*NI.tq;  ang = round(ang/d)*d;
end
x  = 2*pi*(0:M-1).'/M;
S0 = [state(x, ang), state(mod(x - 2*pi/3, 2*pi), ang), state(mod(x - 4*pi/3, 2*pi), ang)];
nd = round(NI.td/dt);
h  = [0:M/2, -(M/2-1):-1].';
Z  = cfg.R + 1i*h*w1*cfg.L;

S  = S0;  v1 = E2*(1 + NI.eps_dc)*ones(M, 1);  v2 = cfg.Vdc - v1;
V0 = (S0 == 1).*v1 - (S0 == -1).*v2;
I  = zeros(M, 3);
for p = 1:3, I(:, p) = current(V0, Z, p); end
for it = 1:max(1, NI.iters)
    for p = 1:3
        S(:, p) = apply_deadtime(S0(:, p), I(:, p), nd);
    end
    for p = 1:3
        s = S(:, p);  i = I(:, p);  ai = abs(i);
        v = (s == 1).*v1 - (s == -1).*v2;
        if NI.drops
            vT = dev.Vce0 + dev.rce*ai;  vD = dev.Vf0 + dev.rf*ai;
            pos = i > 0;
            drop = (s == 1).*(pos.*2.*vT + ~pos.*2.*vD) + ...
                   (s == 0).*(vT + vD) + ...
                   (s == -1).*(~pos.*2.*vT + pos.*2.*vD);
            v = v - sign(i).*drop;
        end
        V0(:, p) = v;
    end
    for p = 1:3, I(:, p) = current(V0, Z, p); end
    if isfinite(NI.Ctot)
        inp = sum(I .* (S == 0), 2);
        inp = inp - mean(inp);
        rip = cumsum(inp)*dt / NI.Ctot;  rip = rip - mean(rip);
        v1 = E2*(1 + NI.eps_dc) + rip;  v2 = cfg.Vdc - v1;
    end
end

vab = V0(:, 1) - V0(:, 2);
Vh  = abs(fft(vab)) / (M/2);  V1 = Vh(2);
n   = 2:49;  nt = n(mod(n, 3) ~= 0);  no = nt(mod(nt, 2) == 1);
W.THD49     = 100*sqrt(sum(Vh(nt + 1).^2))/V1;
W.THD49_odd = 100*sqrt(sum(Vh(no + 1).^2))/V1;
ev = 2:2:50;  W.even_max = 100*max(Vh(ev + 1))/V1;
W.b1 = V1/(sqrt(3)*E2);
W.b1_err_pct = 100*(W.b1 - Ma)/Ma;
W.np_ripple_pp = max(v1) - min(v1);
F = shm_model();  P = F.defaults(numel(ang), 49);
W.dthd = W.THD49 - 100*F.thd(sort(th(:)).', P);
Ih = abs(fft(I(:, 1)))/(M/2);
W.I1_pk = Ih(2);
W.ITHD49 = 100*sqrt(sum(Ih(nt + 1).^2))/Ih(2);
end

function s = state(x, ang)
q = mod(x, 2*pi);  sgn = ones(size(q));  sgn(q >= pi) = -1;
q = mod(q, pi);  q(q > pi/2) = pi - q(q > pi/2);
s = sgn .* mod(sum(bsxfun(@ge, q, ang), 2), 2);
end

function i = current(V0, Z, p)
vN = V0(:, p) - mean(V0, 2);
Ih = fft(vN) ./ Z;  Ih(1) = 0;
i  = real(ifft(Ih));
end

function s = apply_deadtime(s0, i, nd)
s = s0;
if nd <= 0, return; end
M = numel(s0);
prev = s0([end 1:end-1]);
k = find(s0 ~= prev);
for j = k(:).'
    a = prev(j);  b = s0(j);  ij = i(j);
    need = (a == 0 && b == 1 && ij > 0) || (a == 1 && b == 0 && ij < 0) || ...
           (a == 0 && b == -1 && ij < 0) || (a == -1 && b == 0 && ij > 0);
    if need
        idx = mod((j:j+nd-1) - 1, M) + 1;
        s(idx) = a;
    end
end
end

function c = fillc(c, f, v)
if ~isfield(c, f) || isempty(c.(f)), c.(f) = v; end
end
