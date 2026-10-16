function R = shm_switchcheck(th, P, opts)
%SHM_SWITCHCHECK Check an angle set at the level of the switching pattern.
%   R = SHM_SWITCHCHECK(th, P, opts) builds the three-level waveform of the
%   angles th (rad) over one period, integrates its Fourier coefficients
%   exactly, computes them again by FFT of the sampled waveform, and applies a
%   dead time with the output clamped by the load-current direction.
%
%   Options: f1 (50), td_us (0), phi_deg (30), Fs (20e6), nmax (99), verbose
%   Output:  b_model, b_exact, b_fft, max_rel_dev, THD49_model, THD49_exact,
%            THD49_fft, THD49_dead, WTHD_dead, pulses_lost, narrowest_pulse_deg

if nargin < 3, opts = struct(); end
O.f1 = 50; O.td_us = 0; O.phi_deg = 30; O.Fs = 20e6; O.nmax = 99; O.verbose = true;
fn = fieldnames(opts); for i=1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

th = sort(th(:)).';  N = numel(th);
s  = P.s(:).';  s = s(1:N);
T  = 1/O.f1;  w = 2*pi*O.f1;

lvlq = cumsum(s);
ang  = [th, pi - fliplr(th)];
lvl  = [lvlq, fliplr([0 lvlq(1:end-1)])];
ang  = [0, ang, pi];
lvl  = [0, lvl, 0];
ang  = [ang, ang(2:end) + pi];
lvl  = [lvl, -lvl(2:end)];
[ang, ix] = sort(ang);  lvl = lvl(ix);
keep = [true, diff(lvl) ~= 0];
ang = ang(keep);  lvl = lvl(keep);

R.edges_deg = rad2deg(ang);  R.levels = lvl;
gaps = diff([rad2deg(ang), 360]);
if abs(ang(1)) < 1e-12 && lvl(1) == lvl(end)
    gaps = [gaps(1) + gaps(end), gaps(2:end-1)];
end
R.narrowest_pulse_deg = min(gaps(gaps > 1e-9));

ordr = P.ORD(:).';  ordr = [1, ordr(ordr <= O.nmax)];
R.b_exact = pw_fourier(ang, lvl, ordr);
R.b_model = (4./(ordr*pi)) .* (cos(ordr(:)*th) * s(:)).';
R.max_rel_dev = max(abs(R.b_exact - R.b_model)) / max(abs(R.b_model));

n  = round(O.Fs/O.f1);  t = (0:n-1)/O.Fs;
v  = level_at(ang, lvl, w*t);
R.b_fft = fft_amp(v, ordr);

td  = O.td_us*1e-6;
if td > 0
    ic   = sin(ang - deg2rad(O.phi_deg));
    ang2 = [];  lvl2 = [];
    for k = 1:numel(ang)
        if k > 1
            a = lvl(k-1); b = lvl(k);
            if ic(k) > 0, h = min(a,b); else h = max(a,b); end
            ang2(end+1) = ang(k);            lvl2(end+1) = h; %#ok<AGROW>
            ang2(end+1) = ang(k) + w*td;     lvl2(end+1) = b; %#ok<AGROW>
        else
            ang2(end+1) = ang(k);            lvl2(end+1) = lvl(k); %#ok<AGROW>
        end
    end
    [ang2, ix] = sort(ang2);  lvl2 = lvl2(ix);
    ok = ang2 < 2*pi;  ang2 = ang2(ok);  lvl2 = lvl2(ok);
    swallowed = sum(diff(ang2) <= 1e-12);
    bd = pw_fourier(ang2, lvl2, ordr);
    R.pulses_lost = sum(gaps(gaps > 1e-9) < 360*O.f1*td);
else
    bd = R.b_exact;  R.pulses_lost = 0;
end
R.b_dead = bd;

i49 = ordr <= 49;
R.THD49_model = thd(R.b_model, i49);
R.THD49_exact = thd(R.b_exact, i49);
R.THD49_fft   = thd(R.b_fft,   i49);
R.THD49_dead  = thd(R.b_dead,  i49);
wgt = 1 ./ ordr;
m49 = i49; m49(1) = false;
R.WTHD_dead   = sqrt(sum((bd(m49).*wgt(m49)).^2)) / abs(bd(1));
R.b1_model = R.b_model(1);  R.b1_dead = bd(1);

if O.verbose
    fprintf('  edges %d, narrowest pulse %.3f deg', numel(ang), R.narrowest_pulse_deg);
    if td > 0, fprintf(', dead time %.1f us = %.3f deg', O.td_us, 360*O.f1*td); end
    fprintf('\n');
    fprintf('  THD49  model %.4f%%   edge-exact %.4f%%   sampled FFT %.4f%%   (model-vs-exact dev %.2e)\n', ...
        100*R.THD49_model, 100*R.THD49_exact, 100*R.THD49_fft, R.max_rel_dev);
    if td > 0
        fprintf('  with dead time: THD49 %.4f%%  WTHD49 %.4f%%  b1 %.5f -> %.5f  pulses lost %d\n', ...
            100*R.THD49_dead, 100*R.WTHD_dead, R.b1_model, R.b1_dead, R.pulses_lost);
    end
end
end

function b = pw_fourier(ang, lvl, ordr)
a2 = [ang, 2*pi];
b  = zeros(1, numel(ordr));
for j = 1:numel(ordr)
    nq = ordr(j);
    seg = (cos(nq*a2(1:end-1)) - cos(nq*a2(2:end))) / nq;
    b(j) = sum(lvl .* seg) / pi;
end
end

function v = level_at(ang, lvl, x)
x = mod(x, 2*pi);
k = interp1(ang, 1:numel(ang), x, 'previous', 'extrap');
k(isnan(k)) = 1;
k = max(1, min(numel(ang), k));
v = lvl(k);
end

function a = fft_amp(v, ordr)
n = numel(v);  Vf = fft(v) / n;
a = 2*abs(Vf(ordr + 1));
end

function y = thd(b, i49)
m = i49;  m(1) = false;
y = sqrt(sum(b(m).^2)) / abs(b(1));
end
