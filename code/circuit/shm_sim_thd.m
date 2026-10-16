function S = shm_sim_thd(v, t, opts)
%SHM_SIM_THD THD49 of a simulated line-to-line voltage.
%   S = SHM_SIM_THD(vab) or S = SHM_SIM_THD(v, t, opts). The last NCYC cycles
%   are resampled with a zero-order hold and transformed by an FFT over exactly
%   NCYC cycles. Output: V1, Ma_est, THD49 (orders 2..49 without triplens, %),
%   THD49_all, WTHD49, h (orders 1..99, % of V1), window_s.
%
%   Options: f1 (50), ncyc (4), Vdc (800), fs (1e7)

if nargin < 2, t = []; end
if nargin < 3, opts = struct(); end
O.f1 = 50;  O.ncyc = 4;  O.Vdc = 800;  O.fs = 1e7;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

if isa(v, 'timeseries')
    t = v.Time;  v = v.Data;
elseif isstruct(v) && isfield(v, 'time')
    t = v.time;  v = v.signals.values;
end
v = double(squeeze(v));  t = double(t(:));
if size(v, 1) ~= numel(t), v = v.'; end
v = v(:, 1);

T  = 1 / O.f1;
t1 = t(end);  t0 = t1 - O.ncyc * T;
if t0 < t(1)
    error('shm_sim_thd:short', 'the signal is %.3f s long; %d cycles need %.3f s', t1 - t(1), O.ncyc, O.ncyc * T);
end
M  = round(O.ncyc * T * O.fs);
tu = t0 + (0:M-1).' * (O.ncyc * T / M);
[tt, ia] = unique(t, 'last');
vu = interp1(tt, v(ia), tu, 'previous');
vu(isnan(vu)) = v(1);
X  = fft(vu) / M;
n  = 1:99;
amp = 2 * abs(X(O.ncyc * n + 1));
amp = amp(:).';
V1 = amp(1);
band = 2:49;
nt = band(mod(band, 3) ~= 0);
S.V1        = V1;
S.Ma_est    = V1 / (sqrt(3) * O.Vdc / 2);
S.THD49     = 100 * sqrt(sum(amp(nt).^2)) / V1;
S.THD49_all = 100 * sqrt(sum(amp(band).^2)) / V1;
S.WTHD49    = 100 * sqrt(sum((amp(nt) ./ nt).^2)) / V1;
S.h         = 100 * amp / V1;
S.window_s  = [t0 t1];
end
