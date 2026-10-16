function T = shm_sim_run(model, opts)
%SHM_SIM_RUN Simulate the cases in a Simulink model and evaluate THD49.
%   T = SHM_SIM_RUN(model, opts) runs every case of <DIR>/simulink/cases.mat:
%   writes the angles to the base workspace (VAR, and TSVAR for a From
%   Workspace block), simulates, and evaluates the logged line-to-line voltage
%   (SHM_SIM_THD) against the closed form and against the ideal waveform
%   sampled at the same instants. A case whose spectrum is not that of the
%   commanded angles is flagged (angles_used = 0). With step = true, also the
%   branch-switch cases: the angles change at TSW and the phase-current
%   deviation from the new steady state is reported in % of I1.
%   Writes <DIR>/simulink/simulink_check<tag>.mat (struct SIM) and .xlsx.
%
%   Options: DIR ('output'), var ('alpha'), tsvar ('alpha_ts'), units ('rad'),
%            signal ('vab'), isignal ('ia'), npsignal ('dVnp'), f1 (50),
%            Vdc (800), stop (0.2), ncyc (4), step (false), tsw (0.1),
%            stepstop (0.3), only ({} = all), fresh (false), tag (''), tq (0)

if nargin < 2, opts = struct(); end
O.DIR = 'output';  O.var = 'alpha';  O.tsvar = 'alpha_ts';  O.units = 'rad';
O.signal = 'vab';  O.isignal = 'ia';  O.npsignal = 'dVnp';
O.f1 = 50;  O.Vdc = 800;  O.stop = 0.2;  O.ncyc = 4;
O.step = false;  O.tsw = 0.1;  O.stepstop = 0.3;  O.only = {};  O.fresh = false;
O.tag = '';  O.tq = 0;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end
if ischar(O.only), O.only = {O.only}; end

sdir = fullfile(O.DIR, 'simulink');
L = load(fullfile(sdir, 'cases.mat'));  C = L.C;
ST = [];  if isfield(L, 'ST'), ST = L.ST; end
load_system(model);
fprintf('  model %s: %s\n', model, get_param(model, 'FileName'));
if strcmp(get_param(model, 'Dirty'), 'on')
    fprintf('  NOTE: %s has unsaved changes; sim runs the model as it is in memory, changes included\n', model);
end
check_uses(model, O, C, ST);
w = 2*pi*O.f1;  Tp = 1/O.f1;
thdo = struct('f1', O.f1, 'ncyc', O.ncyc, 'Vdc', O.Vdc);
nn = 1:49;  band = nn(mod(nn, 2) == 1 & mod(nn, 3) ~= 0 & nn > 1);

hdr = {'case', 'Ma', 'b1_of_angles', 'THD49_model_pct', 'THD49_sampled_ideal_pct', 'THD49_circuit_pct', ...
       'circuit_minus_model_pp', 'sampled_minus_model_pp', 'circuit_minus_sampled_pp', 'Ma_circuit', ...
       'THD49_incl_triplens_pct', 'WTHD49_circuit_pct', 'I1_A', 'THD49_current_pct', 'NP_mean_V', ...
       'NP_peak_to_peak_V', 'time_step_us', 'angles_used', 'spectrum_mismatch', 'what'};
rows = cell(0, numel(hdr));
wave = struct('name', {}, 't', {}, 'v', {}, 'i', {});
spec = struct('name', {}, 'n', {}, 'circuit', {}, 'model', {});

for i = 1:numel(C)
    if ~isempty(O.only) && ~any(strcmp(C(i).name, O.only)), continue; end
    a = quantise(C(i).alpha(:).', O);  s = signs(C(i));
    put_angles(O, a, a, O.stop, O.stop);
    fprintf('  %-18s Ma %.4f ... ', C(i).name, C(i).Ma);
    so = run_sim(model, O.stop, O);
    v = need(so, O.signal);
    Sv = shm_sim_thd(v, [], thdo);
    [tv, xv] = tsdata(v);
    k = tv >= Sv.window_s(1) - 1e-12;
    ts_us = 1e6 * median(diff(tv(k)));
    vi = (O.Vdc/2) * (pole(a, s, w*tv) - pole(a, s, w*tv - 2*pi/3));
    Sf = shm_sim_thd(vi, tv, thdo);
    mis = sqrt(sum((Sv.h(band) - Sf.h(band)).^2)) / max(Sf.THD49, 0.5);
    b1a = bfield(C(i));  fe = 0;
    if isfinite(b1a), fe = abs(Sv.Ma_est - b1a) / abs(b1a); end
    used = mis <= 0.5 && fe <= 0.05;
    I1 = NaN;  ithd = NaN;  iw = [];
    u = getsig(so, O.isignal);
    if ~isempty(u)
        Si = shm_sim_thd(u, [], thdo);  I1 = Si.V1;  ithd = Si.THD49;
        [ti, xi] = tsdata(u);
    end
    npm = NaN;  npp = NaN;
    u = getsig(so, O.npsignal);
    if ~isempty(u)
        [tu, xu] = tsdata(u);
        k = tu >= Sv.window_s(1);
        npm = mean(xu(k));  npp = max(xu(k)) - min(xu(k));
    end
    rows(end+1, :) = {C(i).name, C(i).Ma, bfield(C(i)), C(i).THD49_model_pct, Sf.THD49, Sv.THD49, ...
        Sv.THD49 - C(i).THD49_model_pct, Sf.THD49 - C(i).THD49_model_pct, Sv.THD49 - Sf.THD49, ...
        Sv.Ma_est, Sv.THD49_all, Sv.WTHD49, I1, ithd, npm, npp, ts_us, double(used), mis, C(i).what}; %#ok<AGROW>
    tw = (tv(end) - 2*Tp) + (0:round(2*Tp/2e-6) - 1).' * 2e-6;
    [tt, iu] = unique(tv, 'last');
    vw = interp1(tt, xv(iu), tw, 'previous');
    if exist('ti', 'var') && ~isempty(ti)
        [tt, iu] = unique(ti, 'last');  iw = interp1(tt, xi(iu), tw, 'linear');
    end
    wave(end+1) = struct('name', C(i).name, 't', single(tw - tw(1)), 'v', single(vw), 'i', single(iw)); %#ok<AGROW>
    spec(end+1) = struct('name', C(i).name, 'n', nn, 'circuit', Sv.h(nn), 'model', model_h(a, s, nn)); %#ok<AGROW>
    clear ti xi
    fprintf('analytic %.3f%%, sampled %.3f%%, circuit %.3f%%', C(i).THD49_model_pct, Sf.THD49, Sv.THD49);
    if used
        fprintf('   [angles used: spectrum mismatch %.3f, fundamental %.2f%% off]\n', mis, 100 * fe);
    else
        fprintf('\n   *** NOT the commanded angles: spectrum mismatch %.2f (limit 0.5), fundamental %.1f%% off\n', ...
            mis, 100 * fe);
        warning('shm_sim_run:angles', ['%s: the circuit''s voltage is not the waveform of the angles this script ' ...
            'wrote into %s. Check that the block feeding the MATLAB Function block reads %s, the units ' ...
            '(opts.units = ''%s''), Vdc (opts.Vdc = %g) and f1 (opts.f1 = %g). This case is kept in the ' ...
            'workbook with angles_used = 0 and left out of the paper''s table and figures.'], ...
            C(i).name, O.var, O.var, O.units, O.Vdc, O.f1);
    end
end

shdr = {'case', 'Ma_from', 'Ma_to', 'largest_angle_move_deg', 'waveform_differs_pct', 'THD49_from_model_pct', ...
        'THD49_to_model_pct', ...
        'switched', 'I1_A', 'steady_state_difference_pct', 'peak_deviation_after_switch_pct', ...
        'first_cycle_offset_pct', 'what'};
srows = cell(0, numel(shdr));
swave = struct('name', {}, 't', {}, 'v', {}, 'i', {}, 'iref', {}, 'e', {}, 'jump_deg', {});
if O.step && ~isempty(ST)
    tsw = round(O.tsw * O.f1) / O.f1;
    for j = 1:numel(ST)
        if ~isempty(O.only) && ~any(strcmp(ST(j).name, O.only)), continue; end
        a1 = quantise(ST(j).alpha1(:).', O);  a2 = quantise(ST(j).alpha2(:).', O);  s = signs(ST(j));
        put_angles(O, a1, a2, tsw, O.stepstop);
        fprintf('  %-18s Ma %.4f -> %.4f at t = %.3f s ... ', ST(j).name, ST(j).Ma1, ST(j).Ma2, tsw);
        so = run_sim(model, O.stepstop, O);
        v = need(so, O.signal);  u = need(so, O.isignal);
        Sv = shm_sim_thd(v, [], thdo);
        [tv, xv] = tsdata(v);
        S1 = shm_sim_thd((O.Vdc/2) * (pole(a1, s, w*tv) - pole(a1, s, w*tv - 2*pi/3)), tv, thdo);
        S2 = shm_sim_thd((O.Vdc/2) * (pole(a2, s, w*tv) - pole(a2, s, w*tv - 2*pi/3)), tv, thdo);
        m1 = norm(Sv.h(band) - S1.h(band));  m2 = norm(Sv.h(band) - S2.h(band));
        d12 = norm(S1.h(band) - S2.h(band));
        if d12 > 3 * min(m1, m2), sw = double(m2 < m1); else, sw = -1; end
        [ti, xi] = tsdata(u);  [ti, iu] = unique(ti, 'last');  xi = xi(iu);
        Si = shm_sim_thd(u, [], thdo);  I1 = Si.V1;
        t0 = ti(end) - Tp;
        tq = (tsw - Tp : 1e-6 : t0).';
        iq = interp1(ti, xi, tq, 'linear');
        ir = interp1(ti, xi, t0 + mod(tq - t0, Tp), 'linear');
        e  = 100 * (iq - ir) / I1;
        pre  = max(abs(e(tq < tsw)));
        post = max(abs(e(tq >= tsw)));
        off1 = abs(mean(e(tq >= tsw & tq < tsw + Tp)));
        wdp = NaN;  if isfield(ST, 'wave_diff_pct'), wdp = ST(j).wave_diff_pct; end
        srows(end+1, :) = {ST(j).name, ST(j).Ma1, ST(j).Ma2, ST(j).jump_deg, wdp, ST(j).THD1_model_pct, ...
            ST(j).THD2_model_pct, double(sw), I1, pre, post, off1, ST(j).what}; %#ok<AGROW>
        tw = (tsw - 2*Tp : 5e-6 : tsw + 3*Tp).';
        [tv, xv] = tsdata(v);  [tv, iv] = unique(tv, 'last');
        iw = interp1(ti, xi, tw, 'linear');  irw = interp1(ti, xi, t0 + mod(tw - t0, Tp), 'linear');
        swave(end+1) = struct('name', ST(j).name, 't', single(tw - tsw), ...
            'v', single(interp1(tv, xv(iv), tw, 'previous')), 'i', single(iw), 'iref', single(irw), ...
            'e', single(100 * (iw - irw) / I1), 'jump_deg', ST(j).jump_deg); %#ok<AGROW>
        if sw == 1
            fprintf('peak deviation %.2f%% of I1 (steady-state difference %.2f%%)\n', post, pre);
        elseif sw == -1
            fprintf(['peak deviation %.2f%% of I1 (steady-state difference %.2f%%); the two angle sets are too ' ...
                'similar to confirm the switch from the spectrum\n'], post, pre);
        else
            fprintf('the angles did NOT change\n');
            warning('shm_sim_run:noswitch', ['%s: the final spectrum is that of the first angle set. Feed the ' ...
                'MATLAB Function block''s alpha input from a From Workspace block (Data: %s, Interpolate ' ...
                'data: off, Form output after final data value: Holding final value).'], ST(j).name, O.tsvar);
        end
    end
end

meta = struct('model', model, 'Vdc', O.Vdc, 'f1', O.f1, 'ncyc', O.ncyc, 'stop', O.stop, ...
    'units', O.units, 'signal', O.signal, 'date', datestr(now, 'yyyy-mm-dd HH:MM'));
SIM = struct('meta', meta, 'header', {hdr}, 'data', {rows}, 'wave', wave, 'spec', spec, ...
    'step_header', {shdr}, 'step_data', {srows}, 'step_wave', swave);
f = fullfile(sdir, ['simulink_check' O.tag '.mat']);
if ~O.fresh && exist(f, 'file')
    old = load(f);
    if isfield(old, 'SIM'), SIM = merge(old.SIM, SIM); end
end
save(f, 'SIM', '-v7');
SH = struct('name', {'simulink_check'}, 'header', {SIM.header}, 'data', {SIM.data}, 'note', {sprintf( ...
    ['Circuit-level check with model %s: THD49 of the logged %s over the last %d cycles at %g Hz, Vdc = %g V, ' ...
     'against the analytic line-to-line THD49 of the same angles. sampled_ideal: the ideal waveform of the ' ...
     'same angles sampled at the instants the circuit was logged at; sampled - model is the effect of the ' ...
     'time step, circuit - sampled what the circuit adds.'], model, O.signal, O.ncyc, O.f1, O.Vdc)});
if ~isempty(SIM.step_data)
    SH(2) = struct('name', 'branch_switch', 'header', {SIM.step_header}, 'data', {SIM.step_data}, 'note', ...
        ['Angle set switched at a whole number of periods. Deviation of the phase current from its new ' ...
         'steady state (the last cycle, repeated), in % of the fundamental current amplitude. ' ...
         'steady_state_difference: the same measure in the cycle before the switch, i.e. how far apart the ' ...
         'two patterns'' steady-state currents are. switched: 1 the final spectrum is that of the new ' ...
         'angles, 0 the model did not take them, -1 the two angle sets are too similar to tell.']);
end
xf = fullfile(sdir, ['simulink_check' O.tag '.xlsx']);
try
    shm_xlsx(xf, SH);
catch err
    fprintf('  workbook not written (%s)\n', err.message);
end
fprintf('  written %s and %s\n', f, xf);
nb = 0;  if ~isempty(rows), nb = sum(cell2mat(rows(:, strcmp(hdr, 'angles_used'))) == 0); end
if nb > 0
    fprintf(2, '  *** %d of %d cases did NOT run the commanded angles; see the warnings above\n', nb, size(rows, 1));
end
if isempty(O.tag), fprintf('  next: run_paper(''full'', ''only'', {''outputs''}) adds Table VIII, Fig. 10 and Fig. 11\n'); end
T = SIM.data;
end

function so = run_sim(model, stop, O)
evalin('base', sprintf('clear %s %s %s', O.signal, O.isignal, O.npsignal));
so = sim(model, 'StopTime', num2str(stop));
end

function check_uses(model, O, C, ST)
if ~isempty(C), a = C(1).alpha(:).'; else, a = ST(1).alpha1(:).'; end
put_angles(O, a, a, 1, 1);
try
    V = Simulink.findVars(model, 'SourceType', 'base workspace');
catch err
    fprintf('  (could not list the variables %s reads: %s; the per-case spectrum check still applies)\n', ...
        model, err.message);
    return
end
users = {};
for k = 1:numel(V)
    if any(strcmp(V(k).Name, {O.var, O.tsvar})), users = [users; V(k).Users(:)]; end %#ok<AGROW>
end
if isempty(users)
    error('shm_sim_run:alpha', ['No block of %s reads %s or %s from the base workspace, so the angles this ' ...
        'script writes never reach the model and every case would simulate the same thing. Feed the ' ...
        'MATLAB Function block''s angle input from a Constant block with Value %s (or a From Workspace ' ...
        'block with Data %s), and save the model.'], model, O.var, O.tsvar, O.var, O.tsvar);
end
fprintf('  the angles reach the model through: %s\n', strjoin(unique(users(:).'), ', '));
end

function put_angles(O, a1, a2, tsw, stop)
if strcmpi(O.units, 'deg'), a1 = a1 * 180/pi;  a2 = a2 * 180/pi; end
assignin('base', O.var, a1);
if tsw >= stop, tt = [0; stop]; else, tt = [0; tsw]; end
X = struct('time', tt, 'signals', struct('values', [a1; a2], 'dimensions', numel(a1)));
assignin('base', O.tsvar, X);
end

function L = pole(a, s, x)
x = mod(x, 2*pi);
sg = ones(size(x));
k = x >= pi;  x(k) = x(k) - pi;  sg(k) = -1;
k = x > pi/2;  x(k) = pi - x(k);
lv = cumsum(s(:).');
L = zeros(size(x));
for j = 1:numel(a), L(x >= a(j)) = lv(j); end
L = sg .* L;
end

function h = model_h(a, s, n)
b  = (4 ./ (n(:) * pi)) .* (cos(n(:) * a(:).') * s(:));
b1 = (4/pi) * (cos(a(:).') * s(:));
h  = 100 * abs(b(:).') / abs(b1);
h(mod(n, 2) == 0 | mod(n, 3) == 0) = 0;
end

function a = quantise(a, O)
if O.tq > 0, d = 2*pi*O.f1*O.tq;  a = round(a/d)*d; end
end

function s = signs(c)
if isfield(c, 's') && ~isempty(c.s), s = c.s(:).'; return; end
if isfield(c, 'alpha'), N = numel(c.alpha); else, N = numel(c.alpha1); end
s = (-1) .^ (0:N-1);
end

function b = bfield(c)
if isfield(c, 'b1'), b = c.b1; else, b = NaN; end
end

function S = merge(A, B)
S = B;
same = struct('data', isfield(A, 'header') && isequal(A.header, B.header), ...
              'step_data', isfield(A, 'step_header') && isequal(A.step_header, B.step_header));
for part = {{'data', 'wave', 'spec'}, {'step_data', 'step_wave', {}}}
    p = part{1};
    if isempty(A.(p{1})), continue; end
    if ~same.(p{1})
        fprintf('  (earlier %s rows have another layout and are dropped)\n', strrep(p{1}, '_', ' '));
        continue
    end
    namesA = A.(p{1})(:, 1);  namesB = {};
    if ~isempty(B.(p{1})), namesB = B.(p{1})(:, 1); end
    keep = ~ismember(namesA, namesB);
    S.(p{1}) = [A.(p{1})(keep, :); B.(p{1})];
    wa = A.(p{2});
    if ~isempty(wa)
        kw = ~ismember({wa.name}, namesB);
        try, S.(p{2}) = [wa(kw), B.(p{2})]; catch, end
    end
    if ~isempty(p{3}) && ~isempty(A.(p{3}))
        sa = A.(p{3});  ks = ~ismember({sa.name}, namesB);
        try, S.(p{3}) = [sa(ks), B.(p{3})]; catch, end
    end
end
end

function v = need(so, name)
v = getsig(so, name);
if isempty(v)
    error('shm_sim_run:signal', ['no logged signal ''%s'' was found. Log it with a To Workspace block ' ...
        '(Variable name: %s, Save format: Timeseries), or pass its name in opts.'], name, name);
end
end

function v = getsig(so, name)
v = [];
try, v = so.get(name); catch, end
if isempty(v), try, v = so.(name); catch, end, end
if isempty(v), try, v = evalin('base', name); catch, end, end
end

function [t, x] = tsdata(v)
if isa(v, 'timeseries')
    t = v.Time;  x = squeeze(v.Data);
elseif isstruct(v) && isfield(v, 'time')
    t = v.time;  x = v.signals.values;
else
    t = (1:numel(v)).';  x = v(:);
end
t = double(t(:));  x = double(x(:));
end
