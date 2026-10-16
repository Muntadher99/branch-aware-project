function SW = shm_sim_switches(model, R, opts)
%SHM_SIM_SWITCHES Circuit simulation of the branch switch at every switch point.
%   SW = SHM_SIM_SWITCHES(model, R, opts) switches, at every owner change of
%   the stored subset, from the old owner's angles to the new owner's (both
%   solved at the switch point) at the start of a period, simulates it with
%   SHM_SIM_RUN, and compares the circuit's current deviation with the exact
%   R-L result of SHM_RL_CURRENT for the model's load. R is the results struct
%   (loaded from <DIR>/paper_results.mat if empty). Writes
%   <DIR>/simulink/switches.mat (struct SW), switches.xlsx and switches/.
%
%   Options: DIR ('output'), Rload, Lload (read from the model if possible),
%            tsw (0.2), stepstop (0.4), ncyc (4), fresh (false)

if nargin < 2, R = []; end
if nargin < 3, opts = struct(); end
O.DIR = 'output';  O.var = 'alpha';  O.tsvar = 'alpha_ts';  O.units = 'rad';
O.signal = 'vab';  O.isignal = 'ia';  O.npsignal = 'dVnp';  O.Vdc = 800;  O.f1 = 50;
O.Rload = 80;  O.Lload = 0.025;  O.tsw = 0.2;  O.stepstop = 0.4;  O.ncyc = 4;  O.fresh = false;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end
if isempty(R)
    L = load(fullfile(O.DIR, 'paper_results.mat'));  R = L.R;
end
D = R.D;  F = shm_model();  P = D.P;  grid = D.grid(:);

own = D.owns(:);
ch = find(own(1:end-1) ~= own(2:end) & isfinite(own(1:end-1)) & isfinite(own(2:end)));
ST = struct('name', {}, 'Ma1', {}, 'Ma2', {}, 'alpha1', {}, 'alpha2', {}, 'alpha1_deg', {}, 'alpha2_deg', {}, ...
            's', {}, 'jump_deg', {}, 'wave_diff_pct', {}, 'THD1_model_pct', {}, 'THD2_model_pct', {}, 'what', {});
for g = ch(:).'
    a = own(g);  b = own(g + 1);
    [xc, thA, thB] = at_switch(F, P, D, g, a, b);
    if isempty(thA), continue; end
    ST(end+1) = struct('name', sprintf('switch_%02d', numel(ST) + 1), 'Ma1', xc, 'Ma2', xc, 'alpha1', thA, ...
        'alpha2', thB, 'alpha1_deg', thA * 180/pi, 'alpha2_deg', thB * 180/pi, 's', P.s(1:numel(thA)), ...
        'jump_deg', max(abs(thB - thA)) * 180/pi, 'wave_diff_pct', shm_wavediff(thA, thB, P.s), ...
        'THD1_model_pct', 100 * F.thd(thA, P), 'THD2_model_pct', 100 * F.thd(thB, P), ...
        'what', sprintf('table switch point: branch %d -> branch %d at Ma %.5f', a, b, xc)); %#ok<AGROW>
end
C = struct('name', {}, 'Ma', {}, 'alpha', {}, 'alpha_deg', {}, 's', {}, 'b1', {}, 'THD49_model_pct', {}, 'what', {});
sub = fullfile(O.DIR, 'simulink', 'switches');
if ~exist(fullfile(sub, 'simulink'), 'dir'), mkdir(fullfile(sub, 'simulink')); end
save(fullfile(sub, 'simulink', 'cases.mat'), 'C', 'ST', '-v7');
fprintf('  %d switch points of the branch-aware table\n', numel(ST));

f = fullfile(sub, 'simulink', 'simulink_check.mat');
done = {};
if ~O.fresh && exist(f, 'file')
    S0 = load(f);  if ~isempty(S0.SIM.step_data), done = S0.SIM.step_data(:, 1); end
end
todo = setdiff({ST.name}, done);
if ~isempty(todo)
    so = struct('DIR', sub, 'var', O.var, 'tsvar', O.tsvar, 'units', O.units, 'signal', O.signal, ...
        'isignal', O.isignal, 'npsignal', O.npsignal, 'Vdc', O.Vdc, 'f1', O.f1, 'ncyc', O.ncyc, ...
        'step', true, 'tsw', O.tsw, 'stepstop', O.stepstop, 'fresh', O.fresh && isempty(done));
    so.only = todo;
    shm_sim_run(model, so);
end
S = load(f);  S = S.SIM;

[O.Rload, O.Lload] = model_load(model, O.Rload, O.Lload, ~isfield(opts, 'Rload'));
pf = O.Rload / hypot(O.Rload, 2*pi*O.f1*O.Lload);
xs = (0:3599).' * (2*pi / 3600);
col = @(nm) find(strcmp(S.step_header, nm), 1);
hdr = {'case', 'Ma', 'from', 'to', 'largest_angle_move_deg', 'waveform_differs_pct', 'switched', ...
       'circuit_steady_state_difference_pct', 'exact_steady_state_difference_pct', ...
       'circuit_peak_deviation_pct', 'exact_deviation_at_period_start_pct', 'what'};
rows = cell(0, numel(hdr));
for j = 1:numel(ST)
    r = find(strcmp(S.step_data(:, 1), ST(j).name), 1);
    if isempty(r), continue; end
    iA = shm_rl_current(ST(j).alpha1, P.s, pf, 0, xs);  iB = shm_rl_current(ST(j).alpha2, P.s, pf, 0, xs);
    I1 = F.b1(ST(j).alpha2, P);
    d = 100 * abs(iA(:, 1) - iB(:, 1)) / I1;
    tok = sscanf(ST(j).what, 'table switch point: branch %d -> branch %d');
    rows(end+1, :) = {ST(j).name, ST(j).Ma1, tok(1), tok(2), ST(j).jump_deg, ST(j).wave_diff_pct, ...
        S.step_data{r, col('switched')}, S.step_data{r, col('steady_state_difference_pct')}, max(d), ...
        S.step_data{r, col('peak_deviation_after_switch_pct')}, d(1), ST(j).what}; %#ok<AGROW>
end
SW = struct('header', {hdr}, 'data', {rows}, 'meta', struct('model', model, 'Rload', O.Rload, 'Lload', O.Lload, ...
    'PF', pf, 'Vdc', O.Vdc, 'f1', O.f1, 'date', datestr(now, 'yyyy-mm-dd HH:MM')));
save(fullfile(O.DIR, 'simulink', 'switches.mat'), 'SW', '-v7');
try
    shm_xlsx(fullfile(O.DIR, 'simulink', 'switches.xlsx'), struct('name', 'branch_switches', 'header', {hdr}, ...
        'data', {rows}, 'note', sprintf(['Branch switch at every switch point of the branch-aware table, model %s, ' ...
        'load %g ohm + %g mH (PF %.3f); circuit vs exact R-L result (shm_rl_current), phase a, %% of I1.'], ...
        model, O.Rload, 1e3 * O.Lload, pf)));
catch err
    fprintf('  workbook not written (%s)\n', err.message);
end
if ~isempty(rows)
    c = cell2mat(rows(:, 10));  e = cell2mat(rows(:, 11));
    fprintf('  %d switches run: circuit peak deviation %.2f..%.2f%% of I1; exact result %.2f..%.2f%%; largest |difference| %.2f pp\n', ...
        size(rows, 1), min(c), max(c), min(e), max(e), max(abs(c - e)));
end
end

function [R, L] = model_load(model, R, L, look)
if ~look, return; end
try
    b = find_system(model, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'MaskType', 'Series RLC Branch');
    v = [];
    for i = 1:numel(b)
        if strcmp(get_param(b{i}, 'BranchType'), 'RL')
            v(end+1, :) = [evalin('base', get_param(b{i}, 'Resistance')), ...
                           evalin('base', get_param(b{i}, 'Inductance'))]; %#ok<AGROW>
        end
    end
    if size(v, 1) == 3 && all(all(abs(v - v(1, :)) <= 1e-12 * abs(v(1, :))))
        R = v(1, 1);  L = v(1, 2);
        fprintf('  load read from %s: %g ohm + %g mH per phase\n', model, R, 1e3 * L);
    else
        fprintf('  the load of %s could not be read: %g ohm + %g mH is used (opts.Rload, opts.Lload)\n', model, R, 1e3 * L);
    end
catch
    fprintf('  the load of %s could not be read: %g ohm + %g mH is used (opts.Rload, opts.Lload)\n', model, R, 1e3 * L);
end
end

function [xc, thA, thB] = at_switch(F, P, D, g, a, b)
grid = D.grid(:);
da = D.THD(g, a) - D.THD(g, b);  db = D.THD(g + 1, a) - D.THD(g + 1, b);
if ~isfinite(D.THD(g + 1, a)), xc = grid(g);
elseif ~isfinite(D.THD(g, b)), xc = grid(g + 1);
elseif all(isfinite([da db])) && da ~= db, xc = grid(g) + da / (da - db) * (grid(g + 1) - grid(g));
else, xc = grid(g + 1); end
thA = one(F, P, D, a, g, xc);  thB = one(F, P, D, b, g, xc);
if isempty(thA) || isempty(thB), thA = []; thB = []; end
end

function th = one(F, P, D, k, g, x)
grid = D.grid(:);  th = [];
f = [g g+1];  f = f(isfinite(D.THD(f, k)));
if isempty(f), return; end
if numel(f) == 2
    x0 = interp1(grid(f), reshape(D.TH(f, :, k), 2, P.N), x, 'linear');
else
    x0 = D.TH(f, :, k);
end
[y, I] = F.solve(x0, x, P);
if I.ok && max(abs(y - sort(x0))) <= pi/180, th = y; else, th = sort(x0); end
end
