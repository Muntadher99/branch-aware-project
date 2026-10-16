function NISIM = shm_sim_nonideal(model, opts)
%SHM_SIM_NONIDEAL Circuit simulations with datasheet switch models (Table 4).
%   NISIM = SHM_SIM_NONIDEAL(model, opts) runs the cases of SHM_SIM_CASES
%   through the IGBT model (NPC_SHM) under each condition of
%   SHM_NONIDEAL_CONDITIONS, by SHM_SIM_RUN, and compares each run with the
%   switching-function model SHM_NONIDEAL_WAVE at the neutral-point offset
%   measured in the circuit. The model reads Ron, Vf, Rond, Vfd, Rs, Cs, Td,
%   C1, C2, VC1_0, VC2_0, Rload, Lload and theta_deg from the base workspace.
%   Writes <DIR>/simulink/simulink_check_ni_<tag>.mat and nonideal.mat/.xlsx.
%
%   Options: DIR ('output'), R (64), L (0.18), Cnom (470e-6), Vdc (800),
%            f1 (50), stop (0.3), ncyc (10), Tvj (25), only_cond, only,
%            step (false), rebuild (false: true rebuilds nonideal.mat from
%            the per-condition files without simulating)

if nargin < 2, opts = struct(); end
O.DIR = 'output';  O.var = 'theta_deg';  O.tsvar = 'theta_ts';  O.units = 'deg';
O.signal = 'vab';  O.isignal = 'ia';  O.npsignal = 'dVnp';
O.R = 64;  O.L = 0.18;  O.Cnom = 470e-6;  O.Vdc = 800;  O.f1 = 50;
O.stop = 0.3;  O.ncyc = 10;  O.Tvj = 25;  O.only_cond = {};  O.only = {};  O.step = false;
O.rebuild = false;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end
if ischar(O.only_cond), O.only_cond = {O.only_cond}; end

sdir = fullfile(O.DIR, 'simulink');
L = load(fullfile(sdir, 'cases.mat'));  C = L.C;
dev = shm_device(O.Tvj);
K   = shm_nonideal_conditions();
cfg = struct('Vdc', O.Vdc, 'f1', O.f1, 'R', O.R, 'L', O.L);

names = {'Ron', 'Vf', 'Rond', 'Vfd', 'Tf', 'Tt', 'Rs', 'Cs', 'Td', 'C1', 'C2', 'VC1_0', 'VC2_0', 'Rload', 'Lload'};
reads = true(size(names));
f = fullfile(sdir, 'nonideal.mat');
if O.rebuild && exist(f, 'file')
    old = load(f);
    if isfield(old, 'NISIM') && isfield(old.NISIM.meta, 'reads')
        reads = ismember(names, old.NISIM.meta.reads);
    end
end
if ~O.rebuild
load_system(model);
set_vars(K(1), dev, O);
a = C(1).alpha(:).' * 180/pi;
assignin('base', O.var, a);
assignin('base', O.tsvar, struct('time', [0; 1], 'signals', struct('values', [a; a], 'dimensions', numel(a))));
try
    V = Simulink.findVars(model, 'SourceType', 'base workspace');
    reads = ismember(names, {V.Name});
    fprintf('  %s reads: %s\n', model, strjoin(names(reads), ', '));
    if any(~reads)
        fprintf(2, '  %s does NOT read: %s -- these are not part of the simulation\n', model, ...
            strjoin(names(~reads), ', '));
    end
    if ~all(ismember({'Rload', 'Lload'}, {V.Name}))
        fprintf(2, ['  the load is fixed inside %s: it must be R = %g ohm and L = %g mH per phase, or the circuit ' ...
            'and the switching-function model describe different converters (type Rload and Lload into the ' ...
            'load block to make this automatic)\n'], model, O.R, 1e3 * O.L);
    end
    if ~ismember('Td', {V.Name})
        fprintf(2, '  %s does not read Td: the dead-time conditions would simulate no dead time\n', model);
    end
catch err
    fprintf('  (could not list the variables %s reads: %s)\n', model, err.message);
end
end

hdr = {'condition', 'tag', 'case', 'Ma', 'THD49_model_pct', 'THD49_circuit_pct', 'circuit_minus_model_pp', ...
       'THD49_switching_function_pct', 'circuit_minus_switching_function_pp', 'fundamental_err_circuit_pct', ...
       'fundamental_err_switching_function_pct', 'THD49_current_pct', 'NP_peak_to_peak_V', ...
       'NP_ripple_switching_function_V', 'NP_offset_measured_pct', 'THD49_switching_function_commanded_pct', ...
       'time_step_us', 'angles_used'};
rows = cell(0, numel(hdr));
steps = {};
for k = 1:numel(K)
    if ~isempty(O.only_cond) && ~any(strcmp(K(k).tag, O.only_cond)), continue; end
    fprintf('\n  === %s ===\n', strrep(strrep(K(k).name, '$\mu$', 'u'), '\,\%', '%'));
    fk = fullfile(sdir, ['simulink_check_ni_' K(k).tag '.mat']);
    if O.rebuild && ~exist(fk, 'file'), continue; end
    if ~O.rebuild
    set_vars(K(k), dev, O);
    so = struct('DIR', O.DIR, 'var', O.var, 'tsvar', O.tsvar, 'units', O.units, 'signal', O.signal, ...
        'isignal', O.isignal, 'npsignal', O.npsignal, 'Vdc', O.Vdc, 'f1', O.f1, 'stop', O.stop, ...
        'ncyc', O.ncyc, 'tag', ['_ni_' K(k).tag], 'tq', K(k).tq, 'fresh', isempty(O.only));
    so.only = O.only;
    shm_sim_run(model, so);
    if O.step && strcmp(K(k).tag, 'all')
        so.step = true;  so.fresh = false;  so.only = {'step_jump', 'step_stay'};
        so.stepstop = O.stop + 0.1;  so.tsw = O.stop - 0.1;
        try
            shm_sim_run(model, so);
        catch err
            fprintf(2, ['  branch switch not run (%s). It needs the angles from a From Workspace block ' ...
                '(Data: %s) and the phase current logged as %s.\n'], err.message, O.tsvar, O.isignal);
        end
    end
    end
    S = load(fullfile(sdir, ['simulink_check_ni_' K(k).tag '.mat']));  S = S.SIM;
    col = @(nm) find(strcmp(S.header, nm), 1);
    NI = struct('drops', K(k).drops, 'td', K(k).td, 'eps_dc', K(k).eps_dc, 'Ctot', 2*O.Cnom, 'tq', K(k).tq);
    for r = 1:size(S.data, 1)
        i = find(strcmp({C.name}, S.data{r, 1}), 1);
        if isempty(i), continue; end
        Wc = shm_nonideal_wave(C(i).alpha, C(i).Ma, NI, cfg, dev);
        npm = S.data{r, col('NP_mean_V')};  Nm = NI;
        if isfinite(npm), Nm.eps_dc = npm / O.Vdc; end
        W = shm_nonideal_wave(C(i).alpha, C(i).Ma, Nm, cfg, dev);
        thc = S.data{r, col('THD49_circuit_pct')};
        mac = S.data{r, col('Ma_circuit')};
        rows(end+1, :) = {K(k).name, K(k).tag, C(i).name, C(i).Ma, C(i).THD49_model_pct, thc, ...
            thc - C(i).THD49_model_pct, W.THD49, thc - W.THD49, 100 * (mac - C(i).b1) / C(i).b1, ...
            100 * (W.b1 - C(i).b1) / C(i).b1, S.data{r, col('THD49_current_pct')}, ...
            S.data{r, col('NP_peak_to_peak_V')}, W.np_ripple_pp, 100 * Nm.eps_dc, Wc.THD49, S.data{r, col('time_step_us')}, ...
            S.data{r, col('angles_used')}}; %#ok<AGROW>
        fprintf('    %-16s circuit %.3f%%, switching-function model %.3f%% (difference %+.3f pp; NP offset %+.2f%%)\n', ...
            C(i).name, thc, W.THD49, thc - W.THD49, 100 * Nm.eps_dc);
    end
    if ~isempty(S.step_data)
        steps = [steps; [repmat({K(k).tag}, size(S.step_data, 1), 1), S.step_data]]; %#ok<AGROW>
        shdr = [{'tag'}, S.step_header];
    end
end

NISIM = struct('header', {hdr}, 'data', {rows}, 'step_header', {{}}, 'step_data', {steps}, ...
    'conditions', K, 'meta', struct('model', model, 'device', dev.name, 'Tvj', O.Tvj, 'R', O.R, 'L', O.L, ...
    'Cnom', O.Cnom, 'Vdc', O.Vdc, 'f1', O.f1, 'ncyc', O.ncyc, 'reads', {names(reads)}, ...
    'not_read', {names(~reads)}, 'date', datestr(now, 'yyyy-mm-dd HH:MM')));
if ~isempty(steps), NISIM.step_header = shdr; end
f = fullfile(sdir, 'nonideal.mat');
if exist(f, 'file')
    old = load(f);
    if isfield(old, 'NISIM') && isequal(old.NISIM.header, hdr) && ~isempty(old.NISIM.data)
        kold = strcat(old.NISIM.data(:, 2), '|', old.NISIM.data(:, 3));
        knew = {};  if ~isempty(rows), knew = strcat(rows(:, 2), '|', rows(:, 3)); end
        keep = ~ismember(kold, knew);
        NISIM.data = [old.NISIM.data(keep, :); rows];
        if ~isempty(old.NISIM.step_data) && isempty(steps)
            NISIM.step_data = old.NISIM.step_data;  NISIM.step_header = old.NISIM.step_header;
        end
    end
end
save(f, 'NISIM', '-v7');
SH = struct('name', {'nonideal'}, 'header', {NISIM.header}, 'data', {NISIM.data}, 'note', {sprintf( ...
    ['Circuit-level check with datasheet switches (%s), model %s, load %g ohm + %g mH per phase, C1 = C2 = %g uF. ' ...
     'THD49 of the line-to-line voltage over the last %d cycles; switching_function: shm_nonideal_wave for the ' ...
     'same condition. Parameters the model does not read: %s.'], dev.name, model, O.R, 1e3 * O.L, 1e6 * O.Cnom, ...
     O.ncyc, nonestr(names(~reads)))});
if ~isempty(NISIM.step_data)
    SH(2) = struct('name', 'nonideal_branch_switch', 'header', {NISIM.step_header}, 'data', {NISIM.step_data}, ...
        'note', 'Branch switch under all non-idealities combined; as in shm_sim_run.');
end
try
    shm_xlsx(fullfile(sdir, 'nonideal.xlsx'), SH);
catch err
    fprintf('  workbook not written (%s)\n', err.message);
end
if ~isempty(NISIM.data)
    d = cell2mat(NISIM.data(:, 9));
    fprintf('\n  %d runs; circuit minus switching-function model: largest |difference| %.3f pp\n', ...
        size(NISIM.data, 1), max(abs(d)));
end
fprintf('  written %s\n  next: run_paper(''full'', ''only'', {''outputs''}) adds the non-ideal table (tab9_nonideal.tex) and the sim.ni.* keys\n', f);
end

function set_vars(k, dev, O)
if k.drops
    v = struct('Ron', dev.rce, 'Vf', dev.Vce0, 'Rond', dev.rf, 'Vfd', dev.Vf0, 'Tf', dev.Tf, 'Tt', dev.Tt);
else
    v = struct('Ron', 1e-3, 'Vf', 0, 'Rond', 1e-3, 'Vfd', 0, 'Tf', 1e-9, 'Tt', 1e-9);
end
v.Rs = 1e5;  v.Cs = Inf;  v.Td = k.td;
v.C1 = O.Cnom * (1 - k.Cmis);  v.C2 = O.Cnom * (1 + k.Cmis);
v.VC1_0 = O.Vdc/2 * (1 + k.eps_dc);  v.VC2_0 = O.Vdc/2 * (1 - k.eps_dc);
v.Rload = O.R;  v.Lload = O.L;
f = fieldnames(v);
for i = 1:numel(f), assignin('base', f{i}, v.(f{i})); end
end

function s = nonestr(c)
if isempty(c), s = 'none'; else, s = strjoin(c, ', '); end
end
