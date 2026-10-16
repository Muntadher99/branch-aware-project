function reproduce_circuit(varargin)
%REPRODUCE_CIRCUIT Rerun the circuit simulations of Sections 6.4 to 6.6.
%   REPRODUCE_CIRCUIT runs three sets of Simulink simulations into
%   output/simulink and then rewrites the figures, tables and numbers:
%     'ideal'        models/test.slx, ideal switches, 64 ohm + 180 mH per phase:
%                    the nine steady-state cases, the branch switches of
%                    Fig. 7, and the branch switch at every switch point
%     'datasheet'    models/NPC_SHM.slx, datasheet IGBT and diode models: the
%                    nine cases under the six conditions of Table 4
%     'second_load'  models/test_80ohm_25mH.slx, the 'ideal' runs with a load of
%                    80 ohm + 25 mH per phase
%   The angle sets come from output/paper_results.mat: run REPRODUCE_PAPER
%   first. About one hour in all on the reference PC.
%
%   REPRODUCE_CIRCUIT('runs', {'ideal'}) runs only the named sets.
%   REPRODUCE_CIRCUIT(..., 'out', folder) uses another output folder.
%
%   Requires Simulink and Simscape Electrical (Specialized Power Systems).
%
%   See also REPRODUCE_PAPER, SHM_SIM_RUN, SHM_SIM_SWITCHES, SHM_SIM_NONIDEAL.

A.runs = {'ideal', 'datasheet', 'second_load'};  A.out = '';
for i = 1:2:numel(varargin), A.(lower(varargin{i})) = varargin{i+1}; end
if ischar(A.runs), A.runs = {A.runs}; end
bad = setdiff(A.runs, {'ideal', 'datasheet', 'second_load'});
if ~isempty(bad), error('reproduce_circuit:runs', 'unknown run: %s', strjoin(bad, ', ')); end

root = setup_paths();
out = A.out;
if isempty(out), out = fullfile(root, 'output'); end
f = fullfile(out, 'paper_results.mat');
if ~exist(f, 'file')
    error('reproduce_circuit:noresults', '%s not found: run reproduce_paper first', f);
end
L = load(f);  R = L.R;
sd = fullfile(out, 'simulink');
diary(fullfile(out, 'circuit_log.txt'));  diary on;
cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
t0 = tic;
fprintf('\n=== reproduce_circuit %s ===\n', datestr(now, 'yyyy-mm-dd HH:MM'));

shm_sim_cases(R, struct('DIR', out));
if any(strcmp('ideal', A.runs))
    fprintf('\n--- ideal switches, model test\n');
    shm_sim_run('test', struct('DIR', out, 'step', true, 'fresh', true));
    shm_sim_switches('test', [], struct('DIR', out, 'fresh', true));
end
if any(strcmp('datasheet', A.runs))
    fprintf('\n--- datasheet switch models, model NPC_SHM\n');
    shm_sim_nonideal('NPC_SHM', struct('DIR', out));
end
if any(strcmp('second_load', A.runs))
    fprintf('\n--- second load, model test_80ohm_25mH\n');
    second_load(R, out, 'test_80ohm_25mH', '_80ohm_25mH');
end
two_loads(sd);
fprintf('\n=== circuit simulations finished in %.1f min ===\n', toc(t0) / 60);
diary off;

run_paper('full', 'dir', out, 'only', {'outputs'});
shm_numbers_supplement([], out);
end

function second_load(R, out, model, tag)
sd = fullfile(out, 'simulink');
tmp = fullfile(sd, ['run' tag]);
ts = fullfile(tmp, 'simulink');
if exist(tmp, 'dir'), rmdir(tmp, 's'); end
mkdir(ts);
copyfile(fullfile(sd, 'cases.mat'), fullfile(ts, 'cases.mat'));
shm_sim_run(model, struct('DIR', tmp, 'step', true, 'fresh', true));
shm_sim_switches(model, R, struct('DIR', tmp, 'fresh', true));
for f = {'simulink_check.mat', 'simulink_check.xlsx', 'switches.mat', 'switches.xlsx'}
    [~, nm, ex] = fileparts(f{1});
    if exist(fullfile(ts, f{1}), 'file'), copyfile(fullfile(ts, f{1}), fullfile(sd, [nm tag ex])); end
end
d = fullfile(sd, ['switches' tag]);
if exist(d, 'dir'), rmdir(d, 's'); end
copyfile(fullfile(ts, 'switches'), d);
rmdir(tmp, 's');
end

function two_loads(sd)
q = dir(fullfile(sd, 'switches_*ohm_*mH.mat'));
f = fullfile(sd, 'switches.mat');
if isempty(q) || ~exist(f, 'file'), return; end
S = {load(fullfile(sd, q(1).name)), load(f)};
S = cellfun(@(x) x.SW, S, 'UniformOutput', false);
out = fullfile(sd, 'switches_two_loads.txt');
fid = fopen(out, 'w');
for k = 1:2
    W = S{k};  h = W.header;  c = @(nm) find(strcmp(h, nm), 1);
    fprintf(fid, 'load %g ohm + %g mH per phase (PF %.3f), model %s\n', W.meta.Rload, 1e3 * W.meta.Lload, ...
        W.meta.PF, W.meta.model);
    fprintf(fid, '%-10s %8s %6s %6s %12s %12s %14s %14s\n', 'case', 'Ma', 'from', 'to', 'circ_peak', ...
        'exact_start', 'circ_ss_diff', 'exact_ss_diff');
    d = zeros(size(W.data, 1), 1);
    for i = 1:size(W.data, 1)
        r = W.data(i, :);
        fprintf(fid, '%-10s %8.4f %6d %6d %12.3f %12.3f %14.3f %14.3f\n', r{1}, r{c('Ma')}, r{c('from')}, ...
            r{c('to')}, r{c('circuit_peak_deviation_pct')}, r{c('exact_deviation_at_period_start_pct')}, ...
            r{c('circuit_steady_state_difference_pct')}, r{c('exact_steady_state_difference_pct')});
        d(i) = abs(r{c('circuit_peak_deviation_pct')} - r{c('exact_deviation_at_period_start_pct')});
    end
    fprintf(fid, 'largest |circuit peak - exact|: %.3f pp of I1 over %d switch points\n\n', max(d), numel(d));
end
fclose(fid);
fprintf('\nevery switch point at both loads: %s\n', out);
type(out);
end
