function R = run_paper(mode, varargin)
%RUN_PAPER Run the computational stages of the study and write all outputs.
%   R = RUN_PAPER(MODE, NAME, VALUE, ...) runs every stage that is not yet in
%   <dir>/cache, loads the others, and writes figures, tables and number files
%   into <dir>. MODE is 'full' (the paper's settings, default) or 'quick'
%   (small grids, for testing the chain). Normally called by REPRODUCE_PAPER.
%
%   Name-value options:
%     'dir'        output folder (default: ./output, or ./output_quick)
%     'force'      recompute every stage (default false)
%     'redo'       cell array of stages to recompute; stages that depend on a
%                  recomputed stage are recomputed as well
%     'only'       cell array of stages to run; the others are loaded from the
%                  cache. {'outputs'} writes figures and tables only.
%     'selection'  include stage E4 (default true)
%     'revision'   include stages IND, CMP, STO, RGS, BSW, TIM (default true)
%
%   Stages (cache file <dir>/cache/<stage>.mat; random stream rng(1000 + k)):
%     D    branch set: multistart, continuation, probes, fold ends, subset
%     S    stored branches extended over (Ma, delta)
%     V    stored branches extended over (Ma, delta, theta_gap)
%     E1   multiplicity and per-branch structure
%     E1b  interior ends of the owning branches
%     E2   lookup table on the grid, corrector, continuity-constrained path
%     E2b  table comparison at the 200 midpoints (Table 1, Figs. 4-5)
%     E3   regression inside a branch, d = 1 (Table 2)
%     E3b  regression inside a branch, d = 2 (Table 2)
%     E3c  regression inside a branch, d = 3 (Table 2)
%     E4   branch selection, d = 1, 2, 3 (Section 4.4)
%     E5   storage, scaling and on-line time
%     E6   independent multistart, spectrum above the band, limits
%     E7   gate pattern, dead time and pulse loss
%     E9   tables under non-ideal switching (Table 4)
%     IND  independent search: interior-point, Halton starts (Section 6.1)
%     CMP  per-order compliance (Section 6.2)
%     STO  Table 1 at equal storage, with published baselines
%     RGS  regression statistics: seeds, model selection, ranking
%     BSW  branch switching under R-L loads, exact (Table 3)
%     TIM  operation count and time of a table read (Section 6.5)
%
%   See also REPRODUCE_PAPER, SHM_FIGURES, SHM_TABLES.

if nargin < 1 || isempty(mode), mode = 'full'; end
A.force = false;  A.redo = {};  A.only = {};  A.dir = '';  A.selection = true;  A.revision = true;
for i = 1:2:numel(varargin), A.(lower(varargin{i})) = varargin{i+1}; end
if ischar(A.redo), A.redo = {A.redo}; end
if ischar(A.only), A.only = {A.only}; end
mode = lower(mode);
if ~any(strcmp(mode, {'full', 'quick'}))
    error('run_paper:mode', 'mode must be ''full'' or ''quick''');
end

here = fileparts(mfilename('fullpath'));
addpath(here);
if isempty(A.dir)
    if strcmp(mode, 'quick'), A.dir = fullfile(pwd, 'output_quick'); else, A.dir = fullfile(pwd, 'output'); end
end
cdir = fullfile(A.dir, 'cache');
if ~exist(cdir, 'dir'), mkdir(cdir); end
diary(fullfile(A.dir, 'run_log.txt'));  diary on;
cleanup = onCleanup(@() diary('off')); %#ok<NASGU>

C  = run_config(mode);
F  = shm_model();
t0 = tic;
fprintf('\n=== run_paper (%s) %s ===\n', mode, datestr(now, 'yyyy-mm-dd HH:MM'));
fprintf('output: %s\n', A.dir);
fprintf('solver backend: %s\n', F.backend());
if strcmp(F.backend(), 'sqp')
    fprintf(['NOTE: fmincon was not found, Octave''s sqp is used instead. The numbers in the paper ' ...
             'must come from a MATLAB run with the Optimization Toolbox.\n']);
end

stages = {'D', 'S', 'V', 'E1', 'E1b', 'E2', 'E2b', 'E3', 'E3b', 'E3c', 'E4', 'E5', 'E6', 'E7', 'E9', ...
          'IND', 'CMP', 'STO', 'RGS', 'BSW', 'TIM'};
if ~A.selection, stages = stages(~strcmp(stages, 'E4')); end
if ~A.revision, stages = stages(~ismember(stages, {'IND', 'CMP', 'STO', 'RGS', 'BSW', 'TIM'})); end
deps = struct('D', {{}}, 'S', {{'D'}}, 'V', {{'D'}}, 'E1', {{'D'}}, 'E1b', {{'D'}}, 'E2', {{'D'}}, ...
    'E2b', {{'D', 'E2'}}, 'E3', {{'D'}}, 'E3b', {{'D', 'S'}}, 'E3c', {{'D', 'V'}}, 'E4', {{'D', 'S', 'V'}}, ...
    'E5', {{'D'}}, 'E6', {{'D'}}, 'E7', {{'D'}}, 'E9', {{'D', 'E2b'}}, ...
    'IND', {{'D', 'E2b'}}, 'CMP', {{'D', 'E2b'}}, 'STO', {{'D', 'E2', 'E2b'}}, 'RGS', {{'D', 'S', 'V'}}, ...
    'BSW', {{'D'}}, 'TIM', {{'D', 'E2b'}});

seedno = struct('D', 1, 'S', 2, 'V', 3, 'E1', 4, 'E2', 5, 'E3', 6, 'E3b', 7, 'E3c', 8, 'E4', 9, ...
                'E5', 10, 'E6', 11, 'E7', 12, 'E9', 13, 'E1b', 14, 'E2b', 15, ...
                'IND', 16, 'CMP', 17, 'STO', 18, 'RGS', 19, 'BSW', 20, 'TIM', 21);
R = struct();
status = struct();  fresh = {};  errs = {};
for i = 1:numel(stages)
    st = stages{i};
    status.(st) = 'skipped';
    if ~isempty(A.only) && ~any(strcmp(st, A.only))
        X = load_stage(cdir, st, stage_opts(C, st));
        if ~isempty(X), R.(st) = X; status.(st) = 'loaded'; end
        continue
    end
    missing = deps.(st);  missing = missing(~cellfun(@(d) isfield(R, d), missing));
    if ~isempty(missing)
        fprintf('\n--- %s skipped: needs %s\n', st, strjoin(missing, ', '));
        continue
    end
    recompute = A.force || any(strcmp(st, A.redo)) || any(ismember(deps.(st), fresh));
    if ~recompute
        X = load_stage(cdir, st, stage_opts(C, st));
        if ~isempty(X)
            R.(st) = X;  status.(st) = 'cached';
            fprintf('\n--- %s: loaded from cache\n', st);
            continue
        end
    end
    fprintf('\n--- %s: %s\n', st, describe(st));
    ts = tic;
    try
        rng(1000 + seedno.(st));
        X = run_stage(st, R, C, F, seedno);
        R.(st) = X;
        save_stage(cdir, st, X, stage_opts(C, st));
        fresh{end+1} = st; %#ok<AGROW>
        status.(st) = sprintf('done in %.1f min', toc(ts) / 60);
        fprintf('--- %s done in %.1f min\n', st, toc(ts) / 60);
    catch err
        status.(st) = 'FAILED';
        errs{end+1} = sprintf('%s: %s', st, err.message); %#ok<AGROW>
        report_error(st, err);
    end
end

R.meta = run_meta(mode, C, F, toc(t0), status, errs);
R.meta.seeds = structfun(@(k) 1000 + k, seedno, 'UniformOutput', false);
sf = fullfile(A.dir, 'simulink', 'simulink_check.mat');
if exist(sf, 'file')
    try
        L = load(sf);  R.SIM = L.SIM;  R.SIM.data_all = R.SIM.data;
        j = find(strcmp(R.SIM.header, 'angles_used'), 1);
        if isempty(j)
            fprintf('\n%s was written by an older shm_sim_run that did not check the angles: not used\n', sf);
            R = rmfield(R, 'SIM');
        else
            ok = cell2mat(R.SIM.data(:, j)) == 1;
            R.SIM.data = R.SIM.data(ok, :);  R.SIM.n_invalid = sum(~ok);
            good = R.SIM.data(:, 1);
            if ~isempty(R.SIM.wave), R.SIM.wave = R.SIM.wave(ismember({R.SIM.wave.name}, good)); end
            if ~isempty(R.SIM.spec), R.SIM.spec = R.SIM.spec(ismember({R.SIM.spec.name}, good)); end
            fprintf('\ncircuit-level results found: %d cases, %d left out because the circuit did not run the commanded angles\n', ...
                sum(ok), sum(~ok));
        end
    catch err
        fprintf('\ncould not read %s: %s\n', sf, err.message);
    end
end
nf = fullfile(A.dir, 'simulink', 'nonideal.mat');
if exist(nf, 'file')
    try
        L = load(nf);  R.SIMNI = L.NISIM;  R.SIMNI.data_all = R.SIMNI.data;
        if ~isempty(R.SIMNI.data)
            ok = cell2mat(R.SIMNI.data(:, strcmp(R.SIMNI.header, 'angles_used'))) == 1;
            R.SIMNI.data = R.SIMNI.data(ok, :);  R.SIMNI.n_invalid = sum(~ok);
            fprintf('non-ideal circuit results found: %d runs, %d left out (commanded angles not run)\n', ...
                sum(ok), sum(~ok));
        end
    catch err
        fprintf('\ncould not read %s: %s\n', nf, err.message);
    end
end
wf = fullfile(A.dir, 'simulink', 'switches.mat');
if exist(wf, 'file')
    try
        L = load(wf);  R.SIMSW = L.SW;
        fprintf('circuit-level branch switches found: %d switch points\n', size(R.SIMSW.data, 1));
    catch err
        fprintf('\ncould not read %s: %s\n', wf, err.message);
    end
end
post = {'figures', 'tables'};
for i = 1:numel(post)
    st = post{i};
    status.(st) = 'skipped';
    if ~isempty(A.only) && ~any(strcmp(st, A.only)) && ~any(strcmp('outputs', A.only)), continue; end
    fprintf('\n--- %s: %s\n', st, describe(st));
    try
        if strcmp(st, 'figures')
            shm_figures(R, struct('DIR', fullfile(A.dir, 'figures')));
        else
            shm_tables(R, struct('DIR', A.dir));
            shm_numbers_ijies(R, A.dir);
            shm_numbers_revision(R, A.dir);
            shm_export_artifact(R, A.dir);
        end
        status.(st) = 'written';
    catch err
        status.(st) = 'FAILED';
        errs{end+1} = sprintf('%s: %s', st, err.message); %#ok<AGROW>
        report_error(st, err);
    end
end
R.meta.status = status;  R.meta.errors = errs;
try
    save(fullfile(A.dir, 'paper_results.mat'), 'R', '-v7');
catch err
    fprintf('    could not save paper_results.mat: %s\n', err.message);
end

fprintf('\n=== summary (%.1f min) ===\n', toc(t0) / 60);
names = [stages, post];
for i = 1:numel(names)
    fprintf('  %-8s %s\n', names{i}, status.(names{i}));
end
if isempty(errs)
    fprintf('\nall stages finished. Start with %s\n', fullfile(A.dir, 'summary.txt'));
else
    fprintf(2, '\n%d stage(s) failed:\n', numel(errs));
    for i = 1:numel(errs), fprintf(2, '  %s\n', errs{i}); end
    fprintf('the full log of this run is %s\n', fullfile(A.dir, 'run_log.txt'));
end
if nargout == 0, clear R; end
end

function C = run_config(mode)
switch mode
    case 'full'
        C.D   = struct('GRID', 0.70:0.0025:1.20, 'ANCHORS', 0.70:0.05:1.20, 'NSTART', 150, 'KEEP', 10, ...
                       'PROBE', 0.705:0.01:1.195, 'PROBE_NSTART', 100, ...
                       'FOLD_ROUNDS', 2, 'FOLD_LOCAL', 30, 'FOLD_NSTART', 100, ...
                       'BIG', [0.7625 0.8625 0.9625 1.0625 1.1625], 'BIG_NSTART', 5000);
        C.S   = struct('MA', 0.70:0.01:1.20, 'R', linspace(1, 1.3, 11), 'GAP', 0.25, 'MAXB', 8);
        C.V   = struct('MA', 0.70:0.025:1.20, 'R', linspace(1, 1.3, 11), 'GAP', linspace(0.25, 2, 11), 'MAXB', 6);
        C.E2  = struct('R', [10 25 50 100], 'R_BASIN', 25, 'CONT_DEG', 5);
        C.E2b = struct('R', [10 25 50 100], 'R_BASIN', 25, 'CONT_DEG', 5, 'R_TEST', 100);
        C.E1b = struct();
        C.E3  = struct('NB', 3, 'K', [6 9 15 25 41], 'SEEDS', 1:5, ...
                       'METHODS', {{'linear', 'spline', 'pchip', 'RBF', 'GP', 'tanhLM', 'MLP', 'ANFIS'}});
        C.E3b = struct('NB', 2, 'Q', [3 4 5 6 8], 'SEEDS', 1:3, ...
                       'METHODS', {{'linear', 'spline', 'RBF', 'GP', 'tanhLM', 'MLP', 'ANFIS'}});
        C.E3c = struct('NB', 2, 'Q', [3 4 5 6], 'SEEDS', 1:3, ...
                       'METHODS', {{'linear', 'spline', 'RBF', 'GP', 'tanhLM', 'MLP', 'ANFIS'}});
        C.E4  = struct('Q1', [5 9 17 33], 'Q2', [3 4 5 6 8], 'Q3', [3 4 5 6], 'SEEDS', 1:5, 'NF', 5, ...
                       'NBOOT', 2000);
        C.E5  = struct();
        C.E6  = struct('NV', 61, 'NSTART', 200, 'BOUND_NSTART', 150);
        C.E7  = struct();
        C.E9  = struct('R', 100, 'STEP', 1);
        C.IND = struct('NSTART', 250, 'ALG', 'interior-point');
        C.CMP = struct();
        C.STO = struct('R', [10 25 50 100], 'SEEDS', 1:5, 'NN_WIDTHS', [2 4 8 16 32 64], 'NN_RESTARTS', 2);
        C.RGS = struct('EXT', struct('WIDTHS', [2 4 8 16 32], 'RESTARTS', 5, 'CV_RESTARTS', 1));
        C.BSW = struct('PF', [0.995 0.95 0.85 0.75 0.6 0.5], 'DELTA', [0 0.02 0.05], 'NTS', 3600, 'R', 100, ...
                       'RAMP_PERIODS', [250 50 10], 'DITHER_PERIODS', 100);
        C.TIM = struct('R', 100, 'NREP', 20000);
    case 'quick'
        mq = struct('EPOCHS', 1500, 'PATIENCE', 300);
        C.D   = struct('GRID', 0.70:0.01:1.20, 'ANCHORS', [0.70 0.80 0.95 1.10 1.20], 'NSTART', 30, 'KEEP', 5, ...
                       'PROBE', [0.75 0.85 1.05 1.15], 'PROBE_NSTART', 15, ...
                       'FOLD_ROUNDS', 1, 'FOLD_LOCAL', 8, 'FOLD_NSTART', 10, ...
                       'BIG', [0.85 1.05], 'BIG_NSTART', 100);
        C.S   = struct('MA', 0.70:0.02:1.20, 'R', linspace(1, 1.3, 4), 'GAP', 0.25, 'MAXB', 4);
        C.V   = struct('MA', 0.70:0.05:1.20, 'R', linspace(1, 1.3, 3), 'GAP', linspace(0.25, 2, 3), 'MAXB', 3);
        C.E2  = struct('R', [5 10 20], 'R_BASIN', 10, 'CONT_DEG', 5);
        C.E2b = struct('R', [5 10 20], 'R_BASIN', 10, 'CONT_DEG', 5, 'R_TEST', 20);
        C.E1b = struct();
        C.E3  = struct('NB', 2, 'K', [6 9 15], 'SEEDS', 1:2, 'mlp', mq, ...
                       'METHODS', {{'linear', 'spline', 'pchip', 'RBF', 'GP', 'tanhLM', 'MLP', 'ANFIS'}});
        C.E3b = struct('NB', 1, 'Q', [3 4], 'SEEDS', 1, 'mlp', mq, ...
                       'METHODS', {{'linear', 'spline', 'RBF', 'GP', 'tanhLM', 'MLP', 'ANFIS'}});
        C.E3c = struct('NB', 1, 'Q', 3, 'SEEDS', 1, 'mlp', mq, ...
                       'METHODS', {{'linear', 'RBF', 'GP', 'tanhLM', 'MLP', 'ANFIS'}});
        C.E4  = struct('Q1', [5 9 17], 'Q2', [3 4], 'Q3', [3], 'SEEDS', 1:2, 'NF', 3, 'NBOOT', 200);
        C.E5  = struct('NREP', 20, 'clsnet', struct('HID', 64, 'NFEAT', 64, 'EPOCHS', 200));
        C.E6  = struct('NV', 7, 'NSTART', 25, 'BOUND_MA', [0.90 1.10], 'BOUND_UB', [89 85], 'BOUND_NSTART', 20);
        C.E7  = struct('TD_US', [0 1 4 16]);
        C.E9  = struct('R', 20, 'STEP', 4);
        C.IND = struct('NSTART', 20, 'ALG', 'interior-point', 'POINTS', 1:5:50);
        C.CMP = struct();
        C.STO = struct('R', [5 10 20], 'SEEDS', 1:2, 'NN_WIDTHS', [2 4 8], 'NN_RESTARTS', 1);
        C.RGS = struct('EXT', struct('WIDTHS', [2 4 8 16], 'RESTARTS', 2, 'CV_RESTARTS', 1));
        C.BSW = struct('PF', [0.85 0.6], 'DELTA', [0 0.05], 'NTS', 720, 'R', 20, ...
                       'RAMP_PERIODS', [50 10], 'DITHER_PERIODS', 20);
        C.TIM = struct('R', 20, 'NREP', 500);
end
end

function o = stage_opts(C, st)
if isfield(C, st), o = C.(st); else, o = struct(); end
if strcmp(st, 'RGS'), o = struct('RGS', o, 'E3', C.E3, 'E3b', C.E3b, 'E3c', C.E3c); end
end

function s = describe(st)
switch st
    case 'D',   s = 'branch set (multistart + continuation + probes + fold-end search + merge + stored subset)';
    case 'S',   s = 'sheets over (Ma, delta)';
    case 'V',   s = 'volumes over (Ma, delta, theta_gap)';
    case 'E1',  s = 'structure of the solution set';
    case 'E1b', s = 'the interior ends of the owning branches: folds, refined continuation';
    case 'E2',  s = 'lookup table, its floor, the corrector, the continuity cost';
    case 'E2b', s = 'the table comparison at the 200 midpoints (Table 1)';
    case 'E3',  s = 'matched-budget regression inside a branch, d = 1';
    case 'E3b', s = 'matched-budget regression inside a branch, d = 2';
    case 'E3c', s = 'matched-budget regression inside a branch, d = 3';
    case 'E4',  s = 'branch selection, d = 1, 2, 3 (corrected protocol)';
    case 'E5',  s = 'storage, scaling and on-line time';
    case 'E6',  s = 'independent multistart, spectrum beyond the band, grid code, angle bound';
    case 'E7',  s = 'gate pattern, dead time and pulse loss';
    case 'E9',  s = 'the tables under non-ideal switching (datasheet IGBTs, dead time, NP offset, timer)';
    case 'IND', s = 'independent check: fmincon interior-point from Halton starts at the midpoints';
    case 'CMP', s = 'grid-code compliance, EN 50160 + CIGRE limits, optimum and tables';
    case 'STO', s = 'Table 1 at equal storage, with sweep, continuous, single-branch and learned baselines';
    case 'RGS', s = 'regression statistics: every seed, model selection, ranking';
    case 'BSW', s = 'branch switch, exact R-L: owner changes x load angle x NP offset, ramps, dither';
    case 'TIM', s = 'on-line cost of reading the tables';
    case 'figures', s = 'figures';
    case 'tables',  s = 'workbook, LaTeX tables, numbers.tex, summary.txt, numbers_ijies.txt';
    otherwise, s = '';
end
end

function X = run_stage(st, R, C, F, seedno)
switch st
    case 'D'
        X = shm_branchset(F.defaults(), C.D);
    case 'S'
        X = shm_extend(R.D, C.S);
    case 'V'
        X = shm_extend(R.D, C.V);
    case 'E1'
        X = shm_x_structure(R.D);
        fprintf('  %s\n  %s\n', X.anchors.note, X.branches.note);
    case 'E1b'
        X = shm_x_folds(R.D, C.E1b);
    case 'E2'
        X = shm_x_lut(R.D, C.E2);
    case 'E2b'
        X = shm_x_lut2(struct('D', R.D, 'E2', R.E2), C.E2b);
    case 'E3'
        X = shm_x_regress(R.D, C.E3);
    case 'E3b'
        o = C.E3b;  o.P = R.D.P;
        X = shm_x_regress_nd(R.S, o);
    case 'E3c'
        o = C.E3c;  o.P = R.D.P;
        X = shm_x_regress_nd(R.V, o);
    case 'E4'
        X = struct();
        base = struct('NF', C.E4.NF, 'SEEDS', C.E4.SEEDS, 'NBOOT', C.E4.NBOOT);
        o = base;  o.Q = C.E4.Q1;
        X.d1 = shm_x_select2({R.D.grid}, R.D.THD(:, R.D.sel), o);
        o = base;  o.Q = C.E4.Q2;
        X.d2 = shm_x_select2({R.S.Ma, R.S.delta}, sheet_thd(R.S, 2), o);
        T3 = sheet_thd(R.V, 3);
        o = base;  o.Q = C.E4.Q3;
        X.d3 = shm_x_select2({R.V.Ma, R.V.delta, R.V.gap}, T3, o);
        o.MIN_ALIVE = 3;
        X.d3c = shm_x_select2({R.V.Ma, R.V.delta, R.V.gap}, T3, o);
    case 'E5'
        X = shm_x_cost(R.D, C.E5);
    case 'E6'
        X = shm_x_validate(R.D, C.E6);
    case 'E7'
        X = shm_x_switching(R.D, C.E7);
        fprintf('  %s\n', X.note);
    case 'E9'
        o = C.E9;  o.TEST = R.E2b.test;
        X = shm_x_nonideal(R.D, o);
    case 'IND'
        X = shm_x_independent(R.D, R.E2b, C.IND);
    case 'CMP'
        X = shm_x_compliance(R.D, R.E2b, C.CMP);
    case 'STO'
        X = shm_x_storage(R.D, R.E2, R.E2b, C.STO);
    case 'RGS'
        cfg = C;  cfg.seed = seedno;
        Eold = struct();
        nm = {'E3', 'E3b', 'E3c'};
        for i = 1:3, if isfield(R, nm{i}), Eold.(nm{i}) = R.(nm{i}); end, end
        X = shm_x_regstats(R.D, R.S, R.V, cfg, Eold);
    case 'BSW'
        X = shm_x_branchswitch(R.D, C.BSW);
    case 'TIM'
        X = shm_x_timing(R.D, R.E2b, C.TIM);
end
end

function T = sheet_thd(S, d)
nk = numel(S.vol);
nM = numel(S.Ma);  nR = numel(S.r);  nG = numel(S.gap);
if d == 2
    T = nan(nM, nR, nk);
    for k = 1:nk
        t = reshape(S.vol(k).thd, nM, nR, nG);  ok = reshape(S.vol(k).ok, nM, nR, nG);
        t(~ok) = NaN;  T(:, :, k) = t(:, :, 1);
    end
else
    T = nan(nM, nR, nG, nk);
    for k = 1:nk
        t = reshape(S.vol(k).thd, nM, nR, nG);  ok = reshape(S.vol(k).ok, nM, nR, nG);
        t(~ok) = NaN;  T(:, :, :, k) = t;
    end
end
end

function X = load_stage(cdir, st, opts)
X = [];
f = fullfile(cdir, [st '.mat']);
if ~exist(f, 'file'), return; end
try
    L = load(f);
    if isfield(L, 'opts') && same_opts(L.opts, opts)
        X = L.X;
    else
        fprintf('\n--- %s: cached result was made with other settings; recomputing\n', st);
    end
catch
    X = [];
end
end

function tf = same_opts(a, b)
if isstruct(a) && isstruct(b)
    fa = sort(fieldnames(a));  fb = sort(fieldnames(b));
    tf = isequal(fa, fb) && isequal(size(a), size(b));
    for i = 1:numel(fa)
        for k = 1:numel(a)
            if ~tf, return; end
            tf = same_opts(a(k).(fa{i}), b(k).(fa{i}));
        end
    end
elseif iscell(a) && iscell(b)
    tf = isequal(size(a), size(b));
    for i = 1:numel(a)
        if ~tf, return; end
        tf = same_opts(a{i}, b{i});
    end
elseif isnumeric(a) && isnumeric(b)
    x = double(a(:));  y = double(b(:));
    tf = isequal(size(a), size(b)) && all(abs(x - y) <= 1e-12 * max(1, abs(y)) | (isnan(x) & isnan(y)));
else
    tf = isequal(a, b);
end
end

function save_stage(cdir, st, X, opts) %#ok<INUSD>
f = fullfile(cdir, [st '.mat']);
save(f, 'X', 'opts', '-v7');
end

function M = run_meta(mode, C, F, t, status, errs)
M.config  = C;
M.backend = F.backend();
M.date    = datestr(now, 'yyyy-mm-dd HH:MM');
M.mode    = mode;
M.platform = platform_str();
M.total_time_s = t;
M.status = status;
M.errors = errs;
end

function report_error(st, err)
fprintf(2, '\n*** %s FAILED: %s\n', st, err.message);
for k = 1:min(4, numel(err.stack))
    fprintf(2, '    at %s, line %d\n', err.stack(k).name, err.stack(k).line);
end
fprintf('    the stages that do not need %s continue.\n', st);
end

function s = platform_str()
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    s = ['Octave ' OCTAVE_VERSION];
else
    s = ['MATLAB ' version];
end
s = [s ', ' computer];
end
