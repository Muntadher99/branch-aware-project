function [C, ST] = shm_sim_cases(R, opts)
%SHM_SIM_CASES Angle sets for the circuit simulations (Section 6.4).
%   [C, ST] = SHM_SIM_CASES(R, opts) writes <DIR>/simulink/cases.mat, cases.csv
%   and steps.csv. C: the optimum at MA, and the plain table with ROWS rows at
%   its worst grid point with the optimum and the owning stored branch there.
%   ST: branch switches at the owner change with the largest waveform change
%   and at the one with the largest angle move, against staying on the old
%   branch. Angles in rad, first quarter cycle.
%
%   Options: DIR ('output'), ROWS ([25 100]), MA ([0.80 1.00 1.20])

if nargin < 1 || isempty(R)
    L = load(fullfile('output', 'paper_results.mat'));  R = L.R;
end
if nargin < 2, opts = struct(); end
O.DIR  = 'output';
O.ROWS = [25 100];
O.MA   = [0.80 1.00 1.20];
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

F = shm_model();  D = R.D;  P = D.P;  grid = D.grid;
C = struct('name', {}, 'Ma', {}, 'alpha', {}, 'alpha_deg', {}, 's', {}, 'b1', {}, ...
           'THD49_model_pct', {}, 'what', {});

for m = O.MA
    [dg, gi] = min(abs(grid - m));
    if dg > 1e-9 || ~isfinite(D.env(gi)), continue; end
    C = addcase(C, F, P, sprintf('opt_Ma%03d', round(100 * m)), grid(gi), D.envTH(gi, :), ...
        'per-Ma optimum');
end

fg = find(isfinite(D.env));
for rows = O.ROWS(:).'
    ri = unique(fg(round(linspace(1, numel(fg), rows))));
    te = setdiff(fg, ri);
    if isempty(te)
        fprintf('  a %d-row table stores every grid point here: no test point, skipped\n', rows);
        continue
    end
    Ya = interp1(grid(ri), D.envTH(ri, :), grid(te), 'linear');
    Rt = F.report(sort(Ya, 2), grid(te), P);
    [~, w] = max(100 * (Rt.THD49 - D.env(te)));
    g = te(w);  tag = sprintf('%03d', rows);
    C = addcase(C, F, P, ['lut' tag], grid(g), Ya(w, :), ...
        sprintf('lookup table, %d rows, linear interpolation, at its worst point', rows));
    C = addcase(C, F, P, ['opt_at_lut' tag], grid(g), D.envTH(g, :), 'per-Ma optimum at the same Ma');
    k = D.owns(g);
    if isfinite(k)
        C = addcase(C, F, P, ['branch_at_lut' tag], grid(g), D.TH(g, :, k), ...
            sprintf('stored branch %d, which owns this Ma (branch-aware table)', k));
    end
end

ST = struct('name', {}, 'Ma1', {}, 'Ma2', {}, 'alpha1', {}, 'alpha2', {}, 'alpha1_deg', {}, ...
            'alpha2_deg', {}, 's', {}, 'jump_deg', {}, 'wave_diff_pct', {}, 'THD1_model_pct', {}, ...
            'THD2_model_pct', {}, 'what', {});
own = D.own(:);
ch = find(own(1:end-1) ~= own(2:end) & isfinite(own(1:end-1)) & isfinite(own(2:end)));
cand = zeros(0, 3);
for g = ch(:).'
    if ~isfinite(D.THD(g + 1, own(g))), continue; end
    cand(end+1, :) = [g, max(abs(D.envTH(g + 1, :) - D.envTH(g, :))) * 180/pi, ...
                      shm_wavediff(D.envTH(g, :), D.envTH(g + 1, :), P.s)]; %#ok<AGROW>
end
if ~isempty(cand)
    [~, q] = max(cand(:, 3));  [~, qa] = max(cand(:, 2));
    pick = {cand(q, 1), ''};
    if cand(qa, 1) ~= cand(q, 1), pick(2, :) = {cand(qa, 1), '_angle'}; end
    for r = 1:size(pick, 1)
        g = pick{r, 1};  k1 = own(g);  k2 = own(g + 1);  sfx = pick{r, 2};
        ST = addstep(ST, F, P, ['step_jump' sfx], grid(g), grid(g + 1), D.envTH(g, :), D.envTH(g + 1, :), ...
            sprintf('optimum: branch %d -> branch %d', k1, k2));
        ST = addstep(ST, F, P, ['step_stay' sfx], grid(g), grid(g + 1), D.envTH(g, :), D.TH(g + 1, :, k1), ...
            sprintf('continuous: branch %d -> branch %d', k1, k1));
    end
end

d = fullfile(O.DIR, 'simulink');
if ~exist(d, 'dir'), mkdir(d); end
save(fullfile(d, 'cases.mat'), 'C', 'ST', '-v7');
fid = fopen(fullfile(d, 'cases.csv'), 'w');
fprintf(fid, 'name,Ma,b1_of_angles,THD49_model_pct');
for j = 1:P.N, fprintf(fid, ',alpha%02d_deg', j); end
fprintf(fid, ',what\n');
for i = 1:numel(C)
    fprintf(fid, '%s,%.6f,%.6f,%.6f', C(i).name, C(i).Ma, C(i).b1, C(i).THD49_model_pct);
    fprintf(fid, ',%.6f', C(i).alpha_deg);
    fprintf(fid, ',"%s"\n', C(i).what);
end
fclose(fid);
fid = fopen(fullfile(d, 'steps.csv'), 'w');
fprintf(fid, 'name,Ma_from,Ma_to,jump_deg,waveform_differs_pct,THD49_from_pct,THD49_to_pct');
for j = 1:P.N, fprintf(fid, ',from%02d_deg', j); end
for j = 1:P.N, fprintf(fid, ',to%02d_deg', j); end
fprintf(fid, ',what\n');
for i = 1:numel(ST)
    fprintf(fid, '%s,%.6f,%.6f,%.4f,%.4f,%.6f,%.6f', ST(i).name, ST(i).Ma1, ST(i).Ma2, ST(i).jump_deg, ...
        ST(i).wave_diff_pct, ST(i).THD1_model_pct, ST(i).THD2_model_pct);
    fprintf(fid, ',%.6f', ST(i).alpha1_deg);  fprintf(fid, ',%.6f', ST(i).alpha2_deg);
    fprintf(fid, ',"%s"\n', ST(i).what);
end
fclose(fid);

fprintf('  %d steady-state cases and %d branch-switch cases written to %s\n', numel(C), numel(ST), d);
for i = 1:numel(C)
    fprintf('    %-18s Ma %.4f   THD49 (model, line-to-line) %8.4f%%\n', C(i).name, C(i).Ma, C(i).THD49_model_pct);
end
for i = 1:numel(ST)
    fprintf('    %-18s Ma %.4f -> %.4f   largest angle move %.1f deg, waveform changes over %.1f%% of the period   (%s)\n', ...
        ST(i).name, ST(i).Ma1, ST(i).Ma2, ST(i).jump_deg, ST(i).wave_diff_pct, ST(i).what);
end
end

function C = addcase(C, F, P, name, Ma, th, what)
th = F.feasible(sort(th(:)).', P);
C(end+1) = struct('name', name, 'Ma', Ma, 'alpha', th, 'alpha_deg', th * 180/pi, 's', P.s(1:numel(th)), ...
    'b1', F.b1(th, P), 'THD49_model_pct', 100 * F.thd(th, P), 'what', what);
end

function ST = addstep(ST, F, P, name, Ma1, Ma2, th1, th2, what)
th1 = F.feasible(sort(th1(:)).', P);  th2 = F.feasible(sort(th2(:)).', P);
ST(end+1) = struct('name', name, 'Ma1', Ma1, 'Ma2', Ma2, 'alpha1', th1, 'alpha2', th2, ...
    'alpha1_deg', th1 * 180/pi, 'alpha2_deg', th2 * 180/pi, 's', P.s(1:numel(th1)), ...
    'jump_deg', max(abs(th2 - th1)) * 180/pi, 'wave_diff_pct', shm_wavediff(th1, th2, P.s), ...
    'THD1_model_pct', 100 * F.thd(th1, P), ...
    'THD2_model_pct', 100 * F.thd(th2, P), 'what', what);
end
