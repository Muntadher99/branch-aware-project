function ok = check_setup()
%CHECK_SETUP Check that this computer can reproduce the paper (about a minute).
%   ok = CHECK_SETUP() prints PASS or FAIL for each check and estimates the
%   run time on this computer:
%     1  MATLAB version (R2016b or newer)
%     2  fmincon (Optimization Toolbox)
%     3  analytic gradients against finite differences
%     4  a multistart solve at Ma = 1.00
%     5  a continuation step
%     6  closed-form spectrum against the exact Fourier integral of the pattern
%     7  DC-link imbalance model against an FFT
%     8  writing an .xlsx file
%     9  exporting a figure as PDF and PNG
%     10 the softmax classifier and the MLP regressor
%     11 the current folder is writable

here = fileparts(mfilename('fullpath'));
addpath(here);
setup_paths();
fprintf('\n=== check_setup ===\n');
tests = {@t_version, @t_backend, @t_gradients, @t_solve, @t_warm, @t_exact, @t_imbalance, ...
         @t_xlsx, @t_figure, @t_learners, @t_writable};
names = {'MATLAB version', 'fmincon available', 'analytic gradients', 'multistart solve', ...
         'continuation step', 'spectrum vs gate pattern', 'imbalance model vs FFT', ...
         'xlsx writer', 'figure export', 'learners', 'writable folder'};
res = true(1, numel(tests));
S = struct();
for i = 1:numel(tests)
    try
        [p, msg, S] = tests{i}(S);
    catch err
        p = false;  msg = err.message;
    end
    res(i) = p;
    tag = 'PASS';  if ~p, tag = 'FAIL'; end
    fprintf('  %2d %-26s %s  %s\n', i, names{i}, tag, msg);
end
ok = all(res);

if isfield(S, 't_cold') && isfield(S, 't_warm')
    n_cold = 11 * 150 + 50 * 100 + 6000 + 61 * 200 + 5 * 5 * 150;
    n_warm = 45 * 201 + 51 * 11 * 8 + 21 * 11 * 11 * 6 + 4 * 190 + 2000;
    t_solve = n_cold * S.t_cold + n_warm * S.t_warm;
    t_learn = 0;
    if isfield(S, 't_mlp'), t_learn = t_learn + (75 + 54) * S.t_mlp; end
    if isfield(S, 't_cls'), t_learn = t_learn + 170 * S.t_cls; end
    fprintf('\n  one cold solve %.0f ms, one warm solve %.0f ms', 1e3 * S.t_cold, 1e3 * S.t_warm);
    if isfield(S, 't_mlp'), fprintf(', one MLP fit %.1f s, one classifier fit %.1f s', S.t_mlp, S.t_cls); end
    fprintf(['\n  estimated time: reproduce_paper(''quick'') %.0f min; stages D to E9 %.1f h ' ...
        '(a full recomputation takes about 4.5 times as long)\n'], ...
        (t_solve + t_learn) / 60 / 25 + 2, (t_solve + t_learn) / 3600);
end
if ok
    fprintf('\nall checks passed. Next: reproduce_paper\n');
else
    fprintf(2, '\n%d check(s) failed: fix these before reproduce_paper.\n', sum(~res));
end
if nargout == 0, clear ok; end
end

function [p, msg, S] = t_version(S)
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    p = true;  msg = sprintf('Octave %s (the paper''s numbers must come from MATLAB)', OCTAVE_VERSION);
    return
end
v = sscanf(version, '%d.%d');
p = v(1) > 9 || (v(1) == 9 && v(2) >= 1);
msg = sprintf('MATLAB %s', version);
if ~p, msg = [msg ' -- R2016b (9.1) or newer is needed']; end
end

function [p, msg, S] = t_backend(S)
F = shm_model();
b = F.backend();
p = strcmp(b, 'fmincon') || exist('OCTAVE_VERSION', 'builtin') ~= 0;
msg = sprintf('backend %s', b);
if strcmp(b, 'sqp'), msg = [msg ' (Octave fallback)']; end
end

function [p, msg, S] = t_gradients(S)
F = shm_model();  P = F.set_r(F.defaults(), 1.2);
rng(3);
x = F.starts(1, P).';
[f, g] = F.objgrad(x, P, 1);
h = 1e-6;  gf = zeros(size(x));
[~, ceq0, ~, ge] = F.eqcon(x, 1.0, P);  gef = zeros(size(x));
for k = 1:numel(x)
    e = zeros(size(x));  e(k) = h;
    gf(k) = (F.objgrad(x + e, P, 1) - F.objgrad(x - e, P, 1)) / (2*h);
    [~, cp] = F.eqcon(x + e, 1.0, P);  [~, cm] = F.eqcon(x - e, 1.0, P);
    gef(k) = (cp - cm) / (2*h);
end
e1 = norm(g - gf) / max(norm(g), eps);
e2 = norm(ge - gef) / max(norm(ge), eps);
p = e1 < 1e-5 && e2 < 1e-6 && isequal(size(g), size(x)) && isequal(size(ge), size(x)) && isscalar(ceq0) && isscalar(f);
msg = sprintf('objective rel. err %.1e, constraint rel. err %.1e', e1, e2);
end

function [p, msg, S] = t_solve(S)
F = shm_model();  P = F.defaults();
rng(1);
X0 = F.starts(12, P);
t = tic;
[X, I] = F.solve(X0, 1.0, P);
S.t_cold = toc(t) / 12;
p = any(I.ok);
uf = unique(I.flag(:)).';
fl = strjoin(arrayfun(@(u) sprintf('%d x%d', u, sum(I.flag == u)), uf, 'UniformOutput', false), ', ');
if p
    [best, j] = min(I.thd + 1e9 * ~I.ok);
    S.xbest = X(j, :);
    msg = sprintf('%d of 12 starts converged, best THD49 %.3f%%, |b1 - Ma| %.1e; exit flags: %s', ...
        sum(I.ok), 100 * best, I.b1err(j), fl);
else
    msg = sprintf('no start converged (exit flags: %s; last error: %s)', fl, I.lasterr);
end
end

function [p, msg, S] = t_warm(S)
F = shm_model();  P = F.defaults();
if ~isfield(S, 'xbest'), p = false;  msg = 'needs test 4'; return; end
t = tic;
[x, I] = F.solve(S.xbest, 1.0025, P);
S.t_warm = toc(t);
d = max(abs(x - S.xbest)) * 180/pi;
p = I.ok && d < 0.5;
msg = sprintf('converged %d, largest angle move %.3f deg, %d iterations', I.ok, d, I.iter);
end

function [p, msg, S] = t_exact(S)
F = shm_model();  P = F.defaults();
if ~isfield(S, 'xbest'), p = false;  msg = 'needs test 4'; return; end
R = shm_switchcheck(S.xbest, P, struct('verbose', false));
p = R.max_rel_dev < 1e-9 && abs(R.THD49_exact - R.THD49_model) < 1e-9;
msg = sprintf('closed form vs exact integral: %.1e relative', R.max_rel_dev);
end

function [p, msg, S] = t_imbalance(S)
F = shm_model();  P = F.set_r(F.defaults(), 1.25);
rng(5);
th = F.starts(1, P);
M = 2^20;  t = (0:M-1).' * 2*pi / M;
w = zeros(M, 1);
q = mod(t, pi);  q(q > pi/2) = pi - q(q > pi/2);
for k = 1:P.N, w = w + P.s(k) * (q >= th(k)); end
w(t >= pi) = -w(t >= pi);
v = w .* (1 + P.delta * (w > 0) - P.delta * (w < 0));
X = fft(v);
n = 1:49;
a = (2 * real(X(n + 1)) / M).';  b = (-2 * imag(X(n + 1)) / M).';
[Bo, Be, b1] = F.harmonics(th, P);
dev = max([abs(b(P.ORD) - Bo(:).'), abs(a(P.ORDE) - Be(:).'), abs(b(1) - b1)]);
thd_fft = sqrt(sum(b(P.ORD).^2) + sum(a(P.ORDE).^2)) / abs(b(1));
thd_mod = F.thd(th, P);
p = dev < 1e-4 && abs(thd_fft - thd_mod) < 1e-3 * thd_mod;
msg = sprintf('r = 1.25: worst coefficient deviation %.1e, THD49 FFT %.4f%% vs model %.4f%%', ...
    dev, 100 * thd_fft, 100 * thd_mod);
end

function [p, msg, S] = t_xlsx(S)
f = [tempname '.xlsx'];
shm_xlsx(f, struct('name', 'check', 'header', {{'a', 'b'}}, 'data', [1 2], 'note', 'self-test'));
p = exist(f, 'file') == 2;
msg = 'multi-sheet workbook written';
if p, delete(f); else, msg = 'no file was produced'; end
end

function [p, msg, S] = t_figure(S)
f = figure('Visible', 'off');
plot(1:3, [1 3 2]);
b = tempname;
isoct = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if ~isoct && any(exist('exportgraphics') == [2 5 6]) %#ok<EXIST>
    exportgraphics(f, [b '.pdf'], 'ContentType', 'vector');
    exportgraphics(f, [b '.png'], 'Resolution', 150);
    how = 'exportgraphics';
else
    print(f, [b '.pdf'], '-dpdf');
    print(f, [b '.png'], '-dpng');
    how = 'print';
end
close(f);
p = exist([b '.pdf'], 'file') == 2 && exist([b '.png'], 'file') == 2;
msg = sprintf('PDF and PNG via %s', how);
if exist([b '.pdf'], 'file'), delete([b '.pdf']); end
if exist([b '.png'], 'file'), delete([b '.png']); end
end

function [p, msg, S] = t_learners(S)
X = linspace(0, 1, 40).';  y = 1 + (X > 0.5);
t = tic;
M = shm_clsnet(X, y, struct('EPOCHS', 400));
S.t_cls = toc(t) * 4000 / 400 * 2.5;
acc = 100 * mean(M.predict(X) == y);
Xr = linspace(0, 1, 25).';  Y = [sin(3 * Xr), cos(2 * Xr)];
t = tic;
Mr = shm_mlp(Xr, Y, struct('EPOCHS', 800, 'PATIENCE', 800));
S.t_mlp = toc(t) * 8000 / 800;
err = max(max(abs(Mr.predict(Xr) - Y)));
p = acc >= 95 && err < 0.05;
msg = sprintf('classifier train accuracy %.0f%%, MLP worst fit error %.3f', acc, err);
end

function [p, msg, S] = t_writable(S)
f = fullfile(pwd, sprintf('.check_setup_%d.tmp', round(1e6 * rand)));
fid = fopen(f, 'w');
p = fid > 0;
if p, fclose(fid); delete(f); msg = pwd; else, msg = ['cannot write to ' pwd]; end
end
