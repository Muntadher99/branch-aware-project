function E = shm_x_cost(D, opts)
%SHM_X_COST Stage E5: storage, scaling with dimension and on-line time.
%   E = SHM_X_COST(D, opts) reports the stored numbers of the plain and
%   branch-aware tables, how table, classifier and MLP sizes scale with the
%   number of operating parameters, and the measured time per call of a table
%   read, a pchip evaluation, the classifier, the MLP and one corrector solve
%   (machine dependent).
%
%   Options: R ([25 100]), KAX ([9 17 25]), DMAX (5), NREP (200), clsnet,
%            QSEL (6), verbose

if nargin < 2, opts = struct(); end
O.R      = [25 100];
O.KAX    = [9 17 25];
O.DMAX   = 5;
O.NREP   = 200;
O.clsnet = struct('HID', 64, 'NFEAT', 64);
O.QSEL   = 6;
O.verbose = true;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

F = shm_model();  P = D.P;  N = P.N;  grid = D.grid;
nb = numel(D.sel);
H = O.clsnet.HID;  NF = O.clsnet.NFEAT;
cls_params = @(d, C) (2*NF)*H + H + H*H + H + H*C + C + (d + 1)*NF;

hs = {'option','rows_R','stored_numbers','bytes_float32','note'};
rs = {};
for R = O.R
    rs(end+1, :) = {'lookup table (per-Ma optimum)', R, R*N, 4*R*N, 'rows x angles'}; %#ok<AGROW>
    rs(end+1, :) = {'branch-aware table', R, R*N + nb, 4*(R*N + nb), 'rows shared among branches + switch points'}; %#ok<AGROW>
end
nchg = D.stats.owner_changes_stored;
rs(end+1, :) = {'label lookup (selection), d = 1', NaN, 2*nchg + 1, 4*(2*nchg + 1), 'owner-change breakpoints + labels'};
rs(end+1, :) = {sprintf('label lookup (selection), d = 3, %dx%dx%d lattice', O.QSEL, O.QSEL, O.QSEL), NaN, O.QSEL^3, 4*O.QSEL^3, 'one label per lattice node'};
rs(end+1, :) = {'softmax classifier, d = 1', NaN, cls_params(1, nb), 4*cls_params(1, nb), 'weights incl. fixed random features'};
rs(end+1, :) = {'softmax classifier, d = 3', NaN, cls_params(3, nb), 4*cls_params(3, nb), 'weights incl. fixed random features'};
nmlp = 2*128 + 128*128 + 128 + 128*N + N + 2*N;
rs(end+1, :) = {'MLP angle regressor (one branch)', NaN, nmlp, 4*nmlp, '[128 128] hidden, SIREN'};
nanf = 2*7 + 7*2*N;
rs(end+1, :) = {'ANFIS angle regressor (one branch)', NaN, nanf, 4*nanf, '7 Gaussian MFs, first-order Sugeno'};

hsc = {'d','K_per_axis','table_sets_per_branch','table_numbers_per_branch', ...
       'offline_solves_per_branch','classifier_numbers','label_table_numbers','mlp_numbers_per_branch'};
rsc = {};
for d = 1:O.DMAX
    for K = O.KAX
        rsc(end+1, :) = {d, K, K^d, K^d * N, K^d, cls_params(d, nb), K^d, nmlp_d(d, N)}; %#ok<AGROW>
    end
end
dstar_sel = zeros(size(O.KAX));  dstar_reg = dstar_sel;
for i = 1:numel(O.KAX)
    K = O.KAX(i);
    dstar_sel(i) = log(cls_params(3, nb)) / log(K);
    dstar_reg(i) = log(nmlp_d(3, N) / N) / log(K);
end

k = D.sel(1);
ex = find(isfinite(D.THD(:, k)));
MaB = grid(ex);  YB = D.TH(ex, :, k);
R = min(25, numel(ex));
ri = unique(round(linspace(1, numel(ex), R)));
pp = pchip(MaB(ri).', YB(ri, :).');
m  = 0.5 * (MaB(1) + MaB(end));
t_lut   = timeit_(@() interp1(MaB(ri), YB(ri, :), m, 'linear'), O.NREP);
t_pchip = timeit_(@() ppval(pp, m), O.NREP);

own = D.owns;  okg = isfinite(own);
Xc = grid(okg);  yc = own(okg);
co = O.clsnet;  co.EPOCHS = 20;
if numel(unique(yc)) > 1
    M = shm_clsnet(Xc, yc, co);
    t_cls = timeit_(@() M.predict(m), O.NREP);
else
    t_cls = NaN;
end
Mm = shm_mlp(MaB(ri), YB(ri, :), struct('EPOCHS', 20, 'PATIENCE', 20));
t_mlp = timeit_(@() Mm.predict(m), O.NREP);

tq = ex(round(numel(ex)/2));
x0 = D.TH(tq, :, k) + 0.002 * randn(1, N);
t_sqp = timeit_(@() F.solve(x0, grid(tq), P), max(5, round(O.NREP / 20)));

ht = {'operation','time_per_call_us','note'};
rt = {'table lookup + linear interpolation', 1e6*t_lut, sprintf('%d rows', R);
      'pchip evaluation inside a branch',      1e6*t_pchip, 'precomputed coefficients';
      'softmax classifier forward pass',       1e6*t_cls, sprintf('%d classes', numel(unique(yc)));
      'MLP forward pass',                      1e6*t_mlp, '[128 128] hidden';
      'SQP corrector (fmincon, warm start)',   1e6*t_sqp, 'one solve from a nearby start'};

E.storage = struct('name', 'E5a_storage', 'header', {hs}, 'data', {rs}, ...
    'note', 'Memory of each deployment option at float32 (4 bytes per stored number).');
E.scaling = struct('name', 'E5b_scaling', 'header', {hsc}, 'data', {rsc}, ...
    'note', sprintf(['Analytic scaling. A table over d parameters at K points per axis stores K^d ' ...
    'angle sets and needs K^d offline solves PER BRANCH; the classifier grows by %d numbers per ' ...
    'added parameter. Selection: a K^d label table outgrows the classifier from d* = %s; ' ...
    'regression: a K^d angle table outgrows one MLP per branch from d* = %s (K = %s).'], NF, ...
    mat2str(round(dstar_sel*100)/100), mat2str(round(dstar_reg*100)/100), mat2str(O.KAX)));
E.timing  = struct('name', 'E5c_online_time', 'header', {ht}, 'data', {rt}, ...
    'note', 'Measured on this machine; compare ratios rather than absolutes.');
E.dstar_sel = dstar_sel;
E.dstar_reg = dstar_reg;
E.KAX = O.KAX;
end

function n = nmlp_d(d, N)
n = d*128 + 128 + 128*128 + 128 + 128*N + N + 2*N;
end

function t = timeit_(f, n)
if exist('timeit', 'file') == 2 && exist('OCTAVE_VERSION', 'builtin') == 0
    t = timeit(f);  return
end
f();
t0 = tic;
for i = 1:n, f(); end
t = toc(t0) / n;
end
