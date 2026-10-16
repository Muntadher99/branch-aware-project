function E = shm_x_regress_nd(V, opts)
%SHM_X_REGRESS_ND Stages E3b/E3c: regression inside a branch, d = 2 and 3.
%   E = SHM_X_REGRESS_ND(S, opts)  over (Ma, delta)        (sheets of SHM_EXTEND)
%   E = SHM_X_REGRESS_ND(V, opts)  over (Ma, delta, gap)   (volumes of SHM_EXTEND)
%   For each of NB branches, a q^d training lattice is placed in the largest
%   box where the branch exists; every other grid point of the box is a test
%   point. Scored as in SHM_X_REGRESS, with each test point's own imbalance
%   and minimum pulse width.
%
%   Options: NB (2), Q ([3 4 5 6]), SEEDS (1:3), MINRUN (5), MINLV (3), P,
%            METHODS, mlp, anfis, tanhlm, gpr, KEEP_SEEDS (false), verbose

if nargin < 2, opts = struct(); end
O.NB      = 2;
O.Q       = [3 4 5 6];
O.SEEDS   = 1:3;
O.MINRUN  = 5;
O.MINLV   = 3;
O.P       = [];
O.mlp     = struct();
O.anfis   = struct();
O.tanhlm  = struct();
O.gpr     = struct();
O.METHODS = {'linear', 'spline', 'MLP', 'ANFIS'};
O.verbose = true;
O.KEEP_SEEDS = false;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

F  = shm_model();
P0 = O.P;  if isempty(P0), P0 = F.defaults(); end
N  = P0.N;
MaV = V.Ma(:);  rV = V.r(:);  dV = V.delta(:);  gV = V.gap(:);
nM = numel(MaV);  nR = numel(rV);  nG = numel(gV);
d  = 2 + (nG > 1);
meth = O.METHODS;

nb = numel(V.vol);
box = repmat(struct('i', [], 'jr', 0, 'kg', 0, 'vol', 0), 1, nb);
for v = 1:nb
    OK = V.vol(v).ok;
    if isempty(OK), continue; end
    OK = reshape(OK, nM, nR, nG);
    best = 0;
    kgs = nG:-1:1;  if d == 3, kgs = kgs(kgs >= min(O.MINLV, nG)); end
    for jr = nR:-1:min(O.MINLV, nR)
        for kg = kgs
            okM = all(all(OK(:, 1:jr, 1:kg), 2), 3);
            [i1, i2] = longest_run(okM);
            len = i2 - i1 + 1;
            if isempty(i1) || len < min(O.MINRUN, nM), continue; end
            vol = len * jr * kg;
            if vol > best
                best = vol;
                box(v) = struct('i', i1:i2, 'jr', jr, 'kg', kg, 'vol', vol);
            end
        end
    end
end
[~, order] = sort([box.vol], 'descend');
order = order([box(order).vol] > 0);
order = order(1:min(O.NB, numel(order)));

hdr = {'dimension','method','branch','points_per_axis','n_train','n_test','angle_MAE_deg', ...
       'angle_MAE_sd_deg','angle_max_err_deg','THD_excess_mean_pp','THD_excess_max_pp', ...
       'pct_infeasible','stored_parameters','seed_CV_pct'};
rows = {};
hb = {'branch','Ma_from','Ma_to','Ma_levels','delta_from','delta_to','delta_levels', ...
      'gap_from_deg','gap_to_deg','gap_levels','box_points'};
rb = {};
curve = [];
seedrows = cell(0, 12);

for u = 1:numel(order)
    v = order(u);  B = box(v);
    im = B.i(:);  jr = (1:B.jr).';  kg = (1:B.kg).';
    TH = reshape(V.vol(v).TH, nM, nR, nG, N);
    TD = reshape(V.vol(v).thd, nM, nR, nG);
    rb(end+1, :) = {V.vol(v).branch, MaV(im(1)), MaV(im(end)), numel(im), dV(1), dV(B.jr), B.jr, ...
                    gV(1), gV(B.kg), B.kg, numel(im) * B.jr * B.kg}; %#ok<AGROW>
    if O.verbose
        fprintf('  branch %d: box Ma %.3f..%.3f x %d delta levels x %d gap levels (%d points)\n', ...
            V.vol(v).branch, MaV(im(1)), MaV(im(end)), B.jr, B.kg, numel(im) * B.jr * B.kg);
    end
    [I1, I2, I3] = ndgrid(im, jr, kg);
    I1 = I1(:);  I2 = I2(:);  I3 = I3(:);
    nbx = numel(I1);
    Yall = zeros(nbx, N);
    for t = 1:nbx, Yall(t, :) = reshape(TH(I1(t), I2(t), I3(t), :), 1, N); end
    Tall = TD(sub2ind([nM nR nG], I1, I2, I3));
    Xall = [MaV(I1), dV(I2), gV(I3)];

    for q = O.Q
        lm = im(unique(round(linspace(1, numel(im), min(q, numel(im))))));
        lr = jr(unique(round(linspace(1, numel(jr), min(q, numel(jr))))));
        lg = kg(unique(round(linspace(1, numel(kg), min(q, numel(kg))))));
        tr = ismember(I1, lm) & ismember(I2, lr) & ismember(I3, lg);
        te = ~tr;
        if sum(te) < 3, continue; end
        ntr = sum(tr);
        Ytr = zeros(numel(lm), numel(lr), numel(lg), N);
        for a = 1:numel(lm)
            for b = 1:numel(lr)
                for c = 1:numel(lg)
                    Ytr(a, b, c, :) = TH(lm(a), lr(b), lg(c), :);
                end
            end
        end
        Xte = Xall(te, :);  Yte = Yall(te, :);  Tte = Tall(te);
        rte = rV(I2(te));  gte = gV(I3(te));
        for m = 1:numel(meth)
            switch meth{m}
                case {'linear', 'spline'}
                    lv = [numel(lm), numel(lr), numel(lg)];
                    if d == 2, lv = lv(1:2); end
                    if any(lv < 2) || (strcmp(meth{m}, 'spline') && any(lv < 4)), continue; end
                    Yp = grid_interp(MaV(lm), dV(lr), gV(lg), Ytr, Xte, d, meth{m});
                    if any(~isfinite(Yp(:))), continue; end
                    r = evaluate(F, P0, Yp, Yte, Xte(:, 1), rte, gte, Tte);
                    rows(end+1, :) = {d, meth{m}, V.vol(v).branch, q, ntr, sum(te), r(1), 0, r(2), ...
                                      r(3), r(4), r(5), ntr * N, 0}; %#ok<AGROW>
                    curve(end+1, :) = [q, m, r(1), r(3), ntr]; %#ok<AGROW>
                case {'MLP', 'tanhLM'}
                    Xin = Xall(tr, 1:d);  Yin = Yall(tr, :);
                    rs = zeros(numel(O.SEEDS), 5);  np = 0;
                    for si = 1:numel(O.SEEDS)
                        if strcmp(meth{m}, 'MLP')
                            mo = O.mlp;  mo.SEED = O.SEEDS(si);
                            M = shm_mlp(Xin, Yin, mo);
                        else
                            mo = O.tanhlm;  mo.SEED = O.SEEDS(si);
                            M = shm_tanhlm(Xin, Yin, mo);
                        end
                        Yp = M.predict(Xte(:, 1:d));
                        rs(si, :) = evaluate(F, P0, Yp, Yte, Xte(:, 1), rte, gte, Tte);
                        if isfield(M, 'nparam'), np = M.nparam; end
                        if O.KEEP_SEEDS
                            wd = NaN;  if isfield(M, 'width'), wd = M.width; end
                            seedrows(end+1, :) = {d, meth{m}, V.vol(v).branch, q, ntr, O.SEEDS(si), rs(si, 1), ...
                                rs(si, 2), rs(si, 3), rs(si, 4), rs(si, 5), wd}; %#ok<AGROW>
                        end
                    end
                    cv = 100 * std(rs(:, 1)) / max(mean(rs(:, 1)), eps);
                    rows(end+1, :) = {d, meth{m}, V.vol(v).branch, q, ntr, sum(te), mean(rs(:,1)), ...
                        std(rs(:,1)), max(rs(:,2)), mean(rs(:,3)), max(rs(:,4)), mean(rs(:,5)), np, cv}; %#ok<AGROW>
                    curve(end+1, :) = [q, m, mean(rs(:, 1)), mean(rs(:, 3)), ntr]; %#ok<AGROW>
                case {'GP', 'RBF'}
                    go = O.gpr;  if strcmp(meth{m}, 'RBF'), go.KERNEL = 'phs'; end
                    M = shm_gpr(Xall(tr, 1:d), Yall(tr, :), go);
                    Yp = M.predict(Xte(:, 1:d));
                    r = evaluate(F, P0, Yp, Yte, Xte(:, 1), rte, gte, Tte);
                    rows(end+1, :) = {d, meth{m}, V.vol(v).branch, q, ntr, sum(te), r(1), 0, r(2), ...
                                      r(3), r(4), r(5), M.nparam, 0}; %#ok<AGROW>
                    curve(end+1, :) = [q, m, r(1), r(3), ntr]; %#ok<AGROW>
                case 'ANFIS'
                    M = shm_anfis(Xall(tr, 1:d), Yall(tr, :), O.anfis);
                    Yp = M.predict(Xte(:, 1:d));
                    r = evaluate(F, P0, Yp, Yte, Xte(:, 1), rte, gte, Tte);
                    rows(end+1, :) = {d, 'ANFIS', V.vol(v).branch, q, ntr, sum(te), r(1), 0, r(2), ...
                                      r(3), r(4), r(5), M.nparam, 0}; %#ok<AGROW>
                    curve(end+1, :) = [q, m, r(1), r(3), ntr]; %#ok<AGROW>
            end
            if O.verbose
                r0 = rows(end, :);
                fprintf('    q=%d n=%4d %-7s MAE %8.4f deg  max %7.3f deg  THD excess %8.4f pp  params %6d\n', ...
                    q, ntr, r0{2}, r0{7}, r0{9}, r0{10}, r0{13});
            end
        end
    end
end

Qs = [];  mae = [];  thx = [];  ntrain = [];
if ~isempty(curve)
    Qs = unique(curve(:, 1)).';
    mae = nan(numel(Qs), numel(meth));  thx = mae;  ntrain = nan(numel(Qs), 1);
    for i = 1:numel(Qs)
        si = curve(:, 1) == Qs(i);
        ntrain(i) = mean(curve(si, 5));
        for m = 1:numel(meth)
            s = si & curve(:, 2) == m;
            if any(s), mae(i, m) = mean(curve(s, 3));  thx(i, m) = mean(curve(s, 4)); end
        end
    end
end

E.name    = sprintf('E3%s_regression_%dD', char('a' + d - 1), d);
E.header  = hdr;
E.data    = rows;
E.boxes   = struct('name', sprintf('E3%s_boxes_%dD', char('a' + d - 1), d), 'header', {hb}, 'data', {rb});
if O.KEEP_SEEDS
    E.seed_rows = seedrows;
    E.seed_header = {'d','method','branch','q','n_train','seed','angle_MAE_deg','angle_max_err_deg', ...
        'THD_excess_mean_pp','THD_excess_max_pp','pct_infeasible','hidden_width'};
end
E.d       = d;
E.curveQ  = Qs(:);
E.curveN  = ntrain;
E.curve   = mae;
E.curveTHD = thx;
E.methods = meth;
E.note = sprintf(['Matched-budget regression inside a branch over %d operating parameters: every ' ...
    'method sees the same q^d training lattice of the same branch and is scored on every other grid ' ...
    'point of the branch''s box. linear = multilinear, spline = tensor-product cubic spline ' ...
    '(interpn), MLP = shm_mlp mean over %d seeds. THD excess is evaluated at each test point''s own ' ...
    'imbalance and minimum pulse width.'], d, numel(O.SEEDS));
end

function [i1, i2] = longest_run(tf)
i1 = [];  i2 = [];  best = 0;  s = 0;
tf = tf(:).';
for k = 1:numel(tf)
    if tf(k)
        if s == 0, s = k; end
        if k - s + 1 > best, best = k - s + 1;  i1 = s;  i2 = k; end
    else
        s = 0;
    end
end
end

function Yp = grid_interp(ga, gb, gc, Ytr, Xq, d, method)
N = size(Ytr, 4);  Yp = nan(size(Xq, 1), N);
for n = 1:N
    if d == 2
        Z = Ytr(:, :, 1, n);
        Yp(:, n) = interpn(ga(:), gb(:), Z, Xq(:, 1), Xq(:, 2), method);
    else
        Z = Ytr(:, :, :, n);
        Yp(:, n) = interpn(ga(:), gb(:), gc(:), Z, Xq(:, 1), Xq(:, 2), Xq(:, 3), method);
    end
end
end

function r = evaluate(F, P0, Yp, Yt, Ma, rr, gg, tb)
ae = abs(Yp - Yt) * 180/pi;
n  = size(Yp, 1);
ex = nan(n, 1);  infe = false(n, 1);
[u, ~, ic] = unique([rr(:), gg(:)], 'rows');
for c = 1:size(u, 1)
    s = find(ic == c);
    P = F.set_gap(F.set_r(P0, u(c, 1)), u(c, 2));
    Yc = Yp(s, :);
    infe(s) = any(diff(Yc, 1, 2) < P.minGap - 1e-12, 2) | Yc(:, 1) < P.lb | Yc(:, end) > P.ub;
    for i = 1:numel(s), Yc(i, :) = F.feasible(Yc(i, :), P); end
    R = F.report(Yc, Ma(s), P);
    ex(s) = 100 * (R.THD49 - tb(s));
end
r = [mean(ae(:)), max(ae(:)), mean(ex), max(ex), 100 * mean(infe)];
end
