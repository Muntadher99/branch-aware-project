function D = shm_branchset(P, opts)
%SHM_BRANCHSET Decompose the SHM solution set into continuous branches.
%   D = SHM_BRANCHSET(P, opts)
%     1. multistart at the ANCHORS (SHM_SEED); the KEEP best optima seed
%     2. continuation of every seed over GRID (SHM_TRACE)
%     3. probe multistarts at PROBE and searches just beyond the interior ends
%        of owning branches; every optimum below the envelope is traced
%     4. merging of coincident traces, lower envelope, and greedy selection of
%        the stored subset whose envelope reproduces the full one to GAP_TOL
%        (mean) and GAP_MAX (worst), in pp
%
%   Options: GRID (0.70:0.0025:1.20), ANCHORS (0.70:0.05:1.20), NSTART (150),
%            KEEP (10), PROBE, PROBE_NSTART (100), FOLD_ROUNDS (2),
%            FOLD_LOCAL (30), FOLD_NSTART (100), BIG, BIG_NSTART (5000),
%            GAP_TOL (0.005), GAP_MAX (0.05), verbose
%   Output:  grid, TH (G x N x K), THD (G x K), env, own (envelope and owner
%            over all branches), sel (stored subset), envs, owns, envTH,
%            branch, seeds, probes, stats, time_s

if nargin < 2, opts = struct(); end
O.GRID     = 0.70:0.0025:1.20;
O.ANCHORS  = 0.70:0.05:1.20;
O.NSTART   = 150;
O.KEEP     = 10;
O.SAME_DEG = 0.10;
O.GAP_TOL  = 0.005;
O.GAP_MAX  = 0.05;
O.KMAX     = 16;
O.PROBE        = [];
O.PROBE_NSTART = 100;
O.PROBE_TOL    = 0.01;
O.THD_SAME     = 1e-6;
O.FOLD_ROUNDS  = 2;
O.FOLD_LOCAL   = 30;
O.FOLD_SIGMA   = 2;
O.FOLD_NSTART  = 100;
O.BIG          = [];
O.BIG_NSTART   = 5000;
O.BIG_SEED     = 777;
O.trace    = struct();
O.seed     = struct();
O.verbose  = true;
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

F = shm_model();
grid = O.GRID(:);  G = numel(grid);  N = P.N;
same = O.SAME_DEG * pi/180;
t0 = tic;  F.nfe_reset();

na = numel(O.ANCHORS);
seeds = repmat(struct('Ma',0,'X',[],'thd',[],'count',[],'nstart',0,'nok',0,'fev',0,'time_s',0), 1, na);
say(O, '  [seed] %d anchors x %d random starts\n', na, O.NSTART);
for a = 1:na
    so = O.seed;  so.NSTART = O.NSTART;
    seeds(a) = shm_seed(O.ANCHORS(a), P, so);
    if isempty(seeds(a).thd)
        say(O, '    Ma %.3f: NO start converged -- this anchor seeds nothing (check the solver exit flags)\n', O.ANCHORS(a));
        continue
    end
    say(O, '    Ma %.3f: %3d converged, %3d distinct optima, best THD49 %.4f%%, best found by %d starts\n', ...
        O.ANCHORS(a), seeds(a).nok, numel(seeds(a).thd), 100*seeds(a).thd(1), seeds(a).count(1));
end
if all(arrayfun(@(s) isempty(s.thd), seeds))
    error('shm_branchset:noconv', ['no random start converged at any anchor. Run check_setup; ' ...
        'if fmincon returns exit flag 0 for most starts, raise the iteration limit in shm_model.']);
end

TH = nan(G, N, 0);  THD = nan(G, 0);  info = struct('lo',{},'hi',{},'ends',{},'steps',{},'sub',{},'src',{});
ntr = 0;  nskip = 0;
for a = 1:na
    [~, ia] = min(abs(grid - O.ANCHORS(a)));
    for j = 1:min(O.KEEP, numel(seeds(a).thd))
        x = seeds(a).X(j, :);
        if onbranch(TH, ia, x, same), nskip = nskip + 1; continue; end
        B = shm_trace(x, O.ANCHORS(a), grid, P, O.trace);
        if isnan(B.lo), continue; end
        ntr = ntr + 1;
        TH(:, :, end+1) = B.TH; %#ok<AGROW>
        THD(:, end+1)   = B.thd; %#ok<AGROW>
        info(end+1) = struct('lo',B.lo,'hi',B.hi,'ends',{B.ends},'steps',B.steps,'sub',B.sub,'src',1); %#ok<AGROW>
    end
    say(O, '    after anchor %.3f: %d traces, %d seeds skipped as already traced\n', ...
        O.ANCHORS(a), ntr, nskip);
end

[TH, THD, info] = merge(TH, THD, info, same, O.THD_SAME);
say(O, '  [merge] %d traces -> %d distinct branches\n', ntr, size(THD, 2));

probes = struct('Ma', {}, 'distinct', {}, 'best_pct', {}, 'env_pct', {}, 'traced', {});
ntr_seed = ntr;
if ~isempty(O.PROBE)
    say(O, '  [probe] %d modulation indices x %d random starts\n', numel(O.PROBE), O.PROBE_NSTART);
    envp = min(THD, [], 2);
    for m = O.PROBE(:).'
        [~, ip] = min(abs(grid - m));
        so = O.seed;  so.NSTART = O.PROBE_NSTART;
        Sp = shm_seed(grid(ip), P, so);
        e = envp(ip);  if ~isfinite(e), e = inf; end
        n0 = ntr;
        for j = 1:numel(Sp.thd)
            if ntr - n0 >= O.KEEP || ~(100 * (e - Sp.thd(j)) > O.PROBE_TOL), break; end
            x = Sp.X(j, :);
            if onbranch(TH, ip, x, same), continue; end
            B = shm_trace(x, grid(ip), grid, P, O.trace);
            if isnan(B.lo), continue; end
            ntr = ntr + 1;
            TH(:, :, end+1) = B.TH; %#ok<AGROW>
            THD(:, end+1)   = B.thd; %#ok<AGROW>
            info(end+1) = struct('lo',B.lo,'hi',B.hi,'ends',{B.ends},'steps',B.steps,'sub',B.sub,'src',2); %#ok<AGROW>
            envp = min(envp, B.thd);
        end
        best = NaN;  if ~isempty(Sp.thd), best = 100 * Sp.thd(1); end
        probes(end+1) = struct('Ma', grid(ip), 'distinct', numel(Sp.thd), 'best_pct', best, ...
                               'env_pct', 100 * e, 'traced', ntr - n0); %#ok<AGROW>
        if ntr > n0
            say(O, '    Ma %.4f: %d new trace(s); envelope there %.4f%% -> %.4f%%\n', ...
                grid(ip), ntr - n0, 100 * e, 100 * envp(ip));
        end
    end
    say(O, '    the probes added %d traces at %d of %d modulation indices\n', ntr - ntr_seed, ...
        sum([probes.traced] > 0), numel(probes));
    [TH, THD, info] = merge(TH, THD, info, same, O.THD_SAME);
    say(O, '  [merge] %d traces -> %d distinct branches\n', ntr, size(THD, 2));
end
ntr_probe = ntr;

big = struct('Ma', {}, 'nstart', {}, 'nok', {}, 'distinct', {}, 'F1', {}, 'F2', {}, 'chao1', {}, ...
             'best_pct', {}, 'best_hits', {}, 'env_pct', {}, 'gap_pp', {}, 'traced', {}, 'time_s', {});
if ~isempty(O.BIG)
    say(O, '  [large multistart] %d modulation indices x %d random starts\n', numel(O.BIG), O.BIG_NSTART);
    rs0 = rng;  rng(O.BIG_SEED);
    envp = min(THD, [], 2);
    for m = O.BIG(:).'
        [~, ip] = min(abs(grid - m));
        so = O.seed;  so.NSTART = O.BIG_NSTART;
        Sp = shm_seed(grid(ip), P, so);
        e = envp(ip);  if ~isfinite(e), e = inf; end
        n0 = ntr;
        for j = 1:numel(Sp.thd)
            if ntr - n0 >= O.KEEP || ~(100 * (e - Sp.thd(j)) > O.PROBE_TOL), break; end
            x = Sp.X(j, :);
            if onbranch(TH, ip, x, same), continue; end
            B = shm_trace(x, grid(ip), grid, P, O.trace);
            if isnan(B.lo), continue; end
            ntr = ntr + 1;
            TH(:, :, end+1) = B.TH; %#ok<AGROW>
            THD(:, end+1)   = B.thd; %#ok<AGROW>
            info(end+1) = struct('lo',B.lo,'hi',B.hi,'ends',{B.ends},'steps',B.steps,'sub',B.sub,'src',3); %#ok<AGROW>
            envp = min(envp, B.thd);
        end
        n = Sp.nok;  S = numel(Sp.thd);  F1 = sum(Sp.count == 1);  F2 = sum(Sp.count == 2);
        ch = S + (n - 1) / max(n, 1) * F1 * (F1 - 1) / (2 * (F2 + 1));
        best = NaN;  hits = 0;  if S > 0, best = 100 * Sp.thd(1);  hits = Sp.count(1); end
        big(end+1) = struct('Ma', grid(ip), 'nstart', O.BIG_NSTART, 'nok', n, 'distinct', S, 'F1', F1, ...
            'F2', F2, 'chao1', ch, 'best_pct', best, 'best_hits', hits, 'env_pct', 100 * e, ...
            'gap_pp', best - 100 * e, 'traced', ntr - n0, 'time_s', Sp.time_s); %#ok<AGROW>
        say(O, ['    Ma %.4f: %d/%d converged, %d distinct (%d found once, Chao1 %.0f); best - envelope %+.4f pp, ' ...
                'best found by %d starts; %d new trace(s) (%.0f s)\n'], grid(ip), n, O.BIG_NSTART, S, F1, ch, ...
                best - 100 * e, hits, ntr - n0, Sp.time_s);
    end
    rng(rs0);
    [TH, THD, info] = merge(TH, THD, info, same, O.THD_SAME);
    say(O, '  [merge] %d traces -> %d distinct branches\n', ntr, size(THD, 2));
end
ntr_big = ntr;

done = false(G, 1);  nfoldpts = 0;
for rnd = 1:O.FOLD_ROUNDS
    if O.FOLD_NSTART + O.FOLD_LOCAL <= 0, break; end
    [envp, ownp] = min(THD, [], 2);  ownp(~isfinite(envp)) = NaN;
    tg = zeros(0, 3);
    for k = unique(ownp(isfinite(ownp))).'
        f = find(isfinite(THD(:, k)));
        if f(1) > 1,   tg(end+1, :) = [f(1) - 1, k, f(1)];   end %#ok<AGROW>
        if f(end) < G, tg(end+1, :) = [f(end) + 1, k, f(end)]; end %#ok<AGROW>
    end
    tg = tg(~done(tg(:, 1)), :);
    if isempty(tg), break; end
    pts = unique(tg(:, 1)).';
    say(O, '  [fold ends] round %d: %d points just beyond the ends of owning branches\n', rnd, numel(pts));
    n0 = ntr;
    for gi = pts
        done(gi) = true;  nfoldpts = nfoldpts + 1;
        rows = tg(tg(:, 1) == gi, :);
        X0 = zeros(0, N);
        for r = 1:size(rows, 1)
            xe = TH(rows(r, 3), :, rows(r, 2));
            X0 = [X0; repmat(xe, O.FOLD_LOCAL, 1) + (O.FOLD_SIGMA*pi/180) * randn(O.FOLD_LOCAL, N)]; %#ok<AGROW>
        end
        so = O.seed;  so.NSTART = O.FOLD_NSTART;  so.X0 = X0;
        Sp = shm_seed(grid(gi), P, so);
        e = envp(gi);  if ~isfinite(e), e = inf; end
        nk = 0;
        for j = 1:numel(Sp.thd)
            if nk >= O.KEEP || ~(100 * (e - Sp.thd(j)) > O.PROBE_TOL), break; end
            x = Sp.X(j, :);
            if onbranch(TH, gi, x, same), continue; end
            B = shm_trace(x, grid(gi), grid, P, O.trace);
            if isnan(B.lo), continue; end
            ntr = ntr + 1;  nk = nk + 1;
            TH(:, :, end+1) = B.TH; %#ok<AGROW>
            THD(:, end+1)   = B.thd; %#ok<AGROW>
            info(end+1) = struct('lo',B.lo,'hi',B.hi,'ends',{B.ends},'steps',B.steps,'sub',B.sub,'src',4); %#ok<AGROW>
        end
        if nk > 0
            say(O, '    Ma %.4f: %d new trace(s); envelope there %.4f%% -> %.4f%%\n', ...
                grid(gi), nk, 100 * e, 100 * min(THD(gi, :)));
        end
    end
    say(O, '    %d new traces\n', ntr - n0);
    if ntr == n0, break; end
    [TH, THD, info] = merge(TH, THD, info, same, O.THD_SAME);
    say(O, '  [merge] %d traces -> %d distinct branches\n', ntr, size(THD, 2));
end

[env0, own0] = min(THD, [], 2);
npt  = sum(isfinite(THD), 1);
keep = npt >= 3 | ismember(1:size(THD, 2), own0(isfinite(env0)).');
if any(~keep)
    TH = TH(:, :, keep);  THD = THD(:, keep);  info = info(keep);
    say(O, '  [clean] %d short traces that own nothing dropped -> %d branches\n', sum(~keep), size(THD, 2));
end
K = size(THD, 2);

[env, own] = min(THD, [], 2);
own(~isfinite(env)) = NaN;
envTH = nan(G, N);
for g = 1:G
    if isfinite(own(g)), envTH(g, :) = TH(g, :, own(g)); end
end

sel = greedy(THD, env, O);
[envs, ow] = min(THD(:, sel), [], 2);
owns = sel(ow);  owns = owns(:);  owns(~isfinite(envs)) = NaN;
gapsel = 100 * (envs - env);
say(O, '  [select] %d of %d branches stored: mean gap %.4f pp, worst %.4f pp\n', ...
    numel(sel), K, mean(gapsel(isfinite(gapsel))), max(gapsel(isfinite(gapsel))));

br = branch_diag(TH, THD, own, info, grid, P);
st = headline(seeds, TH, THD, env, own, envTH, owns, sel, grid, ntr, gapsel);
st.n_probe_traces = ntr_probe - ntr_seed;
st.n_probes       = numel(probes);
st.n_probes_hit   = sum([probes.traced] > 0);
st.n_big_traces   = ntr_big - ntr_probe;
st.n_fold_traces  = ntr - ntr_big;
st.n_fold_points  = nfoldpts;
st.n_dropped      = sum(~keep);
src = [info.src];
nb  = src ~= 3;
st.n_big_branches = sum(~nb);
st.big_owners     = intersect(find(~nb), unique(own(isfinite(own))).');
st.n_big_owners   = numel(st.big_owners);
if any(~nb) && any(nb)
    [env0, own0] = min(THD(:, nb), [], 2);
    idx0 = find(nb);  own0 = idx0(own0);  own0(~isfinite(env0)) = NaN;  own0 = own0(:);
    v0 = isfinite(own0(1:end-1)) & isfinite(own0(2:end));
    st.owners_without_big        = numel(unique(own0(isfinite(own0))));
    st.owner_changes_without_big = sum(own0(1:end-1) ~= own0(2:end) & v0);
    gain = 100 * (env0 - env);  gain(~isfinite(gain)) = 0;
    st.big_gain_max_pp  = max(gain);
    st.big_gain_mean_pp = mean(gain(isfinite(env)));
    og = find(ismember(own, st.big_owners));
    st.big_owned_points = numel(og);
    if ~isempty(og), st.big_owned_Ma = [grid(og(1)), grid(og(end))]; else, st.big_owned_Ma = [NaN NaN]; end
    if st.n_big_owners > 0
        k1 = st.big_owners(1);  f1 = find(isfinite(THD(:, k1)));
        st.big_owner_range = [grid(f1(1)), grid(f1(end))];
    else
        st.big_owner_range = [NaN NaN];
    end
else
    st.owners_without_big = st.n_owners;  st.owner_changes_without_big = st.owner_changes;
    st.big_gain_max_pp = 0;  st.big_gain_mean_pp = 0;  st.big_owned_points = 0;
    st.big_owned_Ma = [NaN NaN];  st.big_owner_range = [NaN NaN];
end

D.grid   = grid;
D.TH     = TH;      D.THD  = THD;
D.env    = env;     D.own  = own;    D.envTH = envTH;
D.sel    = sel;     D.envs = envs;   D.owns  = owns;
D.branch = br;
D.seeds  = seeds;
D.probes = probes;
D.big    = big;
D.src    = [info.src];
D.stats  = st;
D.P      = P;
D.opts   = O;
D.nfe    = F.nfe();
D.time_s = toc(t0);
say(O, '  [branchset] done: %d branches, %d owners, %d owner changes, %.0f s\n', ...
    K, st.n_owners, st.owner_changes, D.time_s);
end

function tf = onbranch(TH, ia, x, same)
tf = false;
if isempty(TH), return; end
R = squeeze(TH(ia, :, :));
if size(R, 1) ~= numel(x), R = R.'; end
if isvector(R) && size(TH,3) == 1, R = R(:); end
d = max(abs(R - repmat(x(:), 1, size(R, 2))), [], 1);
tf = any(d <= same);
end

function [TH, THD, info] = merge(TH, THD, info, same, thdsame)
K = size(THD, 2);  keep = true(1, K);
for a = 1:K
    if ~keep(a), continue; end
    for b = a+1:K
        if ~keep(b), continue; end
        sh = isfinite(THD(:, a)) & isfinite(THD(:, b));
        if ~any(sh), continue; end
        if sum(sh) >= 3
            d = max(abs(TH(sh, :, a) - TH(sh, :, b)), [], 2);
            if mean(d <= same) >= 0.5
                onlyb = ~isfinite(THD(:, a)) & isfinite(THD(:, b));
                TH(onlyb, :, a) = TH(onlyb, :, b);
                THD(onlyb, a)   = THD(onlyb, b);
                info(a).lo = min(info(a).lo, info(b).lo);
                info(a).hi = max(info(a).hi, info(b).hi);
                keep(b) = false;
            end
            continue
        end
        na = sum(isfinite(THD(:, a)));  nb = sum(isfinite(THD(:, b)));
        if min(na, nb) <= 2 && sum(sh) == min(na, nb) && max(abs(THD(sh, a) - THD(sh, b))) <= thdsame
            if na < nb, keep(a) = false;  break; end
            keep(b) = false;
        end
    end
end
TH = TH(:, :, keep);  THD = THD(:, keep);  info = info(keep);
end

function sel = greedy(THD, env, O)
fin = isfinite(env);
K = size(THD, 2);
cur = inf(size(env));  sel = [];
for t = 1:min(O.KMAX, K)
    bestv = inf;  bk = 0;
    for k = 1:K
        if any(sel == k), continue; end
        m = min(cur, THD(:, k));
        g = m(fin) - env(fin);
        miss = ~isfinite(g);
        gg = g;  gg(miss) = 0;
        v = mean(gg) + 0.25 * max(gg) + 10 * mean(miss);
        if v < bestv, bestv = v; bk = k; end
    end
    if bk == 0, break; end
    sel(end+1) = bk; %#ok<AGROW>
    cur = min(cur, THD(:, bk));
    g = 100 * (cur(fin) - env(fin));
    if all(isfinite(g)) && mean(g) < O.GAP_TOL && max(g) < O.GAP_MAX, break; end
end
end

function br = branch_diag(TH, THD, own, info, grid, P)
K = size(THD, 2);  G = numel(grid);
br = repmat(struct('Ma_lo',NaN,'Ma_hi',NaN,'ends','','n',0,'owned',0,'median_step_deg',NaN, ...
    'max_step_deg',NaN,'collapsed_pairs',NaN,'effective_angles',NaN,'theta1_deg',NaN, ...
    'thetaN_deg',NaN,'thd_min',NaN,'thd_mean',NaN), 1, K);
for k = 1:K
    f = find(isfinite(THD(:, k)));
    if isempty(f), continue; end
    X = TH(f, :, k);
    st = max(abs(diff(X, 1, 1)), [], 2) * 180/pi;
    gp = diff(X, 1, 2);
    col = sum(gp < 2 * P.minGap + 1e-12, 2);
    br(k).Ma_lo = grid(f(1));  br(k).Ma_hi = grid(f(end));
    br(k).ends  = sprintf('%s/%s', info(k).ends{1}, info(k).ends{2});
    br(k).n     = numel(f);
    br(k).owned = sum(own == k) / G;
    if ~isempty(st), br(k).median_step_deg = median(st); br(k).max_step_deg = max(st); end
    br(k).collapsed_pairs  = mean(col);
    br(k).effective_angles = P.N - mean(col);
    br(k).theta1_deg = mean(X(:, 1)) * 180/pi;
    br(k).thetaN_deg = mean(X(:, end)) * 180/pi;
    br(k).thd_min  = 100 * min(THD(f, k));
    br(k).thd_mean = 100 * mean(THD(f, k));
end
end

function st = headline(seeds, TH, THD, env, own, envTH, owns, sel, grid, ntr, gapsel)
st.n_anchors      = numel(seeds);
seeds = seeds(arrayfun(@(s) ~isempty(s.thd), seeds));
st.distinct_per_anchor = arrayfun(@(s) numel(s.thd), seeds);
st.mean_distinct  = mean(st.distinct_per_anchor);
st.max_distinct   = max(st.distinct_per_anchor);
st.best_hit_rate  = mean(arrayfun(@(s) s.count(1) / max(s.nok, 1), seeds));
st.thd_spread_pp  = mean(arrayfun(@(s) 100 * (median(s.thd) - s.thd(1)), seeds));
st.n_traces       = ntr;
st.n_branches     = size(THD, 2);
st.n_owners       = numel(unique(own(isfinite(own))));
st.n_stored       = numel(sel);
v = isfinite(own(1:end-1)) & isfinite(own(2:end));
st.owner_changes  = sum(own(1:end-1) ~= own(2:end) & v);
st.owner_changes_stored = sum(owns(1:end-1) ~= owns(2:end) & isfinite(owns(1:end-1)) & isfinite(owns(2:end)));
dj = max(abs(diff(envTH, 1, 1)), [], 2) * 180/pi;
st.env_step_mean_deg   = mean(dj(isfinite(dj)));
st.env_step_max_deg    = max(dj(isfinite(dj)));
chg = own(1:end-1) ~= own(2:end) & v;
st.env_step_at_change_mean_deg = mean(dj(chg & isfinite(dj)));
st.stored_gap_mean_pp  = mean(gapsel(isfinite(gapsel)));
st.stored_gap_max_pp   = max(gapsel(isfinite(gapsel)));
st.coverage            = mean(isfinite(env));
st.env_thd_mean_pct    = 100 * mean(env(isfinite(env)));
st.grid_step           = median(diff(grid));
end

function say(O, varargin)
if O.verbose, fprintf(varargin{:}); end
end
