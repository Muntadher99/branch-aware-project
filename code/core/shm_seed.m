function S = shm_seed(Ma, P, opts)
%SHM_SEED Multistart at one modulation index: the distinct local optima.
%   S = SHM_SEED(Ma, P, opts) solves NSTART random feasible starts (plus the
%   rows of X0) and clusters the converged solutions: two solutions within
%   TOL_DEG of each other in every angle are the same optimum.
%
%   Options: NSTART (150), X0 (extra starts, one per row), TOL_DEG (0.05),
%            solve (options passed to F.solve)
%   Output:  Ma, X (distinct optima, best first, rad), thd, count, nstart,
%            nok (converged starts), fev (objective evaluations), time_s

if nargin < 3, opts = struct(); end
O.NSTART  = 150;
O.X0      = [];
O.TOL_DEG = 0.05;
O.solve   = struct();
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

F  = shm_model();
t0 = tic;  f0 = F.nfe();
X0 = [O.X0; F.starts(O.NSTART, P)];
[X, I] = F.solve(X0, Ma, P, O.solve);

X = X(I.ok, :);  t = I.thd(I.ok);
[t, o] = sort(t);  X = X(o, :);
tol = O.TOL_DEG * pi/180;
rep = zeros(0, P.N);  thd = zeros(0, 1);  cnt = zeros(0, 1);
for i = 1:size(X, 1)
    if ~isempty(rep)
        d = max(abs(rep - repmat(X(i,:), size(rep,1), 1)), [], 2);
        [dm, j] = min(d);
        if dm <= tol, cnt(j) = cnt(j) + 1; continue; end
    end
    rep(end+1, :) = X(i, :); %#ok<AGROW>
    thd(end+1, 1) = t(i); %#ok<AGROW>
    cnt(end+1, 1) = 1; %#ok<AGROW>
end

S.Ma     = Ma;
S.X      = rep;
S.thd    = thd;
S.count  = cnt;
S.nstart = O.NSTART;
S.nok    = sum(I.ok);
S.fev    = F.nfe() - f0;
S.time_s = toc(t0);
end
