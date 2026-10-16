function B = shm_trace(x0, Ma0, grid, P, opts)
%SHM_TRACE Follow one branch of local optima over the modulation-index grid.
%   B = SHM_TRACE(x0, Ma0, grid, P, opts) continues the optimum x0 at Ma0 up
%   and down the grid by natural-parameter continuation with a secant
%   predictor. A step is accepted only if every angle of the corrected solution
%   stays within JUMP_DEG of the predictor; otherwise the step is halved, down
%   to MIN_FRAC of a grid step, and the branch ends there.
%
%   Options: JUMP_DEG (0.5), MIN_FRAC (1/16), solve (options of F.solve)
%   Output:  TH (G x N angles, NaN where absent), thd, lo, hi, ends
%            ('range' or 'fold'), steps, sub

if nargin < 5, opts = struct(); end
O.JUMP_DEG = 0.5;
O.MIN_FRAC = 1/16;
O.solve    = struct();
fn = fieldnames(opts); for i = 1:numel(fn), O.(fn{i}) = opts.(fn{i}); end

F = shm_model();
grid = grid(:);  G = numel(grid);  N = P.N;
jmp = O.JUMP_DEG * pi/180;

[~, i0] = min(abs(grid - Ma0));
x = x0(:).';
if abs(grid(i0) - Ma0) > 1e-12
    [x, I] = F.solve(x, grid(i0), P, O.solve);
    if ~I.ok || max(abs(x - x0(:).')) > jmp
        B = empty_branch(G, N);  return
    end
end

B = empty_branch(G, N);
B.TH(i0, :) = x;
B.thd(i0)   = F.thd(x, P);
B.lo = i0;  B.hi = i0;
B.ends = {'range', 'range'};

for dir = [+1, -1]
    ci = i0;  cur = x;  prev = [];
    while true
        ni = ci + dir;
        if ni < 1 || ni > G, break; end
        [xn, ok, ns] = advance(F, P, cur, prev, grid(ci), grid(ni), jmp, O, 0);
        B.sub = B.sub + ns;
        if ~ok
            if dir > 0, B.ends{2} = 'fold'; else, B.ends{1} = 'fold'; end
            break
        end
        B.TH(ni, :) = xn;
        B.thd(ni)   = F.thd(xn, P);
        B.steps     = B.steps + 1;
        prev = cur;  cur = xn;  ci = ni;
        if dir > 0, B.hi = ni; else, B.lo = ni; end
    end
end
end

function [x, ok, nsub] = advance(F, P, cur, prev, ma, mb, jmp, O, depth)
nsub = 0;
if isempty(prev)
    pred = cur;
else
    pred = cur + (cur - prev);
end
pred = F.feasible(pred, P);
[x, I] = F.solve(pred, mb, P, O.solve);
if I.ok && max(abs(x - pred)) <= jmp && max(abs(x - cur)) <= 4*jmp
    ok = true;  return
end
[x, I] = F.solve(cur, mb, P, O.solve);
if I.ok && max(abs(x - cur)) <= jmp
    ok = true;  return
end
if depth >= round(log2(1/O.MIN_FRAC))
    ok = false;  return
end
mm = 0.5 * (ma + mb);
[xm, ok1, n1] = advance(F, P, cur, [], ma, mm, jmp, O, depth + 1);
nsub = 1 + n1;
if ~ok1, ok = false; return; end
[x, ok, n2] = advance(F, P, xm, [], mm, mb, jmp, O, depth + 1);
nsub = nsub + n2;
end

function B = empty_branch(G, N)
B.TH = nan(G, N);  B.thd = nan(G, 1);
B.lo = NaN;  B.hi = NaN;  B.ends = {'none', 'none'};
B.steps = 0;  B.sub = 0;
end
