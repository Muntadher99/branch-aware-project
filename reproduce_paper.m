function R = reproduce_paper(mode, varargin)
%REPRODUCE_PAPER Reproduce the results of the paper (MATLAB part).
%   REPRODUCE_PAPER regenerates every figure, table and number file from the
%   stage results in results/cache and the circuit results in
%   results/simulink, and writes them to output/ (a few minutes).
%
%   REPRODUCE_PAPER('stages', {'IND', 'CMP'}) recomputes the named stages, and
%   every stage that depends on them, and loads the others from results/cache.
%
%   REPRODUCE_PAPER('all') recomputes every stage from scratch and writes to
%   output_all/ (about 8 hours; see README.md). The run can be interrupted:
%   the same command continues from the first unfinished stage.
%
%   REPRODUCE_PAPER('quick') runs the whole chain on small grids and writes to
%   output_quick/. It tests the installation; its numbers are not the paper's.
%
%   REPRODUCE_PAPER(..., 'out', folder) writes to another folder.
%   R = REPRODUCE_PAPER(...) also returns the results struct.
%
%   The circuit simulations are rerun by REPRODUCE_CIRCUIT. A run is compared
%   with the published results by VERIFY_REPRODUCTION.
%
%   See also REPRODUCE_CIRCUIT, VERIFY_REPRODUCTION, RUN_PAPER, CHECK_SETUP.

if nargin < 1 || isempty(mode), mode = 'outputs'; end
mode = lower(mode);
stages = {};
if strcmp(mode, 'stages')
    if isempty(varargin)
        error('reproduce_paper:stages', 'name the stages to recompute, e.g. reproduce_paper(''stages'', {''IND''})');
    end
    stages = varargin{1};  varargin(1) = [];
    if ischar(stages), stages = {stages}; end
end
A.out = '';
for i = 1:2:numel(varargin), A.(lower(varargin{i})) = varargin{i+1}; end

root = setup_paths();
ref = fullfile(root, 'results');
switch mode
    case {'outputs', 'stages'}, folder = 'output';
    case 'all',                 folder = 'output_all';
    case 'quick',               folder = 'output_quick';
    otherwise
        error('reproduce_paper:mode', 'mode must be ''outputs'', ''stages'', ''all'' or ''quick''');
end
out = A.out;
if isempty(out), out = fullfile(root, folder); end
if ~exist(out, 'dir'), mkdir(out); end

if strcmp(mode, 'quick')
    R = run_paper('quick', 'dir', out);
else
    if any(strcmp(mode, {'outputs', 'stages'}))
        copy_missing(fullfile(ref, 'cache'), fullfile(out, 'cache'));
    end
    copy_missing(fullfile(ref, 'simulink'), fullfile(out, 'simulink'));
    switch mode
        case 'outputs', R = run_paper('full', 'dir', out, 'only', {'outputs'});
        case 'stages',  R = run_paper('full', 'dir', out, 'redo', stages);
        case 'all',     R = run_paper('full', 'dir', out);
    end
end
shm_numbers_supplement(R, out);
fprintf('\nresults written to %s\n', out);
if nargout == 0, clear R; end
end

function copy_missing(src, dst)
if ~exist(src, 'dir'), return; end
if ~exist(dst, 'dir'), mkdir(dst); end
L = dir(src);
for i = 1:numel(L)
    if any(strcmp(L(i).name, {'.', '..'})), continue; end
    s = fullfile(src, L(i).name);  d = fullfile(dst, L(i).name);
    if L(i).isdir
        copy_missing(s, d);
    elseif ~exist(d, 'file')
        copyfile(s, d);
    end
end
end
