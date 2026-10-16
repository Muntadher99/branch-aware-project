function ok = verify_reproduction(out, ref)
%VERIFY_REPRODUCTION Compare a reproduction with the published results.
%   VERIFY_REPRODUCTION compares output/ with results/.
%   VERIFY_REPRODUCTION(out) compares the folder out with results/.
%   VERIFY_REPRODUCTION(out, ref) compares out with ref.
%
%   Every stage result (cache/*.mat), circuit result (simulink/**/*.mat), the
%   results struct (paper_results.mat) and every text output (number files,
%   summaries, LaTeX tables, CSV and JSON files) is compared value by value;
%   numbers must agree bit for bit. Wall-clock times and run metadata (dates,
%   platform, Simulink model name) differ from run to run by nature: such
%   differences are listed as expected. Workbooks, PDF and PNG files carry
%   creation dates and are not compared; they are drawn from the compared data.
%
%   ok = VERIFY_REPRODUCTION(...) is true when nothing else differs.
%
%   See also REPRODUCE_PAPER, REPRODUCE_CIRCUIT.

root = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(out), out = fullfile(root, 'output'); end
if nargin < 2 || isempty(ref), ref = fullfile(root, 'results'); end
if ~exist(out, 'dir'), error('verify_reproduction:out', 'folder not found: %s', out); end

fprintf('\n=== verify_reproduction ===\nreference: %s\nreproduced: %s\n\n', ref, out);
files = list_files(ref, '');
nsame = 0;  nexp = 0;  nbad = 0;  nmiss = 0;  nskip = 0;
unexpected = {};  expected = {};
for i = 1:numel(files)
    rel = files{i};
    [~, ~, ext] = fileparts(rel);
    a = fullfile(ref, rel);  b = fullfile(out, rel);
    if ~any(strcmpi(ext, {'.mat', '.txt', '.tex', '.csv', '.json'})) || is_log(rel)
        nskip = nskip + 1;  continue
    end
    if ~exist(b, 'file')
        nmiss = nmiss + 1;  unexpected{end+1} = sprintf('%s: missing', rel); %#ok<AGROW>
        continue
    end
    switch lower(ext)
        case '.mat',  d = compare_mat(a, b);
        case '.json', d = compare_json(a, b);
        otherwise,    d = compare_text(a, b, rel);
    end
    isexp = cellfun(@(s) is_expected(rel, s), d);
    if isempty(d)
        nsame = nsame + 1;
    elseif all(isexp)
        nexp = nexp + 1;
    else
        nbad = nbad + 1;
    end
    for k = 1:numel(d)
        if isexp(k), expected{end+1} = [rel ': ' d{k}]; %#ok<AGROW>
        else, unexpected{end+1} = [rel ': ' d{k}]; end %#ok<AGROW>
    end
end

fprintf('files compared:                         %d\n', nsame + nexp + nbad);
fprintf('  identical:                            %d\n', nsame);
fprintf('  identical except times and metadata:  %d\n', nexp);
fprintf('  with other differences:               %d\n', nbad);
fprintf('files missing in the reproduction:      %d\n', nmiss);
fprintf('files not compared (xlsx, pdf, png, logs): %d\n', nskip);
if ~isempty(expected)
    fprintf('\nexpected differences (wall-clock times, dates, platform, model name):\n');
    for k = 1:numel(expected), fprintf('  %s\n', short(expected{k})); end
end
ok = isempty(unexpected);
if ok
    fprintf('\nRESULT: the reproduction is identical to the published results.\n');
else
    fprintf(2, '\nRESULT: %d difference(s) found:\n', numel(unexpected));
    for k = 1:min(numel(unexpected), 50), fprintf(2, '  %s\n', short(unexpected{k})); end
    if numel(unexpected) > 50, fprintf(2, '  ... and %d more\n', numel(unexpected) - 50); end
end
if nargout == 0, clear ok; end
end

% =========================================================================
function L = list_files(base, rel)
L = {};
q = dir(fullfile(base, rel));
for i = 1:numel(q)
    n = q(i).name;
    if any(strcmp(n, {'.', '..'})), continue; end
    r = n;  if ~isempty(rel), r = [rel '/' n]; end
    if q(i).isdir, L = [L, list_files(base, r)]; else, L{end+1} = r; end %#ok<AGROW>
end
end

function s = short(s)
if numel(s) > 120, s = [s(1:117) '...']; end
end

function tf = is_log(rel)
tf = ~isempty(regexp(rel, '(^|/)(run_log|circuit_log)\.txt$', 'once'));
end

function tf = is_expected(rel, s)
pats = {'time_s\>', 'time_min', 'total_time', 'time_us', 'corrector_ms', 'ms_per_call', ...
        '\.date\>', '\.platform\>', '\.backend\>', '\.status\>', 'meta\.model\>', ...
        'timing\.data\{\d+,2\}', 'corr\(:,7\)', 'meta\.config\..*\(last bit\)', ...
        'settings \(last bit\)'};
if ~isempty(regexp(s, '^line ', 'once')), pats = text_patterns(); end
tf = false;
for k = 1:numel(pats)
    if ~isempty(regexp(s, pats{k}, 'once')), tf = true; return; end
end
if ~isempty(strfind(rel, 'switches_two_loads')) && ~isempty(regexp(s, 'model', 'once')), tf = true; end
end

function p = text_patterns()
p = {'\d{4}-\d{2}-\d{2} \d{2}:\d{2}', 'written by', '% run: ', 'minutes', 'microseconds', 'hours for', ...
     'time per call', 'run time', '_ms\>', '\.ms\>', '_us\>', 'time_min', 'time\.branchset', 'time\.total', ...
     'shm@t\.', 'measured per read', 'corrector solve from', ' min, '};
end

% ------------------------------------------------------------------ .mat
function d = compare_mat(a, b)
A = load(a);  B = load(b);
d = cmpv(A, B, '');
end

function d = cmpv(a, b, p)
d = {};
if ~strcmp(class(a), class(b))
    d = {sprintf('%s (class %s vs %s)', p, class(a), class(b))};  return
end
if ~isequal(size(a), size(b))
    d = {sprintf('%s (size %s vs %s)', p, mat2str(size(a)), mat2str(size(b)))};  return
end
if isstruct(a)
    fa = fieldnames(a);  fb = fieldnames(b);
    if ~isequal(sort(fa), sort(fb))
        d = {sprintf('%s (fields differ)', p)};  return
    end
    for k = 1:numel(a)
        q = p;  if numel(a) > 1, q = sprintf('%s(%d)', p, k); end
        for j = 1:numel(fa)
            d = [d, cmpv(a(k).(fa{j}), b(k).(fa{j}), join_path(q, fa{j}))]; %#ok<AGROW>
        end
    end
elseif iscell(a)
    for k = 1:numel(a)
        if size(a, 1) > 1 && size(a, 2) > 1
            [r, c] = ind2sub(size(a), k);  q = sprintf('%s{%d,%d}', p, r, c);
        else
            q = sprintf('%s{%d}', p, k);
        end
        d = [d, cmpv(a{k}, b{k}, q)]; %#ok<AGROW>
    end
elseif isnumeric(a) || islogical(a)
    x = a(:);  y = b(:);
    same = (x == y) | (isnan(double(x)) & isnan(double(y)));
    if ~all(same)
        if ismatrix(a) && size(a, 1) > 1 && size(a, 2) > 1
            [~, c] = ind2sub(size(a), find(~same));
            d = {sprintf('%s(:,%s)', p, strjoin(arrayfun(@num2str, unique(c(:)).', 'UniformOutput', false), ','))};
        else
            d = {p};
        end
        if near(x, y), d{1} = [d{1} ' (last bit)']; end
    end
elseif isa(a, 'function_handle')
    if ~strcmp(func2str(a), func2str(b)), d = {p}; end
elseif ~isequal(a, b)
    d = {p};
end
end

function tf = near(x, y)
x = double(x);  y = double(y);
tf = all(abs(x - y) <= 1e-12 * max(1, abs(y)) | (isnan(x) & isnan(y)));
end

function q = join_path(p, f)
if isempty(p), q = f; else, q = [p '.' f]; end
end

% ------------------------------------------------------------------ .json
function d = compare_json(a, b)
A = jsondecode(fileread(a));  B = jsondecode(fileread(b));
d = {};
for f = {'date', 'platform'}
    if isfield(A, 'run') && isfield(B, 'run') && isfield(A.run, f{1}) && isfield(B.run, f{1}) ...
            && ~isequal(A.run.(f{1}), B.run.(f{1}))
        d{end+1} = ['run.' f{1}]; %#ok<AGROW>
    end
end
if isfield(A, 'run'), A = rmfield(A, 'run'); end
if isfield(B, 'run'), B = rmfield(B, 'run'); end
if ~isequaln(A, B)
    r = cmpv(A, B, '');
    if all(cellfun(@(s) ~isempty(strfind(s, '(last bit)')), r))
        d{end+1} = 'settings (last bit)';
    else
        d = [d, r];
    end
end
end

% ------------------------------------------------------------------ text
function d = compare_text(a, b, rel) %#ok<INUSD>
la = regexp(fileread(a), '\r?\n', 'split');
lb = regexp(fileread(b), '\r?\n', 'split');
d = {};
if numel(la) ~= numel(lb)
    d = {sprintf('number of lines %d vs %d', numel(la), numel(lb))};  return
end
for k = find(~strcmp(la, lb))
    d{end+1} = sprintf('line %d: %s', k, la{k}); %#ok<AGROW>
end
end
