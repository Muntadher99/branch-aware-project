function root = setup_paths()
%SETUP_PATHS Add the package's code and Simulink models to the MATLAB path.
%   root = SETUP_PATHS() returns the package folder. Called by every
%   entry point; call it yourself before using a function from code/ directly.

root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(genpath(fullfile(root, 'code')));
addpath(fullfile(root, 'models'));
if nargout == 0, clear root; end
end
