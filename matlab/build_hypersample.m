function build_hypersample()
%BUILD_HYPERSAMPLE  Compile the hypersample MEX gateway (hs_mex) from the C core.
%   Run once from any directory: build_hypersample. Works in MATLAB (needs a C compiler
%   configured with mex -setup C) and in GNU Octave (mkoctfile). The core uses C11 and
%   __int128, so GCC or Clang (Xcode Clang on macOS); MSVC is not supported.
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
srcs = dir(fullfile(root, 'src', '*.c'));
srcs = fullfile(root, 'src', {srcs.name});
inc = ['-I' fullfile(root, 'include')];
gateway = fullfile(here, 'hs_mex.c');
if exist('OCTAVE_VERSION', 'builtin')
    out = fullfile(here, ['hs_mex.' mexext]);
    [status, msg] = system(sprintf('mkoctfile --mex -O2 -std=c11 %s -o "%s" "%s" %s', inc, out, gateway, ...
        sprintf('"%s" ', srcs{:})));
    if status, error('hypersample:build', '%s', msg); end
else
    mex('-O', '-largeArrayDims', 'CFLAGS=$CFLAGS -std=c11', inc, '-outdir', here, gateway, srcs{:});
end
addpath(here);
fprintf('built %s (hypersample %s)\n', fullfile(here, ['hs_mex.' mexext]), hs_mex('version'));
end
