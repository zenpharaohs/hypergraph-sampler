function [E, trials] = hs_config_degrees(deg, k, seed, maxTrials)
%HS_CONFIG_DEGREES  Exactly uniform simple k-uniform hypergraph with degree sequence deg.
%   E = hs_config_degrees(deg, k): deg(v) is the degree of vertex v (sum divisible by k).
%   See hs_config_regular for the method and outputs.
if nargin < 2 || isempty(k), k = 2; end
if nargin < 3 || isempty(seed), seed = uint64(randi(flintmax - 1)); end
if nargin < 4 || isempty(maxTrials), maxTrials = 1e9; end
[E, trials] = hs_mex('config_degrees', deg, k, seed, maxTrials);
end
