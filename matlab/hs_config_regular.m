function [E, trials] = hs_config_regular(n, k, d, seed, maxTrials)
%HS_CONFIG_REGULAR  Exactly uniform simple d-regular k-uniform hypergraph on n vertices.
%   E = hs_config_regular(n, k, d) returns an m-by-k int32 matrix (m = n*d/k), one sorted
%   edge per row, vertices 1..n; k = 2 gives a d-regular graph.
%   [E, trials] = hs_config_regular(n, k, d, seed, maxTrials) also returns the number of
%   configurations tried. Configuration model conditioned on simplicity, by rejection with
%   early abort. Expected trials 1/P(simple): P(simple) -> exp(-(k-1)(d-1)/2) for k >= 3,
%   and exp(-(d-1)/2 - (d-1)^2/4) for graphs (k = 2).
%   Errors with id 'hypersample:saturation' after maxTrials (default 1e9).
if nargin < 4 || isempty(seed), seed = uint64(randi(flintmax - 1)); end
if nargin < 5 || isempty(maxTrials), maxTrials = 1e9; end
[E, trials] = hs_mex('config_regular', n, k, d, seed, maxTrials);
end
