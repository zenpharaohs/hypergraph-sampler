function E = hs_gnm(n, m, k, seed)
%HS_GNM  Uniform G(n,m) on k-sets: m distinct k-subsets of 1..n, uniformly at random.
if nargin < 3 || isempty(k), k = 2; end
if nargin < 4 || isempty(seed), seed = uint64(randi(flintmax - 1)); end
E = hs_mex('gnm', n, m, k, seed);
end
