function E = hs_gnp(n, p, k, seed)
%HS_GNP  Erdos-Renyi G(n,p) on k-sets: each k-subset of 1..n is an edge with probability p.
%   E = hs_gnp(n, p) is an ordinary random graph; hs_gnp(n, p, k) a random k-uniform
%   hypergraph. Geometric skipping over the C(n,k) k-subsets: expected time proportional to
%   the number of edges.
if nargin < 3 || isempty(k), k = 2; end
if nargin < 4 || isempty(seed), seed = uint64(randi(flintmax - 1)); end
E = hs_mex('gnp', n, p, k, seed);
end
