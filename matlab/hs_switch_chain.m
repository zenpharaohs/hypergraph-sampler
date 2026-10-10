function [E, accepted] = hs_switch_chain(E, n, steps, seed)
%HS_SWITCH_CHAIN  Degree-preserving switch chain (uniform is stationary).
%   [E, accepted] = hs_switch_chain(E, n, steps, seed) runs steps Metropolis steps from the
%   simple hypergraph E (m-by-k, vertices 1..n): swap a vertex between two edges, accept if
%   the result is still simple. Irreducible for graphs (k = 2); not guaranteed for k > 2.
if nargin < 4 || isempty(seed), seed = uint64(randi(flintmax - 1)); end
[E, accepted] = hs_mex('switch_chain', E, n, steps, seed);
end
