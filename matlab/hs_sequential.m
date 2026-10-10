function [E, logw] = hs_sequential(deg, k, seed)
%HS_SEQUENTIAL  Sequential stub matching with an exact importance weight (k <= 4).
%   [E, logw] = hs_sequential(deg, k, seed). Each edge is chosen uniformly among the
%   currently valid stub k-sets (V_t of them); w = prod_t V_t, so for samples E_i with
%   weights w_i = exp(logw_i), sum(w_i f(E_i)) / sum(w_i) estimates the uniform average of f.
%   A stuck run returns logw = -Inf (weight 0) and must stay in the sum.
if nargin < 2 || isempty(k), k = 2; end
if nargin < 3 || isempty(seed), seed = uint64(randi(flintmax - 1)); end
[E, logw] = hs_mex('sequential', deg, k, seed);
end
