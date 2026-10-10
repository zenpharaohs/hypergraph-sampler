%DEMO_HYPERSAMPLE  Classroom demo: exact rejection vs weighted sequential sampling.
% Run build_hypersample once first.

% 1. Acceptance of the configuration model for d-regular k-uniform hypergraphs
%    approaches exp(-(k-1)(d-1)/2) as n grows (k >= 3; graphs also lose
%    repeated edges, giving exp(-(d-1)/2 - (d-1)^2/4)).
k = 3; ds = [2 4 6 8]; n = 600; reps = 20;
fprintf('d   mean trials   exp((k-1)(d-1)/2)\n');
for d = ds
    t = zeros(reps, 1);
    for r = 1:reps, [~, t(r)] = hs_config_regular(n, k, d, r); end
    fprintf('%d   %10.1f   %10.1f\n', d, mean(t), exp((k-1)*(d-1)/2));
end

% 2. Importance weights: four hubs of degree 12 among 44 vertices of degree 2.
%    Statistic: number of hub-hub edges. The raw sequential sampler is biased;
%    the self-normalized weighted estimate matches the exact sampler.
deg = [12 12 12 12, 2*ones(1, 44)];
hub = @(E) sum(all(E <= 4, 2));
R = 5000;
exact = zeros(R, 1);
for r = 1:R, exact(r) = hub(hs_config_degrees(deg, 2, r)); end
f = zeros(R, 1); lw = zeros(R, 1);
for r = 1:R, [E, lw(r)] = hs_sequential(deg, 2, 1e6 + r); f(r) = hub(E); end
ok = isfinite(lw); w = zeros(R, 1); w(ok) = exp(lw(ok) - max(lw(ok)));
fprintf('exact %.3f   raw sequential %.3f   weighted %.3f   (stuck %.0f%%)\n', ...
    mean(exact), mean(f(ok)), sum(w .* f) / sum(w), 100 * mean(~ok));
