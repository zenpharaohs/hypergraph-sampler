function test_hypersample()
%TEST_HYPERSAMPLE  Checks for the MATLAB/Octave binding (run after build_hypersample).
%   Statistical exactness is tested on the Python side against full enumeration; here we
%   check the binding: 1-based vertices, shapes, degrees, simplicity, seeds and errors.
E = hs_config_regular(200, 3, 12, 1);
assert(isequal(size(E), [800 3]) && isa(E, 'int32'));
assert(min(E(:)) >= 1 && max(E(:)) <= 200);
assert(hs_is_simple(E, 200));
assert(all(hs_degrees(E, 200) == 12));
assert(isequal(E, hs_config_regular(200, 3, 12, 1)));             % reproducible
assert(isequal(E, hs_config_regular(200, 3, 12, uint64(1))));     % uint64 seed
[~, t] = hs_config_regular(200, 3, 12, 2); assert(t >= 1);

G = hs_config_degrees([3 3 2 2 2 2], 2, 5);
assert(isequal(sort(double(hs_degrees(G, 6)))', [2 2 2 2 3 3]));

[F, acc] = hs_switch_chain(E, 200, 20000, 3);
assert(acc > 0 && hs_is_simple(F, 200) && all(hs_degrees(F, 200) == 12));
[F2, acc2] = hs_switch_chain(double(E), 200, 20000, 3);          % double input accepted
assert(isequal(F, F2) && acc == acc2);

[S, lw] = hs_sequential([12 12 12 12, 2 * ones(1, 44)], 2, 7);
assert(size(S, 2) == 2 && (isfinite(lw) || lw == -Inf));

P = hs_gnp(30, 0.2, 3, 11);
assert(size(P, 2) == 3 && hs_is_simple(P, 30) && isequal(P, hs_gnp(30, 0.2, 3, 11)));
assert(isempty(hs_gnp(10, 0, 2, 1)) && size(hs_gnp(10, 1, 2, 1), 1) == 45);
M = hs_gnm(20, 50, 2, 4);
assert(isequal(size(M), [50 2]) && hs_is_simple(M, 20));
assert(~hs_is_simple(int32([1 2; 2 1]), 3));

threw = false;
try, hs_config_regular(8, 3, 12, 1, 1e4); catch err, threw = strcmp(err.identifier, 'hypersample:saturation'); end
assert(threw, 'saturation should raise hypersample:saturation');
fprintf('test_hypersample: all checks passed (hypersample %s)\n', hs_mex('version'));
end
