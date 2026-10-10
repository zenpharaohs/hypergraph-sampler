function tf = hs_is_simple(E, n)
%HS_IS_SIMPLE  True if no edge repeats a vertex and no edge repeats (vertices 1..n).
tf = hs_mex('is_simple', E, n);
end
