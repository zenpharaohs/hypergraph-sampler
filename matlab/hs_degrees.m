function deg = hs_degrees(E, n)
%HS_DEGREES  Vertex degrees (n-by-1 int32) of the hypergraph E on vertices 1..n.
deg = hs_mex('degrees', E, n);
end
