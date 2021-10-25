function problem = PrepareExchangeMatrix(A2,problem)

[v,c,rs,re] = convertToCSR(A2);
problem.exch_nval = int32(numel(v));
problem.exch_nrow = int32(numel(rs));
problem.exch_val  = single(v);
problem.exch_rows = int32(rs);
problem.exch_rowe = int32(re);
problem.exch_col  = int32(c);