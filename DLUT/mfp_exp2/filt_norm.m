% [B,m]=filt_norm(A,eg,id,m0);
% This function filters out the near-zero elements in A to return a sparse
% matrix B by using the iterative algorithm 4.1 proposed in "High-performance computing of large
% sparse matrix exponential". 
% The matrix B satisfies norm(A-B,id)<=eg
% id = 1, 'inf', or 'fro'.
% m=norm(A-B,id)/ef, where ef is greater than the largest element of abs(A-B).
% m0 is an initial guess value to m.

% Wu Feng writed on 2021 11.29.