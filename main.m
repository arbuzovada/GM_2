H = make_ldpc_mex(12, 16, 3);
% H = [1 1 1 1 0 ; 0 0 1 1 0; 1 0 0 1 1];
[G, ind] = ldpc_gen_matrix(H);
H * G