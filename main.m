n = 16;
k = 4;
m = n - k;

H = make_ldpc_mex(m, n, 3);
% H = [1 1 1 1 0 ; 0 0 1 1 0; 1 0 0 1 1];
rank(full(H))
assert(rank(full(H)) == m)
[G, ind] = ldpc_gen_matrix(H);
H * G