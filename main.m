n = 16;
k = 4;
m = n - k;

H = make_ldpc_mex(m, n, 3);
% H = [1 1 1 1 0 ; 0 0 1 1 0; 1 0 0 1 1];
rank(full(H))
assert(rank(full(H)) == m)
[G, ind] = ldpc_gen_matrix(H);
H * G

% gen messages:
nmes = 1;
U = randi(2, k, nmes) - 1;
% encode:
V = G * U;
% make errors:
W = V;
q = 0.1;
E = binornd(1, q, size(W));
W = xor(V, E);
% decode:
V1 = W;
success = 0;
fail = 0;
for i = 1 : nmes
    E(:, i)'
    [e, status] = ldpc_decoding(mod(H * W(:, i), 2), H, q, 'display', 'true');
    V1(:, i) = xor(V1(:, i), e);
    if (V1(:, i) == V(:, i))
        fprintf('%d message: success!', i);
        success = success + 1;
    else
        fprintf('%d message: fail! Status: %d\n', i, status);
        [E(:, i)'; e']
        fail = fail + 1;
    end
end