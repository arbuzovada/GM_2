% n = 16;
% k = 4;
k = 10;
n = 50;
q = 0.01;
m = n - k;

H = make_ldpc_mex(m, n, 3);
% H = [1 1 1 1 0 ; 0 0 1 1 0; 1 0 0 1 1];
rank(full(H));
assert(rank(full(H)) == m)
[G, ind] = ldpc_gen_matrix(H);
% mod(H * G, 2)

% gen messages:
nmes = 1;
U = randi(2, k, nmes) - 1;
% encode:
V = mod(G * U, 2);
% make errors:
W = V;
q = 0.01;
E = binornd(1, q, size(W));
% W = xor(V, E);
% decode:
V1 = W;
success = 0;
fail = 0;
for i = 1 : nmes
    [e, status] = ldpc_decoding(mod(H * W(:, i), 2), H, q, 'display', ...
        'false', 'lambda', 1);
    V1(:, i) = xor(V1(:, i), e);
    if (V1(:, i) == V(:, i))
        fprintf('%d message: success!\n', i);
        success = success + 1;
    else
        fprintf('%d message: fail! Status: %d\n', i, status);
        fail = fail + 1;
    end
end
success
fail