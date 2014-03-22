function [err_bit, err_block, diver] = ldpc_mc(H, q, num_points)
% This function performs Monte-Carlo
% INPUT:
%    H: m-by-n binary array, parity-check matrix
%    q: double in (0; 0.5), the probability of bit inversion
%    num_points: integer, number of experiments
%
% OUTPUT:
%    err_bit: probability of bit error
%    err_block: probability of block error
%    diver: probability of divergence

    [m, n] = size(H);
    k = n - m;
    
    G = ldpc_gen_matrix(H);
    % gen messages
    U = randi(2, k, num_points) - 1;
    % encode:
    V = mod(G * U, 2);
    % make errors:
    E = binornd(1, q, size(V));
    W = xor(V, E);
    % decode:
    err_bit = 0;
    err_block = 0;
    diver = 0;
    for i = 1 : num_points
        [e, status] = ldpc_decoding(mod(H * W(:, i), 2), H, q, ...
            'display', 'false', 'damping', 7/8, 'schedule', 'parallel', ...
            'max_iter', 200);
        err_bit = err_bit + sum(e ~= E(:, i));
        err_block = err_block + any(e ~= E(:, i));
        if status == 2
            diver = diver + 1;
        end
    end
    err_bit = err_bit / (n * num_points);
    err_block = err_block / num_points;
    diver = diver / num_points;
end