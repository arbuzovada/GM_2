function [e, status] = ldpc_decoding(s, H, q, varargin)
% This function finds errors in message
% INPUT:
%    s: m-by-1 binary array, syndrom
%    H: m-by-n binary array, parity-check matrix
%    q: double in (0; 0.5), the probability of bit inversion
%
% OUTPUT:
%    e: n-by-1 binary array, error vector
%    status: result of decoding
%               0, if OK;
%               1, if beliefs didn't converge;
%               2, if max_iter reached

    % s = He
    MAX_ITER = 1000;
    
    % initialization
    mu_vf = q * ones(n, m, 2); % messages from vertice to factor
    mu_vf(:, :, 1) = 1 - mu_vf(:, :, 1);
    mu_fv = zeros(m, n); % messages from factor to vertice for e_i = 1
    b = zeros(n, 2); % beliefs
    b_old = zeros(n, 2); % beliefs
    
    for t = 1 : MAX_ITER
        % factor to vertice recalculation
        for j = 1 : m
            N_j = find(H(j, :));
            for i = 1 : n
                inds = setdiff(N_j, i);
                p_k = mu_vf(inds(1), j);
                delta = 1 - 2 * mu_vf(inds(1), j);
                for k = inds(2 : end)
                    p_k = p_k * mu_vf(k, j) + ...
                        (1 - p_k) * (1 - mu_vf(k, j));
                    delta = delta * (1 - 2 * mu_vf(inds(1), j));
                end
                if ~s(j)
                    mu_fv(j, i) = (1 - delta) / 2;
                else
                    mu_fv(j, i) = (1 + delta) / 2;
                end
            end
        end
        % vertice to factor and beliefs recalculation
        for i = 1 : n
            N_i = find(H(:, i));
            for j = 1 : m
                mu_vf(i, j, 1) = (1 - q) * ...
                    prod(1 - mu_fv(setdiff(N_i, j), i));
                mu_vf(i, j, 2) = q * prod(mu_fv(setdiff(N_i, j), i));
            end
            b(i, 1) = (1 - q) * prod(1 - mu_fv(N_i, i));
            b(i, 2) = q * prod(mu_fv(N_i, i));
        end
        [~, e] = max(b, [], 2);
        e = e - 1;
        % check stopping criteria
        if (H * e == s)
            status = 0;
            return;
        end
        if (max(b - b_old) < EPS)
            status = 1;
            return;
        end
        b_old = b;
    end
    status = 2;
end