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
%               1, if beliefs converged;
%               2, if max_iter reached

    [m, n] = size(H);
    
    PARALLEL = true;
    LAMBDA = 1;
    MAX_ITER = 200;
    EPS = 1e-4;
    DISPLAY = false;
    for i = 1 : length(varargin)
        if strcmp(varargin{i}, 'schedule')
            PARALLEL = strcmp(varargin{i + 1}, 'parallel');
        end
        if strcmp(varargin{i}, 'damping')
            LAMBDA = varargin{i + 1};
        end
        if strcmp(varargin{i}, 'max_iter')
            MAX_ITER = varargin{i + 1};
        end
        if strcmp(varargin{i}, 'eps')
            EPS = varargin{i + 1};
        end
        if strcmp(varargin{i}, 'display')
            DISPLAY = strcmp(varargin{i + 1}, 'true');
        end
    end
    
    % initialization
    mu_vf = q * ones(n, m, 2); % messages from vertice to factor
    mu_vf(:, :, 1) = 1 - mu_vf(:, :, 1);
    mu_fv = zeros(m, n); % messages from factor to vertice for e_i = 1
    b = zeros(n, 2); % beliefs
    b_old = zeros(n, 2);
    
    for t = 1 : MAX_ITER
        if DISPLAY
            fprintf('Iteration %d\n', t);
        end

        % factor to vertice recalculation
        for j = 1 : m
            N_j = find(H(j, :));
            mask = (eye(length(N_j)) == 0);
            cur = repmat(N_j, length(N_j), 1)';
            inds = cur(mask);
            delta = prod(1 - 2 * reshape(mu_vf(inds, j, 2), ...
                length(N_j) - 1, length(N_j)), 1);
            mu_fv(j, N_j) = LAMBDA * (1 + (2 * s(j) - 1) * delta) / 2 + ...
                (1 - LAMBDA) * mu_fv(j, N_j);
        end
        
        % vertice to factor and beliefs recalculation
        for i = 1 : n
            N_i = find(H(:, i));
            mask = (eye(length(N_i)) == 0);
            cur = repmat(N_i, 1, length(N_i));
            inds = cur(mask);
            mu_vf(i, N_i, 1) = LAMBDA * (1 - q) * ...
                prod(1 - reshape(mu_fv(inds, i), length(N_i) - 1, ...
                length(N_i)), 1) + (1 - LAMBDA) * mu_vf(i, N_i, 1);
            mu_vf(i, N_i, 2) = LAMBDA * q * ...
                prod(reshape(mu_fv(inds, i), length(N_i) - 1, ...
                length(N_i)), 1) + (1 - LAMBDA) * mu_vf(i, N_i, 2);
            b(i, 1) = (1 - q) * prod(1 - mu_fv(N_i, i));
            b(i, 2) = q * prod(mu_fv(N_i, i));
        end
        
        % normalize
        mu_vf = mu_vf ./ (repmat(mu_vf(:, :, 1) + mu_vf(:, :, 2), ...
            [1, 1, 2]));
        
        [~, e] = max(b, [], 2);
        e = e - 1;
        % check stopping criteria
        if (mod(H * e, 2) == s)
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