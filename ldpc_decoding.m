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

    % s = He
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
    b_old = zeros(n, 2); % beliefs
    
    for t = 1 : MAX_ITER
        if DISPLAY
            fprintf('Iteration %d\n', t);
        end
%         fprintf('before: NaN sums %d, %d\n', sum(sum(isnan(mu_fv))), ...
%             sum(sum(sum(isnan(mu_vf)))));

        % factor to vertice recalculation
        for j = 1 : m
            N_j = find(H(j, :));
            for i = N_j
                inds = setdiff(N_j, i);
                if isempty(inds)
                    continue;
                end
                delta = prod(1 - 2 * mu_vf(inds, j, 2));
%                 p_k = mu_vf(inds(1), j);
%                 for k = inds
%                     p_k = p_k * mu_vf(k, j) + ...
%                         (1 - p_k) * (1 - mu_vf(k, j));
%                     delta = delta * (1 - 2 * mu_vf(k, j));
%                 end
                if ~s(j)
                    mu_fv(j, i) = LAMBDA * (1 - delta) / 2 + ...
                        (1 - LAMBDA) * mu_fv(j, i);
                else
                    mu_fv(j, i) = LAMBDA * (1 + delta) / 2 + ...
                        (1 - LAMBDA) * mu_fv(j, i);
                end
            end
        end
        % vertice to factor and beliefs recalculation
        for i = 1 : n
            N_i = find(H(:, i));
            for j = N_i'
                inds = setdiff(N_i, j);
                if isempty(inds)
                    continue;
                end
                mu_vf(i, j, 1) = LAMBDA * (1 - q) * ...
                    prod(1 - mu_fv(inds, i)) + ...
                    (1 - LAMBDA) * mu_vf(i, j, 1);
                mu_vf(i, j, 2) = LAMBDA * q * ...
                    prod(mu_fv(inds, i)) + ...
                    (1 - LAMBDA) * mu_vf(i, j, 2);
            end
            b(i, 1) = (1 - q) * prod(1 - mu_fv(N_i, i));
            b(i, 2) = q * prod(mu_fv(N_i, i));
        end
        % normalize
        mu_vf = mu_vf ./ (repmat(mu_vf(:, :, 1) + mu_vf(:, :, 2), ...
            [1, 1, 2]));
        [~, e] = max(b, [], 2);
        e = e - 1;
%         [mod(H * e, 2)'; s']
        % check stopping criteria
        if (mod(H * e, 2) == s)
            status = 0;
            return;
        end
        if (max(b - b_old) < EPS)
%             [b_old, b]
            status = 1;
            return;
        end
        b_old = b;
    end
    status = 2;
end