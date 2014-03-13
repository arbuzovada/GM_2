function [G, ind] = ldpc_gen_matrix(H)
% This function evaluates distribution p(b | a, d)
% INPUT:
%    H: (n - k)-by-n binary array, parity-check matrix
%
% OUTPUT:
%    G: n-by-k binary array, generator matrix
%    ind: indices for I_k

    [m, n] = size(H);
    k = n - m;
    ind = [];
    cur_e = 1;
    for i = 1 : n
        if (H(cur_e, i) ~= 1)
            % find a lower row with a 1 in desired position
            found = 0;
            for j = (cur_e + 1) : m
                if (H(j, i) == 1)
                    % swap rows
                    tmp = H(j, :);
                    H(j, :) = H(cur_e, :);
                    H(cur_e, :) = tmp;
                    found = 1;
                    break;
                end
            end
            if ~found
                % skip this column
                continue;
            end
        end
        % zero everything higher and lower
        for j = 1 : m
            if (j ~= cur_e) && (H(j, i) == 1)
                H(j, :) = xor(H(j, :), H(cur_e, :));
            end
        end
        ind = [ind, i]; %#ok<AGROW>
        cur_e = cur_e + 1;
        if (cur_e > m)
            break;
        end
    end
    G = zeros(n, k);
    G(setdiff([1 : n], ind), :) = eye(k);
    G(ind, :) = H(:, setdiff([1 : n], ind));
    ind = setdiff([1 : n], ind);
end