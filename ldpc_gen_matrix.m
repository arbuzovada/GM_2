function [G, ind] = ldpc_gen_matrix(H)
% This function produces generator matrix G by given parity-check matrix H.
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
        % find a lower row with a 1 in desired position
        nonzero_row = find(H(cur_e : end, i), 1);
        if isempty(nonzero_row)
            % skip this column
            continue;
        end
        % true index
        nonzero_row = nonzero_row + cur_e - 1;
        % swap rows
        tmp = H(nonzero_row, :);
        H(nonzero_row, :) = H(cur_e, :);
        H(cur_e, :) = tmp;
        % find other 1s
        nonzero_rows = find(H(:, i));
        nonzero_rows = setdiff(nonzero_rows, cur_e);
        % zero everything higher and lower
        H(nonzero_rows, :) = xor(H(nonzero_rows, :), ...
            repmat(H(cur_e, :), size(nonzero_rows, 1), 1));
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