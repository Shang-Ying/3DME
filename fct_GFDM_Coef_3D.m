function [Wcoef, idx] = fct_GFDM_Coef_3D(order, x, y, z, ns)
    n = length(x);    
    xyz = [x, y, z];

    [idx, d] = knnsearch(xyz, xyz, 'k', ns + 1);
    Wcoef = zeros(9, ns + 1, n);

    parfor i = 1:n
        tpx = (x(idx(i, 2:ns+1)) - x(i))';
        tpy = (y(idx(i, 2:ns+1)) - y(i))';
        tpz = (z(idx(i, 2:ns+1)) - z(i))';
        dd = d(i, 2:ns+1) / d(i, ns+1);
        W_func = 1 - 6 * (dd.^2) + 8 * (dd.^3) - 3 * (dd.^4);
        W2 = diag(W_func.^2);

        if order == 2
            P = [tpx; tpy; tpz; (tpx.^2) / 2; tpx .* tpy; tpx .* tpz; (tpy.^2) / 2; tpy .* tpz; (tpz.^2) / 2];
        else
            % Handle other orders if necessary
            continue;
        end

        % Calculate Wcoef for current node
        Wcoef_tp = pinv(P*W2*(P)')*([-sum(P*W2,2),P*W2]);

        % Store result in Wcoef
        Wcoef(:, :, i) = Wcoef_tp;
    end
end
