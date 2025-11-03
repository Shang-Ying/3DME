function [Y_cov] = fct_UnCovMat_3D(pos_est, sill, range, nugget, model)
    if nargin < 5
        model = 'exp';
    end
    Duu = pdist2(pos_est, pos_est);
    switch lower(model)
        case 'exp'
            Y_cov = nugget + sill * exp(-Duu / range);
        case 'sph'
            h = Duu / range;
            Y_cov = sill * (1 - 1.5*h + 0.5*h.^3);
            Y_cov(h > 1) = 0;
            Y_cov = Y_cov + nugget * eye(size(Duu));
        % 可加更多模型
        otherwise
            error('Unknown model type.');
    end
end
