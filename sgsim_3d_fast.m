function [sim_mul] = sgsim_3d_fast(pos_known, val_known, pos_sim, V, options)
% Fast 3D Sequential Gaussian Simulation with knnsearch, vectorized kriging, and parfor
% Supports multiple variogram types
% Call:
%   [sim_mul] = sgsim_3d_fast(pos_known, val_known, pos_sim, V, options);

if ischar(V)
    V = deformat_variogram(V);
end

if ~isfield(options, 'nsim'); options.nsim = 1; end
if ~isfield(options, 'max'); options.max = 10; end

if isempty(val_known)
    default_mean = 0;
else
    default_mean = mean(val_known(:,1));
end

if ~isfield(options, 'mean')
    options.mean = default_mean;
end

n_sim = size(pos_sim, 1);
n_known = size(pos_known, 1);

% Preallocate output
sim_mul = NaN(n_sim, options.nsim);

parfor j = 1:options.nsim
    fprintf('Realization %d/%d\n', j, options.nsim);

    path = randperm(n_sim);
    d_sim = NaN(n_sim, 1);
    pos_done = zeros(n_sim, 3);
    val_done = zeros(n_sim, 1);
    i_use_cond = true(n_known, 1);

    for idx = 1:n_sim
        i = path(idx);
        sim_point = pos_sim(i, :);

        % Gather known + simulated data
        if n_known == 0
            pos_all = pos_done(1:idx-1, :);
            val_all = val_done(1:idx-1);
        else
            pos_all = [pos_done(1:idx-1, :); pos_known(i_use_cond, :)];
            val_all = [val_done(1:idx-1); val_known(i_use_cond, 1)];
        end

        % Select nearest neighbors
        if size(pos_all, 1) > options.max
            [nn_idx, ~] = knnsearch(pos_all, sim_point, 'K', options.max);
            pos_used = pos_all(nn_idx, :);
            val_used = val_all(nn_idx);
        else
            pos_used = pos_all;
            val_used = val_all;
        end

        % Skip to next if no data yet
        if isempty(val_used)
            mean_est = options.mean;
            var_est = sum([V.par1]);
            d_sim(i) = norminv(rand, mean_est, sqrt(var_est));
            pos_done(idx,:) = sim_point;
            val_done(idx) = d_sim(i);
            continue
        end

        % Compute distances
        dists = vecnorm(pos_used - sim_point, 2, 2);
        gamma = zeros(size(dists));
        for v = 1:length(V)
            range = V(v).par2; c0 = V(v).par1; h = dists;
            switch V(v).itype
                case 0, gamma = gamma + c0;
                case 1
                    idx_r = h <= range;
                    g = zeros(size(h));
                    g(idx_r) = c0 * (1.5 * h(idx_r)/range - 0.5 * (h(idx_r)/range).^3);
                    g(~idx_r) = c0;
                    gamma = gamma + g;
                case 2, gamma = gamma + c0 * (1 - exp(-3 * h / range));
                case 3, gamma = gamma + c0 * (1 - exp(-3 * (h / range).^2));
                case 4, gamma = gamma + c0 * h.^range;
                case 6, gamma = gamma + c0 * h;
                otherwise, gamma = gamma + c0;
            end
        end

        % Kriging matrix
        D = pdist2(pos_used, pos_used);
        G = zeros(size(D));
        for v = 1:length(V)
            range = V(v).par2; c0 = V(v).par1; h = D;
            switch V(v).itype
                case 0, G = G + c0;
                case 1
                    idx_r = h <= range;
                    g = zeros(size(h));
                    g(idx_r) = c0 * (1.5 * h(idx_r)/range - 0.5 * (h(idx_r)/range).^3);
                    g(~idx_r) = c0;
                    G = G + g;
                case 2, G = G + c0 * (1 - exp(-3 * h / range));
                case 3, G = G + c0 * (1 - exp(-3 * (h / range).^2));
                case 4, G = G + c0 * h.^range;
                case 6, G = G + c0 * h;
                otherwise, G = G + c0;
            end
        end

        % Solve kriging system
        n = size(G, 1);
        G_full = [G, ones(n, 1); ones(1, n), 0];
        gamma_full = [gamma; 1];
        weights = G_full \ gamma_full;

        mean_est = weights(1:n)' * val_used;
        var_est = sum([V.par1]) - weights(1:n)' * gamma;
        std_est = sqrt(max(var_est, 0));
        d_sim(i) = norminv(rand, mean_est, std_est);

        % Update known set
        pos_done(idx,:) = sim_point;
        val_done(idx) = d_sim(i);

        if n_known > 0
            match_idx = find(all(abs(pos_known - sim_point) < 1e-6, 2));
            i_use_cond(match_idx) = false;
        end
    end

    sim_mul(:, j) = d_sim;
end
end
