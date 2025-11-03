clear; close all; clc;
load("BotBound.dat")
load("TopBound.dat")
load("BotH.dat")
load("TopH.dat")
IndexB = inpolygon(BotH(:,1),BotH(:,2),BotBound(:,1),BotBound(:,2));
BotNodes = BotH(IndexB,:);
IndexT = inpolygon(TopH(:,1),TopH(:,2),TopBound(:,1),TopBound(:,2));
TopNodes = TopH(IndexT,:);
Bds = [BotNodes; TopNodes];
s = 1;
k = boundary(Bds,s);

box = [0,100,0,100,5,28];
ninit = [201,201];
dotmax = 5e6;
radius = 2.4;
tic
xyz = node_drop_3d_GPT(box, ninit, dotmax, radius);
t_np = toc;
Index = intriangulation(Bds,k,xyz);
xyz = xyz(Index,:);

% figure('Position',[0 0 1300 800]);
% scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 15, xyz(:,3), 'filled');
% hold on
% trisurf(k,Bds(:,1),Bds(:,2),Bds(:,3),'Facecolor','k', ...
%     'FaceAlpha',0.1,'LineStyle',":",'Edgecolor', [.5 .5 .5]);
% axis equal
% axis([0 100 0 100 0 30]);
% set(gca,'clim',[5 28]);
% view(8,35);
% fontsize(20,"points")
% xlabel('{\it x} (m)', 'FontSize', 20);
% ylabel('{\it y} (m)', 'FontSize', 20);
% zlabel('{\it z} (m)', 'FontSize', 20);
% ax = gca; ax.FontSize = 20;
% chb = colorbar;
% ylabel(chb, 'Height (m)','FontSize',20);
% hold off

% Triangulation vertex normal
s_xyz = 0.9; % Compact surface boundary
k_xyz = boundary(xyz, s_xyz);  % triangulation on full xyz
TR = triangulation(k_xyz, xyz);
V = vertexNormal(TR);         % 邊界法向量 (unit normal vectors)
nearB = any(V ~= 0, 2);       % logical index for boundary-adjacent points
Bdps_V = [xyz(nearB,:), V(nearB,:)];  % [x, y, z, nx, ny, nz]

% figure('Position',[0 0 1300 800]);
% trisurf(k,Bds(:,1),Bds(:,2),Bds(:,3),'Facecolor','k', ...
%     'FaceAlpha',0.1,'LineStyle',":",'Edgecolor', [.5 .5 .5]);
% axis equal
% axis([0 100 0 100 0 30]);
% set(gca,'clim',[5 28]);
% hold on
% scatter3(Bdps_V(:,1), Bdps_V(:,2), Bdps_V(:,3), 15, 'w', 'filled');
% % scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 15, 'filled');
% quiver3(Bdps_V(:,1),Bdps_V(:,2),Bdps_V(:,3),Bdps_V(:,4),Bdps_V(:,5),Bdps_V(:,6), ...
%     'AutoScale','on','Color',"b","LineWidth",1);
% view(8,35);
% fontsize(20,"points")
% xlabel('{\it x} (m)', 'FontSize', 20);
% ylabel('{\it y} (m)', 'FontSize', 20);
% zlabel('{\it z} (m)', 'FontSize', 20);
% ax = gca; ax.FontSize = 20;
% hold off

% Move boundary nodes upward to the top of the variable matrix (checked)
BdpsNum = ismember(xyz(:,1:3), Bdps_V(:,1:3), 'rows');
xyz_i = xyz;
xyz_i(BdpsNum,:) = [];
tp = [Bdps_V(:,1:3); xyz_i(:,1:3)];

% Right boundary nods set
R_pick = (Bdps_V(:,1) > 67) & (Bdps_V(:,2) < 35) & (Bdps_V(:,3) < 23) & (Bdps_V(:,3) > 11);
R_Bdps_V = Bdps_V(R_pick,:);

% Left boundary nods set
L_pick = (Bdps_V(:,1) < 50) & (Bdps_V(:,2) > 70) & (Bdps_V(:,3) < 23) & (Bdps_V(:,3) > 11);
L_Bdps_V = Bdps_V(L_pick,:);

% Rearrange boundary nodes
Bdps_V(L_pick | R_pick,:) = [];
tp = [R_Bdps_V(:,1:3); L_Bdps_V(:,1:3); Bdps_V(:,1:3); xyz_i(:,:)];
Bd_V = [R_Bdps_V(:,:); L_Bdps_V(:,:); Bdps_V(:,:)];


RBC = 325;
LBC = 320;
order = 2;
ns = 21;
tic
[Wcoef,idx] = fct_GFDM_Coef_3D(order,tp(:,1),tp(:,2),tp(:,3),ns);
nrp = size(R_Bdps_V,1);
nlp = size(L_Bdps_V,1);
nnfp = size(Bdps_V,1);
nip = size(xyz_i,1);
ntp = size(tp,1);

%% Monte Carlo Simulations
nreal = 10000;
u_all = cell(nreal,1);
tp_all = cell(nreal,1);
time_sgs = zeros(nreal,1);
time_solve = zeros(nreal,1);
time_total = zeros(nreal,1);
nreal_list = (1:nreal)';
t_sgs = tic;
V = '1 Sph(20)';
options.nsim = nreal;
options.max = 50;
options.mean = -5;
pos_known = []; val_known = []; pos_sim = tp(:,1:3);
Y_all = sgsim_3d_fast(pos_known, val_known, pos_sim, V, options);  % [ntp × nreal]
time_sgs = toc(t_sgs);
fprintf('SGS time for %d realizations: %.2f s\n', nreal, time_sgs);
parfor ireal = 1:nreal
    fprintf('Solving realization %d \n', ireal);
    tp_local = tp;
    tp_local(:,4) = Y_all(:, ireal);  % 每個 worker 對應一筆資料
    t_solve = tic;
    u = fct_Steady_DeterHead_3D_Acc(Wcoef, idx, ns, tp_local, ...
                                    nrp, nlp, nnfp, nip, ntp, RBC, LBC, Bd_V);
    time_solve(ireal) = toc(t_solve);
    time_total(ireal) = time_sgs / nreal + time_solve(ireal);
    u_all{ireal} = u;
    tp_all{ireal} = tp_local;
end
% %% Plot result
% figure;
% scatter3(tp(:,1), tp(:,2), tp(:,3), 15, u, 'filled');
% axis equal;
% view(8,35);
% colorbar; title('Head distribution (Heterogeneous)');
% xlabel('x'); ylabel('y'); zlabel('z');

%% Node distribution
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
hold on;
scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 15, xyz(:,3), 'filled');
trisurf(k, Bds(:,1), Bds(:,2), Bds(:,3), ...
    'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
axis equal tight;
xlabel('x (m)', 'FontSize', 28);
ylabel('y (m)', 'FontSize', 28);
zlabel('z (m)', 'FontSize', 28);
grid on; view(3);
camlight headlight; lighting gouraud;
colormap;
clim([5 28]);
cb = colorbar;
ylabel(cb, 'Height (m)', 'FontSize', 28);
cb.FontSize = 24;
ax = gca;
ax.FontSize = 28;
ax.FontName = 'Times';
set(findall(gcf,'type','text'), 'FontName', 'Times');

%% Boundary conditions
[~, idxTop] = ismember(round(TopNodes(:,1:3),6), round(Bds,6), 'rows');
[~, idxBot] = ismember(round(BotNodes(:,1:3),6), round(Bds,6), 'rows');
idxTop(idxTop == 0) = [];
idxBot(idxBot == 0) = [];
isTopFace = all(ismember(k, idxTop), 2);
isBotFace = all(ismember(k, idxBot), 2);
isSideFace = ~(isTopFace | isBotFace); 
faceCenter = mean(Bds(k,:), 2); 
v1 = Bds(k(:,1), :); 
v2 = Bds(k(:,2), :);
v3 = Bds(k(:,3), :);
faceCenter = (v1 + v2 + v3) / 3;  % [N_faces × 3]
R_face_mask = (faceCenter(:,1) > 67) & (faceCenter(:,2) < 35) & ...
              (faceCenter(:,3) < 23) & (faceCenter(:,3) > 11);
L_face_mask = (faceCenter(:,1) < 50) & (faceCenter(:,2) > 70) & ...
              (faceCenter(:,3) < 23) & (faceCenter(:,3) > 11);
R_faces = R_face_mask & isSideFace;
L_faces = L_face_mask & isSideFace;
gray_faces = ~(R_faces | L_faces);
% ==================
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
hold on;
trisurf(k(L_faces,:), Bds(:,1), Bds(:,2), Bds(:,3), ...
    'FaceColor', [0.25 0.5 0.95], ...
    'FaceAlpha', 0.9, ...
    'EdgeColor', 'none', ...
    'FaceLighting','gouraud', ...
    'SpecularStrength', 0.3, ...
    'DisplayName', 'Left BC (320 m)');
trisurf(k(R_faces,:), Bds(:,1), Bds(:,2), Bds(:,3), ...
    'FaceColor', [0.95 0.4 0.4], ...
    'FaceAlpha', 0.9, ...
    'EdgeColor', 'none', ...
    'FaceLighting','gouraud', ...
    'SpecularStrength', 0.3, ...
    'DisplayName', 'Right BC (325 m)');
trisurf(k(gray_faces,:), Bds(:,1), Bds(:,2), Bds(:,3), ...
    'FaceColor', [0.8 0.8 0.8], ...
    'FaceAlpha', 0.25, ...
    'EdgeColor', 'none', ...
    'FaceLighting','gouraud', ...
    'SpecularStrength', 0.1, ...
    'DisplayName', 'No Flow');
axis equal tight;
axis([0 100 0 100]);
view(3); grid on;
camlight headlight; lighting gouraud;
xlabel('x (m)', 'FontSize', 28);
ylabel('y (m)', 'FontSize', 28);
zlabel('z (m)', 'FontSize', 28);
legend({'Left BC (320 m)', 'Right BC (325 m)', 'No Flow'}, ...
       'FontSize', 24, 'Location', 'northeast');
ax = gca;
ax.FontSize = 28;
ax.FontName = 'Times';
set(findall(gcf,'type','text'), 'FontName', 'Times');


%% Plot isosurfaces
F = scatteredInterpolant(tp(:,1), tp(:,2), tp(:,3), u_all{3}, 'natural', 'none');
[xq, yq, zq] = meshgrid(...
    linspace(min(tp(:,1)), max(tp(:,1)), 100), ...
    linspace(min(tp(:,2)), max(tp(:,2)), 100), ...
    linspace(min(tp(:,3)), max(tp(:,3)), 80));
vq = F(xq, yq, zq);
shp = alphaShape(tp(:,1), tp(:,2), tp(:,3), 2.5); % radius
in_alpha = inShape(shp, xq, yq, zq);
vq(~in_alpha) = NaN;
u_min = 320; u_max = 325;
iso_vals = linspace(u_min, u_max, 11);
cmap = othercolor('GnBu7', length(iso_vals)-2);

scrsz = get(0,'ScreenSize');
current = figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
hold on;

for i = 2:(length(iso_vals)-1)
    val = iso_vals(i);
    try
        p = patch(isosurface(xq, yq, zq, vq, val));
        isonormals(xq, yq, zq, vq, p);
        set(p, 'FaceColor', cmap(i-1,:), ...
               'EdgeColor', 'none', ...
               'FaceAlpha', 0.5);
    catch ME
        warning("Isosurface at value %.2f failed: %s", val, ME.message);
    end
end
trisurf(k, Bds(:,1), Bds(:,2), Bds(:,3), ...
    'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.25, ...
    'EdgeColor', 'none');
hold off;
axis equal tight;
xlabel('x (m)', 'FontSize', 28);
ylabel('y (m)', 'FontSize', 28);
zlabel('z (m)', 'FontSize', 28);
% title('Head Isosurfaces with Boundary', 'FontSize', 30);
grid on;
view(3); 
camlight headlight; lighting gouraud;
colormap(othercolor('GnBu7'));
clim([u_min, u_max]);
ch = colorbar;
ylabel(ch, 'Head (m)', 'FontSize', 28);
ch.Ticks = u_min:u_max;
ax = gca;
ax.FontSize = 28;
ax.FontName = 'Times';
set(findall(gcf,'type','text'), 'FontName', 'Times');

%% Ensemble Mean and Variance
u_mat = cell2mat(cellfun(@(x) x(:), u_all, 'UniformOutput', false));
u_mat = reshape(u_mat, [], nreal);
u_mean = mean(u_mat, 2);
u_var  = var(u_mat, 0, 2);

%% Ensemble Mean Isosurfaces
Fmean = scatteredInterpolant(tp(:,1), tp(:,2), tp(:,3), u_mean, 'natural', 'none');
[xq, yq, zq] = meshgrid(...
    linspace(min(tp(:,1)), max(tp(:,1)), 100), ...
    linspace(min(tp(:,2)), max(tp(:,2)), 100), ...
    linspace(min(tp(:,3)), max(tp(:,3)), 80));
vq = Fmean(xq, yq, zq);
shp = alphaShape(tp(:,1), tp(:,2), tp(:,3), 2.5);
in_alpha = inShape(shp, xq, yq, zq);
vq(~in_alpha) = NaN;

u_min = 320; u_max = 325;
iso_vals = linspace(u_min, u_max, 11);
cmap = othercolor('GnBu7', length(iso_vals)-2);

scrsz = get(0,'ScreenSize');
figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
hold on;
for i = 2:(length(iso_vals)-1)
    val = iso_vals(i);
    try
        p = patch(isosurface(xq, yq, zq, vq, val));
        isonormals(xq, yq, zq, vq, p);
        set(p, 'FaceColor', cmap(i-1,:), ...
               'EdgeColor', 'none', ...
               'FaceAlpha', 0.5);
    catch ME
        warning("Mean isosurface at %.2f failed: %s", val, ME.message);
    end
end
trisurf(k, Bds(:,1), Bds(:,2), Bds(:,3), ...
    'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.25, ...
    'EdgeColor', 'none');
hold off;
axis equal tight;
xlabel('x (m)', 'FontSize', 28);
ylabel('y (m)', 'FontSize', 28);
zlabel('z (m)', 'FontSize', 28);
grid on; view(3);
camlight headlight; lighting gouraud;
colormap(othercolor('GnBu7'));
clim([u_min, u_max]);
ch = colorbar;
ylabel(ch, 'Head (m)', 'FontSize', 28);
ch.Ticks = u_min:u_max;
ax = gca; ax.FontSize = 28; ax.FontName = 'Times';
set(findall(gcf,'type','text'), 'FontName', 'Times');


%% Ensemble Variance Isosurfaces
Fvar = scatteredInterpolant(tp(:,1), tp(:,2), tp(:,3), u_var, 'natural', 'none');
vq = Fvar(xq, yq, zq);
vq(~in_alpha) = NaN;

v_min = 0; v_max = 0.10;
iso_vals = linspace(v_min, v_max, 11);
cmap = othercolor('RdPu9', length(iso_vals)-2);

figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
hold on;
for i = 2:(length(iso_vals)-1)
    val = iso_vals(i);
    try
        p = patch(isosurface(xq, yq, zq, vq, val));
        isonormals(xq, yq, zq, vq, p);
        set(p, 'FaceColor', cmap(i-1,:), ...
               'EdgeColor', 'none', ...
               'FaceAlpha', 0.5);
    catch ME
        warning("Variance isosurface at %.3f failed: %s", val, ME.message);
    end
end
trisurf(k, Bds(:,1), Bds(:,2), Bds(:,3), ...
    'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.25, ...
    'EdgeColor', 'none');
hold off;
axis equal tight;
xlabel('x (m)', 'FontSize', 28);
ylabel('y (m)', 'FontSize', 28);
zlabel('z (m)', 'FontSize', 28);
grid on; view(3);
camlight headlight; lighting gouraud;
colormap(othercolor('RdPu9'));
clim([v_min, v_max]);
ch = colorbar;
ylabel(ch, 'Head Variance (m²)', 'FontSize', 28);
ch.Ticks = round(linspace(v_min, v_max, 6), 3);
ax = gca; ax.FontSize = 28; ax.FontName = 'Times';
set(findall(gcf,'type','text'), 'FontName', 'Times');

%% CPU time - Cumulative SGS, PDE Solving, and Total
scrsz = get(0,'ScreenSize');
figure('OuterPosition', [0 0 (scrsz(4)+100) scrsz(4)]);
nreal_list = 1:nreal;
time_sgs_each = repmat(time_sgs / nreal, nreal, 1);
cumulative_sgs_time = cumsum(time_sgs_each);
cumulative_solve_time = cumsum(time_solve);
cumulative_total_time = cumulative_sgs_time + cumulative_solve_time;
plot(nreal_list, cumulative_total_time, '-', ...
    'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Total Time'); hold on;
plot(nreal_list, cumulative_sgs_time, '--', ...
    'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'SGS Time');
plot(nreal_list, cumulative_solve_time, ':', ...
    'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Solving PDE Time');
xlabel('Number of Realizations', 'FontSize', 28);
ylabel('Cumulative Time (s)', 'FontSize', 28);
legend('Location', 'northwest');
grid on;
ax = gca;
ax.FontSize = 28;
ax.FontName = 'Times';
set(findall(gcf,'type','text'), 'FontName', 'Times');

%% Sill check
Y_centered = Y_all - mean(Y_all, 2);
Y_cov_MCS = (Y_centered * Y_centered') / (size(Y_all, 2) - 1);
marginal_var = diag(Y_cov_MCS);
fprintf('Mean variance: %.4f\n', mean(marginal_var));
fprintf('Variance range: %.4f ~ %.4f\n', min(marginal_var), max(marginal_var));

%% Range check
n_sample_pairs = 2000;  % 減少計算成本
max_lag = 50;             % 最遠距離（依據你網格大小調整）
n_lags = 30;              % 分成幾個距離區間

idx_i = randi(size(Y_all,1), n_sample_pairs, 1);
idx_j = randi(size(Y_all,1), n_sample_pairs, 1);

diff_sq = mean((Y_all(idx_i,:) - Y_all(idx_j,:)).^2, 2); 
pos_est = tp(:,1:3);

dist_ij = sqrt(sum((pos_est(idx_i,:) - pos_est(idx_j,:)).^2, 2));

edges = linspace(0, max_lag, n_lags+1);
gamma_emp = zeros(n_lags, 1);
h_centers = zeros(n_lags, 1);

for k = 1:n_lags
    in_bin = (dist_ij >= edges(k)) & (dist_ij < edges(k+1));
    gamma_emp(k) = 0.5 * mean(diff_sq(in_bin));
    h_centers(k) = mean(dist_ij(in_bin));
end

% γ(h) = sill*(1.5 h/a - 0.5 (h/a)^3)
fitfun = @(params, h) params(1)*(1.5*(h/params(2)) - 0.5*(h/params(2)).^3).*(h <= params(2)) + params(1)*(h > params(2));
opts = optimset('Display','off');
params0 = [2.67, 20];  % 初始 guess：sill=1, range=20
valid_idx = ~isnan(gamma_emp);
params_fit = lsqcurvefit(fitfun, params0, h_centers(valid_idx), gamma_emp(valid_idx), [], [], opts);

figure;
scatter(h_centers, gamma_emp, 'bo'); hold on;
h_fit = linspace(0, max_lag, 200);
plot(h_fit, fitfun(params_fit, h_fit), 'r', 'LineWidth', 2);
xlabel('Lag distance h'); ylabel('Empirical variogram \gamma(h)');
legend('Empirical', sprintf('Fitted Sph: sill=%.2f, range=%.2f', params_fit));
title('Empirical variogram vs. fitted spherical model');
grid on;