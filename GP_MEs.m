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
s_xyz = 0.9; % Setting s to 0 gives the convex hull, and setting s to 1 gives a compact boundary that envelops the points.
k_xyz = boundary(xyz,s_xyz);
TR = triangulation(k_xyz,xyz);
V = vertexNormal(TR); % normal vectors
nearB = V ~= 0;
Bdps_V = [xyz(nearB(:,1),:) V(nearB(:,1),:)];

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

%% Moment Equations
t_me_total = tic;
Y_mean = -5 * ones(ntp, 1);
Y_cov = fct_UnCovMat_3D(tp(:,1:3), 1, 20, 0, 'sph');
% load("Y_cov_MCS.mat");
% Y_cov = Y_cov_MCS;
tp_cov = [tp(:,1:3), Y_mean, Y_cov];
tp_tvV = [tp(:,1:3) zeros(size(tp,1),5)];
tp_tvV(1:nrp,4) = 1; tp_tvV(1:nrp,5) = RBC;
tp_tvV(nrp+1:nrp+nlp,4) = 1; tp_tvV(nrp+1:nrp+nlp,5) = LBC;
tp_tvV(nrp+nlp+1:nrp+nlp+nnfp,4) = 2; tp_tvV(nrp+nlp+1:nrp+nlp+nnfp,5) = 0;
tp_tvV(nrp+nlp+nnfp+1:nrp+nlp+nnfp+nip,4) = 0; tp_tvV(nrp+nlp+nnfp+1:nrp+nlp+nnfp+nip,5) = 0;
tp_tvV(1:nrp+nlp+nnfp,6:8) = Bd_V(:,4:6);
G_In = (tp_tvV(:,4) == 0);
D_Bd = (tp_tvV(:,4) == 1);
N_Bd = (tp_tvV(:,4) == 2);
G_loop = find(G_In);
D_loop = find(D_Bd);
N_loop = find(N_Bd);
t_me_1st = tic;
u_mean_ME = fct_Steady_1stMEHead_3D(Wcoef, idx, ns, tp_cov, tp_tvV, nrp, nlp, nnfp, nip, ntp, RBC, LBC);
time_ME_1st = toc(t_me_1st);
t_me_2nd = tic;
[CYY, CYu, Cuu] = fct_Steady_2ndMEHead_3D(Wcoef, idx, ns, tp_cov, tp_tvV, nrp, nlp, nnfp, nip, ntp, u_mean_ME);
time_ME_2nd = toc(t_me_2nd);
time_ME_total = toc(t_me_total);
u_var_ME = diag(Cuu);
fprintf('Moment Equation timing:\n');
fprintf('  1st moment solve time  : %.2f s\n', time_ME_1st);
fprintf('  2nd moment solve time  : %.2f s\n', time_ME_2nd);
fprintf('  Total ME computation   : %.2f s\n', time_ME_total);

%% Ensemble Mean Isosurfaces
Fmean = scatteredInterpolant(tp(:,1), tp(:,2), tp(:,3), u_mean_ME, 'natural', 'none');
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
Fvar = scatteredInterpolant(tp(:,1), tp(:,2), tp(:,3), u_var_ME, 'natural', 'none');
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

% %% CPU time
% scrsz = get(0,'ScreenSize');
% figure('OuterPosition',[0 0 (scrsz(4)+100) scrsz(4)]);
% time_sgs_each = repmat(time_sgs / nreal, nreal, 1);
% time_total_each = time_sgs_each + time_solve;
% cumulative_time = cumsum(time_total_each);
% plot(nreal_list, cumulative_time, '-', ...
%     'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Improved SGS + Solving PDE');
% 
% xlabel('Number of Realizations', 'FontSize', 28);
% ylabel('Time (s)', 'FontSize', 28);
% % legend('Location', 'northwest');
% grid on;
% ax = gca;
% ax.FontSize = 28;
% ax.FontName = 'Times';
% set(findall(gcf,'type','text'), 'FontName', 'Times');

%% Post assessment
% load("Results_ME.mat")
load("Results_MC10000.mat")

mae_u_mean_ME = mean(abs(u_mean_ME - u_mean));
mae_u_var_ME = mean(abs(u_var_ME - u_var));

num_MC = 206;
u_all_matrix = cell2mat(u_all(1:num_MC));
u_all_matrix = reshape(u_all_matrix, [], num_MC)';  % reshape [206 x 5075]
u_mean_MC206 = mean(u_all_matrix, 1)'; 
u_var_MC206  = var(u_all_matrix, 0, 1)';  

mae_u_mean_MC206 = mean(abs(u_mean_MC206 - u_mean));
mae_u_var_MC206 = mean(abs(u_var_MC206 - u_var));

fprintf('MAE between u_mean_ME and u_mean: %.6f\n', mae_u_mean_ME);
fprintf('MAE between u_var_ME and u_var: %.6f\n', mae_u_var_ME);
fprintf('MAE between u_mean_MC206 and u_mean: %.6f\n', mae_u_mean_MC206);
fprintf('MAE between u_var_MC206 and u_var: %.6f\n', mae_u_var_MC206);