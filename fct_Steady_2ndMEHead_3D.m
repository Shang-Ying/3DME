function [CYY, CYu, Cuu] = fct_Steady_2ndMEHead_3D(Wcoef, idx, ns, tp_cov, tp_tvV, nrp, nlp, nnfp, nip, ntp, u_mean)
% Initialize Outputs
num_pts = size(tp_cov, 1);
CYY = tp_cov(:, 5:end);
CYu = zeros(num_pts, num_pts);
Cuu = zeros(num_pts, num_pts);

% Extract Wcoef slices for efficiency
W1 = squeeze(Wcoef(1,:,:))';
W2 = squeeze(Wcoef(2,:,:))';
W3 = squeeze(Wcoef(3,:,:))';
W4 = squeeze(Wcoef(4,:,:))';
W7 = squeeze(Wcoef(7,:,:))';
W9 = squeeze(Wcoef(9,:,:))';

%% --- Progress bar setup for CYu ---
dq1 = parallel.pool.DataQueue;
p1 = waitbar(0, 'Computing CYu...');
afterEach(dq1, @(~) updateWaitbar(dq1, num_pts, p1));

parfor dx = 1:num_pts
    CYu(dx,:) = computeCYu_3D(dx, idx, ns, tp_cov, tp_tvV, u_mean, CYY, W1, W2, W3, W4, W7, W9, nrp, nlp, nnfp, nip, ntp);
    send(dq1, []);
end
close(p1);

%% --- Progress bar setup for Cuu ---
dq2 = parallel.pool.DataQueue;
p2 = waitbar(0, 'Computing Cuu...');
afterEach(dq2, @(~) updateWaitbar(dq2, num_pts, p2));

parfor dchi = 1:num_pts
    Cuu(:,dchi) = computeCuu_3D(dchi, idx, ns, tp_cov, tp_tvV, u_mean, CYu, W1, W2, W3, W4, W7, W9, nrp, nlp, nnfp, nip, ntp);
    send(dq2, []);
end
close(p2);
end

%% --- Update waitbar helper function ---
function updateWaitbar(~, total, pHandle)
    persistent count
    if isempty(count)
        count = 1;
    else
        count = count + 1;
    end
    waitbar(count / total, pHandle);
    if count == total
        count = [];
    end
end

%% function of CYu
function f_out = computeCYu_3D(dx, idx, ns, tp_cov, tp_tvV, u_mean, CYY, W1, W2, W3, W4, W7, W9, nrp, nlp, nnfp, nip, ntp)
num_eqs = nrp + nlp + nnfp + ntp;
num_vars = ntp;
nzmax = num_eqs * (ns + 1);
rows = zeros(1, nzmax);
cols = zeros(1, nzmax);
vals = zeros(1, nzmax);
counter = 0;
f_1 = zeros(num_eqs, 1);
for ith = 1:(nrp + nlp)
    counter = counter + 1;
    rows(counter) = ith;
    cols(counter) = ith;
    vals(counter) = 1;
end
for ith = (nrp + nlp + 1):(nrp + nlp + nnfp)
    idxNl = idx(ith, 1:ns);
    len = length(idxNl);
    rows(counter+1:counter+len) = ith;
    cols(counter+1:counter+len) = idxNl;
    vals(counter+1:counter+len) = ...
        tp_tvV(ith, 6) .* W1(ith,1:ns) + ...
        tp_tvV(ith, 7) .* W2(ith,1:ns) + ...
        tp_tvV(ith, 8) .* W3(ith,1:ns);
    counter = counter + len;
end
for ith = 1:ntp
    eq_idx = ith + nrp + nlp + nnfp;
    idxGl = idx(ith, 1:ns);
    W_1 = W1(ith,1:ns);
    W_2 = W2(ith,1:ns);
    W_3 = W3(ith,1:ns);
    W_4 = W4(ith,1:ns);
    W_7 = W7(ith,1:ns);
    W_9 = W9(ith,1:ns);
    CYYdchi = CYY(dx, idxGl);
    umidx = u_mean(idxGl, 1);
    dCYYdx = W_1 * CYYdchi';
    dCYYdy = W_2 * CYYdchi';
    dCYYdz = W_3 * CYYdchi';
    dudx_1 = W_1 * umidx;
    dudy_1 = W_2 * umidx;
    dudz_1 = W_3 * umidx;
    f_1(eq_idx) = -dudx_1 * dCYYdx - dudy_1 * dCYYdy - dudz_1 * dCYYdz;
    len = length(idxGl);
    rows(counter+1:counter+len) = eq_idx;
    cols(counter+1:counter+len) = idxGl;
    vals(counter+1:counter+len) = W_4 + W_7 + W_9;
    counter = counter + len;
end
C_1 = sparse(rows(1:counter), cols(1:counter), vals(1:counter), num_eqs, num_vars);
f_out = sparse(C_1 \ f_1);
end

%% function of Cuu
function f_out = computeCuu_3D(dchi, idx, ns, tp_cov, tp_tvV, u_mean, CYu, W1, W2, W3, W4, W7, W9, nrp, nlp, nnfp, nip, ntp)
num_eqs = nrp + nlp + nnfp + ntp;
num_vars = ntp;
nzmax = num_eqs * (ns + 1);
rows = zeros(1, nzmax);
cols = zeros(1, nzmax);
vals = zeros(1, nzmax);
counter = 0;
f_2 = zeros(num_eqs, 1);
for ith = 1:(nrp + nlp)
    counter = counter + 1;
    rows(counter) = ith;
    cols(counter) = ith;
    vals(counter) = 1;
end
for ith = (nrp + nlp + 1):(nrp + nlp + nnfp)
    idxNl = idx(ith, 1:ns);
    len = length(idxNl);
    rows(counter+1:counter+len) = ith;
    cols(counter+1:counter+len) = idxNl;
    vals(counter+1:counter+len) = ...
        tp_tvV(ith, 6) .* W1(ith,1:ns) + ...
        tp_tvV(ith, 7) .* W2(ith,1:ns) + ...
        tp_tvV(ith, 8) .* W3(ith,1:ns);
    counter = counter + len;
end
for ith = 1:ntp
    eq_idx = ith + nrp + nlp + nnfp;
    idxGl = idx(ith, 1:ns);
    W_1 = W1(ith,1:ns);
    W_2 = W2(ith,1:ns);
    W_3 = W3(ith,1:ns);
    W_4 = W4(ith,1:ns);
    W_7 = W7(ith,1:ns);
    W_9 = W9(ith,1:ns);
    CYudchi = CYu(idxGl, dchi);
    umidx = u_mean(idxGl, 1);
    dCYudx = W_1 * CYudchi;
    dCYudy = W_2 * CYudchi;
    dCYudz = W_3 * CYudchi;
    dudx_2 = W_1 * umidx;
    dudy_2 = W_2 * umidx;
    dudz_2 = W_3 * umidx;
    f_2(eq_idx) = -dudx_2 * dCYudx - dudy_2 * dCYudy - dudz_2 * dCYudz;
    len = length(idxGl);
    rows(counter+1:counter+len) = eq_idx;
    cols(counter+1:counter+len) = idxGl;
    vals(counter+1:counter+len) = W_4 + W_7 + W_9;
    counter = counter + len;
end
C_2 = sparse(rows(1:counter), cols(1:counter), vals(1:counter), num_eqs, num_vars);
f_out = sparse(C_2 \ f_2);
end
