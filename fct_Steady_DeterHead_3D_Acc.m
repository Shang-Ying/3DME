function [u] = fct_Steady_DeterHead_3D_Acc(Wcoef, idx, ns, tp, nrp, nlp, nnfp, nip, ntp, RBC, LBC, Bd_V)
num_eqs = nrp + nlp + nnfp + ntp;
num_vars = ntp;

nzmax = num_eqs * (ns + 1); % 預估非零元數量
rows = zeros(1, nzmax);
cols = zeros(1, nzmax);
vals = zeros(1, nzmax);
counter = 0;

% Dirichlet Boundary Conditions
for ith = 1:(nrp + nlp)
    counter = counter + 1;
    rows(counter) = ith;
    cols(counter) = ith;
    vals(counter) = 1;
end

% No-Flow Boundaries
for ith = (nrp + nlp + 1):(nrp + nlp + nnfp)
    idxNl = idx(ith, 1:ns);
    len = length(idxNl);
    rows(counter+1:counter+len) = ith;
    cols(counter+1:counter+len) = idxNl;
    vals(counter+1:counter+len) = ...
        Bd_V(ith, 4) .* Wcoef(1, 1:ns, ith) + ...
        Bd_V(ith, 5) .* Wcoef(2, 1:ns, ith) + ...
        Bd_V(ith, 6) .* Wcoef(3, 1:ns, ith);
    counter = counter + len;
end

% Governing Equations
for ith = 1:ntp
    eq_idx = ith + nrp + nlp + nnfp;
    idxGl = idx(ith, 1:ns);
    Yval = tp(idxGl, 4);

    dYdx = Wcoef(1, 1:ns, ith) * Yval;
    dYdy = Wcoef(2, 1:ns, ith) * Yval;
    dYdz = Wcoef(3, 1:ns, ith) * Yval;

    Cia = Wcoef(4, 1:ns, ith) + Wcoef(7, 1:ns, ith) + Wcoef(9, 1:ns, ith);
    Cib = dYdx * Wcoef(1, 1:ns, ith);
    Cic = dYdy * Wcoef(2, 1:ns, ith);
    Cid = dYdz * Wcoef(3, 1:ns, ith);

    len = length(idxGl);
    rows(counter+1:counter+len) = eq_idx;
    cols(counter+1:counter+len) = idxGl;
    vals(counter+1:counter+len) = Cia + Cib + Cic + Cid;
    counter = counter + len;
end

% 組裝 sparse 矩陣
C = sparse(rows(1:counter), cols(1:counter), vals(1:counter), num_eqs, num_vars);

% RHS 向量
f = zeros(num_eqs, 1);
f(1:nrp) = RBC;
f((nrp+1):(nrp+nlp)) = LBC;
% no-flow 與 governing eq 預設已為 0

% 解系統
u = C \ f;
end
