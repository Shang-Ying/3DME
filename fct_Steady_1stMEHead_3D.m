function [u_mean] = fct_Steady_1stMEHead_3D(Wcoef, idx, ns, tp_cov, tp_tvV, nrp, nlp, nnfp, nip, ntp, RBC, LBC)
num_eqs = nrp + nlp + nnfp + ntp;
num_vars = ntp;

nzmax = num_eqs * (ns + 1);
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
        tp_tvV(ith, 6) .* Wcoef(1, 1:ns, ith) + ...
        tp_tvV(ith, 7) .* Wcoef(2, 1:ns, ith) + ...
        tp_tvV(ith, 8) .* Wcoef(3, 1:ns, ith);
    counter = counter + len;
end

% Governing Equations
for ith = 1:ntp
    eq_idx = ith + nrp + nlp + nnfp;
    idxGl = idx(ith, 1:ns);
    Yval = tp_cov(idxGl, 4);

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

C = sparse(rows(1:counter), cols(1:counter), vals(1:counter), num_eqs, num_vars);

% RHS
f = zeros(num_eqs, 1);
f(1:nrp) = RBC;
f((nrp+1):(nrp+nlp)) = LBC;

u_mean = C \ f;
end
