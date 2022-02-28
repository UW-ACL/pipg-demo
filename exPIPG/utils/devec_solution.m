function sol = devec_solution(pbm,Z,varargin)
    sol = struct;
    N = pbm.N;
    dt = pbm.dt;
    q = zeros(pbm.nx,N);
    u = zeros(pbm.nu,N-1);
    nxnu = pbm.nx + pbm.nu;
    for j = 1:N-1
        q(:,j) = Z((j-1)*nxnu+1:(j-1)*nxnu+pbm.nx);
        u(:,j) = Z((j-1)*nxnu+pbm.nx+1:j*nxnu);
    end
    q(:,N) = Z(end-pbm.nx+1:end);
    sol.t = 0:dt:(N-1)*dt;
    sol.q = q;
    sol.u = u;
    if nargin == 2
        sol.cost = 0.5*Z'*pbm.vec_pbm.P*Z + pbm.vec_pbm.p'*Z; 
    end
    sol.Z = Z;
end