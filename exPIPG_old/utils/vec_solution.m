function [Z,cost] = vec_solution(pbm,q,u)
    Z = zeros(pbm.vec_pbm.len_Z,1);
    nx = pbm.nx;
    nu = pbm.nu;
    nxnu = nx+nu;
    for j = 1:pbm.N-1
        Z((j-1)*nxnu+1:(j-1)*nxnu+nx) = q(:,j);
        Z((j-1)*nxnu+nx+1:j*nxnu) = u(:,j);
    end
    Z((pbm.N-1)*nxnu+1:(pbm.N-1)*nxnu+nx) = q(:,pbm.N);
    cost = 0.5*Z'*pbm.vec_pbm.P*Z + pbm.vec_pbm.p'*Z;
end