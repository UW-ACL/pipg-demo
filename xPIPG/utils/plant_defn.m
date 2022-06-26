function pp = plant_defn(pp,flag)
%{
05/10/2022
Purnanand Elango

Define the discrete-time plant model 
Input:
    flag:
        0 - Random neutrally stable system
        1 - System of oscillating masses
%}

    switch flag
        case 0
            % Neutrally stable random system
            pp.A = randn(pp.n,pp.n);
            if max(abs(eig(pp.A)))>=1
                pp.A = pp.A/max(abs(eig(pp.A)));
            end
            pp.B = randn(pp.n,pp.m);
        case 1
            assert(abs(pp.m-pp.n/2)<=1e-8,"In the oscillating masses system the number of control inputs should be half the number of states.")
            nb2 = int64(pp.n/2);
            trid_pattern = zeros(1,nb2); trid_pattern(1:2) = [-2 1];
            % Continuous-time system
            Ac = [zeros(nb2)              eye(nb2);
                  toeplitz(trid_pattern)  zeros(nb2)];
            Bc = [zeros(nb2);
                  eye(nb2)];
            % Discrete-time system with ZOH
            pp.A = expm(Ac*pp.T);
            pp.B = Ac\((pp.A - eye(pp.n))*Bc);
        otherwise
            error("Invalid flag for plant definition; should be 0 or 1.")
    end

end