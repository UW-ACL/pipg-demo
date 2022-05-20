function [Ad,Bd] = get_discrete(Ac,Bc,dt,disc_type)
    sys = ss(Ac,Bc,[],[]);
    % zoh: zero-order-hold discretization
    % foh: first-order-hold
    sys_d = c2d(sys,dt,disc_type);
    Ad = sys_d.A;
    Bd = sys_d.B;    
end