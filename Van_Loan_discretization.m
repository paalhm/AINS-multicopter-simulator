function [F_d, Q_d] = Van_Loan_discretization(F,Q,G,dt)
    n   = size(F,1);
    O   = zeros(n);
    M   = [-F,  G*Q*G';
            O,    F'  ];
    E   = expm(M*dt);
    F_d = E(n+1:end,n+1:end)';
    Q_d = F_d'*E(1:n,n+1:end);
end