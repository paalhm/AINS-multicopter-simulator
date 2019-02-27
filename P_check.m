function P = P_check(P, nSteps)
    %Make P Symmetric again!
    steps = 0;
    while(sum(sum(P~=P')) && steps < nSteps)
        P = (P+P')/2;
        steps = steps + 1;
    end
    if ~isreal(P) || min(eig(P)) < 0
        error('P is not real and positive definite');
    end
end