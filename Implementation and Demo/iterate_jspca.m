function [D1, D2, P_bar, Q, C] = iterate_jspca(X, D1, D2, P_bar, lambda, epsilon)

    xxt = X*X';
    rD1 = sqrt(D1);
    rD2 = sqrt(D2);
    m = size(X,1);
    diagIndD = (eye(m)==1);

    % Compute Q_bar according to Eqn. 12
    tmp = (lambda * rD1 * rD2 * rD1 + rD1 * xxt * rD1);
    Q_bar = pinv(tmp) * rD1 * xxt * rD1 * P_bar;

    % Compute Q according to Eqn. 13
    tmp = lambda * D2 + xxt;
    Q = pinv(tmp) * xxt * rD1 * P_bar;

    % Compute P_bar with svd according to Eqn. 19
%     [E,~,U] = svd(rD1 * xxt * rD1 * Q_bar, "econ");
    
    % Old Code - we did full SVD and had to pad U with 0's
       [E,~,U] = svd(rD1 * xxt * rD1 * Q_bar);
       pad = zeros(size(U,1), size(E,1)-size(U,1));
       U = [U, pad];
    
    P_bar = E * U';

    % Compute P according to Eqn. 20
    P = pinv(rD1) * E * U';

    % Compute D1 according to Eqn. 5
    tmp = (X - P * Q' * X);
    D1 = zeros(m);
    nrm = vecnorm(tmp,2,2);
    nrm(nrm==0) = epsilon;   % Adjust for norm==0
    D1(diagIndD) = 1./(2*nrm);
    
    % Compute D2 according to Eqn. 6
    D2 = zeros(m);
    nrm = vecnorm(Q,2,2);
    nrm(nrm==0) = epsilon;   % Adjust for norm==0
    D2(diagIndD) = 1./(2*nrm);

    C = iCalculateCost(X,P,Q,rD1,rD2,lambda);

end

function C = iCalculateCost(X,P,Q, rD1,rD2,lambda)
    tmp =  (X - P*Q'*X);
    tmp = rD1 * tmp;
    
    C = norm(tmp, "fro").^2;
    
    C = C + lambda .* norm(rD2 * Q, "fro").^2;
end