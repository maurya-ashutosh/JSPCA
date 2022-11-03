function  [Q,costs] = jspca(X, lambda, d, numIter, epsilon, rseed)
% Performs Joint-Sparse Principle Component Analysis on the input dataset
%
% 'X' 

if(nargin<6)
    rseed = 492;
end

if(nargin<5)
    epsilon = eps;
end

if(nargin<4)
    numIter = 50;
end

rng(rseed);

[m,n] = size(X,1,2);

D1 = eye(m);
D2 = eye(m);

% Initialize P_bar as rand(m,d) but it should be a unitary matrix
P_bar = iInitializePbar2(m,d);

costs = [];
for i=1:numIter
    [D1, D2, P_bar, Q, C] = iterate_jspca(X, D1, D2, P_bar, lambda, epsilon);
    costs(end+1) = C;

    % Stopping Condition

end

% Normalize each column 'i' of Q so that norm(Q(:,i),2)==1
cnorms = vecnorm(Q,2,1);
Q = Q./(cnorms);

end

%% Functions to initialize P_bar as a unitary matrix

function out = iInitializePbar1(m,d)
    tmp = eye(m);
    out = tmp(:, 1:d);
end

function out = iInitializePbar2(m,d)
    r = rand(m,m);
    [uty,~,~] = svd(r);
    out = uty(:, 1:d);
end