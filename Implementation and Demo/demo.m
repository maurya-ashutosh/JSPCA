
% Run Joint-Sparse PCA on Breast-Cancer (Wisconsin) dataset
% Dataset source: https://archive.ics.uci.edu/ml/datasets/breast+cancer+wisconsin+(diagnostic)

load data;  % Variable 'W' contains the dataset

% Parameters to use
d = 6;
lambda = 3;
numIter = 50;

% Function call
[Q,costs] = jspca(W',lambda,d,numIter);

% Threshold on Q to visualize the sparsity
th = 0.15;
Q(abs(Q)<th) = 0;

sparsityRatio = sum(Q(:)==0)/numel(Q);
jointSparsityRatio = sum(all(Q==0,2)) ./ size(Q,1);
varexp = sum(var(W*Q)) ./ sum(var(W));

% Print Results
fprintf("Out of %d loadings, %d are zero\n", numel(Q), sum(Q(:)==0));
fprintf("Sparsity: \t %2.1f %%\n", sparsityRatio*100);

fprintf("\nOut of %d features, we have removed %d.\n", size(Q,1), sum(all(Q==0,2)));
fprintf("Joint-sparsity ratio: \t %2.1f %%\n", jointSparsityRatio*100);

fprintf("\nVariance Explained: \t %2.1f %%\n", varexp*100);