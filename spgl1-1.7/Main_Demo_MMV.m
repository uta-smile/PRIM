
    
    % -----------------------------------------------------------
    % Solve the multiple measurement vector (MMV) problem
    %
    %    minimize ||Y||_1,2 subject to AW^{-1}Y = B
    %    
    % and the weighted MMV problem (weights on the rows of X):
    %
    %    minimize ||WX||_1,2 subject to AX = B
    %
    % followed by setting Y = WX.
    % -----------------------------------------------------------
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Solve the multiple measurement vector (MMV) problem      \n');
    fprintf('%%                                                          \n');
    fprintf('%% (1) minimize ||Y||_1,2 subject to AW^{-1}Y = B           \n');
    fprintf('%%                                                          \n');
    fprintf('%% and the weighted MMV problem (weights on the rows of X): \n');
    fprintf('%%                                                          \n');
    fprintf('%% (2) minimize ||WX||_1,2 subject to AX = B                \n');
    fprintf('%%                                                          \n');
    fprintf('%% followed by setting Y = WX.                              \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);

    

    % Initialize random number generator
    randn('state',0); rand('state',0);
    
    % Create problem
    m = 100; n = 150; k = 12; l = 6;
    A = randn(m,n);
    p = randperm(n); p = p(1:k);
    X0= zeros(n,l); X0(p,:) = randn(k,l);
    
    weights = 3 * rand(n,1) + 0.1;
    W = spdiags(1./weights,0,n,n);
    
    B = A*W*X0;
    
    % Solve unweighted version
    opts = spgSetParms('verbosity',1);
    x    = spg_mmv(A*W,B,0,opts);
    x1   = x;
    
    % Solve weighted version
    opts = spgSetParms('verbosity',1,'weights',weights);
    x    = spg_mmv(A,B,0,opts);
    x2   = spdiags(weights,0,n,n) * x;
    
    % Plot results
    figure(1); subplot(2,4,6);
    plot(x1(:,1),'b-'); hold on;
    plot(x2(:,1),'b.');
    plot(X0,'ro');
    plot(x1(:,2:end),'-');
    plot(x2(:,2:end),'b.');
    legend('Coefficients (1)','Coefficients (2)','Original coefficients');
    title('(f) Weighted Basis Pursuit with Multiple Measurement Vectors');

    fprintf('\n');
    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf('See figure 1(f).\n');
    fprintf([repmat('-',1,80), '\n']);
    fprintf('\n\n');

    
    % -----------------------------------------------------------
    % Solve the group-sparse Basis Pursuit problem
    %
    %    minimize    sum_i ||y(group == i)||_2
    %    subject to  AW^{-1}y = b,
    %    
    % with W(i,i) = w(group(i)), and the weighted group-sparse
    % problem
    %
    %    minimize    sum_i w(i)*||x(group == i)||_2
    %    subject to  Ax = b,
    %
    % followed by setting y = Wx.
    % -----------------------------------------------------------
    fprintf(['%% ', repmat('-',1,78), '\n']);
    fprintf('%% Solve the group-sparse Basis Pursuit problem            \n');
    fprintf('%%                                                         \n');
    fprintf('%% (1) minimize    sum_i ||y(group == i)||_2               \n');
    fprintf('%%     subject to  AW^{-1}y = b,                           \n');
    fprintf('%%                                                         \n');
    fprintf('%% with W(i,i) = w(group(i)), and the weighted group-sparse\n');
    fprintf('%% problem                                                 \n');
    fprintf('%%                                                         \n');
    fprintf('%% (2) minimize    sum_i w(i)*||x(group == i)||_2          \n');
    fprintf('%%     subject to  Ax = b,                                 \n');
    fprintf('%%                                                         \n');
    fprintf('%% followed by setting y = Wx.                             \n');
    fprintf(['%% ', repmat('-',1,78), '\n']);

   

    % Initialize random number generator
    randn('state',0); rand('state',2); % 2
    
    % Set problem size and number of groups
    m = 100; n = 150; nGroups = 25; groups = [];
    
    % Generate groups with desired number of unique groups
    while (length(unique(groups)) ~= nGroups)
       groups  = sort(ceil(rand(n,1) * nGroups)); % Sort for display purpose
    end

    % Determine weight for each group
    weights = 3*rand(nGroups,1) + 0.1;
    W       = spdiags(1./weights(groups),0,n,n);

    % Create sparse vector x0 and observation vector b
    p   = randperm(nGroups); p = p(1:3);
    idx = ismember(groups,p);
    x0  = zeros(n,1); x0(idx) = randn(sum(idx),1);
    b   = A*W*x0;
    
    % Solve unweighted version
    opts = spgSetParms('verbosity',1);
    x    = spg_group(A*W,b,groups,0,opts);
    x1   = x;

    % Solve weighted version
    opts = spgSetParms('verbosity',1,'weights',weights);
    x    = spg_group(A,b,groups,0,opts);
    x2   = spdiags(weights(groups),0,n,n) * x;
    
    % Plot results
    figure(1); subplot(2,4,7);
    plot(x1); hold on;
    plot(x2,'b+');
    plot(x0,'ro'); hold off;
    legend('Coefficients (1)','Coefficients (2)','Original coefficients');
    title('(g) Weighted Group-sparse Basis Pursuit');

    fprintf('\n');
    fprintf([repmat('-',1,35), ' Solution ', repmat('-',1,35), '\n']);
    fprintf('See figure 1(g).\n');
    fprintf([repmat('-',1,80), '\n']);
   