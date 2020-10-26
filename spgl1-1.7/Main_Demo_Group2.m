clc;clear;

    % Initialize random number generator
    randn('state',0);     
    % Set problem size and number of groups
    m = 100; n = 150; nGroups = 25; groups = [];
    A = randn(m,n);
    % Generate groups with desired number of unique groups
    while (length(unique(groups)) ~= nGroups)
       groups  = sort(ceil(rand(n,1) * nGroups)); % Sort for display purpose
    end

    % Create sparse vector x0 and observation vector b
    p   = randperm(nGroups); p = p(1:3);
    idx = ismember(groups,p);
    x0  = zeros(n,1); x0(idx) = randn(sum(idx),1);
    b   = A*x0;
    
    % Solve unweighted version
    opts = spgSetParms('verbosity',1);
    
    x    = spg_group(A,b,groups,0,opts);
    x1   = x;
    
    % Plot results
    figure(1); 
    plot(x1); hold on;
    plot(x0,'ro'); hold off;
    legend('Coefficients (1)','Coefficients (2)','Original coefficients');
    title('(g) Weighted Group-sparse Basis Pursuit');  