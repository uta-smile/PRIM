function [X, psnr_vals, energy_vals, time_vals, relative_err, err] = TVLR(A, B, pars)
%
%   Arguments
%       A: The downsampling Fourier Transform Operator  
%       B: Multi-channel data : m * n * T
%       pars: pars.lambda_1, pars.lambda_2, pars.gamma, pars.max_iter, pars.tol
%
%   Output 
%       x: denoised image
%       psnr_vals: psnr vals for every iter
%       energy_vals: energy vals for every iter
%       time_vals: time vals for every iter
%

%% Set Default Parameter
if ~isfield(pars, 'verbose')
    pars.verbose = 0;
end
if ~isfield(pars, 'debug_output')
    pars.debug_output = 0;
end
if ~isfield(pars, 'max_iter')
    pars.max_iter = 500;
end
if ~isfield(pars, 'tol')
    pars.tol = 0;
end

%% Data Preprocessing

[m, n, T] = size(B);
N = m*n;
% [Q1, Q2] = SparseGradientMatrix(m, n);
% K = [Q1; Q2];
% if pars.verbose
%     gt = pars.gt;
% end
% warning('off', 'all');

%% Parameter Initialization
lambda_1 = pars.lambda_1; lambda_2 = pars.lambda_2;
Kn = sqrt(8);
tau = 1/Kn; sigma = 1/Kn;
tol = pars.tol;
gamma = 1;
theta = 1; % 1 for , 0 for Arrow-Hurwicz
max_iter = pars.max_iter;

%% Variable Initialization

%x_n = zeros(size_K(2), 1); 
%X_n = zeros(size(B));
X_n = zeros(size(B));  % size of X should be exactly same as B
X_bar = X_n;
Y_n = Ltrans(X_n);
%y_n = K*b; 
relative_err = zeros(max_iter, 1);
if pars.verbose
    energy_vals = zeros(max_iter+1, 1);
    err = zeros(max_iter+1, 1);
    err(1) = norm(X_n(:) - pars.gt(:));
    psnr_vals = zeros(max_iter, 1);
    time_vals = zeros(max_iter, 1);
    % energy_vals(1) = computeEnergy(X_n, A, B, lambda_1, lambda_2);
    t0 = tic();
end


%% Iteration
for i = 1:max_iter
    %
    % Iteration 
    %
    % Iterate rule (Arrow-Hurwicz and so on)
    
    % Ascent in Y
    
    % Ascent in Z
    
    % Update X
    X_bar = X_n_next + theta*(X_n_next-X_n);
    
    % Extrapolation
    
    relative_err(i) = norm(X_n_next(:) - X_n(:))/norm(X_n(:));
    
    % Output Information
    if pars.verbose
        % calc energy value 
        err(i+1) = norm(X_n_next(:) - pars.gt(:));
        % energy_vals(i+1) = computeEnergy(X_n_next, A, B, lambda_1, lambda_2);
        % psnr_vals(i) = UtilPSNR(X_n_next, gt);
        time_vals(i) = toc(t0);
        
        if pars.debug_output
%             fprintf('Iter %d, Err: %.5f, Energy: %.5f, SNR: %.5f\n', ...
%                 i, norm(X_n_next - X_n), energy, psnr_vals(i));
            fprintf('[Iter %d] RMSE: %.5f t=%.2s\n', i, RMSE(X_n_next, pars.gt(:)), time_vals(i));
        end
    end
    
    if relative_err(i) < tol
        break
    end
    
    X_n = X_bar; Y_n = Y_n_next; % tau = theta*tau; sigma = sigma/theta;
end

X = X_n;


end

%% Compute Energy Function
function val = computeEnergy(X, A, B, lambda_1, lambda_2)
% TODO: Modify the energy function
    val = 1/2*(norm(A*X - B)^2) + ... 
        lambda_1 * TVNorm(Ltrans(X)) + ... 
        lambda_2 * NuclearNorm(X);
end

function Y_out = dualNormBallProjection(Y)
    %% Anisotropic TV: Project onto the $\ell^{\infty}$ unit ball
    % Y_out = cellfun(@(X) sign(X) .* min(abs(X), 1), Y, 'UniformOutput', 0);
    
    %% Isotropic TV: Project onto the $\ell^{2, \infty}$ unit ball
    [P, Q] = Y{1:2}; [m, n, T] = size(P); m = m + 1;
    N = [sum(P.^2, 3);zeros(1,n)] + [sum(Q.^2, 3), zeros(m, 1)];
    Nr = sqrt(max(abs(N), 1)); Nr = repmat(Nr, [1,1,T]);
    P = P./Nr(1:m-1, :, :); Q = Q./Nr(:, 1:n-1, :);
    Y_out = {P, Q};
end

%% Gradient of \frac{1}{2} \|AX-B\|^2
%
%   Compute Gradient of \frac{1}{2}\|A.*X-B\|^2 
%

function G = computeGradient(A, X, B)
    G = A'*(A*X-B);
end



%%  TV transformation
%
%   We reference the Lforward and Ltrans function from FISTA-TV 
%
function X=Lforward(P)
    % TODO: Modify TV calculation
    [m, n, T] = size(P{1}); m = m+1;
    X=zeros(m,n,T);
    X(1:m-1,:,:)=P{1};
    X(:,1:n-1,:)=X(:,1:n-1,:)+P{2};
    X(2:m,:,:)=X(2:m,:,:)-P{1};
    X(:,2:n,:)=X(:,2:n,:)-P{2};
end

function P=Ltrans(X)
    % TODO: Modify TV transpose calculation
    [m,n,T]=size(X);
    P{1}=X(1:m-1,:,:)-X(2:m,:,:);
    P{2}=X(:,1:n-1,:)-X(:,2:n,:);
end

%% 
%
%   Tool functions for nuclear norm and nuclear soft-thresholding
%
function val = TVNorm(P)
    val = 0;
    for i = 1:length(P)
        p = P{i};
        val = val + norm(p(:), 1);
    end
end

function val = NuclearNorm(X)
    [~, S, ~] = svd(X);
    val = norm(diag(S), 1);
end

function I_tilde = MatrixShrinkageOperator(I, lambda)
    if lambda ~= 0
        [U, S, V] = svd(I, 'econ');
        s = diag(S);
        s_tilde = sign(s).*max(0, abs(s)-lambda); 
        I_tilde = U*diag(s_tilde)*V';
    else
        I_tilde = I;
    end
end

%
%   Engineering Trick: RangeProjection
%

function X_p = RangeProjection(X, l, u)
    if((l==-Inf)&&(u==Inf))
        project=@(x)x;
    elseif (isfinite(l)&&(u==Inf))
        project=@(x)(((l<x).*x)+(l*(x<=l)));
    elseif (isfinite(u)&&(l==-Inf))
        project=@(x)(((x<u).*x)+((x>=u)*u));
    elseif ((isfinite(u)&&isfinite(l))&&(l<u))
        project=@(x)(((l<x)&(x<u)).*x)+((x>=u)*u)+(l*(x<=l));
    else
        error('lower and upper bound l,u should satisfy l<u');
    end

    X_p = project(X);
end

