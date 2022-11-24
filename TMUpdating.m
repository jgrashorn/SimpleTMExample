
data = [2.9007,24.5910,12.0916]; % "measured" natural frequencies
k0 = [1500,750,1200,1000]'; % nominal values for stiffnesses in physical space

num_threads = 8;
KokkosInitialize(num_threads);

% choice of uncertain stiffnesses, [1,2,3,4]' means all stiffnesses are
% uncertain, [1,2]' would just estimate the first two
Deg = [1,2,3,4]';
Ndeg = length(Deg);
N = 100;
z = randn(Ndeg,N);

opts = MapOptions;
% without setting the bounds MATLAB crashes
% opts.basisLB = -3;
% opts.basisUB = 3;
total_order = 2; % note: 1st order maps work fine

tri_map = CreateTriangular(Ndeg,Ndeg,total_order,opts);

log_target = @(x) wrapper(x,data,Deg,k0);

obj = @(w) objective(log_target,w,tri_map,z);
w0 = tri_map.Coeffs();

options = optimoptions('fminunc','SpecifyObjectiveGradient', true,...
    'Display','iter-detailed','StepTolerance',1e-6,'MaxIterations',400);
    %'CheckGradients',true,'FiniteDifferenceStepSize',1e-5,'FiniteDifferenceType','central');

[~] = fminunc(obj, w0,options);

function [L,dwL]=objective(log_target,coeffs,transport_map,x)
    % KL divergence objective
    num_points = size(x,2);
    transport_map.SetCoeffs(coeffs);
    map_of_x = transport_map.Evaluate(x);
    %logpdf = log_posterior(std_noise,std_prior1,std_prior2,list_t,list_yobs,map_of_x(1,:),map_of_x(2,:));
    [logpdf,gradlogpdf] = log_target(map_of_x);
    logpdf = logpdf';
    log_det = transport_map.LogDeterminant(x);
    L = - sum(logpdf + log_det)/num_points;
    % Gradient of KL divergence objective
    if (nargout > 1)
        sens_vecs = gradlogpdf;
        grad_logpdf = transport_map.CoeffGrad(x,sens_vecs);
        grad_log_det = transport_map.LogDeterminantCoeffGrad(x);
        dwL = - sum(grad_logpdf + grad_log_det,2)/num_points;
    end
end