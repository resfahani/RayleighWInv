function [RVS,RDC,  res] = ADMM(cv, Forw, G, gamma, mu, sigma, Maxiter, verbose, pltverbose )
%%%
% nonlinear dispersion curve inversion using Dix-type global linear approximation
%by Reza Esfahani,  
% An inexact augmented Lagrangian method for nonlinear dispersion-curve inversion 
% using Dix-type global linear approximation Reza Dokht Dolatabadi Esfahani, Ali Gholami, and Matthias Ohrnberger,
% GEOPHYSICS 2020 85:5, EN77-EN85, https://doi.org/10.1190/geo2019-0717.1"

% Input:
% cv :  estimated dispersion curve
% Forw: Forward modeling operator
% G: Dix-type operator
% gamma: Total variation parameter
% sigma: stopping criteria
% Maxiter: Maximum iteration
% verbose: verbose residual
% pltverbose: plot result for each iteration

% return 
% RVS: Estimated shear wave velocity 
% res: residual 


% Global model parameters for plotting

global model

d       = cv.^2;

[m,n]   = size(G);

% differentiating matrix for Total Variation (TV) regularization
D       = get_l_rough(n, 1);
D(n,1)  = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADMM Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LHS       = inv( 1 * D' * D + gamma * G' * G);

y  = d;

% Auxiliary variables
y1 = zeros(n, 1);
lamnda1 = zeros(n, 1);
lamnda  = zeros(m, 1);


res = [];

k = 1
while 1

    x = real(LHS * (  1 * D'* (y1 + lamnda1 ) + gamma * G'* (lamnda + y)));
    
    u = D * x  - lamnda1 ;
    y1 = max(1 - mu *max(abs(u))*1./abs(u),0) .* (u);
    
    lamnda1 = lamnda1 + y1 - D * x ;
    
    RVS = sqrt(x)';

    RDC = Forw(RVS)';

    lamnda = lamnda + 1 * ( y - RDC.^2);
    
    res(k) = norm(RDC-cv)/length(cv);
    
    
    if verbose == 1
        fprintf('Iteration = %g, Residual = %g \n', k, res(k)); 
    end
    
    if pltverbose == 1 
        
    subplot(221); 
    p = plot(res); 
    p.LineWidth = 1.5;
    xlabel('Iteration');
    ylabel('Residual');

    subplot(222); 
    p = plot(model.vsv, model.hzcum, 'k--');
    p.LineWidth = 1.5;
    hold on;
    plot(RVS,model.hzcum);
    axis ij; ylim([0 80]);
    xlabel('Velocity');
    ylabel('Depth');
    hold off;    

    subplot(212);
    p = plot(model.fks,cv, 'k--');
    p.LineWidth = 1.5;
    hold on;
    plot(model.fks,RDC,'r');
    xlabel('Frequency'); 
    ylabel('Phase velocity');
    hold off;
    
    drawnow
    
    end
    
    if res(k)<sigma || k > Maxiter
        
        return
        
    end
    
    k = k + 1;
    
end
