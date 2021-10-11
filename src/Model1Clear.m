clc;clear; close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Velocity Model and Parameter 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters for Dix-type relation (Dix_inv_flag);
% Inv_flag = 1 Rayleigh wave phase velocity, homogeneous form
% 2 = Rayleigh wave phase velocity, power-law form w/0.25 Poisson's ratio 
% 3 = Love wave phase velocity, power-law form 
% 4 = Rayleigh wave group velocity, homogeneous form
% 5 = Rayleigh wave phase velocity, power-law form w/0.3 Poisson's ratio

Dix_inv_flag = 1;

% Total variation regularization parameter
gamma = 5
mu = 0.3
% Frequency range and DC sampling 
fmin = 3;
fmax = 50;
% Frequency samples
Nf = 50;


% Noise standard deviation
Noise_sd = 0
% Stopping criteria 
sigma =0.03;
% Maximum iteration
Maxiter =50;
% Verbose mode
verbose =1;
% Plotting
pltverbose = 1 ;




freq =  linspace(fmin, fmax, Nf);

% Stopping criteria;
sigma = 0.01;
% Noise standard deviation
NS = 0;

% Model construction
vs      = [80 120 180 360 450];  %  S-velocity
vp      = round(poisfun(vs,.25,3));
% vp      = [300 1000 1400 1400 1400] *1;  %  P-velocity
rho     = [1.8 1.8 1.8 1.8 1.8]*1000;    %  Density
z       = [.2 .2 .2 .2 1];   %  grid spacing
nn       = [5 10 40 80 500];  %  number of gr

global model

model   = model_gen(vs,vp,rho,z,nn);  %  model Generator
model.fks =freq;  %  Frequency

Vs      = model.vsv;
fr      = model.fks;

%%%% Forward operator 
Forw    = @(vsv)Raylee_Forward(vsv, model.vpv, model.rhov, model.h, model.fks, model.Nn);

% DATA 
cvorg      = Forw(Vs)';

cv = cvorg + randn(length(cvorg),1)*min(cvorg) * Noise_sd /100 ;


%%%% Dix-type operator
[GI]    = @(cv) Dix_Function(cv',model.fks,model.hzcum,.3, Dix_inv_flag);
G       = GI(cv);

% ADMM inversion
[RVS, RDC, res] = ADMM(cv, Forw, G, gamma, mu, sigma, Maxiter, verbose, pltverbose) 

%% plotting
figure()
FNT = 8
xa=0.6
xi =1
axes('unit','centimeter','position',[.5+xa +xi 4 8])
    plot(RVS,model.hzcum,'k','LineWidth',1.5); 
      hold on; 

    plot(model2.vsv,model2.hzcum,'r--','LineWidth',2);

        
        axis ij;
    ylim([0 50])
    xlim([50,600])
    hold off    
set(gca,'XAxisLocation','top')
xlabel({'Velocity (m/s)'},'fontsize',FNT,'FontUnits','points','interpreter','latex');
ylabel({'Depth (m)'},'fontsize',FNT,'FontUnits','points','interpreter','latex'),
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',FNT)
leg = legend({'Estimated model','True model'},'FontSize',6,'LineWidth',2,'Location','northeast','NumColumns',1)
legend('boxoff')
leg.ItemTokenSize = [15,16];

axes('unit','centimeter','position',[6+xa 4.5+xi 7 3.5]) 
plot(freq,cv,'r--','LineWidth',2);

;hold on;
plot(freq,RDC,'k','LineWidth',1.)
hold off;
ylabel({'Phase velocity (m/s)'},'fontsize',FNT,'FontUnits','points','interpreter','latex');
xlim([min(freq),max(freq)])
ylim([min(cv)-20,max(cv)+10])
xlabel({'Frequency (Hz)'},'fontsize',FNT,'FontUnits','points','interpreter','latex'),
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',FNT)
leg= legend({'observed data','Reconstructed data'},'FontSize',6,'Location','northeast','NumColumns',1)
legend('boxoff')
leg.ItemTokenSize = [16,15];

axes('unit','centimeter','position',[6+xa 0+xi 7 3.5])
plot(res,'k','LineWidth',1.5)
xlim([0,length(res)])
ylim([0,max(res)+.1])

xlabel({'Iteration (k)'},'fontsize',FNT,'FontUnits','points','interpreter','latex');
ylabel({'Misfit function'},'fontsize',FNT,'FontUnits','points','interpreter','latex'),
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',FNT)

%%

set(gcf,'paperpositionmode','auto')
print('-painters','-depsc','-r1000','Model1C')


