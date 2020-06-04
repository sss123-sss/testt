function SugarCreek_Part1_MeatTempering_SP_02032018_1748hrs
close all
clear all
clc
global h Ti Tinf X_water X_b T_f tnodes rnodes Latent X_pr X_lipid X_ash

%=========== INPUTS ===========%
Diameter=0.10; % meters diameter of cylinder
l=0.30; % meters length of cylinder
m = 1; % m can be slab = 0, cylindrical = 1, or spherical = 2.
Processtime_hrs=5; 
Ti=(-18); %C initial temperature
Tinf=(-1); % C air temperature
v=10; % velocity of air in m/s
h= 10.45-v+10*v^1/2; % W/m2 C convective heat transfer coefficient of air (https://www.engineeringtoolbox.com/convective-heat-transfer-d_430.html)
rnodes=20;
tnodes=10;
%__________________________________%

%====== BIOT NUMBER CHECK ======%
% b=Diameter/2; % meters radius of cylinder
% Tref=Ti+Tinf/2;
% kref=0.617; % W/m C
% rhoref= 996; % kg/m^3
% Cpref= 4178; % J/kg C
% V=pi*b^2*l;
% SA=2*pi*(b^2+b*l); 
% Lc=V/SA
% Biot_Number= h*Lc/kref
('Remember to Check Biot Number and only proceed if it is >0.1, otherwise it is a lumped system')
%__________________________________%

%===== PROPERTIES OF THE MEAT ====%
% Add composition and then change the property expressions to reflect
T_f = 272.6942-273.15; % C initial freezing temperature calculated from equation in Freezing and Thawing Calculations (Pham,XXXX)
Latent=18515; % J/mol latent heat of fusion for water   
% compositions from USDA for 3% fat raw meat
% X <- mass franction
X_water = 0.7475;
X_pr = 0.2198; % mass fraction of protein
X_lipid = 0.03;
X_ash = 1.98*10^(-4);
X_b = 0.06156; % unfreezable water mass franction calculated from equation in Freezing and Thawing Calculations (Pham,XXXX) which has an empirical constant for meat along with composition data
%__________________________________%

%==========%==========%==========%==========%==========%==========%
%==========%==========% BELOW IS AUTOMATIC %==========%==========%
%==========%==========%==========%==========%==========%==========%
b=Diameter/2; % meters radius of cylinder
tf=Processtime_hrs*60*60; % seconds total time
r_mesh = linspace(0,b,rnodes); % linspace (a,b,rnodes)
t_span = linspace(0,tf,tnodes); % linspace (t0, tf,tnodes)

sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,r_mesh,t_span);

u = sol(:,:,1); % Extract the first solution component as u.
size(u) % tnodes x rnodes
%ui = sol(j,:,i) approximates component i of the solution at time tspan(j) and mesh points xmesh(:)

X_ice = (X_water-X_b).*(1-(T_f./u)); % 0.68594*(1-272.6942/T);
X_uw = X_water-X_ice; % unfrozen water mass fraction
SIZE_X_ice=size(X_ice)

    for i=1:tnodes
        for j = 1:rnodes
            [rho_AT,k_AT,cp_AT]=SugarCreek_Part2_thermophys_Calculator_SP_02032018_1748hrs(u(i,j),X_uw(i,j),X_ice(i,j),tnodes,rnodes,T_f,Latent,X_water,X_b,X_pr,X_lipid,X_ash);
            rho(i,j)=rho_AT;
            k(i,j)=k_AT;
            cp(i,j)=cp_AT;
        end
    end
figure()
% subplot(1,4,1)
% HH=heatmap(rho); %Note use heatmap instead of HeatMap if there is an error here
% title('Density');
% SIZE_rho=size(rho)

% subplot(1,4,2)
% HH=heatmap(k); %Note use heatmap instead of HeatMap if there is an error here
title('Thermal conductivity at center (r=0)'); hold on;
 SIZE_k=size(k)
% xlabel('Center -----(nodes in radial distance, meters)------ Surface')
% ylabel('Final time ----(nodes in time, seconds )------ Initial time')
set(gca,'fontsize',14)
plot(t_span./(60*60),k(:,1),'LineWidth',3); hold on;
% xlabel('time (hrs)'); ylabel('k');
plot(t_span./(60*60),X_uw(:,1),'LineWidth',3); hold on;
plot(t_span./(60*60),X_ice(:,1),'LineWidth',3); xlabel('time (hrs)');
legend('k','Unfrozen water','ice');
ylabel('k and Water mass fraction');


% subplot(1,4,3)
% HH=heatmap(cp); %Note use heatmap instead of HeatMap if there is an error here
% title('Specific heat capacity');
% SIZE_cp=size(cp)
% 
% subplot(1,4,4)
% HH=heatmap(u);%Note use heatmap instead of HeatMap if there is an error here
% title('TEMPERATURE');
% figure()
% HH=heatmap(X_ice.*100)
% xlabel('Center ---------------(nodes in radial distance, meters)---------------- Surface')
% ylabel('Final time ---------(nodes in time, seconds )----------- Initial time')
% set(gca,'fontsize',20)
% 

% 
% % figure()
% % surf(r_mesh,t_span,u) % A surface plot is often a good way to study a solution.
% % title('Numerical solution computed with 20 mesh points.')
% % xlabel('Distance x (m)')
% % ylabel('Time t (s)')
% % hold off
% 
figure()
% A solution profile can also be illuminating.
subplot(1,2,1)
plot(r_mesh,u(1,:),'LineWidth',3); hold on;
plot(r_mesh,u(round(length(t_span)/4),:),'LineWidth',3); hold on;
plot(r_mesh,u(round(length(t_span)/2),:),'LineWidth',3); hold on;
plot(r_mesh,u(round(length(t_span)*3/4),:),'LineWidth',3); hold on;
plot(r_mesh,u(end,:),'LineWidth',3); hold on;
legend({'initial t= 0','t=1/4 tf','half-time t= tf/2','t=3/4 tf','final t= tf'},'Location','southeast','FontSize',18)
title('Temperature-distance profile','FontSize',18)
xlabel('Distance from center (m)','FontSize',14)
ylabel('Temperature (C)','FontSize',14)
set(gca,'fontsize',20)

subplot(1,2,2)
plot(t_span./(60*60),u(:,1),'LineWidth',3); hold on;
plot(t_span./(60*60),u(:,round(length(r_mesh)/4)),'LineWidth',3); hold on;
plot(t_span./(60*60),u(:,round(length(r_mesh)/2)),'LineWidth',3); hold on;
plot(t_span./(60*60),u(:,round(length(r_mesh)*3/4)),'LineWidth',3); hold on;
plot(t_span./(60*60),u(:,end),'LineWidth',3); hold on;
legend({'center r=0','r=1/4 R','half-way r=R/2','r= 3/4 R','surface r=R'},'Location','southeast','FontSize',18)
title('Temperature-time profile','FontSize',18)
xlabel('Time (hours)','FontSize',14)
ylabel('Temperature (C)','FontSize',14)
set(gca,'fontsize',20)
% --------------------------------------------------------------

function [c,f,s] = pdex1pde(x,t,u,DuDx)
global h Ti Tinf X_water X_b T_f tnodes rnodes Latent X_pr X_lipid X_ash

% Cp= 2.0082*10^3+(1.2089*u)-(1.3129*10^-3)*u.^2; % J/kg K specific heat capacity of protein Choi-Okos http://b.web.umkc.edu/beckerb/publications/journals/thermophysical.pdf

X_ice = (X_water-X_b).*(1-(T_f./u)); % 0.68594*(1-272.6942/T);
X_uw = X_water-X_ice; % unfrozen water mass fraction
[rho,k,cp]=SugarCreek_Part2_thermophys_Calculator_SP_02032018_1748hrs(u,X_uw,X_ice,tnodes,rnodes,T_f,Latent,X_water,X_b,X_pr,X_lipid,X_ash);           

alpha=k/(rho*cp);

c = 1/alpha;
f = DuDx;
s = 0;
% --------------------------------------------------------------

function u0 = pdex1ic(x)
global h Ti Tinf X_water X_b T_f tnodes rnodes Latent X_pr X_lipid X_ash
u0 = Ti; %Initial temperature throughout the meat
% --------------------------------------------------------------

function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
global h Ti Tinf X_water X_b T_f tnodes rnodes Latent X_pr X_lipid X_ash

X_ice = (X_water-X_b).*(1-(T_f./ur)); % 0.68594*(1-272.6942/T);
X_uw = X_water-X_ice; % unfrozen water mass fraction
[rho,k,cp]=SugarCreek_Part2_thermophys_Calculator_SP_02032018_1748hrs(ur,X_uw,X_ice,tnodes,rnodes,T_f,Latent,X_water,X_b,X_pr,X_lipid,X_ash);           

pl = 0; %for m>0 and x=0, pdepe ignore this bound. See matlab documentation
ql = 1;  %for m>0 and x=0, pdepe ignore this bound. See matlab documentation
pr=h*(Tinf-ur);
qr=(-k);
