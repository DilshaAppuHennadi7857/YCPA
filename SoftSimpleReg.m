%% ELEC 4700 - PA 10, EM and the Yee Cell
%
% This program is simulates an EM wave as it propagates through a waveguide
% and encounters an inclusion.
%

winstyle = 'docked';
% winstyle = 'normal';

set(0,'DefaultFigureWindowStyle',winstyle)
set(0,'defaultaxesfontsize',18)
set(0,'defaultaxesfontname','Times New Roman')
% set(0,'defaultfigurecolor',[1 1 1])

% clear VARIABLES;
clear
global spatialFactor;
global c_eps_0 c_mu_0 c_c c_eta_0
global simulationStopTimes;
global AsymForcing
global dels
global SurfHxLeft SurfHyLeft SurfEzLeft SurfHxRight SurfHyRight SurfEzRight



dels = 0.75;
spatialFactor = 1;

c_c = 299792458;                  % speed of light
c_eps_0 = 8.8542149e-12;          % vacuum permittivity
c_mu_0 = 1.2566370614e-6;         % vacuum permeability
c_eta_0 = sqrt(c_mu_0/c_eps_0);


tSim = 200e-15 % simulation time
f = 90e12;%20e12;%800e12;%230e12; % frequency
lambda = c_c/f; % wave length

xMax{1} = 20e-6;
nx{1} = 200;
ny{1} = 0.75*nx{1};


Reg.n = 1;

mu{1} = ones(nx{1},ny{1})*c_mu_0;

epi{1} = ones(nx{1},ny{1})*c_eps_0; % set epsilon for the region
% add inclusions
% epi{1}(125:150,5:15)= c_eps_0*11.3;
% epi{1}(125:150,25:35)= c_eps_0*11.3;
% epi{1}(125:150,45:55)= c_eps_0*11.3;
% epi{1}(125:150,65:75)= c_eps_0*11.3;
% epi{1}(125:150,85:95)= c_eps_0*11.3;
epi{1}(125:150,120:130)= c_eps_0*11.3; %eye
epi{1}(50:75,120:130)= c_eps_0*11.3; %eye
epi{1}(95:105,80:90)= c_eps_0*11.3; %nose
epi{1}(125:150,50:60)= c_eps_0*11.3; %mouth
epi{1}(50:75,50:60)= c_eps_0*11.3; %mouth
epi{1}(50:150,40:50)= c_eps_0*11.3; %mouth
% epi{1}(125:150,:)= c_eps_0*11.3; % A wall

sigma{1} = zeros(nx{1},ny{1});
sigmaH{1} = zeros(nx{1},ny{1});

dx = xMax{1}/nx{1};
dt = 0.25*dx/c_c;
nSteps = round(tSim/dt*2);
yMax = ny{1}*dx;
nsteps_lamda = lambda/dx

movie = 1;
Plot.off = 0;
Plot.pl = 0;
Plot.ori = '13';
Plot.N = 100;
Plot.MaxEz = 1.1;
Plot.MaxH = Plot.MaxEz/c_eta_0;
Plot.pv = [0 0 90];
Plot.reglim = [0 xMax{1} 0 yMax];

%%
%
% This bc structure is used to hold characteristics of the electric and
% megnetice components of the wave. bc{1}.s sets up the source(s) while
% there are also boundary conditions for each xm, xp, ym, and yp (further
% down).
%
bc{1}.NumS = 3; %1;

% Source 1
bc{1}.s(1).xpos = nx{1}/(4) + 1;
bc{1}.s(1).type = 'ss';
bc{1}.s(1).fct = @PlaneWaveBC;

% Source 2
bc{1}.s(2).xpos = nx{1}/(1) + 1;
bc{1}.s(2).type = 'ss';
bc{1}.s(2).fct = @PlaneWaveBC;

% Source 3
bc{1}.s(3).xpos = nx{1}/(10) + 15;
bc{1}.s(3).type = 'ss';
bc{1}.s(3).fct = @PlaneWaveBC;

% mag = -1/c_eta_0;
mag = 1;
phi = 0;
omega = f*2*pi;
betap = 0;
t0 = 30e-15;
st = 15e-15;%-0.05;%
s = 0;
y0 = yMax/2;
sty = 1.5*lambda;
bc{1}.s(1).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'}; % Source 1
bc{1}.s(2).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'}; % Source 2
bc{1}.s(3).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'}; % Source 3



Plot.y0 = round(y0/dx);

% set boundary conditions
bc{1}.xm.type = 'e';
bc{1}.xp.type = 'e';
bc{1}.ym.type = 'e';
bc{1}.yp.type = 'e';

pml.width = 20 * spatialFactor;
pml.m = 3.5;

Reg.n  = 1;
Reg.xoff{1} = 0;
Reg.yoff{1} = 0;

RunYeeReg






