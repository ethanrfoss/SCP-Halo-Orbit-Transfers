
%% CR3BPOptTransfers
% Description: Studies the CR3BP and performs interior Halo Orbit Transfers
% using Seqeuential Convex Programming for the Earth-Moon System
% Author: Ethan Foss
% Email: erfoss@ucsd.edu
% Notes: A convex solver called CVX must be configured in Matlab for the
% code to run. It can be easily installed for free at http://cvxr.com/cvx/
% All necessary functions are included in this file.
% The code takes about 5 minutes to run.
function CR3BPOptTransfers
%% System Parameters(Earth-Moon):
m1 = 5.9722*10^24; %Earth Mass(kg)
m2 = 7.34767*10^22; %Moon Mass(kg)
Re = 6370*1000; % Earth Radius (m)
Rm = 1740*1000; % Moon Radius (m)
D = 382500*1000; % Earth Moon Distance (m)
TU = 2*pi*3.747*10^5; % Time Unit (s)
DU = D; % Distance Unit (m)

mu = (m1/m2+1)^(-1); %Jacobi Mu Non-Dimensional Constant

e = .0549; % Eccentricity of Earth-Moon System
psi = 0;

% Lagrange Points of Earth Moon System:
L1 = [.83692,0];
L2 = [1.15568,0];
L3 = [-1.00506,0];
L4 = [.48785,.86603];
L5 = [.48785,-.86603];

%% Load in Earth and Moon Images
load topo;
topoMoon = imread('moon.jpg');
[Xsphere,Ysphere,Zsphere] = sphere(50);

%% Scales:
ScaleEarth = 5;
ScaleMoon = 5;
ScaleU = .02;

%% Jacobi Synodic Frame Plot:
p0 = [.6;0;.2]; v0 =[0;.4;0];
x0 = [p0;v0]; x0 = [.795, 0, 0, 0, 0, 0]'; P = 20;
tspan = [0 P];

[t,x] = ode45(@(t,x)JacobiEOM(t,x,mu),tspan,x0,odeset('AbsTol',1e-12,'RelTol',1e-9));

C = @(x) 2*((1-mu)./sqrt((x(:,1)+mu).^2+x(:,2).^2+x(:,3).^2)+mu./sqrt((x(:,1)-1+mu).^2+x(:,2).^2+x(:,3).^2))+x(:,1).^2+x(:,2).^2-(x(:,4).^2+x(:,5).^2+x(:,6).^2);
E = -C(x);

[X,Y,CC] = HillCurve(mu,[-1.4 1.4],[-1.4 1.4]);

% for i = 1:size(CC,1)
%     for j = 1:size(CC,2)
%         if CC(i,j) > 4
%             CC(i,j) = NaN;
%         end
%     end
% end
% 
% figure(1); hold on;axis equal; axis([-.2 1.4 -.8 .8 -1 0]);
% surf(X,Y,(CC-4)); shading interp; colorbar;
% plot([-1.4 1.4],[0 0],'--k');
% plot([0 0],[-1.4 1.4],'--k');
% earth = surf(-mu+Xsphere*ScaleEarth*Re/D,Ysphere*ScaleEarth*Re/D,Zsphere*ScaleEarth*Re/D);
% set(earth,'facecolor','texturemap','cdata',topo,'edgecolor','none');
% moon = surf(1-mu+Xsphere*ScaleMoon*Rm/D,Ysphere*ScaleMoon*Rm/D,Zsphere*ScaleMoon*Rm/D);
% set(moon,'facecolor','texturemap','cdata',topoMoon,'edgecolor','none');

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

figure(1); %sgtitle('Planar Circular Restricted Three Body Problem Plot');
subplot(2,2,1); hold on; axis equal; grid on; axis([-1.4 1.4 -1.4 1.4]);
xlabel('$\xi$'); ylabel('$\eta$'); %%title('Circular-Planar Synodic Frame Three-Body Problem','interpreter','latex','fontsize',20);
plot([-1.4 1.4],[0 0],'--k');
plot([0 0],[-1.4 1.4],'--k');
earth = surf(-mu+Xsphere*ScaleEarth*Re/D,Ysphere*ScaleEarth*Re/D,Zsphere*ScaleEarth*Re/D);
set(earth,'facecolor','texturemap','cdata',topo,'edgecolor','none');
moon = surf(1-mu+Xsphere*ScaleMoon*Rm/D,Ysphere*ScaleMoon*Rm/D,Zsphere*ScaleMoon*Rm/D);
set(moon,'facecolor','texturemap','cdata',topoMoon,'edgecolor','none');
plot(0,0,'k.');

plot(L1(1),L1(2),'b*');
plot(L2(1),L2(2),'b*');
plot(L3(1),L3(2),'b*');
plot(L4(1),L4(2),'b*');
plot(L5(1),L5(2),'b*');

xp = plot3(x(:,1),x(:,2),x(:,3));

contour(X,Y,CC,[-E(1) -E(1)]);

subplot(2,2,2); hold on; axis equal; axis([-1.4*D 1.4*D -1.4*D 1.4*D]); %grid on;
xI = RotatingToInertial(t,x',psi,mu,D);
plot3(xI(1,:),xI(2,:),xI(3,:));
earth = surf(-mu+Xsphere*ScaleEarth*Re,Ysphere*ScaleEarth*Re,Zsphere*ScaleEarth*Re);
set(earth,'facecolor','texturemap','cdata',topo,'edgecolor','none');
plot3(D*cos(2*pi*t+psi),D*sin(2*pi*t+psi),zeros(1,length(t)),'k--');

subplot(2,2,[3 4]); hold on;
xlabel('Time'); ylabel('Jacobi Integral');%%title('Jacobi Integral','interpreter','latex','fontsize',35);
xE = plot(t,E);

%% Generate Periodic Non-Planar Orbit About L1:
Az = .08; % Approximate Z-Amplitude of Periodic Orbit
L = 1;
Lxi = L1(1);
NS = -1; % 1 for North, -1 for South
[x0,P] = RichardsonICs(Az,mu,L,Lxi,NS);

[XT,P,tT,~,phif,xfdot] = ShootingMethod(x0,P,mu,[0;0;0],[3,5],[2,4,6]);
x0 = XT{end}(1,1:6)';

tspan = [0 P];
[t,x] = ode45(@(t,x)JacobiEOM(t,x,mu),tspan,x0,odeset('AbsTol',1e-12,'RelTol',1e-9));

E = -C(x);

figure(2); hold on; axis equal; grid on; axis([-.2 1.4 -.8 .8]);
xlabel('$\xi$','interpreter','latex'); ylabel('$\eta$','interpreter','latex'); %%title('Circular-Planar Synodic Frame Three-Body Problem','interpreter','latex','fontsize',20);
plot([-1.4 1.4],[0 0],'--k');
plot([0 0],[-1.4 1.4],'--k');
earth = surf(-mu+Xsphere*ScaleEarth*Re/D,Ysphere*ScaleEarth*Re/D,Zsphere*ScaleEarth*Re/D);
set(earth,'facecolor','texturemap','cdata',topo,'edgecolor','none');
moon = surf(1-mu+Xsphere*ScaleMoon*Rm/D,Ysphere*ScaleMoon*Rm/D,Zsphere*ScaleMoon*Rm/D);
set(moon,'facecolor','texturemap','cdata',topoMoon,'edgecolor','none');
plot(0,0,'k.');

plot(L1(1),L1(2),'b*');
plot(L2(1),L2(2),'b*');
plot(L3(1),L3(2),'b*');
plot(L4(1),L4(2),'b*');
plot(L5(1),L5(2),'b*');

for i = 1:length(XT)
    plot3(XT{i}(:,1),XT{i}(:,2),XT{i}(:,3));
end

plot3(x(:,1),x(:,2),x(:,3));

[X,Y,CC] = HillCurve(mu,[-1.4 1.4],[-1.4 1.4]);
contour(X,Y,CC,[-E(1) -E(1)]);

%% Generate Halo Families Using L1 Halo Orbit:
dS = .001;
Nm = 700;
Np = 200;
x0S = [1,3,5];
xfS = [2,4,6];
Df0 = [phif(xfS,x0S) xfdot(xfS)];
delX = null(Df0);
XHALOL1m{1} = XT{end};
PHALOL1m(1) = P;
EHALOL1m(1) = -C(XHALOL1m{1}(1:1:6));
XHALOL1p{1} = XT{end};
PHALOL1p(1) = P;
EHALOL1p(1) = -C(XHALOL1p{1}(1:1:6));
tHALOL1m{1} = tT;

for i = 2:Nm
    
    [XHALOL1m{i},PHALOL1m(i),tHALOL1m{i},Converged,delX] = PseudoArclengthContinuation(XHALOL1m{i-1}(1,1:6)',-delX,PHALOL1m(i-1),mu,[0;0;0;dS],x0S,xfS);
    
    if Converged == false
        fprintf('Pseudo Arclength Continuation Algorithm Failed to Converge at %dth step\n',i);
    else
        fprintf('Pseudo Arclength Continuation Algorithm Successfully Converged at %dth step\n',i);
    end
    
    EHALOL1m(i) = -C(XHALOL1m{i}(1:1:6));
    
end

for i = 2:Np
    
    [XHALOL1p{i},PHALOL1p(i),tHALOL1p{i},Converged,delX] = PseudoArclengthContinuation(XHALOL1p{i-1}(1,1:6)',delX,PHALOL1p(i-1),mu,[0;0;0;dS],x0S,xfS);
    
    if Converged == false
        fprintf('Pseudo Arclength Continuation Algorithm Failed to Converge at %dth step\n',i);
    else
        fprintf('Pseudo Arclength Continuation Algorithm Successfully Converfed at %dth step\n',i);
    end
    
    EHALOL1p(i) = -C(XHALOL1p{i}(1:1:6));
    
end

XHALOL1 = [XHALOL1p{end:-1:2} XHALOL1m];
PHALOL1 = [PHALOL1p(end:-1:2) PHALOL1m];
EHALOL1 = [EHALOL1p(end:-1:2) EHALOL1m];
tHALOL1 = [tHALOL1p{end:-1:2} tHALOL1m];

%% Plot L1 Halo Orbits:

f = figure(3); hold on; axis equal; grid on; axis([.7 1.1 -.2 .2]);
xlabel('$\xi$','interpreter','latex'); ylabel('$\eta$','interpreter','latex'); %%title('Circular-Planar Synodic Frame Three-Body Problem','interpreter','latex','fontsize',20);
plot([-1.4 1.4],[0 0],'--k');
plot([0 0],[-1.4 1.4],'--k');
earth = surf(-mu+Xsphere*ScaleEarth*Re/D,Ysphere*ScaleEarth*Re/D,Zsphere*ScaleEarth*Re/D);
set(earth,'facecolor','texturemap','cdata',topo,'edgecolor','none');
moon = surf(1-mu+Xsphere*ScaleMoon*Rm/D,Ysphere*ScaleMoon*Rm/D,Zsphere*ScaleMoon*Rm/D);
set(moon,'facecolor','texturemap','cdata',topoMoon,'edgecolor','none');
plot(0,0,'k.');

plot(L1(1),L1(2),'b*');
plot(L2(1),L2(2),'b*');
plot(L3(1),L3(2),'b*');
plot(L4(1),L4(2),'b*');
plot(L5(1),L5(2),'b*');

for i = 1:length(XHALOL1)
    color = interp1([min(EHALOL1) max(EHALOL1)],autumn(2),EHALOL1(i));
    plot3(XHALOL1{i}(:,1),XHALOL1{i}(:,2),XHALOL1{i}(:,3),'Color',color);
    plot3(XHALOL1{i}(:,1),-XHALOL1{i}(:,2),XHALOL1{i}(:,3),'Color',color);
end

plot3(XHALOL1{Np}(:,1),XHALOL1{Np}(:,2),XHALOL1{Np}(:,3),'Color','g');

% Rotating Halo Plot:
set(f,'WindowState','maximized'); set(f,'Color','w'); axis vis3d;
% for i = 1:359
%     view(i,20);
%     frame = getframe(f);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if i == 1
%         imwrite(imind,cm,[cd '\' 'RotatingFamilyL1' '.gif'],'gif','DelayTime',0, 'Loopcount',inf); 
%     else
%         imwrite(imind,cm,[cd '\' 'RotatingFamilyL1' '.gif'],'gif','DelayTime',0,'WriteMode','append'); 
%     end
% end

%% Generate Periodic Non-Planar Orbit About L2:
Az = .08; % Approximate Z-Amplitude of Periodic Orbit
L = 2;
Lxi = L2(1);
NS = -1; % 1 for North, -1 for South
[x0,P] = RichardsonICs(Az,mu,L,Lxi,NS);

[XT,P,tT,~,phif,xfdot] = ShootingMethod(x0,P,mu,[0;0;0],[3,5],[2,4,6]);
x0 = XT{end}(1,1:6)';

tspan = [0 P];
[t,x] = ode45(@(t,x)JacobiEOM(t,x,mu),tspan,x0,odeset('AbsTol',1e-12,'RelTol',1e-9));

E = -C(x);

f = figure(2); set(f,'color','w'); hold on; axis vis3d
plot3(x(:,1),x(:,2),x(:,3));

% Rotating Halo Plot:
% set(f,'WindowState','maximized');
% for i = 1:359
%     view(i,20);
%     frame = getframe(f);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if i == 1
%         imwrite(imind,cm,[cd '\' 'RotatingHalo' '.gif'],'gif','DelayTime',0, 'Loopcount',inf); 
%     else
%         imwrite(imind,cm,[cd '\' 'RotatingHalo' '.gif'],'gif','DelayTime',0,'WriteMode','append'); 
%     end
% end

%% Generate Halo Families Using L2 Halo Orbit:
dS = .001;
Nm = 700;
Np = 200;
x0S = [1,3,5];
xfS = [2,4,6];
Df0 = [phif(xfS,x0S) xfdot(xfS)];
delX = null(Df0);
XHALOL2m{1} = XT{end};
PHALOL2m(1) = P;
EHALOL2m(1) = -C(XHALOL2m{1}(1:1:6));
XHALOL2p{1} = XT{end};
PHALOL2p(1) = P;
EHALOL2p(1) = -C(XHALOL2p{1}(1:1:6));
tHALOL2m{1} = tT;

for i = 2:Nm
    
    [XHALOL2m{i},PHALOL2m(i),tHALOL2m{i},Converged,delX] = PseudoArclengthContinuation(XHALOL2m{i-1}(1,1:6)',-delX,PHALOL2m(i-1),mu,[0;0;0;dS],x0S,xfS);
    
    if Converged == false
        fprintf('Pseudo Arclength Continuation Algorithm Failed to Converge at %dth step\n',i);
    else
        fprintf('Pseudo Arclength Continuation Algorithm Successfully Converged at %dth step\n',i);
    end
    
    EHALOL2m(i) = -C(XHALOL2m{i}(1:1:6));
    
end

for i = 2:Np
    
    [XHALOL2p{i},PHALOL2p(i),tHALOL2p{i},Converged,delX] = PseudoArclengthContinuation(XHALOL2p{i-1}(1,1:6)',delX,PHALOL2p(i-1),mu,[0;0;0;dS],x0S,xfS);
    
    if Converged == false
        fprintf('Pseudo Arclength Continuation Algorithm Failed to Converge at %dth step\n',i);
    else
        fprintf('Pseudo Arclength Continuation Algorithm Successfully Converfed at %dth step\n',i);
    end
    
    EHALOL2p(i) = -C(XHALOL2p{i}(1:1:6));
    
end

XHALOL2 = [XHALOL2p{end:-1:2} XHALOL2m];
PHALOL2 = [PHALOL2p(end:-1:2) PHALOL2m];
EHALOL2 = [EHALOL2p(end:-1:2) EHALOL2m];
tHALOL2 = [tHALOL2p{end:-1:2} tHALOL2m];

%% Plot L2 Halo Orbits:

f = figure(5); hold on; axis equal; grid on; axis([.9 1.3 -.2 .2]);
xlabel('$\xi$','interpreter','latex'); ylabel('$\eta$','interpreter','latex'); %%title('Circular-Planar Synodic Frame Three-Body Problem','interpreter','latex','fontsize',20);
plot([-1.4 1.4],[0 0],'--k');
plot([0 0],[-1.4 1.4],'--k');
earth = surf(-mu+Xsphere*ScaleEarth*Re/D,Ysphere*ScaleEarth*Re/D,Zsphere*ScaleEarth*Re/D);
set(earth,'facecolor','texturemap','cdata',topo,'edgecolor','none');
moon = surf(1-mu+Xsphere*ScaleMoon*Rm/D,Ysphere*ScaleMoon*Rm/D,Zsphere*ScaleMoon*Rm/D);
set(moon,'facecolor','texturemap','cdata',topoMoon,'edgecolor','none');
plot(0,0,'k.');

plot(L1(1),L1(2),'b*');
plot(L2(1),L2(2),'b*');
plot(L3(1),L3(2),'b*');
plot(L4(1),L4(2),'b*');
plot(L5(1),L5(2),'b*');

for i = 1:length(XHALOL2)
    color = interp1([min(EHALOL2) max(EHALOL2)],autumn(2),EHALOL2(i));
    plot3(XHALOL2{i}(:,1),XHALOL2{i}(:,2),XHALOL2{i}(:,3),'Color',color);
end

plot3(XHALOL2{Np}(:,1),XHALOL2{Np}(:,2),XHALOL2{Np}(:,3),'Color','g');

% Rotating Halo Plot:
% set(f,'WindowState','maximized'); set(f,'Color','w'); axis vis3d;
% for i = 1:359
%     view(i,20);
%     frame = getframe(f);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if i == 1
%         imwrite(imind,cm,[cd '\' 'RotatingFamilyL2' '.gif'],'gif','DelayTime',0, 'Loopcount',inf); 
%     else
%         imwrite(imind,cm,[cd '\' 'RotatingFamilyL2' '.gif'],'gif','DelayTime',0,'WriteMode','append'); 
%     end
% end

%% Plot L1 and L2 Halo Orbits Together:

figure(6); hold on; axis equal; grid on; axis([.7 1.3 -.3 .3]);
xlabel('$\xi$','interpreter','latex'); ylabel('$\eta$','interpreter','latex'); %%title('Circular-Planar Synodic Frame Three-Body Problem','interpreter','latex','fontsize',20);
plot([-1.4 1.4],[0 0],'--k');
plot([0 0],[-1.4 1.4],'--k');
earth = surf(-mu+Xsphere*ScaleEarth*Re/D,Ysphere*ScaleEarth*Re/D,Zsphere*ScaleEarth*Re/D);
set(earth,'facecolor','texturemap','cdata',topo,'edgecolor','none');
moon = surf(1-mu+Xsphere*ScaleMoon*Rm/D,Ysphere*ScaleMoon*Rm/D,Zsphere*ScaleMoon*Rm/D);
set(moon,'facecolor','texturemap','cdata',topoMoon,'edgecolor','none');
plot(0,0,'k.');

plot(L1(1),L1(2),'b*');
plot(L2(1),L2(2),'b*');
plot(L3(1),L3(2),'b*');
plot(L4(1),L4(2),'b*');
plot(L5(1),L5(2),'b*');

for i = 1:length(XHALOL1)
    color = interp1([min([EHALOL1 EHALOL2]) max([EHALOL1 EHALOL2])],autumn(2),EHALOL2(i));
    plot3(XHALOL1{i}(:,1),XHALOL1{i}(:,2),XHALOL1{i}(:,3),'Color',color);
end

for i = 1:length(XHALOL2)
    color = interp1([min([EHALOL1 EHALOL2]) max([EHALOL1 EHALOL2])],autumn(2),EHALOL2(i));
    plot3(XHALOL2{i}(:,1),XHALOL2{i}(:,2),XHALOL2{i}(:,3),'Color',color);
end

%% Poincare Maps from Orbits:

% L1
for i = 1:length(XHALOL1) 
    L1PState(i,:) = XHALOL1{i}(1,:);  
end
% J vs. Zeta
figure(7); 
subplot(1,3,1); hold on;grid on;
plot(EHALOL1,L1PState(:,3),'Color',[0 0 1],'Linewidth',2);
ylabel('$\zeta$','interpreter','latex'); xlabel('J','interpreter','latex');
subplot(1,3,2); hold on; grid on;
plot(EHALOL1,L1PState(:,5),'Color',[0 1 1],'Linewidth',2);
ylabel('$\dot{\eta}$','interpreter','latex'); xlabel('J','interpreter','latex');
subplot(1,3,3); hold on; grid on;
plot(EHALOL1,L1PState(:,1),'Color',[1 0 1],'Linewidth',2);
ylabel('$\xi$','interpreter','latex'); xlabel('J','interpreter','latex');

% L2
for i = 1:length(XHALOL2) 
    L2PState(i,:) = XHALOL2{i}(1,:);  
end
% J vs. Zeta
figure(8);
subplot(1,3,1); hold on; grid on;
plot(EHALOL2,L2PState(:,3),'Color',[0 0 1],'Linewidth',2);
ylabel('$\zeta$','interpreter','latex'); xlabel('J','interpreter','latex');
subplot(1,3,2); hold on; grid on;
plot(EHALOL2,L2PState(:,5),'Color',[0 1 1],'Linewidth',2);
ylabel('$\dot{\eta}$','interpreter','latex'); xlabel('J','interpreter','latex');
subplot(1,3,3); hold on; grid on;
plot(EHALOL2,L2PState(:,1),'Color',[1 0 1],'Linewidth',2);
ylabel('$\xi$','interpreter','latex'); xlabel('J','interpreter','latex');


%% SCP Paramters:
% Parameters:
SCP.lambda = 100000; % Slack Variable Penalty
SCP.lambda0 = 100000; % Initial Slack Variable Penalty
SCP.lambdaf = 100000; % Final Slack Variable Penalty
SCP.Nsub = 10; % Discretization Integration Nodes
SCP.eta = .01; % State Trust Region
SCP.etap = t(end)*.05; % Time Scaling Trust Region
SCP.etau = .5; % Input Trust Region
SCP.eps = .01; % Convergence Value
SCP.pmin = .8; % Min Time scaling
SCP.pmax = 1.2; % Max time Scaling
SCP.deltaMax = 1*10^-5; % Maximum Slack
SCP.iterMax = 7;

%% Spacecraft Paramaters:
S.Tmax = 8*10^-5; % Maximum Thrust Capability [N]
S.m0 = 14; % Initial Mass [kg]
S.mu = mu;
MU = S.m0;

% Scale Spacecraft Parameters by Canonical Units:
S.Tmax = S.Tmax*TU^2/DU/MU;
S.m0 = S.m0/MU;

%% Generate Transfers between L1 Orbits
i1 = 250; % Index of Initial Orbit
i2 = 650; % Index of Final Orbit

% Initial Guess Parameters
No = 5; % Number of Interior Orbits in Transfer
N = 500; % Number of Nodes in Transfer (about 100 per orbit recomended)

% Perform Transfer:
[tL1,xL1,uL1,t0L1,x0L1] = OrbitTransfer(SCP,S,N,No,i1,i2,XHALOL1,EHALOL1,tHALOL1);

% Plot Transfer:
PlotTransfer(9,[.7 1.1 -.2 .2 -.2 .2],mu,ScaleEarth,ScaleMoon,ScaleU,Re,Rm,topo,topoMoon,D,L1,L2,L3,L4,L5,XHALOL1{i1},XHALOL1{i2},x0L1,xL1,tL1,uL1);

% Animate Orbit:
% saveName = 'L1Transfer';
% AnimateTransfer(10,[.7 1.1 -.2 .2 -.2 .2],mu,ScaleEarth,ScaleMoon,ScaleU,Re,Rm,topo,topoMoon,D,L1,L2,L3,L4,L5,XHALOL1{i1},XHALOL1{i2},xL1,tL1,uL1,500,saveName);

%% Generate Transfers between L2 Orbits
i1 = 450; % Index of Initial Orbit
i2 = 250; % Index of Final Orbit

% Initial Guess Parameters
No = 5; % Number of Orbits in Transfer
N = 500; % Number of Nodes in Transfer

% Perform Transfer:
[tL2,xL2,uL2,t0L2,x0L2] = OrbitTransfer(SCP,S,N,No,i1,i2,XHALOL2,EHALOL2,tHALOL2);

% Plot Transfer:
PlotTransfer(11,[.9 1.5 -.3 .3 -.3 .3],mu,ScaleEarth,ScaleMoon,ScaleU,Re,Rm,topo,topoMoon,D,L1,L2,L3,L4,L5,XHALOL2{i1},XHALOL2{i2},x0L2,xL2,tL2,uL2);

% Animate Orbit:
% saveName = 'L2Transfer';
% AnimateTransfer(12,[.9 1.5 -.3 .3 -.3 .3],mu,ScaleEarth,ScaleMoon,ScaleU,Re,Rm,topo,topoMoon,D,L1,L2,L3,L4,L5,XHALOL2{i1},XHALOL2{i2},xL2,tL2,uL2,500,saveName);


end

%% Rotating to Inertial Frame Conversion:
% Inputs:
%  t - time vector
%  x - state vector in synodic frame
%  phi - Phase shift
%  mu - Nondimensional constant of CR3BBP
%  D - Earth Radius
% Outputs:
%  xI - Inertial Coordinates
function [xI] = RotatingToInertial(t,x,phi,mu,D)

xI = zeros(size(x));

for i = 1:length(x)
   
    R = [cos(t(i)+phi) -sin(t(i)+phi) 0;
         sin(t(i)+phi) cos(t(i)+phi) 0;
         0 0 1];
     
    xI(:,i) = [R zeros(3,3);zeros(3,3) R]*(x(:,i)-[mu;0;0;0;0;0])*D;
    
end

end

%% Equations of Motion for Circular Restricted Three Body Problem:
% Inputs:
%  t - time vector
%  x - state vector in synodic frame
%  mu - Nondimensional constant of CR3BBP
% Outputs:
%  xdot - Change in state
function xdot = JacobiEOM(t,x,mu)

rho13 = ((x(1)+mu)^2+x(2)^2+x(3)^2)^(3/2);
rho23 = ((x(1)+mu-1)^2+x(2)^2+x(3)^2)^(3/2);
xdot = [x(4)
        x(5)
        x(6)
        -(1-mu)*(x(1)+mu)/rho13-mu*(x(1)+mu-1)/rho23+2*x(5)+x(1)     
        -((1-mu)/rho13+mu/rho23)*x(2)-2*x(4)+x(2)
        -((1-mu)/rho13+mu/rho23)*x(3)];

end

%% Hill Limiting Curve Generation:
function [X,Y,C] = HillCurve(mu,xlim,ylim)

[X,Y] = meshgrid(linspace(xlim(1),xlim(2),1000),linspace(ylim(1),ylim(2),1000));

r1 = sqrt((X+mu).^2+Y.^2);
r2 = sqrt((X-1+mu).^2+Y.^2);
C = 2*((1-mu)*(.5*r1.^2+r1.^(-1))+mu*(.5*r2.^2+r2.^(-1)));

end

%% Initial Conditions for Determining out of Plane Periodic Orbits, from Richardson's Paper:
% Inputs:
%  Az - Z-Amplitude
%  mu - Nondimensional constant of CR3BBP
%  L - Lagrange Point Choice
%  Lxi - Lagrange Point Coords
%  NS - North/SOuth Direction
% Outputs:
%  x0 - Initial Conditions
%  P - Initial Period
function [x0,P] = RichardsonICs(Az,mu,L,Lxi,NS)

if L == 1
    gamL = (1-mu)-Lxi;
    n = 1:4;
    c = 1/gamL^3*(mu+(-1).^(n).*(1-mu).*gamL.^(1+n)./(1-gamL).^(n+1));
    n = 2-NS;
elseif L == 2
    gamL = Lxi-(1-mu);
    n = 1:4;
    c = 1/gamL^3*((-1).^(n)*mu+(-1).^(n).*(1-mu).*gamL.^(1+n)./(1+gamL).^(n+1));
    n = 2+NS;
elseif L == 3
    gamL = (1-mu)-Lxi;
    n = 1:4;
    c = 1/gamL^3*((1-mu)+mu*gamL.^(1+n)./(1+gamL).^(n+1));
    n = 2-NS;
else
    display('Invalid Libration Point');
end

lam = sqrt((-(c(2)-2)+sqrt((c(2)-2)^2+4*(c(2)-1)*(1+2*c(2))))/2);
k = 1/(2*lam)*(lam^2+1+2*c(2));

d1 = 3*lam^2/k*(k*(6*lam^2-1)-2*lam);
d2 = 8*lam^2/k*(k*(11*lam^2-1)-2*lam);

b21 = -3*c(3)*lam/(2*d1)*(3*k*lam-4);
b22 = 3*c(3)*lam/d1;

d21 = -c(3)/(2*lam^2);

a21 = 3*c(3)*(k^2-2)/(4*(1+2*c(2)));
a22 = 3*c(3)/(4*(1+2*c(2)));
a23 = -3*c(3)*lam/(4*k*d1)*(3*k^3*lam-6*k*(k-lam)+4);
a24 = -3*c(3)*lam/(4*k*d1)*(2+3*k*lam);
a31 = -9*lam/(4*d2)*(4*c(3)*(k*a23-b21)+k*c(4)*(4+k^2))+(9*lam^2+1-c(2))/(2*d2)*(3*c(3)*(2*a23-k*b21)+c(4)*(2+3*k^2));
a32 = -1/d2*(9*lam/4*(4*c(3)*(k*a24-b22)+k*c(4))+3/2*(9*lam^2+1-c(2))*(c(3)*(k*b22+d21-2*a24)-c(4)));

b31 = 3/(8*d2)*(8*lam*(3*c(3)*(k*b21-2*a23)-c(4)*(2+3*k^2))+(9*lam^2+1+2*c(2))*(4*c(3)*(k*a23-b21)+k*c(4)*(4+k^2)));
b32 = 1/d2*(9*lam*(c(3)*(k*b22+d21-2*a24)-c(4))+3/8*(9*lam^2+1+2*c(2))*(4*c(3)*(k*a24-b22)+k*c(4)));


d31 = 3/(64*lam^2)*(4*c(3)*a24+c(4));
d32 = 3/(64*lam^2)*(4*c(3)*(a23-d21)+c(4)*(4+k^2));

a1 = -3/2*c(3)*(2*a21+a23+5*d21)-3/8*c(4)*(12-k^2);
a2 = 3/2*c(3)*(a24-2*a22)+9/8*c(4);

s1 = 1/(2*lam*(lam*(1+k^2)-2*k))*(3/2*c(3)*(2*a21*(k^2-2)-a23*(k^2+2)-2*k*b21)-3/8*c(4)*(3*k^4-8*k^2+8));
s2 = 1/(2*lam*(lam*(1+k^2)-2*k))*(3/2*c(3)*(2*a22*(k^2-2)+a24*(k^2+2)+2*k*b22+5*d21)+3/8*c(4)*(12-k^2));

l1 = a1 + 2*lam^2*s1;
l2 = a2+2*lam^2*s2;

del = lam^2-c(2);
deln = (2-n); % 1 or 3

Az = Az/gamL;
Ax = sqrt((-del-l2*Az^2)/l1);

Om = 1+s1*Ax^2+s2*Az^2;

if L == 2
    xi0 = a21*Ax^2+a22*Az^2+Ax+(a23*Ax^2-a24*Az^2)-(a31*Ax^3-a32*Ax*Az^2);
    eta0 = 0;
    zeta0 = -deln*Az+deln*d21*Ax*Az*(-2)-deln*(d32*Az*Ax^2-d31*Az^3);
    
    xidot0 = 0;
    etadot0 = lam*Om*(-k*Ax+2*(b21*Ax^2-b22*Az^2)-3*(b31*Ax^3-b32*Ax*Az^2));
    zetadot0 = 0;
else
    xi0 = a21*Ax^2+a22*Az^2-Ax+(a23*Ax^2-a24*Az^2)+(a31*Ax^3-a32*Ax*Az^2);
    eta0 = 0;
    zeta0 = deln*Az+deln*d21*Ax*Az*(-2)+deln*(d32*Az*Ax^2-d31*Az^3);
    
    xidot0 = 0;
    etadot0 = lam*Om*(k*Ax+2*(b21*Ax^2-b22*Az^2)+3*(b31*Ax^3-b32*Ax*Az^2));
    zetadot0 = 0;
end

P = 2*pi/(lam*Om);

x0 = [xi0;eta0;zeta0;xidot0;etadot0;zetadot0]*gamL + Lxi*[1;zeros(5,1)];

end

%% Equations of Motion for CR3BP and State Transition Matrix Propogation:
% Inputs:
%  t - time vector
%  xphi - state vector with flat STM appended
%  mu - Nondimensional constant of CR3BBP
% Outputs:
%  dxphi - Change in State
function dxphi = XSTMEOM(t,xphi,mu)

x = xphi(1:6);
phi = reshape(xphi(7:end),6,6);

rho13 = ((x(1)+mu)^2+x(2)^2+x(3)^2)^(3/2);
rho23 = ((x(1)+mu-1)^2+x(2)^2+x(3)^2)^(3/2);
xdot = [x(4)
        x(5)
        x(6)
        -(1-mu)*(x(1)+mu)/rho13-mu*(x(1)+mu-1)/rho23+2*x(5)+x(1)     
        -((1-mu)/rho13+mu/rho23)*x(2)-2*x(4)+x(2)
        -((1-mu)/rho13+mu/rho23)*x(3)];

rho1 = (x(1)+mu).^2+x(2).^2+x(3).^2;
rho2 = (x(1)-1+mu).^2+x(2).^2+x(3).^2;

Uxx = [(mu - 1)/rho1^(3/2) - mu/rho2^(3/2) + (3*mu*(2*mu + 2*x(1) - 2)*(mu + x(1) - 1))/(2*rho2^(5/2)) - (3*(2*mu + 2*x(1))*(mu + x(1))*(mu - 1))/(2*rho1^(5/2)) + 1,(3*x(2)*mu*(mu + x(1) - 1))/rho2^(5/2) - (3*x(2)*(mu + x(1))*(mu - 1))/rho1^(5/2),(3*mu*x(3)*(mu + x(1) - 1))/rho2^(5/2) - (3*x(3)*(mu + x(1))*(mu - 1))/rho1^(5/2)
       -x(2)*((3*(2*mu + 2*x(1))*(mu - 1))/(2*rho1^(5/2)) - (3*mu*(2*mu + 2*x(1) - 2))/(2*rho2^(5/2))),(mu - 1)/rho1^(3/2) - mu/rho2^(3/2) + x(2)*((3*x(2)*mu)/rho2^(5/2) - (3*x(2)*(mu - 1))/rho1^(5/2)) + 1,x(2)*((3*mu*x(3))/rho2^(5/2) - (3*x(3)*(mu - 1))/rho1^(5/2))
       -x(3)*((3*(2*mu + 2*x(1))*(mu - 1))/(2*rho1^(5/2)) - (3*mu*(2*mu + 2*x(1) - 2))/(2*rho2^(5/2))),x(3)*((3*x(2)*mu)/rho2^(5/2) - (3*x(2)*(mu - 1))/rho1^(5/2)), x(3)*((3*mu*x(3))/rho2^(5/2) - (3*x(3)*(mu - 1))/rho1^(5/2)) + (mu - 1)/rho1^(3/2) - mu/rho2^(3/2)];
   
A = [zeros(3,3),eye(3,3)
     Uxx,[0 2 0;-2 0 0;0 0 0]];

phidot = A*phi;

dxphi = [xdot; reshape(phidot,36,1)];

end

%% Shooting Method For Determining Periodic Orbits:
% Inputs:
%  x0 - initial guess
%  P - initial Period
%  mu - Nondimensional constant of CR3BBP
%  xd - Desired Final Conditions
%  x0s - Initial Conditions Indices
%  xfS - Final Conditions Indices
% Outputs:
%  XT - Structure of trajectories
%  P - Final Period
%  t - Time vector
%  Converged - Convergence boolean
%  phif - Final STM
%  xfdot - Final change in state
function [XT,P,t,Converged,phif,xfdot] = ShootingMethod(x0,P,mu,xd,x0S,xfS)

tf = P/2;
maxiter = 100;
tol = 10^(-8);
i = 0;

%figure(1); hold on;

while 1
    
    tspan = linspace(0,tf,1000);
    xphi0 = [x0; reshape(eye(6,6),36,1)];
    [t,xphi] = ode45(@(t,xphi)XSTMEOM(t,xphi,mu),tspan,xphi0,odeset('AbsTol',1e-12,'RelTol',1e-9));
    
    xf = xphi(end,1:6)';
    phif = reshape(xphi(end,7:end),6,6);
    
    xfdot = JacobiEOM(t,xf,mu);
    xe = xd-xf(xfS);
    
    if norm(xe)<tol || i>maxiter
        break;
    end
    
    K = [phif(xfS,x0S) xfdot(xfS)];
    Delta = K'*inv(K*K')*xe;
    
    x0(x0S) = x0(x0S)+Delta(1:end-1);
    tf = tf + Delta(end);
    
    i = i + 1;
    XT{i} = xphi(:,1:6);
    
end

if i>maxiter
%     disp('Shoting Algorithm Did Not Converge Before Maximum Iteration Limit');
    Converged = false;
else
%     fprintf('Shoting Algorithm Converged in %d iterations',i);
    Converged = true;
end

P = tf*2;

end

%% Pseudo-Arclength Method for Determining Orbit Families
% Inputs:
%  x0 - initial guess
%  delXprev - Previous tangent correction
%  P - initial Period
%  mu - Nondimensional constant of CR3BBP
%  xd - Desired Final Conditions
%  x0s - Initial Conditions Indices
%  xfS - Final Conditions Indices
% Outputs:
%  XT - Structure of trajectories
%  P - Final Period
%  t - Time vector
%  Converged - Convergence boolean
%  delX - New tangent correction
function [XT,P,t,Converged,delX] = PseudoArclengthContinuation(x0,delXprev,P,mu,xd,x0S,xfS)

tf = P/2;
maxiter = 100;
tol = 10^-8;
i = 0;

x0prev = x0;
tfprev = tf;

while 1
    
    tspan = linspace(0,tf,1000);
    xphi0 = [x0; reshape(eye(6,6),36,1)];
    [t,xphi] = ode45(@(t,xphi)XSTMEOM(t,xphi,mu),tspan,xphi0,odeset('AbsTol',1e-12,'RelTol',1e-9));
    
    xf = xphi(end,1:6)';
    phif = reshape(xphi(end,7:end),6,6);
    
    xfdot = JacobiEOM(t,xf,mu);
    
    xe = xd - [xf(xfS);([x0(x0S); tf] - [x0prev(x0S); tfprev])'*delXprev];
    Df = [phif(xfS,x0S) xfdot(xfS)];
    
    if norm(xe)<tol || i>maxiter
        break;
    end
    
    K = [Df;delXprev'];
    
    Delta = K'*inv(K*K')*xe;
    
    x0(x0S) = x0(x0S)+Delta(1:end-1);
    tf = tf + Delta(end);
    
    i = i + 1;
    
    XT = xphi(:,1:6);
    
end

if i>maxiter
%     disp('Pseudo Arclength Did Not Converge Before Maximum Iteration Limit');
    Converged = false;
else
%     disp(sprintf('Pseudo Arclength Converged in %d iterations',i));
    Converged = true;
end

delX = null(Df);
P = tf*2;
XT = [XT; [XT(end-1:-1:1,1) -XT(end-1:-1:1,2) XT(end-1:-1:1,3) -XT(end-1:-1:1,4) XT(end-1:-1:1,5) -XT(end-1:-1:1,6)]];
t = [t; 2*t(end)-t(end-1:-1:1)];

end




%% Orbit Transfer
% Performs the Orbit Transfer using SCP and given a set of Halo Orbit
% Families, the indices of the start and final orbit, the number of nodes,
% and the number of intermediate orbits.
function [t,x,u,t0,x0] = OrbitTransfer(SCP,S,N,No,i1,i2,X,E,T)

% Generate Initial Guess between L1 Orbits
C1 = E(i1);
C2 = E(i2);

Cs = linspace(C1,C2,No);
is = interp1(E,1:length(E),Cs,'nearest');

x0 = X{is(1)}(1,1:6);
t0 = 0;

% Append Orbits:
for i = is
    
    xi = X{i}(:,1:6);
    ti = T{i};
    
    x0 = [x0;xi(2:end,1:6)];
    t0 = [t0; t0(end)+ti(2:end)];
    
end

% Interpolate Nodes:
t = linspace(t0(1),t0(end),N);
x = interp1(t0,x0,t);
u = ([0;0;0]*zeros(1,length(t)))';

% Call Sequential Program
[x,u,p,deltaMax] = SequentialAlgorithm(S,@DynamicsContinuous,SCP,x',u',t0(end),x0(1,1:6)',x0(end,1:6)',N);

if deltaMax < SCP.deltaMax
    fprintf('Sequential Convex Program Converged!\n');
else
    fprintf('Sequential Convex Program Uncoverged, Infeasibility value: %d \n Try Increasing Orbit Count\n',deltaMax);
end

t = t*p;
x = x';
u = u';

end

%% Sequential Algorithm
% SCP Algorithm that Solves the Trajectory Optimziation Problem given some
% initial and final conditions, and an initial guess.
function [x,u,p,deltaMax] = SequentialAlgorithm(S,DynamicsContinuous,SCP,xk0,uk0,p0,xic,xtc,N)

tk = linspace(0,1,N);
xbar = xk0;
ubar = uk0;
pbar = p0;

Nx = size(xbar,1);
Nu = size(ubar,1);

pmin = SCP.pmin*p0;
pmax = SCP.pmax*p0;

i = 1;

f = figure; hold on; view(90,0);
plot3(xk0(1,:),xk0(2,:),xk0(3,:));

fprintf('Starting SCP Algorithm...\n');

while 1
    
    % Convexify Problem, obtain discretized dynamics
    fprintf('Convexifying Problem...\n');
    [delta,DynamicsDiscrete] = Convexify(S,DynamicsContinuous,SCP,xbar,ubar,pbar,tk);
    fprintf('Convexification Complete\n');
    
    % Solve Problem as a SOCP
    fprintf('Calling Convex Solver...\n');
    [L,x,u,p,deltaMax] = ConvexSolver(DynamicsDiscrete,xbar,ubar,pbar,SCP,xic,xtc,pmin,pmax,N,Nx,Nu);
%     Jbar = NonLinearCost(x,u,xic,xtc,delta,SCP.lambda,SCP.lambda0,SCP.lambdaf);
    fprintf('Convex Problem Solved\n');

    plot3(x(1,:),x(2,:),x(3,:));
    drawnow;
    
    % Check Convergence:
    if i == 1
        Lbar = L;
        xbar = x;
        ubar = u;
        pbar = p;
        i = i + 1;
    else
        e = (Lbar - L)/L;
        if e <= SCP.eps && deltaMax < SCP.deltaMax
            break;
        else
            xbar = x;
            ubar = u;
            pbar = p;
            Lbar = L;
            fprintf('Current Iteration of SCP Algorithm: %d \n Current Error: %d\n',i,e);
            i = i + 1;
        end
    end
    
    if i > SCP.iterMax
        fprintf('Iteration Limit of SCP Algorithm Reached, Exiting\n');
        break;
    end
        
    % Check Convergence Criterion
%     e = abs(Jbar - L)/abs(Jbar);
%     if e <= SCP.eps
%         fprintf('SCP Algorithm Converged');
%         break;
%     else
%         xbar = x;
%         ubar = u;
%         pbar = p;
%         
%         fprintf('Current Iteration of SCP Algorithm: %d \n Current Error: %d',i,e);
%         i = i + 1;
%         
%     end
    
end

close(f);

end

%% Convex Solver
% Solves the Convex Subproblem given Discrete Matrices, a reference
% trajectory, and SCP parameters
function [L,x,u,p,deltaMax] = ConvexSolver(DynamicsDiscrete,xbar,ubar,pbar,SCP,xic,xtc,pmin,pmax,N,Nx,Nu)

l = [zeros(6,1); -1];
E = eye(Nx,Nx);
H0 = eye(6);
Hf = eye(6);

eta = SCP.eta;
etap = SCP.etap;
etau = SCP.etau;
lambda0 = SCP.lambda0;
lambdaf = SCP.lambdaf;
lambda = SCP.lambda;

Nx0 = size(xic);
Nxf = size(xtc);

% Begin Optimization
cvx_begin

variables x(Nx,N) u(Nu,N) p v(Nx,N) s(N) vic(Nx0) sic vtc(Nxf) stc su(N)

minimize(sum(su) + lambda*sum(s) + lambda0*sic + lambdaf*stc)

subject to

% Boundary Conditions:
xic == H0*x(:,1) + vic;
xtc == Hf*x(:,N) + vtc;

% Slack Variable Constraints:
norm(vic,1) <= sic;
norm(vtc,1) <= stc;

for k = 1:N
    
    % Dynamics:
    if k ~= N
        x(:,k+1) == DynamicsDiscrete.Ad{k}*x(:,k) + DynamicsDiscrete.Bdm{k}*u(:,k) + DynamicsDiscrete.Bdp{k}*u(:,k+1) + DynamicsDiscrete.Sd{k}*p + DynamicsDiscrete.Rd{k} + E*v(:,k);
    end
    
    % Slack Variable Constraints:
    norm(v(:,k),1) <= s(k);
    norm(u(1:3,k)) <= su(k);
    
    % Parameter Constraints:
    pmin <= p <= pmax;
    
    % Thrust Constraint
    norm(u(1:3,k)) <= 1;
    
    % Parameter Trust Region:
    norm(x(:,k)-xbar(:,k),2) <= eta;
    norm(u(:,k)-ubar(:,k),2) <= etau;
    norm(p-pbar,2) <= etap;
end

cvx_end

L = cvx_optval;
deltaMax = max(s);

end


% function J = NonLinearCost(xk,uk,xic,xtc,delta,lambda,lambda0,lambdaf)
% 
% H0 = eye(6);
% Hf = eye(6);
% J =  sum(vecnorm(uk,2)) + lambda*sum(abs(delta)) + lambda0*norm((H0*xk(:,1)-xic),1) + lambdaf*norm((Hf*xk(:,end)-xtc),1);
% 
% end

%% Convexify
% Convexifies the Dynamics of the problem
function [delta,DynamicsMatrices] = Convexify(S,Dynamics,SCP,xk,uk,p,tk)

[Ad,Bdm,Bdp,Sd,Rd,delta] = Propagate(S,@Dynamics,SCP,xk,uk,tk,p);

DynamicsMatrices.Ad = Ad;
DynamicsMatrices.Bdm = Bdm;
DynamicsMatrices.Bdp = Bdp;
DynamicsMatrices.Sd = Sd;
DynamicsMatrices.Rd = Rd;

end

%% Propogate
% Propogates the Discrete matrices forward
function [Ad,Bdm,Bdp,Sd,Rd,delta] = Propagate(S,Dynamics,SCP,xk,uk,tk,p)

Nx = size(xk,1);
Nu = size(uk,1);
N = size(xk,2);

Feasible = true;

for k = 1:N-1

    phi0 = eye(Nx);
    Bm0 = zeros(Nx,Nu);
    Bp0 = zeros(Nx,Nu);
    S0 = zeros(Nx,1);
    R0 = zeros(Nx,1);
    P0 = Flatten(xk(:,k),phi0,Bm0,Bp0,S0,R0,Nx,Nu);
    
    tspan = linspace(tk(k),tk(k+1),SCP.Nsub);
    P = RK4(@(t,x)Derivatives(t,x,p,S,@Dynamics,tk(k+1),tk(k),uk(:,k+1),uk(:,k),Nx,Nu),tspan,P0);
    
    [xd,Pphi,PBm,PBp,Ps,Pr] = Unflatten(P,Nx,Nu);
    
    delta(k) = norm(xd-xk(:,k+1));
    
%     if delta > SCP.FeasibilityTolerance
%         Feasible = false;
%     end
    
    Ad{k} = Pphi;
    Bdm{k} = Pphi*PBm;
    Bdp{k} = Pphi*PBp;
    Sd{k} = Pphi*Ps;
    Rd{k} = Pphi*Pr;
    

end

end

%% Unflatten
% Unflattens the discrete vector
function [x,Pphi,PBm,PBp,PS,PR] = Unflatten(P,Nx,Nu)

x = P(1:Nx);
Pphi = reshape(P(Nx+1:Nx+Nx^2),[Nx,Nx]);
PBm = reshape(P(Nx+Nx^2+1:Nx+Nx^2+Nx*Nu),[Nx,Nu]);
PBp = reshape(P(Nx+Nx^2+Nx*Nu+1:Nx+Nx^2+2*Nx*Nu),[Nx,Nu]);
PS = P(Nx+Nx^2+2*Nx*Nu+1:2*Nx+Nx^2+2*Nx*Nu);
PR = P(2*Nx+Nx^2+2*Nx*Nu+1:3*Nx+Nx^2+2*Nx*Nu);

end

%% Flatten
% Flattens the discrete matrices
function P = Flatten(x,Phi,PBm,PBp,PS,PR,Nx,Nu)

P = [x;
     reshape(Phi,[Nx^2,1]);
     reshape(PBm,[Nx*Nu,1]);
     reshape(PBp,[Nx*Nu,1]);
     PS;
     PR];

end

%% RK4
% RK4 Numerical Integration
function x = RK4(xdot,tspan,x0)

x = x0;
dt = tspan(2)-tspan(1);

for i = 1:length(tspan)-1
    
    f1 = xdot(tspan(i),x);
    f2 = xdot(tspan(i)+dt/2,x+f1*dt/2);
    f3 = xdot(tspan(i)+dt/2,x+f2*dt/2);
    f4 = xdot(tspan(i)+dt,x+f3*dt);

    x = x + dt*(f1/6+(f2+f3)/3+f4/6);

end

end

%% Derivatives
% Computes the Derivatives for Discretization
function [Pdot] = Derivatives(t,P,p,S,Dynamics,tk1,tk,uk1,uk,Nx,Nu)

[x,Pphi,PBm,PBp,PS,PR] = Unflatten(P,Nx,Nu);

lambdakm = (tk1-t)/(tk1-tk);
lambdakp = 1-lambdakm;

u = lambdakm*uk + lambdakp*uk1;

[f,A,B,S,R] = Dynamics(S,x,u,p);

psi = Pphi^-1;

Pxdot = f;
Pphidot = A*Pphi;
PBmdot = psi*lambdakm*B;
PBpdot = psi*lambdakp*B;
PSdot = psi*S;
PRdot = psi*R;

Pdot = Flatten(Pxdot,Pphidot,PBmdot,PBpdot,PSdot,PRdot,Nx,Nu);

end

%% Dynamics
% Linearized Dynamics of CR3BP
function [f,A,B,S,R] = Dynamics(S,x,u,p)

Tmax = S.Tmax;
m = S.m0;
mu = S.mu;

rho13 = ((x(1)+mu)^2+x(2)^2+x(3)^2)^(3/2);
rho23 = ((x(1)+mu-1)^2+x(2)^2+x(3)^2)^(3/2);

rho1 = (x(1)+mu).^2+x(2).^2+x(3).^2;
rho2 = (x(1)-1+mu).^2+x(2).^2+x(3).^2;

Uxx = [(mu - 1)/rho1^(3/2) - mu/rho2^(3/2) + (3*mu*(2*mu + 2*x(1) - 2)*(mu + x(1) - 1))/(2*rho2^(5/2)) - (3*(2*mu + 2*x(1))*(mu + x(1))*(mu - 1))/(2*rho1^(5/2)) + 1,(3*x(2)*mu*(mu + x(1) - 1))/rho2^(5/2) - (3*x(2)*(mu + x(1))*(mu - 1))/rho1^(5/2),(3*mu*x(3)*(mu + x(1) - 1))/rho2^(5/2) - (3*x(3)*(mu + x(1))*(mu - 1))/rho1^(5/2)
       -x(2)*((3*(2*mu + 2*x(1))*(mu - 1))/(2*rho1^(5/2)) - (3*mu*(2*mu + 2*x(1) - 2))/(2*rho2^(5/2))),(mu - 1)/rho1^(3/2) - mu/rho2^(3/2) + x(2)*((3*x(2)*mu)/rho2^(5/2) - (3*x(2)*(mu - 1))/rho1^(5/2)) + 1,x(2)*((3*mu*x(3))/rho2^(5/2) - (3*x(3)*(mu - 1))/rho1^(5/2))
       -x(3)*((3*(2*mu + 2*x(1))*(mu - 1))/(2*rho1^(5/2)) - (3*mu*(2*mu + 2*x(1) - 2))/(2*rho2^(5/2))),x(3)*((3*x(2)*mu)/rho2^(5/2) - (3*x(2)*(mu - 1))/rho1^(5/2)), x(3)*((3*mu*x(3))/rho2^(5/2) - (3*x(3)*(mu - 1))/rho1^(5/2)) + (mu - 1)/rho1^(3/2) - mu/rho2^(3/2)];
   

f = p*[x(4)
        x(5)
        x(6)
        -(1-mu)*(x(1)+mu)/rho13-mu*(x(1)+mu-1)/rho23+2*x(5)+x(1)     
        -((1-mu)/rho13+mu/rho23)*x(2)-2*x(4)+x(2)
        -((1-mu)/rho13+mu/rho23)*x(3)];

A = p*[zeros(3,3),eye(3,3)
       Uxx,[0 2 0;-2 0 0;0 0 0]];

B = p*[ 0,    0,    0;
        0,    0,    0;
        0,    0,    0;
     Tmax/m,    0,    0;
        0, Tmax/m,    0;
        0,    0, Tmax/m];
    

S = f/p;

R = -A*x-B*u;

end

%% PlotTransfer
% Plots the Transfer
function PlotTransfer(figNum,ax,mu,ScaleEarth,ScaleMoon,ScaleU,Re,Rm,topo,topoMoon,D,L1,L2,L3,L4,L5,X0,Xf,xg,x,t,u)

figure(figNum); hold on; axis equal; grid on; axis(ax);
xlabel('$\xi$','interpreter','latex'); ylabel('$\eta$','interpreter','latex'); %%title('Circular-Planar Synodic Frame Three-Body Problem','interpreter','latex','fontsize',20);
plot([-1.4 1.4],[0 0],'--k');
plot([0 0],[-1.4 1.4],'--k');
[Xsphere,Ysphere,Zsphere] = sphere(50);
earth = surf(-mu+Xsphere*ScaleEarth*Re/D,Ysphere*ScaleEarth*Re/D,Zsphere*ScaleEarth*Re/D);
set(earth,'facecolor','texturemap','cdata',topo,'edgecolor','none');
moon = surf(1-mu+Xsphere*ScaleMoon*Rm/D,Ysphere*ScaleMoon*Rm/D,Zsphere*ScaleMoon*Rm/D);
set(moon,'facecolor','texturemap','cdata',topoMoon,'edgecolor','none');
plot(0,0,'k.');

plot(L1(1),L1(2),'b*');
plot(L2(1),L2(2),'b*');
plot(L3(1),L3(2),'b*');
plot(L4(1),L4(2),'b*');
plot(L5(1),L5(2),'b*');

plot3(xg(:,1),xg(:,2),xg(:,3),'Linewidth',1,'Color','k','Linestyle','--');
plot3(X0(:,1),X0(:,2),X0(:,3),'g','Linewidth',2);
plot3(Xf(:,1),Xf(:,2),Xf(:,3),'r','Linewidth',2);
plot3(x(:,1),x(:,2),x(:,3),'m','Linewidth',2);
% for i = 1:length(u)
%     quiver3(x(i,1),x(i,2),x(i,3),ScaleU*u(i,1),ScaleU*u(i,2),ScaleU*u(i,3),'b');
% end


end

%% AnimateTransfer
%  Animates the transfer
function AnimateTransfer(figNum,ax,mu,ScaleEarth,ScaleMoon,ScaleU,Re,Rm,topo,topoMoon,D,L1,L2,L3,L4,L5,X0,Xf,x,t,u,Na,saveName)

% Animate Orbit:

f = figure(figNum); view(-35,20); hold on; axis equal; grid on; axis(ax);
set(f,'Color','w');
xlabel('$\xi$','interpreter','latex'); ylabel('$\eta$','interpreter','latex'); %%title('Circular-Planar Synodic Frame Three-Body Problem','interpreter','latex','fontsize',20);
plot([-1.4 1.4],[0 0],'--k');
plot([0 0],[-1.4 1.4],'--k');
[Xsphere,Ysphere,Zsphere] = sphere(50);
earth = surf(-mu+Xsphere*ScaleEarth*Re/D,Ysphere*ScaleEarth*Re/D,Zsphere*ScaleEarth*Re/D);
set(earth,'facecolor','texturemap','cdata',topo,'edgecolor','none');
moon = surf(1-mu+Xsphere*ScaleMoon*Rm/D,Ysphere*ScaleMoon*Rm/D,Zsphere*ScaleMoon*Rm/D);
set(moon,'facecolor','texturemap','cdata',topoMoon,'edgecolor','none');
plot(0,0,'k.');

plot(L1(1),L1(2),'b*');
plot(L2(1),L2(2),'b*');
plot(L3(1),L3(2),'b*');
plot(L4(1),L4(2),'b*');
plot(L5(1),L5(2),'b*');

plot3(X0(:,1),X0(:,2),X0(:,3),'g');
plot3(Xf(:,1),Xf(:,2),Xf(:,3),'r');

ta = linspace(t(1),t(end),Na);
xa = interp1(t,x,ta);
ua = interp1(t,u,ta);

pa = scatter3(xa(1,1),xa(1,2),xa(1,3),60,'k','filled');
pl = plot3(xa(1,1),xa(1,2),xa(1,3),'m');
qa = quiver3(xa(1,1),xa(1,2),xa(1,3),ScaleU*ua(1,1),ScaleU*ua(1,2),ScaleU*ua(1,3),'b');

for i = 2:length(ta)
    set(pa,'XData',xa(i,1),'YData',xa(i,2),'ZData',xa(i,3));
    set(pl,'XData',xa(1:i,1),'YData',xa(1:i,2),'ZData',xa(1:i,3));
    set(qa,'XData',xa(i,1),'YData',xa(i,2),'ZData',xa(i,3),'UData',ScaleU*ua(i,1),'VData',ScaleU*ua(i,2),'WData',ScaleU*ua(i,3));
    if ~strcmp(saveName,'noSave')
        frame = getframe(f);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 2
            imwrite(imind,cm,[cd '\' saveName '.gif'],'gif','DelayTime',0, 'Loopcount',inf); 
        else
            imwrite(imind,cm,[cd '\' saveName '.gif'],'gif','DelayTime',0,'WriteMode','append'); 
        end
    end
end

end