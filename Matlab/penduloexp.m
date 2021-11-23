%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Análise Dinâmica de um Pêndulo Forçado - Expoente de Lyapunov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,"defaultaxesfontsize",15)
set(0,"defaulttextfontname","arial")
set(0,"defaulttextfontsize",15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parâmetros de Simulação

l = 9.8;  % comprimento do fio [m]
g = 9.8;  % gravidade [m/s^2]
c = 0.5;  % amortecimento [N.s/m]

N  = 15000; % número de mostras
dt = 0.1; % período de amostragem
t = 0:dt:(N-1)*dt; % vetor tempo de simulação

% Força
A = 1.2; % condição que dá caos
%A = 0.5;    % amplitude da força de excitação
Omega = 2/3;  % frequência de excitação em rad/s
F = A*sin(Omega*t); % sinal de força de excitação

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotando as Linhas do Campo

y1 = linspace(-6,9,20);
y2 = linspace(-4,3,20);

% creates two matrices one for all the x-values on the grid, and one for
% all the y-values on the grid. Note that x and y are matrices of the same
% size and shape, in this case 20 rows and 20 columns
[x,y] = meshgrid(y1,y2);
u = zeros(size(x));
v = zeros(size(x));
% we can use a single loop over each element to compute the derivatives at
% each point (y1, y2)
%t=0; % we want the derivatives at each point at t=0, i.e. the starting time

for i = 1:numel(x)
    Yprime(1) =y(i);
    Yprime(2) = F(1)-(g/l)*sin(x(i));
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end

%figure(1)
%quiver(x,y,u,v,'r'); hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pêndulo Livre
% condições iniciais  

theta(1) = pi/4;
dottheta(1) = 0;
ddottheta(1) = F(1)-(g/l)*sin(theta(1)) - c*dottheta(1);
  
for k = 2:N
  % Aproximando usando o método de Runge-Kutta de 4.º ordem          
  % Parâmetros k para aproximar theta(i) e dottheta(i)
	% vetor de estados é z = [theta dottheta]'
	% dz = F(z); dz(1) = dotz(1) = dottheta
	% dz(2) = -(g/l)*sin(theta) = ddottheta

  % Cálculo das constantes
	k1a = dottheta(k-1);
	k1b = F(k-1)-(g/l)*sin(theta(k-1)) - c*dottheta(k-1); % ddotheta(i-1)

	k2a = dottheta(k-1) + dt*0.5*k1b; 
	k2b = F(k-1)-(g/l)*sin(theta(k-1) + dt*0.5*k1a) - c*k2a;

	k3a = dottheta(k-1) + dt*0.5*k2b;
	k3b =  F(k-1)-(g/l)*sin(theta(k-1) +dt*0.5*k2a)-c*k3a;

	k4a = dottheta(k-1) + dt*k3b;
	k4b =  F(k-1)-(g/l)*sin(theta(k-1) + dt*k3a)-c*k4a;

	% aproximação de theta(i)
	theta(k) = theta(k-1) + dt*(k1a+2*k2a+2*k3a+k4a)/6;
	% aproximação de dottheta(i)
 	dottheta(k) = dottheta(k-1) + dt*(k1b+2*k2b+2*k3b+k4b)/6;
  % aproximação de ddotheta(i)
  ddottheta(k) = F(k)-(g/l)*sin(theta(k)) - c*dottheta(k);

end

theta1 = theta;
  
figure(1)
plot(t,theta1,'linewidth',6,'-b'); hold on
xlim([0 100]);
ylabel('$\theta$ [rad]')
xlabel('$Tempo [s]')
%saveas(1,'theta.tex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesma condição inicial com leve diferença

theta(1) = pi/4 + 0.001;
dottheta(1) = 0;
ddottheta(1) = F(1)-(g/l)*sin(theta(1)) - c*dottheta(1);
  
for k = 2:N
  % Aproximando usando o método de Runge-Kutta de 4.º ordem          
  % Parâmetros k para aproximar theta(i) e dottheta(i)
	% vetor de estados é z = [theta dottheta]'
	% dz = F(z); dz(1) = dotz(1) = dottheta
	% dz(2) = -(g/l)*sin(theta) = ddottheta

  % Cálculo das constantes
	k1a = dottheta(k-1);
	k1b = F(k-1)-(g/l)*sin(theta(k-1)) - c*dottheta(k-1); % ddotheta(i-1)

	k2a = dottheta(k-1) + dt*0.5*k1b; 
	k2b = F(k-1)-(g/l)*sin(theta(k-1) + dt*0.5*k1a) - c*k2a;

	k3a = dottheta(k-1) + dt*0.5*k2b;
	k3b =  F(k-1)-(g/l)*sin(theta(k-1) +dt*0.5*k2a)-c*k3a;

	k4a = dottheta(k-1) + dt*k3b;
	k4b =  F(k-1)-(g/l)*sin(theta(k-1) + dt*k3a)-c*k4a;

	% aproximação de theta(i)
	theta(k) = theta(k-1) + dt*(k1a+2*k2a+2*k3a+k4a)/6;
	% aproximação de dottheta(i)
 	dottheta(k) = dottheta(k-1) + dt*(k1b+2*k2b+2*k3b+k4b)/6;
  % aproximação de ddotheta(i)
  ddottheta(k) = F(k)-(g/l)*sin(theta(k)) - c*dottheta(k);

end
  
theta2 = theta;
  
figure(1)
plot(t,theta2,'linewidth',10,'--r');
xlim([0 100]);
ylim([-30 30])
ylabel('$\theta$ [rad]')
xlabel('Tempo [s]')
saveas(1,'thetacc.tex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expoente de Lyapunov - Estimativa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltaTheta = theta1 - theta2;
% Delta Theta = exp(lambda t);
% ln(Delta Theta) = lambda t
lambda = log(deltaTheta)./t;
%[11.2 -6.8]
%[40.8 -0.26]
k = (+6.8-0.26)/(40.8-11.2);
y = k*t-10;
% Visualizando a sensibilidade as condições iniciais
figure(2)
plot(t,log(deltaTheta),'linewidth',8); hold on;
plot(t,y,'linewidth',10,'--r')
xlim([0 100]);
ylim([-15 5])
ylabel('$\ln(\Delta \theta)$ [rad]')
xlabel('Tempo [s]')
saveas(2,'deltathetacc.tex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%