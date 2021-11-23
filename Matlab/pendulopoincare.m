%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Análise Dinâmica de um Pêndulo Forçado - Secções de Poincaré
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set(0,"defaultaxesfontsize",15)
%set(0,"defaulttextfontname","arial")
%set(0,"defaulttextfontsize",15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parâmetros de Simulação

l = 9.8;  % comprimento do fio [m]
g = 9.8;  % gravidade [m/s^2]
c = 0.5;  % amortecimento [N.s/m]

N  = 150000; % número de mostras
dt = 0.1; % período de amostragem
t = 0:dt:(N-1)*dt; % vetor tempo de simulação

% Força
A = 1.2; % condição que dá caos
%A = 0.5;    % amplitude da força de excitação
Omega = 2/3;  % frequência de excitação em rad/s
F = A*sin(Omega*t); % sinal de força de excitação

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Only plot omega and theta point when omega is in phase with the driving force Omega_

theta1 = mod(theta+pi, 2*pi)-pi; % convertando tudo no intervalo -pi a pi
%theta1 = theta;
I=find(abs(rem(t,2*pi/Omega)) > dt/2);
theta1(I)=NaN;
dottheta(I)=NaN;

figure(1)
scatter (theta1(1500:end),dottheta(1500:end),2); hold on
%plots the numerical solution 
plot(theta1(1500:end),dottheta(1500:end),'.')
xlabel('$\theta$ [rad]')
ylabel('$\dot\theta$ [rad/s]')
xlim([-4 4])
ylim([-2 1])
saveas(1,'poincareccaos.tex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%