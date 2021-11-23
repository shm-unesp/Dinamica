%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Análise Dinâmica de um Pêndulo Forçado - Diagramas de bifurcação
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'defaultaxesfontname','Euclid');
set(0,'defaultaxesfontsize',12)
set(0,'defaulttextfontsize',12)
set(0, 'defaultTextInterpreter', 'latex'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parâmetros de Simulação

l = 9.8;  % comprimento do fio [m]
g = 9.8;  % gravidade [m/s^2]
c = 0.5;  % amortecimento [N.s/m]

% Força
%A = 1.2; % condição que dá caos
%A = 0.5;    % amplitude da força de excitação
Omega = 2/3;  % frequência de excitação em rad/s
N  = 100000; % número de mostras
dt = 0.1; % período de amostragem
t = 0:dt:(N-1)*dt; % vetor tempo de simulação

A = 1:0.001:1.4;    % amplitudes testadas

for i = 1:length(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = A(i)*sin(Omega*t); % sinal de força de excitação

time = [];
Theta = [];
time = t;

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

% Ajustando para entre -pi e pi
Theta = mod(theta+pi,2*pi)-pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tirando o transiente
I=find (time< 3*pi*300);
time(I)=NaN;
Theta(I)=NaN;

% Tirando o valores que não estão em fase
Z=find(abs(rem(time, 2*pi/Omega)) > dt/2);
time(Z)=NaN;
Theta(Z)=NaN;
% Remove all NaN values from the array to reduce dataset size
time(isnan(time)) = [];
Theta(isnan(Theta)) = [];

figure(1)
plot(A(i),Theta,'.b'); hold on
xlabel('$A$','fontsize',12)
ylabel('$\theta$','fontsize',12)

end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%