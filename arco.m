%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Anlise Dinmica de uma Partcula em um Mecanismo na Forma de Arco com
%  duas rotaes consecutivas em torno de eixos diferentes
%
% Autor: Samuel da Silva, 
% Data: 08/setembro/09
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;            % limpa a tela
clear;          % deleta todas as variveis na memria

%--------------------------------------------------------------------------
% Parmetros Geomtricos e Fsicos do Problema
%--------------------------------------------------------------------------

R = 0.5;        % Raio do Arco [m]
g = 9.8;        % Gravidade [m/s^2]
Omega = 5;      % Rotao da base mvel B1 [rad/s]

N = 1000;                   % Numero de amostras temporais
T = 8;                      % Tempo da simulao [s]
deltaT = T/(N-1);           % Passo de integrao da equao do movimento
t = 0:deltaT:(N-1)*deltaT;  % Vetor tempo

Theta = Omega*t;       % Deslocamento angular Theta1 [rad]

% Vetor posio da partcula A na base B2
B2r_OA  = [0 R 0]';

%--------------------------------------------------------------------------
% Equao do Movimento No-Linear: 
% ddotalpha = sin(alpha)*(Omega^2*cos(alpha)-g/R)
%--------------------------------------------------------------------------

% Soluo da EDO usando aproximao em srie de Taylor de 1. ordem (Mtodo
% de Euler)

% Inicializando vetores de deslocamento, velocidade e acelerao angular
alpha = zeros(N,1);
dotalpha= zeros(N,1);
ddotalpha = zeros(N,1);

% Condies iniciais
alpha(1) = pi/2;
dotalpha(1) =0;
ddotalpha(1) = sin(alpha(1))*(Omega*Omega*cos(alpha(1))-g/R);

for i=2:N
    
    % Aproximao de alpha(i) e dotalpha(i) usando valores em i-1
    alpha(i) = alpha(i-1)+deltaT*dotalpha(i-1);
    dotalpha(i) = dotalpha(i-1)+deltaT*ddotalpha(i-1);
    
    % Aproximaa de ddotalpha(i)
    ddotalpha(i) = sin(alpha(i))*(Omega*Omega*cos(alpha(i))-g/R);

end

%--------------------------------------------------------------------------
% Representao da trajetria na base Inercial

for i=1:N
    
    % Matriz de Transformao de I para B1
    Ttheta = [cos(Theta(i)) 0 -sin(Theta(i));
           0            1      0       ;
          sin(Theta(i)) 0      cos(Theta(i))];
          
   % Matriz de Transformao de B1 para B2
   Talpha = [cos(alpha(i))  sin(alpha(i))  0;
         -sin(alpha(i))   cos(alpha(i)) 0;
         0                 0             1];
     
   Ir_OA(:,i) = Ttheta'*Talpha'*B2r_OA; 

end

%--------------------------------------------------------------------------
% Anlise grfica dos resultados

figure(1)
subplot(311)
plot(t,alpha,'linewidth',1)
ylabel('Desl. Angular \alpha [rad]');
xlabel('Tempo [s]');

subplot(312)
plot(t,dotalpha,'linewidth',1)
ylabel('Vel. Angular [rad/s]');
xlabel('Tempo [s]');

subplot(313)
plot(t,ddotalpha,'linewidth',1)
ylabel('Acel. Angular [rad/s^2]');
xlabel('Tempo [s]');

saveas(1,'nomefig1.eps')

figure(2)
title('Orbitas')
subplot(121)
plot(Ir_OA(1,:),Ir_OA(3,:))
xlabel('X')
ylabel('Z')

subplot(122)
plot(Ir_OA(1,:),Ir_OA(2,:))
xlabel('X')
ylabel('Y')

figure(3)
plot3(Ir_OA(1,:),Ir_OA(3,:),-Ir_OA(2,:)); grid on
xlabel('X')
ylabel('Z')
zlabel('Y')
AZ = -32; EL=54;
view(AZ,EL)
title('Orbita 3D do movimento da partcula no sistema inercial')

%--------------------------------------------------------------------------
% Animao da Orbita realizada pela Parcula A
%--------------------------------------------------------------------------

% Gravando filme da animao
%mov = avifile('bolaarco.avi')

%numframe = 40;  % numero de frames (quadros) usados na animao

%figure(4)
%Xmin = -R; Xmax = R;
%Ymin = -R;  Ymax = R;
%Zmin = -R; Zmax = 0;

% criando o arco no espao
%angulo = 0:0.1:pi;
%r_arco = -R*exp(j*angulo);
%R_arco = [real(r_arco); imag(r_arco); zeros(1,length(r_arco))];

%set(gca,'nextplot','replacechildren');

%for i=1:numframe
   
%   num = 20*i;   % para controlar o tempo de cada quadro
   
   % plota a trajetria
 %  plot3(Ir_OA(3,:),Ir_OA(1,:),-Ir_OA(2,:),'--b','linewidth',.5); hold on; grid on 
   
   % plota arco girando na base B2
    % Matriz de Transformao de I para B1
  %  Ttheta = [cos(Theta(num)) 0 -sin(Theta(num));
   %        0            1      0       ;
    %      sin(Theta(num)) 0      cos(Theta(num))];
     % 
   %IR_arco = Ttheta'*R_arco;
  % plot3(IR_arco(3,:),IR_arco(1,:),IR_arco(2,:),'k');
   
   % plota partcula no instante num
   %%plot3(Ir_OA(3,num),Ir_OA(1,num),-Ir_OA(2,num),'ro','linewidth',2); 
   %axis([Xmin Xmax Ymin Ymax Zmin Zmax]);
   %xlabel('z')
   %ylabel('x')
   %zlabel('y')
   %view(AZ,EL)
   %hold off
   %Traj(i) = getframe;
   %pause(.1)
   %mov = addframe(mov,Traj(i)); 
   
%end

%mov = close(mov);

%--------------------------------------------------------------------------






