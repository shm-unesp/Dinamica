%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Anlise Dinmica de um sistema com duas rotaes consecutivas em torno
%  de eixos diferentes
%
% Autor: Samuel da Silva
% Data: 09/outubro/2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;            % limpa a tela
clear;          % deleta todas as variveis na memria

%--------------------------------------------------------------------------
% Parmetros Geomtricos e Fsicos do Problema
%--------------------------------------------------------------------------

Xmin = -0.3;
Xmax = 0.3;
Ymin = -0.3;
Ymax = 0.3;
Zmin = -0.2;
Zmax = 0.2;

r = 0.1;        % Barra OA [m]
l = 0.2;        % Comprimento do fio AB [m] - sistema  simtrico
g = 9.8;        % Gravidade [m/s^2]
omega = -1;     % (negativo para ficar igual a P1) Rotao da base mvel B1 [rad/s]

N = 1000;                   % Numero de amostras temporais
T = 3;                      % Tempo da simulao [s]
deltaT = T/(N-1);           % Passo de integrao da equao do movimento
t = 0:deltaT:(N-1)*deltaT;  % Vetor tempo

Alpha = omega*t;       % Deslocamento angular Alpha [rad]

% Vetor posio OA na base mvel B1
B1r_OA  = [r 0 0]';

% Vetor posio AB na base mvel B2
B2r_AB = [0 -l 0]';

%--------------------------------------------------------------------------
% Equao do Movimento No-Linear na base mvel B2
%--------------------------------------------------------------------------

% Soluo da EDO usando aproximao em srie de Taylor de 1. ordem (Mtodo
% de Euler)

% Inicializando vetores de deslocamento, velocidade e acelerao angular
theta = zeros(N,1);
dottheta= zeros(N,1);
ddottheta = zeros(N,1);

% Condies iniciais
% Partcula B
theta(1) = pi/8;
dottheta(1) =0;
ddottheta(1) =  (r/l)*omega*omega*cos(theta(1)) + omega*omega*sin(theta(1))*cos(theta(1))-(g/l)*sin(theta(1));

% Partcula B'
theta1(1) = pi/4;
dottheta1(1) =0;
ddottheta1(1) =  (r/l)*omega*omega*cos(theta1(1)) + omega*omega*sin(theta1(1))*cos(theta1(1))-(g/l)*sin(theta1(1));

for i=2:N
    
    % Equao No-linear para partcula B:
    
    % Aproximao de theta(i) e dottheta(i) usando valores em i-1
    theta(i) = theta(i-1)+deltaT*dottheta(i-1);
    dottheta(i) = dottheta(i-1)+deltaT*ddottheta(i-1);
    
    % Aproximaa de ddottheta(i)
    ddottheta(i) = (r/l)*omega*omega*cos(theta(i-1)) + omega*omega*sin(theta(i-1))*cos(theta(i-1))-(g/l)*sin(theta(i-1));

    % Equao No-linear para partcula B':
    
    % Aproximao de theta(i) e dottheta(i) usando valores em i-1
    theta1(i) = theta1(i-1)+deltaT*dottheta1(i-1);
    dottheta1(i) = dottheta1(i-1)+deltaT*ddottheta1(i-1);
    
    % Aproximaa de ddottheta(i)
    ddottheta1(i) = (r/l)*omega*omega*cos(theta1(i-1)) + omega*omega*sin(theta1(i-1))*cos(theta1(i-1))-(g/l)*sin(theta1(i-1));
    
end

%--------------------------------------------------------------------------
% Representao da trajetria na base inercial

for i=1:N
    
    % Matriz de transformao de I para B1
    T_omega = [cos(Alpha(i))    0      -sin(Alpha(i));
               0                1           0;
               sin(Alpha(i))    0         cos(Alpha(i))];
    
    
    % Matriz de transformanao de B1 para B2
    T_theta = [cos(theta(i))     sin(theta(i))    0;
               -sin(theta(i))     cos(theta(i))   0;
               0                      0           1];

% Trajetria da massa B no sistema inercial
Ir_OB(:,i) = T_omega'*B1r_OA + T_omega'*T_theta'*B2r_AB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Matriz de transformao de I para B1
    T_omega = [cos(Alpha(i))    0      -sin(Alpha(i));
               0                1           0;
               sin(Alpha(i))    0         cos(Alpha(i))];
    
    
    % Matriz de transformanao de B1 para B3
    T_theta1 = [cos(theta1(i))     sin(theta1(i))    0;
               -sin(theta1(i))     cos(theta1(i))   0;
               0                      0           1];

% Trajetria da massa B' no sistema inercial
Ir_OBl(:,i) = T_omega'*B1r_OA + T_omega'*T_theta1'*B2r_AB;

end

%--------------------------------------------------------------------------
% Anlise grfica dos resultados
%--------------------------------------------------------------------------

figure(1)
subplot(311)
plot(t,theta,'linewidth',1)
ylabel('Desl. Angular \theta [rad]');
xlabel('Tempo [s]');

subplot(312)
plot(t,dottheta,'linewidth',1)
ylabel('Vel. Angular [rad/s]');
xlabel('Tempo [s]');

subplot(313)
plot(t,ddottheta,'linewidth',1)
ylabel('Acel. Angular [rad/s^2]');
xlabel('Tempo [s]');

figure(2)
plot3(Ir_OB(1,:),Ir_OB(3,:),Ir_OB(2,:)); grid on; hold on
plot3(-Ir_OBl(1,:),-Ir_OBl(3,:),Ir_OBl(2,:),'r'); grid on
xlabel('X')
ylabel('Z')
zlabel('Y')
legend('B','Blinha')
AZ = 32; EL=54;
view(AZ,EL)
title('Orbita 3D do movimento da partcula B e Blinha no sistema inercial')

%--------------------------------------------------------------------------
% Animao da Orbita realizada pelas Parculas B e B'
%--------------------------------------------------------------------------

% Gravando filme da animao
%mov = avifile('centrifugo.avi')

%numframe = 40;  % numero de frames (quadros) usados na animao

%figure(4)

%for i=1:numframe 
   
 %   num = 20*i; % para controlar o tempo de cada quadro
    
    % Plota trajetria a ser executada por B:
  %  plot3(Ir_OB(1,:),Ir_OB(3,:),Ir_OB(2,:),'--k','linewidth',.5); hold on;
    
    % Plota trajetria a ser executada por B':
   % plot3(-Ir_OBl(1,:),-Ir_OBl(3,:),Ir_OBl(2,:),'--g','linewidth',.5); hold on;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plota a barra OA girando na base B1
    % Matriz de Transformao de I para B1
    %T1 = [cos(Alpha(num))    0      -sin(Alpha(num));
     %          0                1           0;
      %         sin(Alpha(num))    0         cos(Alpha(num))];
           
    % Matriz de transformanao de B1 para B2
    %T2 = [cos(theta(num))     sin(theta(num))    0;
     %          -sin(theta(num))     cos(theta(num))   0;
      %         0                      0           1];
           
    % Matriz de transformanao de B2 para B3
    %T2l = [cos(theta1(num))     sin(theta1(num))    0;
     %          -sin(theta1(num))     cos(theta1(num))   0;
      %         0                      0           1];   
                   
    % Coordenadas dos ns da Barra
    %index(1,:) = T1'*[0 0 0]';          % N 1
    %index(2,:) = T1'*B1r_OA;          % N 2
    
    %elemento(1,:) = [1 2];         % Elemento 1 entre n 1 e 2 (OA)
  %  elemento(2,:) = [2 3];         % Elemento 2 entre n 2 e 3 (AB)
   
   %     aux1 = elemento(1,1);
    %    aux2 = elemento(1,2);
     %   x1 = [index(aux1,1),index(aux1,2),index(aux1,3)];
      %  x2 = [index(aux2,1),index(aux2,2),index(aux2,3)];
       % OA = [x1 ;x2];
        
   % Coordenadas dos ns do fio
        %index(2,:) = T1'*B1r_OA;
        %index(3,:) = Ir_OB(:,num);

        %aux3 = elemento(2,1);
        %aux4 = elemento(2,2);
        %y1 = [index(aux3,1),index(aux3,2),index(aux3,3)];
        %y2 = [index(aux4,1),index(aux4,2),index(aux4,3)];
        %AB = [y1 ;y2];
       
        %plot3(OA(:,1),OA(:,3),OA(:,2),'LineStyle','-','Color','r'); hold on; grid on
        %plot3(AB(:,1),AB(:,3),AB(:,2),'LineStyle','-','Color','r'); hold on; grid on
        
    %%%%%%%%%%%%%%%%%%%
        % Coordenadas dos ns da Barra
    %index(1,:) = T1'*[0 0 0]';          % N 1
    %index(2,:) = T1'*B1r_OA;          % N 2
    
    %elemento(1,:) = [1 2];         % Elemento 1 entre n 1 e 2 (OA)
    %elemento(2,:) = [2 3];         % Elemento 2 entre n 2 e 3 (AB)
   
     %   aux1 = elemento(1,1);
      %  aux2 = elemento(1,2);
       % x1 = [index(aux1,1),index(aux1,2),index(aux1,3)];
        %x2 = [index(aux2,1),index(aux2,2),index(aux2,3)];
        %OA = [x1 ;x2];
        
   % Coordenadas dos ns do fio
        %index(2,:) = T1'*B1r_OA;
        %index(3,:) = Ir_OBl(:,num);

        %aux3 = elemento(2,1);
        %aux4 = elemento(2,2);
        %y1 = [index(aux3,1),index(aux3,2),index(aux3,3)];
        %y2 = [index(aux4,1),index(aux4,2),index(aux4,3)];
        %AB = [y1 ;y2];
        
        %plot3(-OA(:,1),-OA(:,3),OA(:,2),'LineStyle','-','Color','r'); hold on; grid on
        %plot3(-AB(:,1),-AB(:,3),AB(:,2),'LineStyle','-','Color','r'); hold on; grid on
 
    % Plota movimento da partcula B:
    %plot3(Ir_OB(1,num),Ir_OB(3,num),Ir_OB(2,num),'o','linewidth',2); 
    
    % Plota movimento da partcula B':
   % plot3(-Ir_OBl(1,num),-Ir_OBl(3,num),Ir_OBl(2,num),'o','linewidth',2);
    
    %axis([Xmin Xmax Ymin Ymax Zmin Zmax]);
    %xlabel('X')
    %ylabel('Z')
    %zlabel('Y')
    %AZ = 32; EL=54;
    %view(AZ,EL)
	%hold off
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   % Traj(i) = getframe;
    %pause(.1)
    %mov = addframe(mov,Traj(i));         % Salva a animao em um arquivo AVI
    
%end

%mov = close(mov);


%--------------------------------------------------------------------------






