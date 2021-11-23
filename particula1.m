%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Solução da Eq. do Movimento do Exemplo 8 do Livro do Prof. Ilmar Santos
% "Dinâmica de Sistemas Mecânicos" e Animação do Movimento da Partícula
%
% Autor: Samuel da Silva
% Data: 05/setembro/09
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;            % limpa a tela
clear;          % deleta todas as variáveis na memória

%--------------------------------------------------------------------------
% Parâmetros Geométricos e Físicos do Problema
%--------------------------------------------------------------------------

h = 1400;       % Altura da haste [m]
g = 9.8;        % Gravidade [m/s^2]
Omega = 1;     % Rotação da base móvel B1 [rad/s]

N = 1000;               % Numero de amostras temporais
T = 8;                  % Tempo da simulação [s]
deltaT = T/(N-1);       % Passo de integração da equação do movimento
t = 0:deltaT:(N-1)*deltaT;  % Vetor tempo

Theta1 = Omega*t;       % Deslocamento angular Theta1 [rad]
Theta2 = pi/4;          % Angulo Theta2 constante [rad]

%--------------------------------------------------------------------------
% Equação do Movimento: ddotx - Omega^2*cos^2(Theta2)x=gsin(Theta2)
%--------------------------------------------------------------------------

% Solução da EDO usando aproximação em série de Taylor de 1.º ordem (Método
% de Euler)

% Inicializando vetores de deslocamento, velocidade e aceleração
x = zeros(N,1);
dotx = zeros(N,1);
ddotx = zeros(N,1);

% Condições iniciais
x(1) = 0;
dotx(1) =0;
ddotx(1) = g*sin(Theta2) + Omega*Omega*(cos(Theta2)*cos(Theta2))*x(1);

xexata(1) = (g*sin(Theta2)/(Omega*Omega*cos(Theta2)*cos(Theta2)))*(cosh(Theta1(1)*cos(Theta2))-1);

for i=2:N       % Abre um laço de cálculo de t=t_0+deltat até t=t_0+Ndeltat
   
    % Aproximação de x(i) e dotx(i) usando valores em i-1
    
    x(i) = x(i-1)+deltaT*dotx(i-1);
    dotx(i) = dotx(i-1)+deltaT*ddotx(i-1);
    
    ddotx(i) = g*sin(Theta2) + Omega*Omega*(cos(Theta2)*cos(Theta2))*x(i);
    
    % Solução exata (Analítica)
    xexata(i) = (g*sin(Theta2)/(Omega*Omega*cos(Theta2)*cos(Theta2)))*(cosh(Theta1(i)*cos(Theta2))-1);
    
end

% Descrevendo os vetores deslocamento, velocidade e aceleração na base
% inercial

B2rBA = [x zeros(N,1) zeros(N,1)]';         % vetor posição aproximada(trajetória na base B2)
B2rBAex = [xexata' zeros(N,1) zeros(N,1)]';  % vetor posição exata (trajetória na base B2)
IrOB = [0 h 0]';        % vetor posição de O até B no sistema inercial

for i =1:N
    
    % Matrizes de Transformação
    
    T1 = [cos(Theta1(i))    0     -sin(Theta1(i));
                   0           1            0       ;
               sin(Theta1(i))   0          cos(Theta1(i))];
           
    T2 = [cos(Theta2)    -sin(Theta2)   0;
                sin(Theta2)    cos(Theta2)  0;
                0                 0         1];
            
    % Transformando trajetória da base móvel B2 para base inercial I
    
    IrOA(:,i) = IrOB + T1'*T2'*B2rBA(:,i);              % Trajetória base inercial aproximada
    IrOAex(:,i) = IrOB + T1'*T2'*B2rBAex(:,i);          % Trajetória base inercial exata
end

% Plotando no Espaço Trajetória Aproximada e Exata

figure(1)
plot3(IrOA(3,:),IrOA(1,:),IrOA(2,:)); grid on
xlabel('z')
ylabel('x')
zlabel('y')
saveas(1,'fig1.eps')

hold on

plot3(IrOAex(3,:),IrOAex(1,:),IrOAex(2,:),'r'); grid on
legend('Trajetória Aproximada','Trajetória Exata')
title('Trajetória realizada pela partícula A no sistema inercial')
AZ = 163;
EL = 32;
view(AZ,EL);

%--------------------------------------------------------------------------
% Visualização no Movimento da Partícula - Animação (Existe inúmeras formas
% de animar o movimento da partícula).
% Consulte help no Matlab dos comandos getframe, addframe, movie para maiores
% detalhes
%--------------------------------------------------------------------------

numframe = 50;    % número de quadros

% Mecanismo na posição inicial
%                        (4)
%                           \
%                            \
%                             \(3)
%                              |
%                              |
%                              |
%                              |
%                              |
%       (1)-----------------(2)
    
% Gravando o filme da animação    
mov = avifile('exemplo8.avi')

figure(2)
Xmax = 2000;
Xmin = -2000;

Ymin = -2000;
Ymax = 2000;

Zmin = 0;
Zmax = h;

set(gca,'nextplot','replacechildren');

for i=1:numframe 
    
    num = 20*i; % para controlar o tempo de cada quadro
    
  
    % Geometria do mecanismo
    d = sin(Theta2)*max(x); 
    l = cos(Theta2)*max(x);
   
    % Plotando mecanismo
    % Matriz de transformação da base inercial para movel B1
     T1 = [cos(Theta1(num))    0     -sin(Theta1(num));
                   0           1            0       ;
               sin(Theta1(num))   0          cos(Theta1(num))];
           
     %X1 = T1'*[h 0 0]';     % Para plotar o mecanismo girando com Omega
     %Z1 = T1'*[0 0 h]';
     %Y1 = T1'*[0 h 0]';

    % Coordenadas dos Nós do Mecanismo
    index(1,:) = T1'*[0 0 0]';          % Nó 1
    index(2,:) = T1'*[l 0 0]';          % Nó 2
    index(3,:) = T1'*[l h-d 0]';        % Nó 3
    index(4,:) = T1'*[0 h 0]';          % Nó 4

    elemento(1,:) = [1 2];         % Elemento 1 entre nó 1 e 2
    elemento(2,:) = [2 3];         % Elemento 2 entre nó 2 e 3
    elemento(3,:) = [3 4];         % Elemento 4 entre nó 3 e 4

    nele = 3;       % o mecanismo é composto por 3 barras
    
    for j=1:nele

        aux1 = elemento(j,1);
        aux2 = elemento(j,2);
        x1 = [index(aux1,1),index(aux1,2),index(aux1,3)];
        x2 = [index(aux2,1),index(aux2,2),index(aux2,3)];
        X = [x1 ;x2];
        plot3(X(:,3),X(:,1),X(:,2),'LineStyle','-','Color','r'); hold on; grid on
    end

    plot3(IrOA(3,:),IrOA(1,:),IrOA(2,:),'--k','linewidth',.5); 
    plot3(IrOA(3,num),IrOA(1,num),IrOA(2,num),'o','linewidth',2); 
    axis([Xmin Xmax Ymin Ymax Zmin Zmax]);
    xlabel('z')
	ylabel('x')
    zlabel('y')
	hold off
    Traj(i) = getframe;
    pause(.1)
    mov = addframe(mov,Traj(i));         % Salva a animação em um arquivo AVI
end

mov = close(mov);

%--------------------------------------------------------------------------

