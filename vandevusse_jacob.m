%% Modelo do espaço de estados para o Reator de Van de Vusse

% COQ 792 - Controle de Processos

% Autor: Gustavo Luís Rodrigues Caldas

% Parâmetros do modelo
%   k10 k20 k30 E1 E2 E3 Ro cp DeltaH1 DeltaH2 DeltaH3 Kw Ar V 

% Variáveis Distúrbios
%  T0 - temperatura na carga
%  Ca0 - concentração de A na entrada do reator

% Variáveis de entrada
%  Vazão de alimentação ou ainda F/V (velocidade espacial - space velocity
%  Tk - temperatura na camisa

% Variáveis de estado
% Ca e Cb - concentração de A e B no reator
% T - temperatura no reator

function [A,B,C,D]=vandevusse_jacob(x_ss,u,p,d)
%Declarando o vetor de parâmetros p
k10 = p(1); 
k20 = p(2);
k30 = p(3);
E1  = p(4);
E2  = p(5);
E3  = p(6);
Ro  = p(7);
cp  = p(8); 
DeltaH1 = p(9);
DeltaH2 = p(10);
DeltaH3 = p(11);
Kw = p(12);
Ar = p(13);
V = p(14);

% Declarando as entradas u
Fv = u(1); % F/V
Tk = u(2);

%Distúrbios
Ca0 = d(1);
T0 = d(2);

% Declarando as variáveis de estado
Ca = x_ss(1,1);
Cb = x_ss(2,1);
T  = x_ss(3,1);

%Modelo
%dx(1,1) = Fv*(Ca0-x(1))-(k10)*exp(-E1/((x(3)+273.15)))*x(1) - (k30)*exp(-E3/((x(3)+273.15)))*(x(1))^2; 
%dx(2,1) = -Fv*x(2) + (k10)*exp(-E1/((x(3)+273.15)))*x(1) - (k20)*exp(-E2/((x(3)+273.15)))*x(2);
%dx(3,1) = (1/(Ro*cp))*((k10)*exp(-E1/((x(3)+273.15)))*x(1)*DeltaH1 + (k20)*exp(-E2/((x(3)+273.15)))*x(2)*DeltaH2 + DeltaH3*(k30)*exp(-E3/((x(3)+273.15)))*(x(1))^2) + Fv*(T0-x(3)) + Kw*Ar/(Ro*cp*V)*(Tk-x(3));

% Elemento A11 = dx1/dCa
A11 = -Fv -(k10)*exp(-E1/((x_ss(3)+273.15)))-2*(k30)*exp(-E3/((x_ss(3)+273.15)))*(x_ss(1));
% Elemento A12 = dx1/dCb
A12 = 0;
% Elemento A13 = dx3/dT
A13= -(k10)*(E1/(x_ss(3)+273.15)^2)*exp(-E1/((x_ss(3)+273.15)))*x_ss(1) - (k30)*(E3/(x_ss(3)+273.15)^2)*exp(-E3/((x_ss(3)+273.15)))*(x_ss(1))^2;

% Elemento A21 = dx2/dCa
A21 = (k10)*exp(-E1/((x_ss(3)+273.15)));
% Elemento A22 = dx2/dCb
A22 = -Fv-(k20)*exp(-E2/((x_ss(3)+273.15)));
% Elemento A23 = dx2/dT
A23 = (k10)*(E1/(x_ss(3)+273.15)^2)*exp(-E1/((x_ss(3)+273.15)))*x_ss(1)-(k20)*(E2/(x_ss(3)+273.15)^2)*exp(-E2/((x_ss(3)+273.15)))*x_ss(2);

% Elemento A31 = dx3/dCa
A31 = (1/(Ro*cp))*((k10)*exp(-E1/((x_ss(3)+273.15)))*DeltaH1 + 2*DeltaH3*(k30)*exp(-E3/((x_ss(3)+273.15)))*(x_ss(1)));
% Elemento A32 = dx3/dCb
A32 = (1/(Ro*cp))*(k20)*exp(-E2/((x_ss(3)+273.15)))*DeltaH2;
% Elemento A33 = dx3/dT
A33 = (1/(Ro*cp))*((k10)*(E1/(x_ss(3)+273.15)^2)*exp(-E1/((x_ss(3)+273.15)))*x_ss(1)*DeltaH1 + (k20)*(E2/(x_ss(3)+273.15)^2)*exp(-E2/((x_ss(3)+273.15)))*x_ss(2)*DeltaH2 + DeltaH3*(k30)*(E3/(x_ss(3)+273.15)^2)*exp(-E3/((x_ss(3)+273.15)))*(x_ss(1))^2)-Fv- Kw*Ar/(Ro*cp*V);

A = [A11 A12 A13;A21 A22 A23;A31 A32 A33]; %matriz de estados


% Elemento B11 = dx1/dFv
B11 = (Ca0-x_ss(1));
% Elemento B12 = dx1/dTk
B12 = 0;

% Elemento B11 = dx2/dFv
B21 = -x_ss(2);

% Elemento B22 = dx2/dTk
B22 = 0;

% Elemento B31 = dx3/dFv
B31 = (T0-x_ss(3));
% Elemento B32 = dx3/dTk
B32 = Kw*Ar/(Ro*cp*V);

B = [B11 B12;B21 B22;B31 B32]; %matriz das entradas

C = eye(3); % matriz das saídas 

%Matriz dos distúrbios
% Elemento Gama11 = dx1/dCa0
%Gama11 = Fv;
% Elemento D12 = dx1/dT0
%Gama12 = 0;
% Elemento Gama21 = dx2/dCa0
%Gama21 = 0;
% Elemento Gama22 = dx2/dT0
%Gama22 = 0;
% Elemento Gama31 = dx3/dCa0
%Gama31 = 0;
% Elemento Gama32 = dx2/dT0
%Gama32 = Fv;


D = [0 0;0 0;0 0];

end

