%% Modelo de simulação dinâmica para o Reator de Van de Vusse

% COQ 792 - Controle de Processos

% Autor: Gustavo Luís Rodrigues Caldas

% Parâmetros do modelo
%   k10 k20 k30 E1 E2 E3 Ro cp DeltaH1 DeltaH2 DeltaH3 Tk Kw Ar V 

% Variáveis Distúrbios
%  T0 - temperatura na carga
%  Ca0 - concentração de A na entrada do reator

% Variáveis de entrada
%  Vazão de alimentação ou ainda F/V (velocidade espacial - space velocity
%  Tk - temperatura na camisa

% Variáveis de estado
% Ca e Cb - concentração de A e B no reator
% T - temperatura no reator


function dx=model_vandevusse(t,x,u,p,d)
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

% Pertubação em degrau

%if t>0.1
%   u(1) = u(1) + 50;
%end

% Declarando as entradas u
Fv = u(1);
Tk = u(2);

%Distúrbios
Ca0 = d(1);
T0 = d(2);

% Declarando as variáveis de estado
Ca = x(1,1);
Cb = x(2,1);
T  = x(3,1);

%Modelo
dx(1,1) = Fv*(Ca0-x(1))-(k10)*exp(-E1/((x(3)+273.15)))*x(1) - (k30)*exp(-E3/((x(3)+273.15)))*(x(1))^2; 
dx(2,1) = -Fv*x(2) + (k10)*exp(-E1/((x(3)+273.15)))*x(1) - (k20)*exp(-E2/((x(3)+273.15)))*x(2);
dx(3,1) = (1/(Ro*cp))*((k10)*exp(-E1/((x(3)+273.15)))*x(1)*DeltaH1 + (k20)*exp(-E2/((x(3)+273.15)))*x(2)*DeltaH2 + DeltaH3*(k30)*exp(-E3/((x(3)+273.15)))*(x(1))^2)+ Fv*(T0-x(3)) + Kw*Ar/(Ro*cp*V)*(Tk-x(3));

end