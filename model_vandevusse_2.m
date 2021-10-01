%% Modelo de simulação dinâmica para o Reator de Van der Vusse

% COQ 792 - Controle de Processos

% Autor: Gustavo Luís Rodrigues Caldas

% Esta é uma versão adaptada do arquivo "vandevusse.m" para o uso de s-function. 

% Parâmetros do modelo
%   k10 k20 k30 E1 E2 E3 Ro cp DeltaH1 DeltaH2 DeltaH3 Tk Kw Ar V 

% Variáveis perturbadoras (Entradas)
%  Vazão de alimentação ou ainda F/V (velocidade espacial - space velocity
%  T0 - temperatura na carga
%  Ca0 - concentração de A na entrada do reator
%  Tk - temperatura na camisa

% Variáveis de estado
% Ca e Cb - concentração de A e B no reator
% T - temperatura no reator


function dx=model_vandevusse_2(t,x,u)
%Declarando os parâmetros p
% Parâmetros
k10=1.287*10^12;
k20=1.287*10^12;
k30=9.043*10^9;
E1=9758.3;
E2=9758.3;
E3=8560;
Ro=0.9342;
cp=3.01;
DeltaH1=-4.20;
DeltaH2=11;
DeltaH3=41.85;  
Kw=4.032;
Ar=0.215;
V=10;

% Declarando as entradas u
Fv = u(1);
Tk = u(2);

%Distúrbios
Ca0 = u(3);
T0 = u(4);

% Declarando as variáveis de estado
Ca = x(1,1);
Cb = x(2,1);
T  = x(3,1);

%Modelo
dx(1,1) = Fv*(Ca0-x(1))-(k10)*exp(-E1/((x(3)+273.15)))*x(1) - (k30)*exp(-E3/((x(3)+273.15)))*(x(1))^2; 
dx(2,1) = -Fv*x(2) + (k10)*exp(-E1/((x(3)+273.15)))*x(1) - (k20)*exp(-E2/((x(3)+273.15)))*x(2);
dx(3,1) = (1/(Ro*cp))*((k10)*exp(-E1/((x(3)+273.15)))*x(1)*DeltaH1 + (k20)*exp(-E2/((x(3)+273.15)))*x(2)*DeltaH2 + DeltaH3*(k30)*exp(-E3/((x(3)+273.15)))*(x(1))^2)+ Fv*(T0-x(3)) + Kw*Ar/(Ro*cp*V)*(Tk-x(3));

end