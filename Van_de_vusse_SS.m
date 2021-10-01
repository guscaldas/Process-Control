clear 
clc

%% Simulação usando FT para o Reator de Van de Vusse

% COQ 792 - Controle de Processos

% Autor: Gustavo Luís Rodrigues Caldas

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

%%Distúrbios
Ca0=5.1;
T0=130;
d = [Ca0;T0];

% Condições iniciais
% Vetor p
p = [k10; k20; k30; E1; E2; E3; Ro; cp; DeltaH1; DeltaH2; DeltaH3; Kw; Ar; V]; 
x0 = [Ca0;0;T0];

%% Obtenção das funções de transferência - Letra b
Fv_loop = [20,45,60,120];
for i=Fv_loop
    %i
    %Entradas
    Fv=i;
    Tk=128.95;
    u = [Fv;Tk];
    
    %Resolução em torno de um estado estacionário
    options = optimoptions('fsolve','Display','off');
    x_ss = fsolve(@(x) model_vandevusse(0,x,u,p,d),x0,options);

    [A,B,C,D]=vandevusse_jacob(x_ss,u,p,d);

    %FTs
    [num1,den1] = ss2tf(A,B,C,D,1);

    tf11 = tf(num1(2,:),den1); %y1/u1
    tf21 = tf(num1(3,:),den1); %y2/u1
    
    [num2,den2] = ss2tf(A,B,C,D,2);
    tf12 = tf(num2(2,:),den2); %y1/u2
    tf22 = tf(num2(3,:),den2); %y2/u2

end
