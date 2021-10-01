clear 
clc

%% Cálculo dos ganhos relativos para o Reator de Van de Vusse

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

%Entradas
Fv=160;
Tk=128.95;
u = [Fv;Tk];

%%Distúrbios
Ca0=5.1;
T0=130;
d = [Ca0;T0];

% Condições iniciais
% Vetor p
p = [k10; k20; k30; E1; E2; E3; Ro; cp; DeltaH1; DeltaH2; DeltaH3; Kw; Ar; V]; 
x0 = [Ca0;0;T0];

% Ganhos relativos

% Alocando memória

Fvspan = 0.5:0.5:160;
lambda11 = zeros(1,length(Fvspan));
lambda12 = zeros(1,length(Fvspan));
lambda21 = zeros(1,length(Fvspan));
lambda22 = zeros(1,length(Fvspan));

ganho11 = zeros(1,length(Fvspan));
ganho12 = zeros(1,length(Fvspan));
ganho21 = zeros(1,length(Fvspan));
ganho22 = zeros(1,length(Fvspan));


Nied = zeros(1,length(Fvspan));

for i=1:length(Fvspan)
    %Resolução do problema estacionário: dx/dt = 0
    u(1) = Fvspan(i);
    options = optimoptions('fsolve','Display','off');
    x_ss = fsolve(@(x) model_vandevusse(0,x,u,p,d),x0,options);
    [A,B,C,D]=vandevusse_jacob(x_ss,u,p,d);
    
    [num1,den1] = ss2tf(A,B,C,D,1);
    tf11= tf(num1(2,:),den1); %y1/u1
    tf21= tf(num1(3,:),den1); %y2/u1
    %avaliando em s=0
    gains11 = evalfr(tf11,0); 
    gains21 = evalfr(tf21,0);
    ganho11(i) = gains11; %colocando em uma nova matriz
    ganho21(i) = gains21;
    
    [num2,den2] = ss2tf(A,B,C,D,2);
    tf12 = tf(num2(2,:),den2); %y1/u2
    tf22 = tf(num2(3,:),den2); %y2/u2
    %avaliando em s=0
    gains12 = evalfr(tf12,0);
    gains22 = evalfr(tf22,0);
    ganho12(i) = gains12; %colocando em uma nova matriz
    ganho22(i) = gains22;
    
    ganhototal = [gains11 gains12;gains21 gains22]; %matriz G(0)
    
    %Cálculo do RGA
    R = transpose(inv(ganhototal));
    Lambda(1,1) = R(1,1)*ganhototal(1,1);
    Lambda(1,2) = R(1,2)*ganhototal(1,2);
    Lambda(2,1) = R(2,1)*ganhototal(2,1);
    Lambda(2,2) = R(2,2)*ganhototal(2,2);
    lambda11(i) = Lambda(1,1);
    lambda12(i) = Lambda(1,2);
    lambda21(i) = Lambda(2,1);
    lambda22(i) = Lambda(2,2);
    
     % Índice de Niederlinski

    %Determinante
    deter = det(ganhototal);

    %Cálculo dos coeficientes
    Nied(i) = (deter)/(gains11*gains22);
end

%Plotando
figure();
plot(Fvspan,lambda11,'b-');
xlabel({'$(F/V)_e$'},'Interpreter','latex')
ylabel({'$\lambda_{11}$'},'Interpreter','latex')

figure();
plot(Fvspan,Nied,'b-');
xlabel({'$(F/V)_e$'},'Interpreter','latex')
ylabel('{Indice de Niederlinski}','Interpreter','latex')

%set(gcf,'Renderer','zbuffer')
%print(gcf,'Niederlinski.jpg','-djpeg','-r300')

figure(4);
plot(Fvspan,lambda12,'b-');
xlabel({'$(F/V)_e$'},'Interpreter','latex')
ylabel({'$\lambda_{12}$'},'Interpreter','latex')

%figure(5);
%plot(Fvspan,lambda21,'b-');
%xlabel({'$(F/V)_e$'},'Interpreter','latex')
%ylabel({'$\lambda_{21}$'},'Interpreter','latex')

%figure(6);
%plot(Fvspan,lambda22,'b-');
%xlabel({'$(F/V)_e$'},'Interpreter','latex')
%ylabel({'$\lambda_{22}$'},'Interpreter','latex')





