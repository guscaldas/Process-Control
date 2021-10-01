clear
clc

%% Simulação dinâmica para o Reator de Van de Vusse

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

%Vetor p
p = [k10; k20; k30; E1; E2; E3; Ro; cp; DeltaH1; DeltaH2; DeltaH3; Kw; Ar; V]; 

% Variáveis de entrada
Fv=60;
Tk=128.95;
u = [Fv;Tk];

%Distúrbios
Ca0=5.1;
T0=130;
d = [Ca0;T0];
%% Resolução do problema estacionário: dx/dt = 0

% Condições iniciais
x0 = [Ca0;0;T0];
options = optimoptions('fsolve','Display','off');
x_ss = fsolve(@(x) model_vandevusse(0,x,u,p,d),x0,options);
%% Mapeamento estacionário de x vs u

Fvspan = 160:-0.5:0.1; % de 0.1 a 160 h^-1
for i=1:length(Fvspan)
    u(1) = Fvspan(i);  
    x = fsolve(@(x) model_vandevusse(0,x,u,p,d),x0,options);
    ca(i)=x(1,1);
    cb(i)=x(2,1);
    T(i)=x(3,1);
end

%Plotando concentração de cb
figure(1)
plot(Fvspan,cb,'b-')
xlabel({'$(F/V)_e$'},'Interpreter','latex')
ylabel({'$C_{b,e}$'},'Interpreter','latex')

%Salvando figuras

%set(gcf,'Renderer','zbuffer')
%print(gcf,'Mapeamento_EE_van_de_vusse.jpg','-djpeg','-r300')

%% Problema dinâmico - caso contínuo

% definição do intervalo de integração
tstep = 0.001; 
tmax = 0.3; %horas
tspan = 0:tstep:tmax;

% Condições iniciais
x0 = [Ca0;0;T0];

%Distúrbios
Ca0=5.1;
T0=130;
d = [Ca0;T0];

x = zeros(3);

% Variáveis de entrada
%Fv_range=[20,45,60,120];
Tk=128.95;

%for i=1:length(Fv_range) %plotando para múltiplos Fv
Fv = 60; %somente um
u = [Fv;Tk]; %vetor de entradas

[t,x] = ode45(@(t,x) model_vandevusse(t,x,u,p,d),tspan, x0);

Ca = x(:,1);
Cb = x(:,2);
T  = x(:,3);


subplot(3,1,1)
plot(t,Ca,'linewidth',1.5,'DisplayName','10')
xlabel({'$t$'},'Interpreter','latex');
ylabel({'$C_{a}$'},'Interpreter','latex')
hold on
plot([t(1);t(end)],[x_ss(1);x_ss(1)],'k-','linewidth',1.5)

subplot(3,1,2)
plot(t,Cb,'linewidth',1.5)
xlabel({'$t$'},'Interpreter','latex');
ylabel({'$C_{b}$'},'Interpreter','latex')
hold on
plot([t(1);t(end)],[x_ss(2);x_ss(2)],'k-','linewidth',1.5)

subplot(3,1,3)
plot(t,T,'linewidth',1.5)
xlabel({'$t$'},'Interpreter','latex');
ylabel({'$T$'},'Interpreter','latex')
hold on
plot([t(1);t(end)],[x_ss(3);x_ss(3)],'k-','linewidth',1.5)

%end

%Salvando figuras
%set(gcf,'Renderer','zbuffer')
%print(gcf,'Degrau dinâmico.jpg','-djpeg','-r300')

