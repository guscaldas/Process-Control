clear 
clc

%% Simulação usando FT para o Reator de Van de Vusse

% COQ 792 - Controle de Processos

% Autor: Gustavo Luís Rodrigues Caldas

% Numeradores e denominadores 
num = [12.8 -18.9 6.6 -19.4];
den = [16.7 1; 21.0 1; 10.9 1;14.4 1];

%% Obtenção das funções de transferência - Letra b
tf11 = tf(num(1),den(1,:),'IODelay', 1);
tf12 = tf(num(2),den(2,:),'IODelay', 3);
tf21 = tf(num(3),den(3,:), 'IODelay', 7); 
tf22 = tf(num(4),den(4,:), 'IODelay', 3); 
Gp = [tf11 tf12;tf21 tf22];
%% Cálculo do RGA - Letra c
%avaliando em s=0
gains11 = evalfr(tf11,0); 
gains12 = evalfr(tf12,0);
gains21 = evalfr(tf21,0);
gains22 = evalfr(tf22,0);
ganhototal = [gains11 gains12;gains21 gains22]; %matriz G(0)

%Cálculo do RGA
Lambda = zeros(2); %alocando memória
R = transpose(inv(ganhototal));
Lambda(1,1) = R(1,1)*ganhototal(1,1);
Lambda(1,2) = R(1,2)*ganhototal(1,2);
Lambda(2,1) = R(2,1)*ganhototal(2,1);
Lambda(2,2) = R(2,2)*ganhototal(2,2);

%Distúrbios
% Numeradores e denominadores 
numd = [3.8;4.9];
dend = [14.9 1; 13.2 1];

tfd1 = tf(numd(1),dend(1,:),'IODelay',8.1);
tfd2 = tf(numd(2),dend(2,:),'IODelay', 3.4);

%% Sintonização ZN para uma FOPDT
%Contrária a Bristol
%Caso u1-y2
tau21 = 10.9;
theta21 = 7;
K21 = 6.6;
[Kc21,ti21] = ZN(tau21,theta21,K21);

%Caso u2-y1
tau12 = 21.0;
theta12 = 3;
K12 = -18.9;
[Kc12,ti12] = ZN(tau12,theta12,K12);

%Favorável a Bristol
%Caso u1-y1
tau11 = 16.7;
theta11 = 1;
K11 = 12.8;
[Kc11,ti11] = ZN(tau11,theta11,K11);

%Caso u2-y2
tau22 = 14.4;
theta22 = 3;
K22 = -19.4;
[Kc22,ti22] = ZN(tau22,theta22,K22);

%McAvoy
Kc11_McAv = Kc11*(Lambda(1,1) - sqrt(Lambda(1,1)^2 - Lambda(1,1)));
Kc22_McAv = Kc22*(Lambda(2,2) - sqrt(Lambda(2,2)^2 - Lambda(2,2)));

%% Desacopladores

%Desacopladores estáticos
D = (ganhototal)\[gains11 0;0 gains22];
Gl1_EE = D(1,2)/D(1,1); %D11 = D22
Gl2_EE = D(2,1)/D(2,2);

%Pode-se conferir que:
%Gl1_EE = - gains12/gains11;
%Gl2_EE = - gains21/gains22;

%Sintonizando
tf1_ef = tf11 + tf12*Gl2_EE;
[Gm,~,Wcg,~] = margin(tf1_ef);
Kc11_DE = Gm/2.2;
ti11_DE = (2*pi()/Wcg)/1.2;

tf2_ef = tf22 + tf21*Gl1_EE;
[Gm,~,Wcg,~] = margin(tf2_ef);
Kc22_DE = -Gm/2.2;
ti22_DE = (2*pi()/Wcg)/1.2;

%Desacoplador dinâmico simplificado

%Gl1 = -tf12/tf11;
Gl1 = -tf(num(2),den(2,:))/tf(num(1),den(1,:))*tf(1, 'IODelay',2);

%Gl2 = -tf21/tf22;
Gl2 = -tf(num(3),den(3,:))/tf(num(4),den(4,:))*tf(1, 'IODelay',4);


%Sintonizando
[Gm,~,Wcg,~] = margin(tf11 + tf12*Gl2);
Kc11_D = Gm/2.2;
ti11_D = (2*pi()/Wcg)/1.2;

[Gm,~,Wcg,~] = margin(tf22 + tf21*Gl1);
Kc22_D =  -Gm/2.2;
ti22_D = (2*pi()/Wcg)/1.2;

%Desacoplador dinâmico ideal
%D = Gp\[tf11 0;0 tf22];


%Salvando figuras
%set(0,'ShowHiddenHandles','On')
%set(gcf,'Units','centimeters','PaperUnits','centimeters')
%pos = get(gcf,'Position');
%set(gcf,'PaperPosition',[0 0 pos(3) pos(4)],'Papersize',[ pos(3),pos(4) ]);
%set(gcf,'InvertHardcopy','off','Renderer','painters')
%xlabel('Tempo(s)')
%ylabel('Variável controlada (%ST)')
%set(gcf,'Renderer','zbuffer')
%print(gcf,'Wood-berry_MA_DesacopladorE_2fechado_y2.jpg','-djpeg','-r300')

function [Kc,ti]=ZN(tau,theta,K)
%FOPDT
% Sintonização em malha aberta fechada para um PI
Kcu = (1 + 2*tau/theta)/K;
wu = 2/theta*sqrt(1 +(theta/tau)); 
%wu = sqrt(2/(theta*tau)*(1+Kcu*K));
Kc = Kcu/2.2;
ti = (2*pi()/wu)/1.2;

%Sintonização em malha aberta
%Kc = (tau/theta)*(0.9/K);
%ti = 3.33*theta;
end