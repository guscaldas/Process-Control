%% Modelo de simulação dinâmica para o Reator de Van de Vusse no Simulink

% COQ 792 - Controle de Processos

% Autor: Gustavo Luís Rodrigues Caldas

function [sys,x0]=sfunction_vandevusse(t,x,u,flag)
if flag==0 % Inicializa��o
   
%  [sys]=[neq(estados continuos),0(estados discretos),nsai,nent,0,0]
   [sys]=[3,0,3,4,0,0];

  % condicoes iniciais de integracao
    Ca0=5.1;
    T0=130;
    x0 = [Ca0;0;T0];

elseif flag==1 % Derivadas
   [sys]=model_vandevusse_2(t,x,u);
   
elseif flag==3 % Saídas
     c=[1 0 0
       0 1 0
       0 0 1];
   [sys]=c*x;
else
   [sys]=[];
end