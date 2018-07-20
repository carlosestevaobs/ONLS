%% Defini��o do Alvo com o softwre CVX
% Trabalho 2 - Otimiza��o com restri��es de igualdade
% Docente: Jos� Daniel de Alencar Santos
% Discentes: Carlos Estev�o Bastos Sousa

clear;
clc;
close all;
tic

load A
load x0;
n = size(A, 2);                 % 100 vari�veis
p = size(A, 1);                 % 30 equa��es
b =  A * x0;     

%% Resolvendo o problema com a biblioteca CVX
cvx_begin
    variable x(n);
    minimize (sum(entr( x )) * - 1)
    subject to
        A * x == b;  
cvx_end
Alvo = cvx_optval

toc