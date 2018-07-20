%% Definição do Alvo com o softwre CVX
% Trabalho 2 - Otimização com restrições de igualdade
% Docente: José Daniel de Alencar Santos
% Discentes: Carlos Estevão Bastos Sousa

clear;
clc;
close all;
tic

load A
load x0;
n = size(A, 2);                 % 100 variáveis
p = size(A, 1);                 % 30 equações
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