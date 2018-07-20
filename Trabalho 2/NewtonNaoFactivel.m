%% Otimização Não Linear de Sistemas
% Trabalho 2 - Otimização com restrições de igualdade
% Docente: José Daniel de Alencar Santos
% Discentes: Carlos Estevão Bastos Sousa

%% Parâmetros iniciais
clear;
clc;
close all;
tic 

load A;
n = size(A, 2);                 % 100 variáveis
p = size(A, 1);                 % 30 equações

x0 = ones(n,1);                 % entradas não factíveis

iteracoes = 100;
alfa = 0.01;
beta = 0.5;
e = 1e-6;
x = x0;
nu = zeros(p,1);
p_otimo = -34.3473;


if (min(x0) <= 0) || (norm(A * x0 - b) > e)
    disp('Entrada não factível');
    
    for iter = 1:iteracoes
        
        %% Vetor do resíduo
        % r(x, v) = (rdual(x,v), rpri(x,v))
        r = [1+log(x)+A'*nu; A*x-b];      
        
        %% Resolução da equação ax = b baseado na matriz KKT        
        grad = 1 + log(x);
        
        hess = diag(1./x);
        
        KKT_A = (-[(hess) A'; A zeros(p,p)]);
        drw =  KKT_A \ r;        

        dx = drw(1 : n);        
        Dnu = drw(n + [1:p]);
        
        %% Verifica o critério de parada
        if (norm(r) < e)
            break; 
        end
        %% Backtracking line search
        t = 1;
        while (min(x + t * dx) <= 0)
            t = beta * t;
        end      
        while norm([1 + log(x + t * dx) + A' *(nu + Dnu); A *(x + dx) - b]) > (1 - alfa * t) * norm(r)
            t = beta * t; 
        end
        %% Atualizar x e v
        ValoresX(:,iter) = x;
        x = x + t * dx; 
        nu = nu + t * Dnu;
        
        %% Cálculo do erro
        f = x' * log(x);        
        erroTolerancia(iter) = f - e;
        erroP_otimo(iter) = f - p_otimo;
        res_dual(iter) = norm(r(1:n));
        res_pri(iter) = norm(r(n+[1:p]));
    end
    figure(1)
    loglog(erroP_otimo,'bo-');
    hold on;
    title("Erro com base no ponto ótimo dado pela biblioteca CVX");
    xlabel('Iteração');
    ylabel('Erro');
    hold off; 
    
    figure(2);
    loglog(erroTolerancia,'bo-');
    hold on;
    title("Erro com base na tolerância (e)");
    xlabel('Iteração');
    ylabel('Erro');
    hold off;
    
    figure(3)
    iter = [1:6];
    loglog(iter,res_dual,'bo-', iter,res_pri, 'ro-');
    legend('Resíduo dual','Resíduo primal');
    
else
    disp('Entrada factível');
end
toc

res_dual = res_dual'
res_pri = res_pri'