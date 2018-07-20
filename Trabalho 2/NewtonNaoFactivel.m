%% Otimiza��o N�o Linear de Sistemas
% Trabalho 2 - Otimiza��o com restri��es de igualdade
% Docente: Jos� Daniel de Alencar Santos
% Discentes: Carlos Estev�o Bastos Sousa

%% Par�metros iniciais
clear;
clc;
close all;
tic 

load A;
n = size(A, 2);                 % 100 vari�veis
p = size(A, 1);                 % 30 equa��es

x0 = ones(n,1);                 % entradas n�o fact�veis

iteracoes = 100;
alfa = 0.01;
beta = 0.5;
e = 1e-6;
x = x0;
nu = zeros(p,1);
p_otimo = -34.3473;


if (min(x0) <= 0) || (norm(A * x0 - b) > e)
    disp('Entrada n�o fact�vel');
    
    for iter = 1:iteracoes
        
        %% Vetor do res�duo
        % r(x, v) = (rdual(x,v), rpri(x,v))
        r = [1+log(x)+A'*nu; A*x-b];      
        
        %% Resolu��o da equa��o ax = b baseado na matriz KKT        
        grad = 1 + log(x);
        
        hess = diag(1./x);
        
        KKT_A = (-[(hess) A'; A zeros(p,p)]);
        drw =  KKT_A \ r;        

        dx = drw(1 : n);        
        Dnu = drw(n + [1:p]);
        
        %% Verifica o crit�rio de parada
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
        
        %% C�lculo do erro
        f = x' * log(x);        
        erroTolerancia(iter) = f - e;
        erroP_otimo(iter) = f - p_otimo;
        res_dual(iter) = norm(r(1:n));
        res_pri(iter) = norm(r(n+[1:p]));
    end
    figure(1)
    loglog(erroP_otimo,'bo-');
    hold on;
    title("Erro com base no ponto �timo dado pela biblioteca CVX");
    xlabel('Itera��o');
    ylabel('Erro');
    hold off; 
    
    figure(2);
    loglog(erroTolerancia,'bo-');
    hold on;
    title("Erro com base na toler�ncia (e)");
    xlabel('Itera��o');
    ylabel('Erro');
    hold off;
    
    figure(3)
    iter = [1:6];
    loglog(iter,res_dual,'bo-', iter,res_pri, 'ro-');
    legend('Res�duo dual','Res�duo primal');
    
else
    disp('Entrada fact�vel');
end
toc

res_dual = res_dual'
res_pri = res_pri'