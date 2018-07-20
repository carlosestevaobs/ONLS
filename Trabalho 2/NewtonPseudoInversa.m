%% Otimiza��o N�o Linear de Sistemas
% Trabalho 2 - Otimiza��o com restri��es de igualdade
% Docente: Jos� Daniel de Alencar Santos
% Discentes: Carlos Estev�o Bastos Sousa

%% Par�metros iniciais
clear;
clc;
close all;
tic

load A
n = size(A, 2);                 % 100 vari�veis
p = size(A, 1);                 % 30 equa��es

load x0
b = A * x0;                     % Ax = b

iteracoes = 100;                % N�mero m�ximo de itera��es
alfa = 0.01;                    % Alfa e beta para o backtracking
beta = 0.5;
e = 1e-6;                       % Toler�ncia
x = x0;                         % Ponto inicial

p_otimo = -34.3473;
lambda = [];

if (min(x0) <= 0) || (norm(A * x0 - b) > e)
    disp('Entrada n�o fact�vel');
else
    disp('Entradas fact�veis');
    for iter = 1:iteracoes
        %% C�lculo do valor da fun��o
        valFun = x'*log(x);
        %% Gradiente da fun��o
        grad = 1+log(x);
        %% Hessiana da fun��o
        hess = diag(x.^-2);        
        %% Resolu��o da equa��o Ax = b (matriz KKT) com pseudo-inversa
        KKT_A = [hess A'; A zeros(p,p)];
        KKT_B = [-grad; zeros(p,1)];
        dxw = pinv(KKT_A) \ KKT_B;
        %% C�lculo da dire��o
        dx = dxw(1:n);
        %% Verifica o crit�rio de parada
        % lambda = (-grad'*dx)^0.5; da forma abaixo ficou mais r�pido  encontrar o m�nimo
        lambda = grad'*dx;
        if ((lambda^2)/2 < e)
            break;
        end
        %% Backtracking line search
        t=1;
        while (min(x + t*dx) <= 0)
            t = beta * t;
        end
        while ((x + t*dx)'*log(x + t*dx) >= valFun + t * alfa * lambda)
            t = beta * t;
        end
         %% Atualizar X
         ValoresX(:,iter) = x;
        x = x + t * dx;
        %% C�lculo do erro
        
        f = x' * log(x);
        AcF(iter,:) =  f;
        valorF(iter,:) = f; 
        erroTolerancia(iter) = f - e;
        erroP_otimo(iter) = f - p_otimo;
    end
    %% Plotagem do erro
    figure(1);
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
end
toc