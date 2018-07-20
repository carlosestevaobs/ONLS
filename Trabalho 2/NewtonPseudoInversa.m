%% Otimização Não Linear de Sistemas
% Trabalho 2 - Otimização com restrições de igualdade
% Docente: José Daniel de Alencar Santos
% Discentes: Carlos Estevão Bastos Sousa

%% Parâmetros iniciais
clear;
clc;
close all;
tic

load A
n = size(A, 2);                 % 100 variáveis
p = size(A, 1);                 % 30 equações

load x0
b = A * x0;                     % Ax = b

iteracoes = 100;                % Número máximo de iterações
alfa = 0.01;                    % Alfa e beta para o backtracking
beta = 0.5;
e = 1e-6;                       % Tolerância
x = x0;                         % Ponto inicial

p_otimo = -34.3473;
lambda = [];

if (min(x0) <= 0) || (norm(A * x0 - b) > e)
    disp('Entrada não factível');
else
    disp('Entradas factíveis');
    for iter = 1:iteracoes
        %% Cálculo do valor da função
        valFun = x'*log(x);
        %% Gradiente da função
        grad = 1+log(x);
        %% Hessiana da função
        hess = diag(x.^-2);        
        %% Resolução da equação Ax = b (matriz KKT) com pseudo-inversa
        KKT_A = [hess A'; A zeros(p,p)];
        KKT_B = [-grad; zeros(p,1)];
        dxw = pinv(KKT_A) \ KKT_B;
        %% Cálculo da direção
        dx = dxw(1:n);
        %% Verifica o critério de parada
        % lambda = (-grad'*dx)^0.5; da forma abaixo ficou mais rápido  encontrar o mínimo
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
        %% Cálculo do erro
        
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
end
toc