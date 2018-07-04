%% Gradiente Descendente com Busca Backtracking
% Otimização Não Linear de Sistemas
% Docente: Daniel Alencar
% Discentes: Carlos Estevão e Ramiro

%% Limpar dados 
    clc
    clear
    format long

%% Inicialização das variáveis e funções
    % Variáveis simbólicas
    syms x1 x2 s;
    fx = @(x1, x2) exp(x1 + 3*x2 - 0.1) + exp(x1 - 3*x2 - 0.1)+ exp(-x1 - 0.1);
    
    
    % Chute inicial    
    x0 = randn(1, 2);

    % Tolerância
    n = 10^(-1);
    
    % Máximo de iterações
    MaxIte = 1000;    
    
    % Verifica em que iteração convergiu
    contaIter = 1;
    
    % x recebe o chute inicial
    x = x0;
        
    %% Cálculo do Gradiente e da Hessiana
    grad = [diff(fx,x1) diff(fx,x2)];
    hess = [diff(grad(1),x1)   diff(grad(1),x2)
             diff(grad(2),x1)   diff(grad(2),x2)];
    
    VGrad = inf;
    DeltaX = [];
    
    
    % Tamanho do Passo
    t = 1;
        
    while ((norm(VGrad) > n) & (contaIter <= MaxIte))
       
        %% Passo 1 - Computa o passo de Newton e decremento
        % Cálculo do gradiente e hessiana
        VGrad = vpa(subs(grad,[x1,x2],x))
        VHess = vpa(subs(hess,[x1,x2],x));        
        
        
        %% Passo 2 - Critério de parada
        DeltaX = (-pinv(VHess) * VGrad')';
        lambda = VGrad * pinv(VHess) * VGrad';
                
        if (lambda/2) <= n
            break;
        end
        
        %% Passo 3 - Busca BackTracking
        t = Backtrack(fx, x, DeltaX, VGrad); 
        passos(contaIter,:) = x;
        
        %% Passo 4 - Atualização ->  X = X + T*(DeltaXnt)
        x = x +  t * DeltaX;       
        contaIter = contaIter + 1;
        
    end
    
    %% Dados e plotagem
    Iteracoes = 1:contaIter-1;
    Tabela = table(Iteracoes', passos);
    disp(Tabela);
    disp('Valor ótimo:');
    disp(x);   
    Plotar(fx, passos);
    