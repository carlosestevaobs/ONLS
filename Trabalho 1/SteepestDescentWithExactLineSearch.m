%% Descida mais Íngreme com Busca Exata
% Otimização Não Linear de Sistemas
% Docente: Daniel Alencar
% Discentes: Carlos Estevão e Ramiro

%% Limpar dados 
    clc
    clear
    format long

%% Inicialização das variáveis e funções
    % Variáveis simbólicas
    syms x1 x2;
    
    % Função
    fx = @(x1, x2) exp(x1 + 3*x2 - 0.1) + exp(x1 - 3*x2 - 0.1)+ exp(-x1 - 0.1);
    
    % Chute inicial    
    x0 = randn(1, 2);

    % Tolerância
    n = 10^(-2);
    
    % Máximo de iterações
    MaxIte = 1000;    
       
    % Verifica em que iteração convergiu
    contaIter = 1;
   
    % x recebe o chute inicial
    x = x0;
    
    % Acumula os passos
    passos = [];
    
    %% Cálculo do Gradiente 
    grad = [diff(fx,x1) diff(fx,x2)];
    
    % Definição da direção
    direcao =  -vpa(subs(grad,[x1,x2],x));

%% Início do algoritmo
    while ((norm(direcao) >= n) & (contaIter <= MaxIte))
        %% Passo 1 - Calcula a direção da descida mais íngreme
        direcao = vpa(subs(grad,[x1,x2],x));        
        [~,i] = max(abs(direcao));
        Delta_x = zeros(size(x));
        Delta_x(i) = -direcao(i);               
        
        %% Passo 2 - Busca Linha Exata
        t = Exact(fx, x, direcao);
        passos = [passos;x];
        
        %% Passo 3 - Atualização ->  X = X + T*(DeltaXsd)  
        x = x + t*Delta_x;
        contaIter = contaIter + 1;
    end
    
    %% Dados e plotagem
    Iteracoes = 1:contaIter-1;
    Tabela = table(Iteracoes', passos);
    disp(Tabela);
    disp('Valor ótimo:');
    disp(x);   
    Plotar(fx, passos);    
    
    