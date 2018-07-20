    % Variáveis simbólicas
    syms x1 x2 x3;
    
    % Função
    fx = @(x1, x2) x1 * log(x1) + x2 * log(x2) ;
    figure(1)
    ezcontourf(fx)
    title('Exemplo da função com N = 2')
    
    figure(2)
    ezsurfc(fx)
    title('Exemplo da função com N = 2')
    