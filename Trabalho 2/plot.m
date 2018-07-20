    % Vari�veis simb�licas
    syms x1 x2 x3;
    
    % Fun��o
    fx = @(x1, x2) x1 * log(x1) + x2 * log(x2) ;
    figure(1)
    ezcontourf(fx)
    title('Exemplo da fun��o com N = 2')
    
    figure(2)
    ezsurfc(fx)
    title('Exemplo da fun��o com N = 2')
    