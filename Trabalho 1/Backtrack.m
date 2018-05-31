%% Gradiente Descendente com Busca Backtracking
% Otimização Não Linear de Sistemas
% Docente: Daniel Alencar
% Discentes: Carlos Estevão e Ramiro

function t = Backtrack(fx, x, Direcao, Gradiente)
    syms x1 x2;
    % t := 1
    t=1;
    
    % 0 < ? < 0.5, 0 < ? < 1
    alpha = 0.2; 
    beta = 0.5;
    
    % while (f(x + t?x) > f(x) + ?t?f(x)T ?x)
    while (vpa(subs(fx,[x1,x2],x + t*Direcao)) > (fx(x(:,1),x(:,2)) + alpha*t*Gradiente*Direcao'))
        % t := ?t
        t = beta*t;
    end
end