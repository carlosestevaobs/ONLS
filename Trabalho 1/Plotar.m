%% Gradiente Descendente com Busca Backtracking
% Otimização Não Linear de Sistemas
% Docente: Daniel Alencar
% Discentes: Carlos Estevão e Ramiro

function Plotar(f, passos)    
    fcontour(f);
    hold on
    plot(passos(:,1),passos(:,2), '*-r')   
    hold off
end