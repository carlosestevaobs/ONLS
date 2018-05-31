%% Gradiente Descendente com Busca Backtracking
% Otimização Não Linear de Sistemas
% Docente: Daniel Alencar
% Discentes: Carlos Estevão e Ramiro

function t = Exact(fx, x, Direcao)
    s = 0.01:0.01:.99;
    
    % t = argmin s?0 -> f(x + s?x).
    ray = repmat(x,size(s'))+s'*Direcao;
    [~,iDirecao] = min(fx(ray(:,1), ray(:,2)));
    t = s(iDirecao);
end