function cmap = blue_white_red(n)
    if nargin < 1
        n = 256;  % numero di colori
    end

    % Numero di colori per ciascuna metà
    n_half = floor(n / 2);

    % Colori da blu scuro ad azzurro → bianco
    r1 = linspace(0, 1, n_half)';
    g1 = linspace(0, 1, n_half)';
    b1 = ones(n_half, 1);

    % Colori da bianco → rosa → rosso scuro
    r2 = ones(n - n_half, 1);
    g2 = linspace(1, 0, n - n_half)';
    b2 = linspace(1, 0, n - n_half)';

    % Metti insieme
    cmap = [r1 g1 b1; r2 g2 b2];
end