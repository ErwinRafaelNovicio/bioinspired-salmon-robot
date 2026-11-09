function animazione_quadrati(x_centri, y_temp, xgrid, xgi, L, Tmax, delay, nome_gif)
% ANIMAZIONE_QUADRATI
% Crea un'animazione con quadrati orientati secondo atan2(y,x)
% e (opzionalmente) salva una GIF animata.
%
% INPUT:
%   x_centri : vettore (1xN) con coordinate X dei centri (fisse)
%   y_temp   : matrice (N x Tmax) con coordinate Y nel tempo
%   xgrid, xgi : vettori di riferimento per il modello
%   L        : lato dei quadrati
%   Tmax     : numero di frame
%   delay    : tempo tra frame (s)
%   nome_gif : (facoltativo) nome file per salvare la GIF (es. 'animazione.gif')
%
% Esempio:
%   animazione_quadrati(x, y_temp, xgrid, xgi, 0.6, Tmax, 0.05, 'test.gif')

    N = numel(x_centri);
    assert(size(y_temp,1) == N, 'y_temp deve avere N righe (una per ogni centro)');
    assert(size(y_temp,2) == Tmax, 'y_temp deve avere Tmax colonne (una per ogni frame)');

    % --- Pre-elaborazione ---
    L_mod = diff(xgrid);
    Id_centri_mod = round(diff(xgi)/2) + xgi(1:end-1);
    x_centri_mod = x_centri(Id_centri_mod);
    y_centri_mod = y_temp(Id_centri_mod,:);
    y_mod = y_temp(xgi,:);

    % --- Errore ---
    err = y_temp - interp1(xgrid, y_mod, x_centri);
    err_in_x = sum(err,2);
    error = norm(err,2);
    error_max = norm(err,'inf');
    pos_max_error_in_x = x_centri(abs(err_in_x)==max(abs(err_in_x)));

    fprintf('Errore totale: %.4g\n', error);
    fprintf('Errore massimo: %.4g\n', error_max);
    fprintf('Posizione max errore: %.3f\n', pos_max_error_in_x);

    % --- Figure e subplot ---
    % --- Figure e subplot ---
figure;
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); % tutto schermo

    x_fun = linspace(min(x_centri)-2, max(x_centri)+2, 400);

    %% --- SUBPLOT 1: salmon data ---
    subplot(3,1,1);
    hold on; grid on;
    xlabel('X'); ylabel('Y'); title('Salmon movement');
    quadrati_salmon = gobjects(1, N);
    for i = 1:N
        [Xq, Yq] = square_oriented(x_centri(i), y_temp(i,1), L);
        quadrati_salmon(i) = plot(Xq, Yq, 'r-', 'LineWidth', 2);
    end
    axis([min(x_fun) max(x_fun) -20 20]);

    %% --- SUBPLOT 2: model fit ---
    subplot(3,1,2);
    hold on; grid on;
    xlabel('X'); ylabel('Y'); title('Model fit');
    xline(xgrid,'--');
    quadrati_model = gobjects(1, numel(x_centri_mod));
    for i = 1:numel(x_centri_mod)
        [Xqm, Yqm] = square_oriented(x_centri_mod(i), y_centri_mod(i,1), L_mod(i));
        quadrati_model(i) = plot(Xqm, Yqm, 'r-', 'LineWidth', 2);
    end
    axis([min(x_fun) max(x_fun) -20 20]);

    %% --- SUBPLOT 3: trajectories ---
    subplot(3,1,3);
    hold on; grid on;
    xlabel('X'); ylabel('Y');
    title('Comparison of trajectories');
    xline(xgrid,'--');
    linea_salmon = plot(nan, nan, 'b-', 'LineWidth', 1.5);
    linea_model = plot(nan, nan, 'k--', 'LineWidth', 1.5);
    hinge_model = plot(nan, nan, 'o', 'LineWidth', 1.5);
    axis([min(x_fun) max(x_fun) -20 20]);
legend('salmon body','robot','hinges',Location='southwestoutside')
    %% --- Animazione ---
    for k = 1:Tmax
        % Subplot 1
        subplot(3,1,1);
        for i = 1:N
            [Xq, Yq] = square_oriented(x_centri(i), y_temp(i,k), L);
            set(quadrati_salmon(i), 'XData', Xq, 'YData', Yq);
        end

        % Subplot 2
        subplot(3,1,2);
        for i = 1:numel(x_centri_mod)
            [Xqm, Yqm] = square_oriented(x_centri_mod(i), y_centri_mod(i,k), L_mod(i));
            set(quadrati_model(i), 'XData', Xqm, 'YData', Yqm);
        end

        % Subplot 3
        subplot(3,1,3);
        set(linea_salmon, 'XData', x_centri, 'YData', y_temp(:,k));
        set(linea_model, 'XData', xgrid, 'YData', y_mod(:,k));
        set(hinge_model, 'XData', xgrid, 'YData', y_mod(:,k));

        drawnow;

        % --- Salva frame in GIF (se richiesto) ---
        if nargin >= 8 && ~isempty(nome_gif)
           
            frame = getframe(gcf);
            im = frame2im(frame);
            [A,map] = rgb2ind(im,256);
            if k == 1
                imwrite(A,map,nome_gif,'gif','LoopCount',inf,'DelayTime',delay);
            else
                imwrite(A,map,nome_gif,'gif','WriteMode','append','DelayTime',delay);
            end
        else
            pause(delay);
        end
    end
end

% --- Funzione ausiliaria per i quadrati ---
function [Xq, Yq] = square_oriented(x, y, L)
    theta = atan2(y, x);
    halfL = L/2;
    halfLy=halfL/2;
    X0 = [-halfL,  halfL,  halfL, -halfL, -halfL];
    Y0 = [-halfLy, -halfLy,  halfLy,  halfLy, -halfLy];
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    P = R * [X0; Y0];
    Xq = P(1,:) + x;
    Yq = P(2,:) + y;
end
