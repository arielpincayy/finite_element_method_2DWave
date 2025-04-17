function [u_aprox] = fem_edp(lx, ly, n, m, c, delta_t, f, f_source, t_end)

    %lx -> tamano x
    %ly -> tamano y
    %n -> n puntos en y
    %m -> m puntos en x
    %c -> contante velocidad
    %delta_t -> espacio entre cada instante de tiempo
    %f -> funcion solucion
    %f_source -> funcion fuente
    %t_end -> tiempo final

    top = zeros(1,n);
    bottom = zeros(1,n);
    left = zeros(m-2,1);
    right = zeros(m-2,1);

    num_steps = floor(t_end/delta_t);

    T = linspace(delta_t,t_end,num_steps);
    errors = zeros(1,num_steps);

    z_min = -1;
    z_max = 1;
    
    [X, Y] = meshgrid(linspace(0, lx, n), linspace(0, ly, m));


    [triangles, points, neighs, adjs, stencil] = triangulation_mesh(lx, ly, top, bottom, left, right, 0);

    %CALCULO DE U^N-1 -----------------------------------------------------------------------------------------

    u_old = zeros(n,m);
    
    for i=1:n
        for j=1:m
            point = stencil(i,j);
            u_old(i,j) = f(points(point, 1), points(point, 2), 0);
        end
    end

    u_old(1, :) = top;            
    u_old(end, :) = bottom;       
    u_old(2:end-1, 1) = left;     
    u_old(2:end-1, end) = right;  

    %CALCULO DE U^N-----------------------------------------------------------------------------

    df_u = zeros(n,m);
 
    for i=1:n
        for j=1:m
            point = stencil(i,j);
            df_u(i, j) = (f(points(point, 1), points(point, 2), 1e-10) - f(points(point, 1), points(point, 2), 0)) / 1e-10;
        end
    end

    rx = (c * delta_t / (lx/n))^2;
    ry = (c * delta_t / (ly/m))^2;

    u = zeros(n, m);  % u^1 (primer paso temporal)

    % Fórmula especial para el primer paso usando la derivada inicial
    for i = 2:n-1
        for j = 2:m-1
            point = stencil(i,j);
            u(i,j) = u_old(i,j) + delta_t * df_u(i,j) + 0.5 * (rx * (u_old(i+1,j) + u_old(i-1,j)) + ry * (u_old(i,j+1) + u_old(i,j-1)) - 2 * (rx + ry) * u_old(i,j)) + 0.5 * delta_t^2 * f_source(points(point, 1),points(point, 2),T(1));
        end
    end

    u(1, :) = top;            % Borde superior
    u(end, :) = bottom;       % Borde inferior
    u(2:end-1, 1) = left;     % Borde izquierdo
    u(2:end-1, end) = right;  % Borde derecho

    [u_exact] = exact(c,lx,ly,T(1),n,m,f,0);
    errors(1,1) = max(max((abs(u_exact - u))));

    figure;
    surf(X, Y, u, 'LineStyle', 'none', 'FaceColor', 'interp');
    colorbar; % Mostrar la barra de colores
    title(['Wave equation on t = ', num2str(delta_t)]);
    xlabel('x');
    ylabel('y');
    zlabel('u(x, y)');
    view(3); % Vista 3D
    colorbar;
    zlim([z_min, z_max]);
    pause(0.01);

    %ITERACION SOBRE LOS TRIANGULOS//////////////////////////////////////////////////////////////////////////////////////////////////////
    [M, K, a_hat] = matrices(triangles, points, n, m);
    %///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
    vect_u = zeros((n-2)*(m-2),1);
    vect_u_old = zeros((n-2)*(m-2),1);
    
    for i=2:n-1
        for j=2:m-1
            point = stencil(i,j);
            vect_u(point) = u(i,j);
            vect_u_old(point) = u_old(i,j);
        end
    end

    for t=2:num_steps

    % CALCULO DE F------------------------------------------------------------------------------------------------------------------------
        F_cpltd = zeros(n*m);
        for tk = 1:size(triangles, 1)
            triangle = triangles(tk, :); % Triángulo actual
            node_indices = triangle(1:3); % Índices de los nodos del triángulo
            area = triangle(4);

            for i=1:3
                i_g = node_indices(i);

                if i == 1
                    f_b = @(x,y)(1-x-y);
                elseif i == 2
                    f_b = @(x,y)(x);
                elseif i == 3
                    f_b = @(x,y)(y);
                end

                F_cpltd(i_g) = F_cpltd(i_g) + ...
                               (area/60)*(3*(f_source(a_hat(1,1),a_hat(1,2),T(t))*f_b(a_hat(1,1),a_hat(1,2)) ...
                                           + f_source(a_hat(2,1),a_hat(2,2),T(t))*f_b(a_hat(2,1),a_hat(2,2)) ...
                                           + f_source(a_hat(3,1),a_hat(3,2),T(t))*f_b(a_hat(3,1),a_hat(3,2))) ... 
                                        +...
                                          8*(f_source(a_hat(4,1),a_hat(4,2),T(t))*f_b(a_hat(4,1),a_hat(4,2)) ...
                                           + f_source(a_hat(5,1),a_hat(5,2),T(t))*f_b(a_hat(5,1),a_hat(5,2)) ...
                                           + f_source(a_hat(6,1),a_hat(6,2),T(t))*f_b(a_hat(6,1),a_hat(6,2))) ...
                                        +...
                                         27*(f_source(a_hat(7,1),a_hat(7,2),T(t))*f_b(a_hat(7,1),a_hat(7,2))));
            end
        end

        F = F_cpltd(1:(n-2)*(m-2));
        
        %RESOLUCION Y CONTRUCCION DE SISTEMA-------------------------------------------------------------------------------------------------------------------------------------------------------
        b = F  + (2*M*vect_u - M*vect_u_old)/(delta_t^2);
        A = (M/(delta_t^2) + K);
    
        [u_aprox_vect, iter, err] = gaussSeidel(A,b,1e-6,100);
        u_aprox = zeros(n,m);
    
    
        for i=2:n-1
            for j=2:m-1
                point = stencil(i,j);
                u_aprox(i,j) = u_aprox_vect(point);
            end
        end

        [u_exact] = exact(c,lx,ly,T(t),n,m,f,0);
        errors(1,t) = max(max((abs(u_exact - u_aprox))));
    
        u_aprox(1, :) = top;            % Borde superior
        u_aprox(end, :) = bottom;       % Borde inferior
        u_aprox(2:end-1, 1) = left;     % Borde izquierdo
        u_aprox(2:end-1, end) = right;  % Borde derecho
    
        surf(X, Y, u_aprox, 'LineStyle', 'none', 'FaceColor', 'interp');
        colorbar; % Mostrar la barra de colores
        title(['Wave equation on t = ', num2str(t * delta_t)]);
        xlabel('x');
        ylabel('y');
        zlabel('u(x, y)');
        view(3); % Vista 3D
        colorbar;
        zlim([z_min, z_max]);
        pause(0.01);

        vect_u_old = vect_u;
        vect_u = u_aprox_vect;
    end

    [u_exact] = exact(c,ly,lx,T(end),n,m,f,1);
    [error_mesh] = error_graph(u_exact, u_aprox, n, m, ly, lx);

    figure;
    plot(T,errors);
    set(gca, 'FontSize', 16);
    title('Error curve');
    xlabel('Time (seconds)');
    ylabel('Error');
    colorbar;
end
