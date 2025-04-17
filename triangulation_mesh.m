function [triangles, points, neighs, adjs, stencil] = triangulation_mesh(lx, ly, top, bottom, left, right, plot)

    n = size(top,2);
    m = size(left,1) + 2;

    stencil = zeros(n,m);

    pos = 1;
    for i=2:n-1
        for j=2:m-1
            stencil(i,j) = pos;
            pos = pos + 1;
        end
    end

    for i=1:m
        stencil(1,i) = pos;
        stencil(end,i) = m + pos;
        pos = pos + 1; 
    end

    pos = m + pos;

    for i=2:n-1
        stencil(i,1) = pos;
        stencil(i,end) = (n-2) + pos;
        pos = pos + 1;
    end

    % Creación de los puntos usando meshgrid
    px = linspace(0, lx, n);
    py = linspace(0, ly, m);

    points = zeros(n*m,2);

    for i = 1:n
        for j = 1:m
            pos = stencil(i, j);
            points(pos, :) = [px(i), py(j)];
        end
    end 

    % Triangulación
    triangles = delaunay(points(:,1), points(:,2));
    num_triangles = size(triangles, 1);
    triangles_matrix = zeros(num_triangles, 4);
    
    % Inicialización de las listas de adyacencia
    adjs = cell(size(points, 1), 1);
    for t = 1:num_triangles
        tri_indices = triangles(t, :);
        adjs{tri_indices(1)} = [adjs{tri_indices(1)}, t];
        adjs{tri_indices(2)} = [adjs{tri_indices(2)}, t];
        adjs{tri_indices(3)} = [adjs{tri_indices(3)}, t];

        % Calcular área usando fórmula determinante
        tri_points = points(tri_indices, :);
        x = tri_points(:, 1);
        y = tri_points(:, 2);
        %area = abs(x(1)*(y(2)-y(3)) + x(2)*(y(3)-y(1)) + x(3)*(y(1)-y(2))) / 2;
        area = 0.5 * det([x(1), y(1), 1; x(2), y(2), 1; x(3), y(3), 1]);
        triangles_matrix(t, :) = [tri_indices, area];
    end
    
    % Cálculo de vecinos usando un mapa de bordes
    edge_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    neighs = zeros(num_triangles, 3);
    for t = 1:num_triangles
        tri_indices = triangles(t, :);
        edges = [tri_indices([1, 2]); tri_indices([2, 3]); tri_indices([3, 1])];
        for e = 1:3
            edge = sort(edges(e, :));
            key = sprintf('%d-%d', edge(1), edge(2));
            if isKey(edge_map, key)
                edge_map(key) = [edge_map(key), t];
            else
                edge_map(key) = t;
            end
        end
    end

    % Asignar vecinos
    for t = 1:num_triangles
        tri_indices = triangles(t, :);
        edges = [tri_indices([1, 2]); tri_indices([2, 3]); tri_indices([3, 1])];
        neighbor_count = 0;
        for e = 1:3
            edge = sort(edges(e, :));
            key = sprintf('%d-%d', edge(1), edge(2));
            neighbors = edge_map(key);
            neighbors(neighbors == t) = []; % Excluir triángulo actual
            if ~isempty(neighbors)
                neighbor_count = neighbor_count + 1;
                neighs(t, neighbor_count) = neighbors(1);
            end
        end
    end

    % Grafica de la malla y los números de los puntos
    if(plot)
        figure;
        triplot(triangles, points(:, 1), points(:, 2), 'k');
        hold on;
        for i = 1:size(points, 1)
            text(points(i, 1), points(i, 2), num2str(i), 'Color', 'red', 'FontSize', 18, 'HorizontalAlignment', 'center');
        end
        hold off;
        title('Triangulación de la Malla con Puntos Numerados');
        xlabel('x');
        ylabel('y');
        axis equal;
    end

    % Salida
    triangles = triangles_matrix;
end
