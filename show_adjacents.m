function adjacent_triangles = show_adjacents(point_id, n, m, lx, ly)
    % Llama a la función `triangulation_mesh` para generar la malla
    [triangles, points, ~, adjs] = triangulation_mesh(n, m, lx, ly, 0);

    % Verifica si el punto solicitado es válido
    num_points = size(points, 1);
    if point_id < 1 || point_id > num_points
        error('El ID del punto debe estar entre 1 y %d.', num_points);
    end

    % Obtén los triángulos adyacentes al punto dado
    adjacent_triangles = adjs{point_id};
    
    % Verifica si hay triángulos adyacentes
    if isempty(adjacent_triangles)
        warning('El punto %d no tiene triángulos adyacentes.', point_id);
        return;
    end

    % Graficar la malla completa como fondo
    figure;
    triplot(triangles(:, 1:3), points(:, 1), points(:, 2), 'k'); % Malla completa en negro
    hold on;
    axis equal;

    % Resalta los triángulos adyacentes
    for i = 1:length(adjacent_triangles)
        % Índices de los vértices del triángulo
        tri_indices = triangles(adjacent_triangles(i), 1:3);
        tri_points = points(tri_indices, :);
        centroid = mean(tri_points, 1);

        % Rellenar triángulo
        fill(tri_points(:, 1), tri_points(:, 2), 'cyan', 'FaceAlpha', 0.5, 'EdgeColor', 'blue');

        % Etiqueta del triángulo
        text(centroid(1), centroid(2), sprintf('%d', adjacent_triangles(i)), ...
             'Color', 'black', 'FontSize', 12, 'HorizontalAlignment', 'center');
    end

    % Resalta el punto especificado
    plot(points(point_id, 1), points(point_id, 2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);

    % Etiqueta el punto
    text(points(point_id, 1), points(point_id, 2), sprintf('P%d', point_id), ...
         'Color', 'red', 'FontSize', 12, 'HorizontalAlignment', 'left');

    % Configuración de visualización
    title(sprintf('Punto %d y sus triángulos adyacentes', point_id));
    hold off;

    % Imprime los resultados
    fprintf('El punto %d es adyacente a los siguientes triángulos:\n', point_id);
    disp(adjacent_triangles);
end
