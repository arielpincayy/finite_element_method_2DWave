function [M, K, a_hat] = matrices(triangles, points, n, m)

    m_st = [2, 1, 1;
            1, 2, 1;
            1, 1, 2];

    a_hat = [0,0;
             1,0;
             0,1;
             0,1/2;
             1/2,1/2;
             1/2,0;
             1/3,1/3];


    % Inicializar matrices globales
    num_nodes = size(points, 1); % Total de nodos
    M = zeros(num_nodes, num_nodes);
    K = zeros(num_nodes, num_nodes);

    % Ensamblar matrices globales
    for tk = 1:size(triangles, 1)
        triangle = triangles(tk, :); % Triángulo actual
        node_indices = triangle(1:3); % Índices de los nodos del triángulo
        area = triangle(4); % Área del triángulo

        % Matriz de masa local
        for i = 1:3
            for j = 1:3
                i_g = node_indices(i); % Índice global del nodo i
                j_g = node_indices(j); % Índice global del nodo j
                M(i_g, j_g) = M(i_g, j_g) + m_st(i, j) * (area / 12);
            end
        end

        % Coordenadas de los nodos del triángulo
        p1 = points(node_indices(1), :);
        p2 = points(node_indices(2), :);
        p3 = points(node_indices(3), :);

        x1 = p1(1); y1 = p1(2);
        x2 = p2(1); y2 = p2(2);
        x3 = p3(1); y3 = p3(2);

        % Coeficientes b_i y c_i
        b1 = (y2 - y3);
        b2 = (y3 - y1);
        b3 = (y1 - y2);

        c1 = (x3 - x2);
        c2 = (x1 - x3);
        c3 = (x2 - x1);

        % Matriz de rigidez local
        K_k = [
            ((b1 * b1) + (c1 * c1)), ((b1 * b2) + (c1 * c2)), ((b1 * b3) + (c1 * c3));
            ((b2 * b1) + (c2 * c1)), ((b2 * b2) + (c2 * c2)), ((b2 * b3) + (c2 * c3));
            ((b3 * b1) + (c3 * c1)), ((b3 * b2) + (c3 * c2)), ((b3 * b3) + (c3 * c3))
        ];

        % Ensamblar en la matriz global de rigidez
        for i = 1:3
            for j = 1:3
                i_g = node_indices(i); % Índice global del nodo i
                j_g = node_indices(j); % Índice global del nodo j
                K(i_g, j_g) = K(i_g, j_g) + (K_k(i, j) * (1 / (4 * area)));
            end
        end

    end

    num_boundary = 2*m + 2*(n-2);
    % Número de nodos internos
    num_internal = size(points, 1) - num_boundary;

    % Reducción de M y K
    M_reduced = M(1:num_internal, 1:num_internal);
    K_reduced = K(1:num_internal, 1:num_internal);

    M = M_reduced;
    K = K_reduced;
end