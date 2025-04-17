function val = f_source_test(x, y, t, x0, y0, t0, r, w0)
    % x, y: Coordenadas espaciales
    % t: Tiempo
    % x0, y0: Centro espacial
    % t0: Centro temporal
    % r: Radio de influencia
    % w0: Amplitud y escala temporal

    % Calcular gamma
    gamma = (x - x0)^2 + (y - y0)^2;

    % Evaluar la función w(t)
    wt = -8 * w0 * (t - t0) * exp(-((t - t0)^2) / (16 * w0^2));

    % Evaluar s(x, y, t)
    if gamma <= r^2
        val = (cos((pi / 2) * gamma / r^2)^2) * wt;
    else
        val = 0;
    end
end
