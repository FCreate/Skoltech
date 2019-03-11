x = rand(2,10);
cvx_begin
    variable A(2,2)
    variable b(2)
    maximize(log_det( A ))
    subject to
        norms( A * x + b * ones( 1, 10 ), 2 ) <= 1;
cvx_end
points  = linspace( 0, 2 * pi, 100);
figure  = A \ [ cos(points) - b(1) ; sin(points) - b(2) ];
plot(figure(1,:), figure(2,:), '--', x(1,:), x(2,:), 'r*');