xx = [0,1];
yy = xx;
zz = xx;

basis = zeros(8,8);
for i_point = 1:8
    [i_x,i_y,i_z] = ind2sub([2 2 2],i_point);
    x = xx(i_x); y = yy(i_y); z = zz(i_z);
    for i_term = 1:8
        [x_power,y_power,z_power] = ind2sub([2 2 2],i_term);
        basis(i_point,i_term) = x^(x_power-1)*y^(y_power-1)*z^(z_power-1);
    end
end
A = inv(basis);


func = @(x,y,z) log(1+2*x^2+3*y^2+4*z^2);
for i_point = 1:8
    [i_x,i_y,i_z] = ind2sub([2 2 2],i_point);
    x = xx(i_x); y = yy(i_y); z = zz(i_z);
    nodesval(i_point) = func(x,y,z);
end

%% Chebyshev

xx = ChebyshevRoots( 5, 'Tn', [-1 1] );
yy = ChebyshevRoots( 5, 'Tn', [-1 1] );
zz = ChebyshevRoots( 5, 'Tn', [-1 1] );
degree = 3;
for i_point = 1:5^3
    [i_x,i_y,i_z] = ind2sub([5 5 5],i_point);
    x = xx(i_x); y = yy(i_y); z = zz(i_z);
    [basis(i_point,:),~] = ChebyshevND(degree,[x y z]);
end
A = (basis'*basis)\(basis')

