function [basis,order_table] = ChebyshevND(degree,input)
% Evaluate value of N dimensional complete product basis for any degree at x
% x is n-by-1 or 1-by-n vector between [-1,1].
% Author: Linxi Chen
x_length = length(input);
n = x_length;
d = degree;
if mod(degree,1) ~= 0
    error('Degree not integer');
elseif degree == 0
    basis = 1;
    order_table = zeros(1,x_length);
else
    [basis,order_table] = ChebyshevND(degree-1,input);
    % Beautiful algorithm from Roger Stafford on MATLAB central
    c = nchoosek(1:d+n-1,n-1);
    m = size(c,1);
    t = ones(m,d+n-1);
    t(repmat((1:m).',1,n-1)+(c-1)*m) = 0;
    u = [zeros(1,m);t.';zeros(1,m)];
    v = cumsum(u,1);
    add_table = diff(reshape(v(u==0),n+1,m),1).';
    %==========================================================
    order_table = [order_table;add_table];
    add_terms = size(add_table,1);
    temp_stuff = zeros(1,n);
    for i_term = 1:add_terms
        for j = 1:x_length
            temp_stuff(j) = chebypoly(add_table(i_term,j),input(j));
        end
        basis = [basis,prod(temp_stuff)];
    end
end
