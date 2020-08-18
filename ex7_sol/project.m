function index = project( in, values )
% takes scalar value in and projects in on grid values
% return the corresponding index of values (which might lie outside of
% values if in is outside values
% ASSUMPTION:   the argument "values" contains a vector of equally spaced
%               points
    min = values(1);
    max = values(end);
    N = length(values);
    index = round(((in-min)/(max-min))*(N-1)) + 1;
end

