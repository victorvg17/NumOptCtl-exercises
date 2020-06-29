function index = project( a, values )
% takes scalar value a and projects it on the grid defined by values
% returns the corresponding index (which might lie outside values if 
% a lies outside the grid values)

% ASSUMPTION:   the argument "values" contains a vector of ascending 
%               equally spaced points

% Input:
%           a:          scalar value
%           values:     grid of equally spaced values

% Output:
%           index: index corresponding to a on grid values 

    min = values(1);
    max = values(end);
    N = length(values);
    index = round(((a-min)/(max-min))*(N-1)) + 1;
end

