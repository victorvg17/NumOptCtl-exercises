function [ dx ] = ballistic_dynamics( x ) 

w = 2;
d = .1;
g = 9.81;

v = x(3:4);
dnorm = norm( v - [w;0]);

dx  = [ v(1);
        v(2);
        -(v(1) - w) * dnorm  * d ;
        -v(2)  * dnorm * d - g];
end