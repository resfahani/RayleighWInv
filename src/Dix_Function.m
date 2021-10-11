function [G] = Inv_Tsai(c,fr,z,v, inv_flag)



z = [z inf];
% v = 0.25;

% Set parameter for inversion
%inv_flag = 5;

% set some factors that are used to define kernel for homogeneous form

% f is the factor that multiplies shear wave velocity to get Rayleigh
% wave velocity
f = sqrt((8/3) - ((-16 + 56*v - 40*(v^2))/...
    (3*(2^(2/3))*(-1 + v)*(-11 + 78*v - 123*(v^2) + 56*(v^3) + ...
    3*sqrt(3)*sqrt(-5 + 36*v - 94*(v^2) + 148*(v^3) - 165*(v^4) + ...
    112*(v^5) - 32*(v^6)))^(1/3))) + ((1/(3*(-1 + v)))*(2^(2/3))*(-11 + ...
    78*v - 123*(v^2) + 56*(v^3) + 3*sqrt(3)*sqrt(-5 + 36*v - 94*(v^2) + ...
    148*(v^3) - 165*(v^4) + 112*(v^5) - 32*(v^6)))^(1/3)));
% make sure there isn't a small imaginary part
f = real(f);
% this is the denominator of the Dix expression in Haney and Tsai (2015), 
% multiplied by 8
g = -((((-2 + f^2)^3)/((1 - f^2)^(3/2))) - ...
    (((2*sqrt(2)*sqrt((-2 + f^2 + 2*v - 2*(f^2)*v)/(-1 + v))*(4 - 4*v + ...
    (f^2)*(-1 + 2*v))))/((2 - 2*v + (f^2)*(-1 + 2*v)))) + ((8*(-2 + ...
    f^2)*(2 - 2*f^2 + sqrt(2 - 2*(f^2))*sqrt((-2 + f^2 + 2*v - ...
    2*(f^2)*v)/(-1 + v))))/((-1 +f^2)*(2*sqrt(1 - f^2) + ...
    sqrt(2)*sqrt((-2 + f^2 + 2*v - 2*(f^2)*v)/(-1 + v))))));

% these are factors and exponential terms used in the homogeneous 
% formulation
fctr1 = -((f^2 - 2)^2)*(8 - 8*(f^2) + f^4)/((1 - f^2)^1.5)/g;
exp1 = 2*sqrt(1 - f^2);
fctr2 = -8*(f^2 - 2)*2*(1+sqrt(1 - f^2)*sqrt(1 + ((f^2)*((2*v - 1)/...
    (2 - 2*v)))))/sqrt(1 - f^2)/g;
exp2 = sqrt(1 - f^2)+sqrt(1 + ((f^2)*((2*v - 1)/(2 - 2*v))));
fctr3 = 4*sqrt(1 + ((f^2)*((2*v - 1)/(2 - 2*v))))*(16*(v-1) + ...
    (f^2)*(f^2 - 8)*(2*v - 1))/(2 - 2*v + (f^2)*(2*v - 1))/g;
exp3 = 2*sqrt(1 + ((f^2)*((2*v - 1)/(2 - 2*v))));
 

w = 2*pi*fr; k = w./c; kz = k'*z;

% Phase velocity, homogeneous form
f_homo = fctr3*exp(-exp3*kz)+fctr2*exp(-exp2*kz)+fctr1*exp(-exp1*kz);
f_homo(:,end) = 0;

% Group velocity, homogeneous form
f_homo_U = fctr3*(1-exp3*kz).*exp(-exp3*kz)+...
           fctr2*(1-exp2*kz).*exp(-exp2*kz)+...
           fctr1*(1-exp1*kz).*exp(-exp1*kz);        
f_homo_U(:,end) = 0;        

% Poisson ratio of 0.3
f_ray2 = -130.253*exp(-1.8362*kz) + 6.88812*exp(-1.7556*kz) + ...
    271.540*exp(-1.7380*kz) - 12.1077*exp(-1.6859*kz) - ...
    183.478*exp(-1.6750*kz) + 1.70238*exp(-1.6574*kz) - ...
    142.361*exp(-1.6398*kz) + 341.126*exp(-1.6053*kz) + ...
    4.76194*exp(-1.5877*kz) - 159.030*exp(-1.5356*kz);
% Poisson ratio of 0.25
f_ray1 = -103.14*exp(-1.8586*kz) + 6.1446*exp(-1.7714*kz) + ...
    217.120*exp(-1.7555*kz) - 10.312*exp(-1.7012*kz) - ...
    160.68*exp(-1.6842*kz) + 1.2856*exp(-1.6683*kz) - ...
    115.00*exp(-1.6524*kz) + 294.66*exp(-1.6140*kz) + ...
    4.0924*exp(-1.5981*kz) - 135.49*exp(-1.5438*kz);
% Love waves
f_love = -(1+0.85^2)*exp(-2*0.85*kz);

% Construct G matrix
if inv_flag == 1
    G = diff(f_homo,1,2);
elseif inv_flag == 2
    G = diff(f_ray1,1,2);
elseif inv_flag == 3
    G = diff(f_love,1,2);
elseif inv_flag == 4
    G = diff(f_homo_U,1,2);
elseif inv_flag == 5
    G = diff(f_ray2,1,2);
elseif inv_flag == 6
    G = ((1-alph)^2)*diff(f_ray1,1,2);
elseif inv_flag == 7
    G = ((1-alph)^2)*diff(f_love,1,2);
else
    G = ((1-alph)^2)*diff(f_ray2,1,2);
end
