function [out] = h(X,T)

%% Description
%
% Function for specific enthalpy calculation 
% of gas with composition X and temperature T,
% using the 9 Term NASA Polynomials.
%
%   Cp/R = a1 T^-2 + a2 T^-1 + a3 + a4 T + a5 T^2 + a6 T^3 + a7 T^4 
%    
%   H/RT = -a1 T^-2 + (a2/T) ln(T) + a3 + (a4/2) T + (a5/3) T^2 + (a6/4)
%          T^3 + (a7/5) T^4 + a8/T  
%
%   S/R = -(a1/2) T^-2 - (a2/T) + a3 ln(T) + a4 T + (a5/2) T^2 + (a6/3)
%         T^3 + (a7/4) T^4 + a9
%    
%% Input:
%
% T - Temperature of mixture                                   [K]
%
% X - Structure array: each element X(i) is a structure of an element
%     containing also the molar fraction, the molar mass and the coefficients
%
% a1, a2, a3, a4, a5, a6, a7, a8, a9 - 200 - 1000 K polynomial coefficients
% b1, b2, b3, b4, b5, b6, b7, b8, b9  - 1000 - 6000 K polynomial coefficients 
%
% R = 8.314472 - universal gas constant                        [J.mol^-1.K^-1]
%
% MM - molar mass                                              [kg.kmol^-1]
%
% Tref = 298.15                                                [K]
%
%% Output
%
% h - Specific enthalpy of gas mixture                         [kJ.kg^-1.K^-1]
%     defined for temperature T =<200;1000> T =<1000;6000>
%
%% The 9 term polynomials actually include 20 constants. 
%  The first set of 9 constants belong to the 200 - 1000 K polynomial, 
%  the second set of 9 constants belong to the 1000 - 6000 K polynomial 
%  and the tenth constant is H298/R ≡ ΔfH298/R.
%
%%

if T < 200 || T > 6000
   % printf('Temperature %d out of range 200K - 6000K \n', T );
end

R=8.314472;  % [J.mol^-1.K^-1]
base = [ -1*T^(-2), 1/T*log(T), 1, 1/2*T, 1/3*T^2, 1/4*T^3, 1/5*T^4, 1/T];
MMi=[];
hT = [];
for i=1:columns(X)
   
  if T<1000
    hT(i) =  sum(X(i).a(1:8).*base)*X(i).fraction;
    MMi(i) =  X(i).fraction * X(i).MM;
  elseif T>1000
    hT(i) =  sum(X(i).b(1:8).*base)*X(i).fraction;    
    MMi(i) =  X(i).fraction * X(i).MM;
  end
  
end

out=R*T*sum(hT)/sum(MMi);

endfunction
