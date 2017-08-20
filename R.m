function [out] = R(X)
%% Description
%
% Function for specific gas costant R of gas with composition X.
%
%% Input:
%
% X - Structure array: each element X(i) is a structure of an element
%     containing also the molar fraction, the molar mass and the coefficients
%
%% Used coefficients and functions:
%
% R = 8.314472 universal gas constant          [J.mol^-1.K^-1]
% MM - molar mass                              [kg.kmol^-1]
%
%% Output
%
% R(X) - specific gas costant                  [J.kg^-1.K^-1]
%
%%
for i=1:columns(X)
     R(i)=X(i).fraction/X(i).MM;
end
out = 8.314472*sum(R);
endfunction
