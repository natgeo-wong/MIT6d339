function [u, Troot] = ThermalFin(mesh, mu)
%
% -------------------------------------------------------------------------
%   Computes the temperature distribution and root temperature for a fin
%   using the Finite Element Method.
% -------------------------------------------------------------------------
%
%   INPUT   mesh    grid lebel (coarse, medium, or fine)
%           mu      thermal conductivities of sections and Biot number 1x5
%
%   OUTPUT  u       temperature disctribution in the fin
%           Troot   root temperature
%
% -------------------------------------------------------------------------


% Parameters setup: rearranges the values in mu
kappa = ones(6,1);
kappa(1:4) = mu(1:4);    % Bi = mu(6)
kappa(6) = mu(5);


% Initialization
A = sparse(mesh.nodes, mesh.nodes);
F = sparse(mesh.nodes,1);


% Domain Interior
for i = 1:5         % interior regions
    for n = 1:length(mesh.theta{i})
        phi = mesh.theta{i}(n,:)';
        ....
        A(phi,phi) = A(phi,phi) + Alocal;    
    end
end


% Boundaries (not root)
i = 6;
    for n = 1:length(mesh.theta{i})
        phi = mesh.theta{i}(n,:)';
        ...
  
        A(phi,phi) = A(phi,phi) + Alocal;
    end


% Root Boundary
i = 7;
    for n = 1:length(mesh.theta{i})
        phi = mesh.theta{i}(n,:)';
        ...
            
        F(phi) = F(phi) + Flocal;
    end


% Solve for the temperature distribution
u = full(A\F);      % use full since A and F are sparse


% Compute the temperature at the root
...
    