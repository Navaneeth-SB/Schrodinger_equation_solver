% Program to calculate the lowest eigenvalues and eigenfunctions of the 1D Schrödinger equation using finite-difference methods.
% Author : Navaneeth_S.B
% Date: January 2025
% Description: This script solves the 1D Schrödinger equation numerically
%              and computes the eigenvalues and eigenfunctions using 
%              finite-difference methods.
% LinkedIn : https://www.linkedin.com/in/navaneeth-shetty-b-8b7536318/

% Constants
m = 9.1e-31;                   % Mass of particle (electron) in Kg
hcut = 1.0546e-34;              % Reduced plank's constant in J.s

% Discretization parameter
L = 1e-9;                      % Interval length (domain)
N = 100;                       % No. of grid points.
x = linspace(0, L, N);         % coordinate vector
dx = L/N ;                     % coordinate setup (grid spacing)

%Two-point finite difference representation of Derivative. 
D = (diag(ones((N-1),1),1) - diag(ones((N-1),1),-1)) / (2*dx);

% Boundary condition for D matrix
% modifying D so that it is consistent with f(0) = f(L) = 0
D(1, 1) = 0; D(1, 2) = 0; D(2, 1) = 0;
D(N, N-1) = 0; D(N, N) = 0;

% Laplacian operator usins three point difference method
Lap = (diag(ones((N-1),1),-1) - 2*diag(ones(N,1),0) ...
    + diag(ones((N-1),1),1)) / (dx*2);

% Kinetic energy operator
T = -((hcut^2)/(2*m))*Lap;

% Potential energy operator (infinite potential well)
U_diag = zeros(N, 1); % Modify this to add potential terms
U = diag(U_diag);

% Hamiltonian operator
H = T + U; 

% Solve eigenvalue problem
[eigenfunction , eigenvalue] = eig(H); %this gives us eigenvalue and function from H
E = diag(eig(H));
Eev = E/(1.6e-19); 
% plot first three eigen functions
for n = 1:3
    subplot(3,1,n);
    plot(x, eigenfunction(:,n));
end