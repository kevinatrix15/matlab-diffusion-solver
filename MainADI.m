clear all, clc, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REQUIREMENTS:
% 1. Make flexible enought to handle a number of solution schemes (e.g.,
% Crank-Nicolson, Gauss-Seidell w/ SOR, ADI, Multigrid, etc.)
% 2. Add functionality to handle time-varying source terms.
% 3. Handle different boundary conditions (Dirichlet, Neumann, mixed).

% FUTURE DEVELOPMENTS:
% 1. Extend to 3D solution.
% 2. Handle non-uniform spacing (eventually).
% 3. Set up for octree implementation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx = 0.01;
dy = 0.01;
Lx = 0.3;
Ly = 0.3;
nx = Lx/dx + 1;
ny = Ly/dy + 1;

alph = 1.1234e-4;
Tinit = 0;
Tb1 = Tinit + 40;
Tb2 = Tinit + 0;
Tb3 = Tinit + 10;
Tb4 = Tinit + 0;

%   +-----------------------+
%   |           3           |
%   |                       |
%   |                       |
%   |                       |
%   |                       |
%   | 4                   2 |   y
%   |                       |
%   |                       |   ^
%   |                       |   |
%   |           1           |   |
%   +-----------------------+   +----> x

dt = 0.5; % TODO: determine this based on stability condition
tol = 0.0001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tfield = zeros(ny,nx);
Tfield = ones(ny,nx)*Tinit;

% set boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tfield(1:ny,1) = Tb4;
Tfield(1:ny,nx) = Tb2;
Tfield(1,1:nx) = Tb1;
Tfield(ny,1:nx) = Tb3;

Tp = Tfield;  % temperatures at half step

% build matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rx = 0.5*alph*dt/(dx^2)
ry = 0.5*alph*dt/(dy^2)

A1 = zeros(nx-2, nx-2);
A1(1,1) = 1+2*rx;
A1(1,2) = -rx;
for i=2:nx-3
  A1(i, i-1) = -rx;
  A1(i, i) = 1+2*rx;
  A1(i, i+1) = -rx;
end
A1(nx-2,nx-3) = -rx;
A1(nx-2,nx-2) = 1+2*rx;

A2 = zeros(ny-2, ny-2);
A2(1,1) = 1+2*ry;
A2(1,2) = -ry;
for i=2:ny-3
  A2(i, i-1) = -ry;
  A2(i, i) = 1+2*ry;
  A2(i, i+1) = -ry;
end
A2(ny-2,ny-3) = -ry;
A2(ny-2,ny-2) = 1+2*ry;



% we wish to solve to steady state

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize loop variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TV = 1;
prevTV = 0;
dTV = 1;
err = 0;
Tprev = Tfield;
time = 0;

invA1 = inv(A1);
invA2 = inv(A2);

tic
while dTV >= tol

  % perform solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % x-sweep
  for i=2:ny-1
    %rhs1 = zeros(nx-2, 1);
    for j=2:nx-1
      rhs1(j-1, 1) = ry*Tfield(i-1,j) + (1-2*ry)*Tfield(i,j) + ry*Tfield(i+1,j);
    end
    Tp(i,2:nx-1) = invA1 * rhs1;
  end

  % y-sweep
  for j=2:nx-1
    %rhs = zeros(ny-2, 1);
    for i=2:ny-1
      rhs(i-1, 1) = rx*Tp(i,j-1) + (1-2*rx)*Tp(i,j) + rx*Tp(i,j+1);
    end
    Tfield(2:ny-1,j) = invA2 * rhs;
  end

  % reset BCs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % TODO: probably not necessary
  %Tfield(1:ny,nx) = Tb2;
  %Tfield(1:ny,1) = Tb4;
  %Tfield(1,1:nx) = Tb1;
  %Tfield(ny,1:nx) = Tb3;

  % check convergence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % err = 0;
  for i=2:ny-1
    for j=2:nx-1
      err = err + abs(Tfield(i,j) - Tprev(i,j));
    end
  end
  TV = 1/(nx*ny)*err
  dTV = TV - prevTV
  prevTV = TV
  time = time + 1

  Tprev = Tfield;

  % plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(1)
  contourf(Tfield(:,:));
  colorbar
  %caxis([150, 340])
  pause(0.05)

end
toc
