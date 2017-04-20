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
mode = 1; % 0-Jacobi, 1-Gauss-Seidel w/ SOR
omeg = 1.0;

dx = 0.01;
dy = 0.01;
Lx = 0.3;
Ly = 0.3;
nx = Lx/dx + 1;
ny = Ly/dy + 1;

alph = 1.1234e-4;
Tinit = 300;
Tb1 = Tinit + 40;
Tb2 = Tinit + 30;
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

dt = 0.18; % TODO: determine this based on stability condition
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

tic
while dTV >= tol

  % call solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (mode == 0)
    [Tfield] = jacobi(Tprev, alph, dx, dy, dt, nx, ny);
  elseif (mode==1)
    %[Tfield] = gaussSeidelSOR(Tprev, alph, dx, dy, dt, nx, ny, omeg);
    rx = alph*dt/(dx^2);
    ry = alph*dt/(dy^2);
    B = dx/dy;
    % loop over domain in column-major order for efficiency
    for col=2:nx-1
      for row=2:ny-1
        %Tfield(row,col) = Tprev(row,col) +...
        %  rx*(Tfield(row, col-1) - 2*Tprev(row,col) + Tprev(row,col+1)) +...
        %  ry*(Tfield(row-1, col) - 2*Tprev(row,col) + Tprev(row+1,col));
        %Tfield(row,col) = (1-omeg)*Tprev(row,col) +...
        %  omeg*alph*dt/(2*(1+B^2))*...
        %  (Tfield(row, col-1) - 2*Tprev(row,col) + Tprev(row,col+1) +...
        %  B^2*(Tfield(row-1, col) - 2*Tprev(row,col) + Tprev(row+1,col)));
        Tfield(row,col) = (1-omeg)*Tprev(row,col) +...
          omeg*alph*dt/(dx^2)*...
          (Tfield(row, col-1) - 2*Tprev(row,col) + Tprev(row,col+1) +...
          B^2*(Tfield(row-1, col) - 2*Tprev(row,col) + Tprev(row+1,col)));
      end
    end
  end
  assert(size(Tfield, 1) == size(Tprev, 1) && ...
      size(Tfield, 2) == size(Tprev,2), 'Updated T field changed its size.');

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
