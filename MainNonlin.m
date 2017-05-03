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

numDim = 3; % 2 for 2D, 3 for 3D

dx = 0.005;
dy = 0.005;
dz = 0.005;
Lx = 0.3;
Ly = 0.3;
Lz = 0.025;
nx = Lx/dx + 1;
ny = Ly/dy + 1;
nz = Lz/dz + 1;

alph = 1.1234e-4;
Tinit = 300;
Tb1 = Tinit+5;
Tb2 = Tinit;
Tb3 = Tinit;
Tb4 = Tinit;
Tb5 = Tinit;  % z-bottom
Tb6 = Tinit;  % z-top

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

dt = 0.01; % TODO: determine this based on stability condition
tol = 0.0001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize field and set boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (numDim == 2)
  Tnp1 = ones(ny,nx)*Tinit;
  Tnp1(1:ny,1) = Tb4;
  Tnp1(1:ny,nx) = Tb2;
  Tnp1(1,1:nx) = Tb1;
  Tnp1(ny,1:nx) = Tb3;
elseif (numDim == 3)
  Tnp1 = ones(ny,nx,nz)*Tinit; % TODO: verify 3D column major access for y vs z
  Tnp1(:, 1, :) = Tb4;  % -yz plane
  Tnp1(:, nx, :) = Tb2; % +yz plane
  Tnp1(:, :, 1) = Tb5;  % -xy plane
  Tnp1(:, :, nz) = Tb6; % +xy plane
  Tnp1(1, :, :) = Tb1;  % -xz plane
  Tnp1(ny, :, :) = Tb3; % +xz plane
end
Tn = Tnp1;  % temperatures at current step
Tnm1 = Tnp1;  % temperatures at previous step

% volumetric source %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
volFlux = 0.0003;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize loop variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TV = 1;
prevTV = 0;
dTV = 1;
err = 0;
time = 0;

% setup laser flux
if (numDim == 2)
  Sp = zeros(ny, nx);
  speed = 0.2;
  angle = 0.5*pi;
  diameter = 0.01;
  power = 100;
  % fluxDens = power/(dx*dy);
  fluxDens = volFlux;
  startX = 0.015;
  startY = 0*Ly;
  startTime = 0;

  laser1 = laserSource2D(nx, ny, dx, dy, 0.7*speed, angle, diameter, ...
      fluxDens, 3*dx, 0, startTime)
  laser2 = laserSource2D(nx, ny, dx, dy, 1*speed, angle, diameter, ...
      fluxDens, 0.333*Lx, 0, startTime)
  laser3 = laserSource2D(nx, ny, dx, dy, 1.2*speed, angle, diameter, ...
      fluxDens, 0.6667*Lx, 0, startTime)
  laser4 = laserSource2D(nx, ny, dx, dy, 1.7*speed, angle, diameter, ...
      fluxDens, Lx-3*dx, 0, startTime)

elseif (numDim == 3)
  Sp = zeros(ny, nx, nz);
  speed = 0.1;
  angle = 0.0;
  diameter = 0.02;
  power = 10;
  fluxDens = power*(dx*dy*dz);
  %fluxDens = 3*volFlux;
  startX = 0;
  startY = 0.5*Ly;
  startTime = 0;
  laser1 = laserSource3D(nx, ny, nz, dx, dy, dz, speed, angle, diameter, ...
      fluxDens, startX, startY, startTime)
end

tic
while dTV >= tol
  if (numDim == 2)
    Sp = zeros(ny, nx);
    Sp = laser1.getFaceFluxAtTime(time*dt, Sp);
    Sp = laser2.getFaceFluxAtTime(time*dt, Sp);
    Sp = laser3.getFaceFluxAtTime(time*dt, Sp);
    Sp = laser4.getFaceFluxAtTime(time*dt, Sp);
  elseif (numDim == 3)
    Sp = zeros(ny, nx, nz);
    Sp = laser1.getFaceFluxAtTime(time*dt, Sp);
  end

  % call solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (numDim == 2)
    [Tnp1] = nonlinJacobi(Tnm1, Tn, Sp, alph, dx, dy, dt, nx, ny);
  elseif (numDim == 3)
    [Tnp1] = nonlinJacobi3D(Tnm1, Tn, Sp, alph, dx, dy, dz, dt, nx, ny, nz);
  end
  assert(size(Tnp1, 1) == size(Tnm1, 1) && ...
      size(Tnp1, 2) == size(Tnm1,2), 'Updated T field changed its size.');

  % check convergence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % err = 0;
  for i=2:ny-1
    for j=2:nx-1
      if (numDim==2)
        err = err + abs(Tnp1(i,j) - Tnm1(i,j));
      elseif (numDim==3)
        for k=2:nz-1
          err = err + abs(Tnp1(i,j,k) - Tnm1(i,j,k));
        end
      end
    end
  end
  if (numDim==2)
    TV = 1/(nx*ny)*err
  elseif (numDim==3)
    TV = 1/(nx*ny*nz)*err
  end
  dTV = TV - prevTV
  prevTV = TV
  time = time + 1

  Tnm1 = Tn;
  Tn = Tnp1;

  % plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(1)
  Tmax = 950;
  %Tmax = max(max(Tnp1));
  if (numDim==2)
    contourf(Tnp1(:,:), [0.99*Tmax 0.98*Tmax 0.96*Tmax 0.93*Tmax 0.90*Tmax 0.88*Tmax 0.84*Tmax ...
      0.8*Tmax 0.75*Tmax 0.7*Tmax 0.65*Tmax 0.6*Tmax 0.55*Tmax 0.50*Tmax ...
      0.45*Tmax  0.40*Tmax 0.34*Tmax min(min(Tnp1))]);
  elseif (numDim==3)
    Tplane = squeeze(Tnp1(:, :, nz));
    contourf(Tplane)
%    contourf(Tplane, [0.99*Tmax 0.98*Tmax 0.96*Tmax 0.93*Tmax 0.90*Tmax 0.88*Tmax 0.84*Tmax ...
%      0.8*Tmax 0.75*Tmax 0.7*Tmax 0.65*Tmax 0.6*Tmax 0.55*Tmax 0.50*Tmax ...
%      0.45*Tmax]);
%    contourf(Tplane, [0.99*Tmax 0.98*Tmax 0.96*Tmax 0.93*Tmax 0.90*Tmax 0.88*Tmax 0.84*Tmax ...
%      0.8*Tmax 0.75*Tmax 0.7*Tmax 0.65*Tmax 0.6*Tmax 0.55*Tmax 0.50*Tmax ...
%      0.45*Tmax  0.40*Tmax 0.34*Tmax min(min(Tnp1))]);
  end
  colorbar
  %caxis([300, 680])
  pause(0.05)

end
toc
