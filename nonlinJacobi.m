function [Tnp1] = nonlinJacobi(Tnm1, Tn, Sp, alph, dx, dy, dt, nx, ny)

  % TODO: use Sp, or switch to face flux??

  Tnp1 = Tn;

  % TODO: expand this for 3D
  Aew = dy;
  Ans = dx;
  Vol = dx*dy;

  % loop over domain in column-major order for efficiency
  for col=2:nx-1
    for row=2:ny-1
      % get neighboring diffusivity
      alphP = getDiffusivity(Tn(row,col));
      alphE = 0.5*(getDiffusivity(Tn(row,col+1)) + alphP);
      alphW = 0.5*(getDiffusivity(Tn(row,col-1)) + alphP);
      alphN = 0.5*(getDiffusivity(Tn(row+1,col)) + alphP);
      alphS = 0.5*(getDiffusivity(Tn(row-1,col)) + alphP);
      
      % determine neighbor coefficients
      aE = alphE*Aew/dx;
      aW = alphW*Aew/dx;
      aN = alphN*Ans/dy;
      aS = alphS*Ans/dy;
      aP = aE + aW + aN + aS - Sp(row,col);

      Tnp1(row,col) = 1/(1+dt/Vol*aP)*(...
        (1-dt/Vol*aP)*Tnm1(row,col) +...
        (aE*Tn(row, col+1) + aW*Tn(row, col-1) +...
         aN*Tn(row+1, col) + aS*Tn(row-1, col))*(2*dt/Vol));
      % TODO: experiment with SOR
      % TODO: add face source term
    end
  end
  
  % maintain boundary conditions
  Tnp1(:,1) = Tn(:,1);
  Tnp1(:,nx) = Tn(:,nx);
  Tnp1(1,:) = Tn(1,:);
  Tnp1(ny,:) = Tn(ny,:);
end
