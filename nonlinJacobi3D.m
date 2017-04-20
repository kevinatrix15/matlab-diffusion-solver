function [Tnp1] = nonlinJacobi(Tnm1, Tn, Sp, alph, dx, dy, dz, dt, nx, ny, nz)

  % TODO: use Sp, or switch to face flux??

  Tnp1 = Tn;

  Aew = dy*dz;
  Ans = dx*dz;
  Aud = dx*dy;
  Vol = dx*dy*dz;

  % loop over domain in column-major order for efficiency
  for lay=2:nz-1
    % get neighboring diffusivity
    alphU = 0.5*(alph + alph); % TODO: add function call for nonlinear diffus.
    alphD = 0.5*(alph + alph); % TODO: add function call for nonlinear diffus.

    % determine neighbor coefficients
    aU = alphU*Aud/dz;
    aD = alphD*Aud/dz;

    for col=2:nx-1
      % get neighboring diffusivity
      alphE = 0.5*(alph + alph); % TODO: add function call for nonlinear diffus.
      alphW = 0.5*(alph + alph); % TODO: add function call for nonlinear diffus.

      % determine neighbor coefficients
      aE = alphE*Aew/dx;
      aW = alphW*Aew/dx;

      for row=2:ny-1
        % get neighboring diffusivity
        alphN = 0.5*(alph + alph); % TODO: add function call for nonlinear diffus.
        alphS = 0.5*(alph + alph); % TODO: add function call for nonlinear diffus.
        
        % determine neighbor coefficients
        aN = alphN*Ans/dy;
        aS = alphS*Ans/dy;
        aP = aE + aW + aN + aS + aU + aD- Sp(row,col);

        Tnp1(row,col,lay) = 1/(1+dt/Vol*aP)*(...
          (1-dt/Vol*aP)*Tnm1(row,col,lay) +...
          (aE*Tn(row, col+1,lay) + aW*Tn(row, col-1,lay) +...
          aU*Tn(row, col,lay+1) + aD*Tn(row, col,lay-1) +...
          aN*Tn(row+1, col,lay) + aS*Tn(row-1, col,lay))*(2*dt/Vol));
        % TODO: experiment with SOR
        % TODO: add face source term for convective flux
      end
    end
  end
  
  % maintain boundary conditions
  Tnp1(:,1, :) = Tn(:,1, :);
  Tnp1(:,nx, :) = Tn(:,nx, :);
  Tnp1(1,:, :) = Tn(1,:, :);
  Tnp1(ny,:, :) = Tn(ny,:, :);
  Tnp1(:,:, 1) = Tn(:,:, 1);
  Tnp1(:,:, nz) = Tn(:,:, nz);
end

