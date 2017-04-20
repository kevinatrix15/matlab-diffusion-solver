function [T] = gaussSeidelSOR(Tprev, alph, dx, dy, dt, nx, ny, omeg)

  T = Tprev;

  rx = alph*dt/(dx^2);
  ry = alph*dt/(dy^2);
  % loop over domain in column-major order for efficiency
  for col=2:nx-1
    for row=2:ny-1
      T(row,col) = (1-omeg)*Tprev(row,col) +...
        omeg*rx * (T(row, col-1) - 2*Tprev(row,col) + Tprev(row,col+1)) +...
        omeg*ry * (T(row-1, col) - 2*Tprev(row,col) + Tprev(row+1,col));
    end
  end
  
  T(:,1) = Tprev(:,1);
  T(:,nx) = Tprev(:,nx);
  T(1,:) = Tprev(1,:);
  T(ny,:) = Tprev(ny,:);
end

