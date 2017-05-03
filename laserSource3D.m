classdef laserSource3D
  properties
    m_speed
    m_angle % in radians
    m_diameter
    m_fluxDens
    m_startX
    m_startY
    m_nx
    m_ny
    m_nz
    m_dx
    m_dy
    m_dz
    m_startTime
  end
  methods
    % constructor
    function obj = laserSource3D(nx, ny, nz, dx, dy, dz, speed, angle, ...
        diameter, fluxDens, startX, startY, startTime)
      obj.m_nx = nx;
      obj.m_ny = ny;
      obj.m_nz = nz;
      obj.m_dx = dx;
      obj.m_dy = dy;
      obj.m_dz = dz;
      obj.m_speed = speed;
      obj.m_angle = angle;
      obj.m_diameter = diameter;
      obj.m_fluxDens = fluxDens;
      obj.m_startX = startX;
      obj.m_startY = startY;
      obj.m_startTime = startTime;
    end
    function modGridFlux = getFaceFluxAtTime(obj, t, gridFlux)
      dTime = t - obj.m_startTime;
      laserPosX = (obj.m_speed*dTime)*cos(obj.m_angle) + obj.m_startX;
      laserPosY = (obj.m_speed*dTime)*sin(obj.m_angle) + obj.m_startY;
      xIdx = min(obj.m_nx, round(laserPosX/obj.m_dx) + 1);
      yIdx = min(obj.m_ny, round(laserPosY/obj.m_dy) + 1);

      % initially, only assign flux at top face
      xStartIdx = max(1, round((laserPosX - 0.5*obj.m_diameter)/obj.m_dx));
      yStartIdx = max(1, round((laserPosY - 0.5*obj.m_diameter)/obj.m_dy));
      xEndIdx = min(obj.m_nx, round((laserPosX + 0.5*obj.m_diameter)/obj.m_dx));
      yEndIdx = min(obj.m_ny, round((laserPosY + 0.5*obj.m_diameter)/obj.m_dy));

      k = obj.m_nz; % TODO: modify to attenuate with z
      modGridFlux = gridFlux;
      for j=yStartIdx:yEndIdx
        for i=xStartIdx:xEndIdx
          modGridFlux(j, i, k) = gridFlux(j, i, k) + obj.m_fluxDens;
        end
      end
    end % getFaceFluxAtTime
  end % methods section
end
