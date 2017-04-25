classdef laserSource2D
  properties
    m_speed
    m_angle % in radians
    m_diameter
    m_fluxDens
    m_startX
    m_startY
    m_nx
    m_ny
    m_dx
    m_dy
    m_startTime
  end
  methods
    % constructor
    function obj = laserSource2D(nx, ny, dx, dy, speed, angle, ...
        diameter, fluxDens, startX, startY, startTime)
      obj.m_nx = nx;
      obj.m_ny = ny;
      obj.m_dx = dx;
      obj.m_dy = dy;
      obj.m_speed = speed;
      obj.m_angle = angle;
      obj.m_diameter = diameter;
      obj.m_fluxDens = fluxDens;
      obj.m_startX = startX;
      obj.m_startY = startY;
      obj.m_startTime = startTime;
    end
%    function endTime = getLineEndTime(obj)
%      endX = obj.m_ 
%    end
    function gridFlux = getFaceFluxAtTime(obj, t, gridFlux)
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

      for j=yStartIdx:yEndIdx
        for i=xStartIdx:xEndIdx
          gridFlux(j, i) = gridFlux(j, i) + obj.m_fluxDens;
        end
      end
    end % getFaceFluxAtTime
  end % methods section
end

