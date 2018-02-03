function [force] = rotateMat(v, angle)

    RMat = [cosd(angle) -sind(angle);
            sind(angle) cosd(angle)]; % rotation matrix 

    force = v * RMat;
end

