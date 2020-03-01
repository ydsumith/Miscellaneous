function [x] = F05_rotate(x1,y1,ang_theta)
    Xcent = mean(x1);
    Ycent = mean(y1);
    shift_x = x1 - Xcent;
    shift_y = y1 - Ycent;
    COORDXY = [shift_x shift_y];
    COORDXY = COORDXY';

    R_theta = [cos(ang_theta) -sin(ang_theta);...
               sin(ang_theta) cos(ang_theta)];
    ROTATED_M = R_theta * COORDXY;
    shift_x = ROTATED_M(1,:) + Xcent;
    shift_y = ROTATED_M(2,:) + Ycent;
    x = shift_x';
end