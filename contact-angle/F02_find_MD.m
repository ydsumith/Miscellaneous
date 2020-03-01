function [xy_COUNT,MODX,MODY,MODZ] = F02_find_MD(x,y,z,MD_cutoffMAX,MD_cutoffMIN,dropheight)
n = length(x);
x1 = x;
x2 = z;
x1m = mean(x1);
x2m = mean(x2);
Xc = [(x1-x1m) (x2-x2m)];
Cx = (Xc' * Xc)/(n-1);
Cinvx = inv(Cx);
MD1 = zeros(n,1);
for p1 = 1:n
    MD1(p1) = sqrt(Xc(p1,:) * Cinvx * Xc(p1,:)');
end
MD = MD1;

ii = 1;
while  ii <= n
    if (MD(ii) > MD_cutoffMAX || MD(ii) < MD_cutoffMIN) && z(ii) > dropheight % remove the vapors
        x(ii) = [];
        if length(y) ~= 0
            y(ii) = [];
        end
        z(ii) = [];
        MD(ii) = [];
        n = n - 1;
    else
        ii = ii + 1;
    end
end

xy_COUNT = n;
MODX = x;
MODY = y;
MODZ = z;
end