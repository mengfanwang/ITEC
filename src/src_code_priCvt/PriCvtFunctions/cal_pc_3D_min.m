function eig3 = cal_pc_3D_min(xx, yy, zz, xy, xz, yz) 


a = -(xx + yy + zz);
a = a(:);
b = -(yz.^2 + xy.^2 + xz.^2 - xx.*zz -yy.*zz -xx.*yy);
b = b(:);
c = -(xx.*yy.*zz + 2*xy.*yz.*xz - xx.*(yz.^2) - zz.*(xy.^2) - yy.*(xz.^2));
c = c(:);

p = b - a.^2/3;
q = 2*(a.^3)/27 - (a.*b)/3 + c;

% result3 = zeros(length(xx(:)),3);
% 
% result3(:,1) = 2/sqrt(3)*(sqrt(-p)).*cos(1/3.*acos(3*q./(2*p).*sqrt(-3./p)) - 2*pi/3) - a/3;
% result3(:,2) = 2/sqrt(3)*(sqrt(-p)).*cos(1/3.*acos(3*q./(2*p).*sqrt(-3./p)) - 4*pi/3) - a/3;
% result3(:,3) = 2/sqrt(3)*(sqrt(-p)).*cos(1/3.*acos(3*q./(2*p).*sqrt(-3./p))) - a/3;
% 
% eig3 = zeros(size(xx));
% eig3(:) = max(result3,[],2);

eig3 = zeros(size(xx));
temp=3*q./(2*p).*sqrt(-3./p);
temp=min(max(temp,-1),1);
eig3(:) = 2/sqrt(3)*(sqrt(-p)).*cos(1/3.*acos(temp) - 4*pi/3) - a/3;
end