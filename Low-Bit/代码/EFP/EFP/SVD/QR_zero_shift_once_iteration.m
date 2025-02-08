function [tri,U,VtT] = QR_zero_shift_once_iteration(tri,U,Vt,i_min,i_max)
n = i_max+1;
oldc = 1;
x = tri(i_min+1,i_min+1);
y = tri(i_min+1,i_min+2);
for i = i_min:n-2
    [c,s,r] = csr(x,y);
    VtT = Vt';
    [VtT(:,i+1),VtT(:,i+2)] = updatecsvv(c,s,VtT(:,i+1),VtT(:,i+2));
    Vt = VtT';
    if i ~= i_min
        tri(i,i+1) = i_hp_mul(-olds,r);
%         tri(i,i+1)=-olds*r;
    end
    x = i_hp_mul(oldc,r);
%     x=oldc*r;
    y = i_hp_mul(-tri(i+2,i+2),s);
%     y=-tri(i+2,i+2)*s;
    h = i_hp_mul(tri(i+2,i+2),c);
%     h=tri(i+2,i+2)*c;
    [c,s,r] = csr(x,y);
    [U(:,i+1),U(:,i+2)] = updatecsvv(c,s,U(:,i+1),U(:,i+2));
    tri(i+1,i+1) = r;
    x = h;
    if i ~= n-2
        y = tri(i+2,i+3);
    end
    oldc = c;
    olds = s;
end
tri(n-1,n) = i_hp_mul(-h,s);
% tri(n-1,n)=-h*s;
tri(n,n) = i_hp_mul(h,c);
% tri(n,n)=h*c;
end