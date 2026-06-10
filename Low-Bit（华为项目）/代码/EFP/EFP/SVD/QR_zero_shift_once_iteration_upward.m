function [tri,U] = QR_zero_shift_once_iteration_upward(tri,U,i_min,i_max)
n = i_max+1;
oldc = 1;
x = tri(n,n);
y = tri(n-1,n);
for i = n-1:-1:i_min
    [c,s,r] = csr(x,y);
    [U(:,i+1),U(:,i)] = updatecsvv(c,s,U(:,i+1),U(:,i));
    if i ~= n-1
        tri(i+1,i+2) = i_hp_mul(-olds,r);
%         tri(i+1,i+2)=-olds*r;
    end
    x = i_hp_mul(oldc,r);
%     x=oldc*r;
    y = i_hp_mul(-tri(i,i),s);
%     y=-tri(i,i)*s;
    h = i_hp_mul(tri(i,i),c);
%     h=tri(i,i)*c;
    [c,s,r] = csr(x,y);
    oldc = c;
    olds = s;
    tri(i+1,i+1) = r;
    x = h;
    if(i ~= 1)
        y = tri(i-1,i);
    end
end
tri(i_min+1,i_min+2) = i_hp_mul(-h,s);
% tri(i_min+1,i_min+2)=-h*s;
tri(i_min+1,i_min+1)=i_hp_mul(h,c);
% tri(i_min+1,i_min+1)=h*c;
end