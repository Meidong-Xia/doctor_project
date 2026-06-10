function [v1,v2] = updatecsvv(c,s,v1,v2)
[n,~] = size(v1);
for i = 1:n
    t = v1(i);
    ct = i_hp_mul(c,t);
    sv2 = i_hp_mul(s,v2(i));
    v1(i) = i_hp_sub(ct,sv2);
%     v1(i) = c*t-s*v2(i);
    st = i_hp_mul(s,t);
    cv2 = i_hp_mul(c,v2(i));
    v2(i) = i_hp_add(st,cv2);
%     v2(i)=s*t+c*v2(i);
end
end