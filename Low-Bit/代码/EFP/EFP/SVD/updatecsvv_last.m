function [v1,v2]=updatecsvv_last(c1,s1,c2,s2,v1,v2)
[n,~] = size(v1);
for i = 1:n
    t = v1(i);
    c1t = i_hp_mul(c1,t);
    s1v2 = i_hp_mul(s1,v2(i));
    v1(i) = i_hp_add(c1t,s1v2);
%     v1(i) = c1*t+s1*v2(i);
    s2t = i_hp_mul(s2,t);
    c2v2 = i_hp_mul(c2,v2(i));
    v2(i) = i_hp_add(s2t,c2v2);
%     v2(i) = s2*t+c2*v2(i);
end
end