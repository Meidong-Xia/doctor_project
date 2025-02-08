function [G,U1,U2,V1,V2] = mul_22_submatrix(G,U1,U2,V1,V2)
G = real(G);
g11 = G(1,1);
g12 = G(1,2);
g22 = G(2,2);
a = real(i_hp_mul(g11,g11));
% a = g11^2;
g122 = real(i_hp_mul(g12,g12));
g222 = real(i_hp_mul(g22,g22));
b = real(i_hp_add(g122,g222));
% b = g12^2+g22^2;
c=real(i_hp_mul(g11,g12));
% c = g11*g12;
c2 = real(i_hp_mul(c,2));
if abs(c2)>1e-3
    bsa = real(i_hp_sub(b,a));
    theta = real(i_hp_div(bsa,c2));
    if abs(theta) > 255
        theta = real(i_hp_div(theta,100));
        theta2 = real(i_hp_mul(theta,theta));
        result = real(i_hp_add(0.0001,theta2));
        result1 = real(i_hp_sqrt(result));
        fenmu = real(i_hp_add(abs(theta),result1));
        t = real(i_hp_div(sign(theta),fenmu));
        t = real(i_hp_div(t,100));
    else
        theta2 = real(i_hp_mul(theta,theta));
        result = real(i_hp_add(1,theta2));
        result1 = real(i_hp_sqrt(result));
        fenmu = real(i_hp_add(abs(theta),result1));
        t = real(i_hp_div(sign(theta),fenmu));
    end
else
    t = 1e-3;
end
% theta = i_hp_div(bsa,c2);
% theta=(b-a)/(2*c);
% theta2 = i_hp_mul(theta,theta);
% optheta2 = i_hp_add(1,theta2);
% sq = i_hp_sqrt(optheta2);
% atpsq = i_hp_add(abs(theta),sq);
% t = i_hp_div(sign(theta),atpsq);
% t=sign(theta)/(abs(theta)+sqrt(1+theta^2));
t2 = real(i_hp_mul(t,t));
opt2 = real(i_hp_add(1,t2));
sqopt = real(i_hp_sqrt(opt2));
cs = real(i_hp_div(1,sqopt));
% cs=1/sqrt(1+t^2);
sn = real(i_hp_mul(cs,t));
% sn=cs*t;
G = i_hp_matrix_mul(complex(G),complex([cs,sn;-sn,cs]));
G = real(G);
% G=G*[cs,sn;-sn,cs];
[V1,V2] = updatecsvv(cs,sn,V1,V2);
G112 = i_hp_mul(G(1,1),G(1,1));
G212 = i_hp_mul(G(2,1),G(2,1));
gag = i_hp_add(G112,G212);
alpha = i_hp_sqrt(gag);
% alpha=sqrt(G(1,1)^2+G(2,1)^2);
G122 = i_hp_mul(G(1,2),G(1,2));
G222 = i_hp_mul(G(2,2),G(2,2));
gag = i_hp_add(G122,G222);
beta = i_hp_sqrt(gag);
% beta=sqrt(G(1,2)^2+G(2,2)^2);
c1 = i_hp_div(G(1,1),alpha);
% c1=G(1,1)/alpha;
c2 = i_hp_div(G(2,2),beta);
% c2=G(2,2)/beta;
s1 = i_hp_div(G(2,1),alpha);
% s1=G(2,1)/alpha;
s2 = i_hp_div(G(1,2),beta);
% s2=G(1,2)/beta;
G = i_hp_matrix_mul(complex([c1,s1;s2,c2]),complex(G));
% G=[c1,s1;s2,c2]*G;
G = real(G);
[U1,U2] = updatecsvv_last(c1,s1,c2,s2,U1,U2);
end