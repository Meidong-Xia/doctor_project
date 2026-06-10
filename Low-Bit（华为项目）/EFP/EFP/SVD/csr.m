function [c,s,rr] = csr(x,y)
if(y == 0)
    c = 1;
    s = 0;
    rr = x;
else
    if(abs(y) > abs(x))
        tao = real(i_hp_div(-x,y));
%         tao=-x/y;
        tao2 = real(i_hp_mul(tao,tao));
        opt2 = real(i_hp_add(1,tao2));
        s = real(i_hp_sqrt(opt2));
%         s=sqrt(1+tao^2);
        rr = real(i_hp_mul(-y,s));
%         rr=-y*s;
        s = real(i_hp_div(1,s));
%         s=1/s;
        c = real(i_hp_mul(s,tao));
%         c=s*tao;
    else
        tao = real(i_hp_div(-y,x));
%         tao=-y/x;
        tao2 = real(i_hp_mul(tao,tao));
        opt2 = real(i_hp_add(1,tao2));
        c = real(i_hp_sqrt(opt2));
%         c=sqrt(1+tao^2);
        rr = real(i_hp_mul(x,c));
%         rr=x*c;
        c = real(i_hp_div(1,c));
%         c=1/c;
        s = real(i_hp_mul(c,tao));
%         s=c*tao;
    end
end
end