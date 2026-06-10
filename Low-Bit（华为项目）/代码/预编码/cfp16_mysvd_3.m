% 半精度下的SVD，输入任意规模实数、复数矩阵Z，如4*4、8*4、4*8均可实现分解Z=U*A*V'
% 但是由于精度限制，可能陷入死循环或者产生NaN，于是在本文件第156-158行设置了迭代次数上限为10000，可以根据实际需求将其更改（矩阵规模越大需要的迭代次数越多）
% 对于生成的SVD分解，需要外部调用（可以使用norm(Z-U*A*V')，正常运转时取值在[0,1]）手动将NaN值和溢出错误的分解排除
% 经测试，在64*8矩阵规模下基本相对误差已经达到50%左右，矩阵规模越大越容易溢出错误
% 这个版本为1.0（2023/11/29）更改，调整迭代次数和迭代判断精度来尽可能减少64*8矩阵的计算速度
% 这个版本为2.0（2023/11/29）更改，尝试调整设置主对角线非零
function [U,A,Vt]=cfp16_mysvd_3(mat)
[size_matm,size_matn]=size(mat);
if size_matm>size_matn
    islong=1;
else
    islong=0;
end
if ~islong
    mat=mat';
end
[U,A,V]=twodiag(mat);
[U,A,Vt]=bi_diag_svd(U,A,V);
[m0,n0]=size(A);
[lenu,~]=size(U);
[lenv,~]=size(Vt);
n=min(m0,n0);
for i=1:n-1
    if real(A(i,i))<0
        A(i,i)=-A(i,i);
        for k=1:lenu
            U(k,i)=-U(k,i);
        end
    end
end
for i=1:n-1
    for j=i:n
        if(abs(A(i,i))<abs(A(j,j)))
            tmp1=A(i,i);
            A(i,i)=A(j,j);
            A(j,j)=tmp1;
            for k=1:lenu
                temp2=U(k,i);
                U(k,i)=U(k,j);
                U(k,j)=temp2;
            end
            for k=1:lenv
                tmp2=Vt(i,k);
                Vt(i,k)=Vt(j,k);
                Vt(j,k)=tmp2;
            end
        end
    end
end
if ~islong
    T=Vt';
    Vt=U';
    U=T;
    A=A';
end
Vt=Vt';
[m,n]=size(A);
ZERO=zeros(m,n);
for i=1:min(m,n)
    ZERO(i)=A(i,i);
end
A=zeros(m,n);
for i=1:min(m,n)
    A(i,i)=ZERO(i);
end
end


function [U,A,V]=twodiag(A0)
% 二对角化,使得U'*A0*Vt'=A,对于复数矩阵而言，生成的U,V为复数矩阵，A为实数二对角矩阵。
A=A0;
[m,n]=size(A);
if m>=n
    is=1;
else
    is=0;
end
if ~is
    A=A';
end
[m,n]=size(A);
U=eye(m);
V=eye(n);

for j=1:n
    UK=eye(m);
    VK=eye(n);
    [UJH,~,~,~]=newht(A(j:m,j));
%   UK(j:m,j:m)=UK(j:m,j:m)*UJH';
    [~,~,UK(j:m,j:m)]=i_hp_matrix_mul(complex(UK(j:m,j:m)),complex(UJH'));
%   U=U*UK;
    [~,~,U]=i_hp_matrix_mul(complex(U),complex(UK));
%   A(j:m,j:n)=UJH*A(j:m,j:n);
    [~,~,A(j:m,j:n)]=i_hp_matrix_mul(complex(UJH),complex(A(j:m,j:n)));
    if j<n-1
        [VJH,~,~,~]=newht(A(j,j+1:n)');
%       VK(j+1:n,j+1:n)=VK(j+1:n,j+1:n)*VJH';
        [~,~,VK(j+1:n,j+1:n)]=i_hp_matrix_mul(complex(VK(j+1:n,j+1:n)),complex(VJH'));
%       V=V*VK;
        [~,~,V]=i_hp_matrix_mul(complex(V),complex(VK));
%       A(j:m,j+1:n)=A(j:m,j+1:n)*VJH';
        [~,~,A(j:m,j+1:n)]=i_hp_matrix_mul(complex(A(j:m,j+1:n)),complex(VJH'));
    end
end
if ~is
    A=A';
    u=U;
    U=V;
    V=u;
end
end

function [sigma_low_bound,lamda,u]=estimation(B)
[m0,n0]=size(B);
n=min(m0,n0);
lamda=zeros(n,1);
lamda(n)=abs(B(n,n));
for j=n-2:-1:0
    [~,~,fenmu]=i_hp_add(lamda(j+2),abs(B(j+1,j+2)));
    [~,~,result0]=i_hp_div(lamda(j+2),fenmu);
    [~,~,lamda(j+1)]=i_hp_mul(abs(B(j+1,j+1)),result0);
%     lamda(j+1)=abs(B(j+1,j+1))*(lamda(j+2)/(lamda(j+2)+abs(B(j+1,j+2))));
end
u=zeros(n,1);
u(1)=abs(B(1,1));
for j=0:n-2
    [~,~,fenmu]=i_hp_add(u(j+1),abs(B(j+1,j+2)));
    [~,~,result0]=i_hp_div(u(j+1),fenmu);
    [~,~,u(j+2)]=i_hp_mul(abs(B(j+2,j+2)),result0);
%     u(j+2)=abs(B(j+2,j+2))*(u(j+1)/(u(j+1)+abs(B(j+1,j+2))));
end
B_Infinity=min(lamda);
B_1=min(u);
sigma_low_bound=min(B_Infinity,B_1);

end

function [U,BI,V]=bi_diag_svd(U,BI,V)
underflow=1e-2;
tol=1e-2;
epcl=1e-6;
iter_num=0;
old_i_low=-1;
old_i_high=-1;
[m,n]=size(BI);
fudge=min(m,n);
% maxiter=4*m*n;
BI=real(BI);
% maxit=n*n;
[~,~,maxit]=r_hp_mul(n,n);% 避免溢出
[low_bound_sigma,lad,~]=estimation(BI);
max_bound_sigma=max(abs(BI),[],'all');
[~,~,tl]=r_hp_mul(tol,low_bound_sigma);
[~,~,mu]=r_hp_mul(maxit,underflow);
thresh=max(tl,mu);
while true
    if(iter_num>maxit)
        break;
    end
    iter_num=iter_num+1;
    i_u=n-1;
    while((i_u>=1)&&(abs(BI(i_u,i_u+1))<=1e-4))
        i_u=i_u-1;
    end
    if i_u==0
        break
    end
    i_l=i_u-1;
    if i_l~=0
        while((abs(BI(i_l,i_l+1))>1e-4)&&(i_u>=1))
            i_l=i_l-1;
            if i_l==0
                break
            end
        end
    end
    if i_u==i_l+1
        [BI(i_l+1:i_u+1,i_l+1:i_u+1),U(:,i_l+1),U(:,i_u+1),V(:,i_l+1),V(:,i_u+1)]=mul_22_submatrix(BI(i_l+1:i_u+1,i_l+1:i_u+1),U(:,i_l+1),U(:,i_u+1),V(:,i_l+1),V(:,i_u+1));
        continue
    end
    if(old_i_low~=i_l||old_i_high~=i_u)
        if(abs(BI(i_l+1,i_l+1))>=abs(BI(i_u+1,i_u+1)))
            direction=1;%'down'改为1
        else
            direction=2;%'up'改为2
        end
    end
    old_i_low=i_l;
    old_i_high=i_u;
    direction=1;

    if direction==1
        [~,~,result0]=i_hp_div(abs(BI(i_u,i_u+1)),lad(i_u+1));
        if(result0<=tol)
            BI(i_u,i_u+1)=0;
        end
    end
    %compute shift
    [~,~,ft]=r_hp_mul(fudge,tol);
    [~,~,ftl]=r_hp_mul(ft,low_bound_sigma);
    [~,~,result1]=r_hp_div(ftl,max_bound_sigma);
    if (result1<=epcl)
        shift=0;
    else
        if direction==1
            s=BI(i_u+1,i_u+1);
            [~,~,r0]=r_hp_mul(BI(i_u,i_u),BI(i_u,i_u));
            [~,~,r1]=r_hp_mul(BI(i_u-1,i_u),BI(i_u-1,i_u));
            [~,~,r2]=r_hp_add(r0,r1);
            [~,~,r3]=r_hp_mul(BI(i_u+1,i_u+1),BI(i_u+1,i_u+1));
            [~,~,r4]=r_hp_mul(BI(i_u,i_u+1),BI(i_u,i_u+1));
            [~,~,r5]=r_hp_add(r3,r4);
            [~,~,r6]=r_hp_sub(r2,r5);
            [~,~,d]=r_hp_div(r6,2);
%             d=((BI(i_u,i_u)*BI(i_u,i_u)+BI(i_u-1,i_u)*BI(i_u-1,i_u))-(BI(i_u+1,i_u+1)*BI(i_u+1,i_u+1)+BI(i_u,i_u+1)*BI(i_u,i_u+1)))/2;
            [~,~,r7]=r_hp_mul(BI(i_u+1,i_u+1),BI(i_u+1,i_u+1));
            [~,~,r8]=r_hp_mul(BI(i_u,i_u+1),BI(i_u,i_u+1));
            [~,~,r9]=r_hp_add(r7,r8);
            [~,~,r10]=r_hp_add(r9,d);
            [~,~,d2]=r_hp_mul(d,d);
            [~,~,r11]=r_hp_mul(BI(i_u,i_u),BI(i_u,i_u));
            [~,~,r12]=r_hp_mul(BI(i_u,i_u+1),BI(i_u,i_u+1));
            [~,~,r13]=r_hp_mul(r11,r12);
            [~,~,r14]=r_hp_add(d2,r13);
            [~,r15]=r_hp_sqrt(r14);
            [~,~,r16]=r_hp_mul(sign(d),r15);
            [~,~,shift]=r_hp_sub(r10,r16);
%             shift=(BI(i_u+1,i_u+1)*BI(i_u+1,i_u+1)+BI(i_u,i_u+1)*BI(i_u,i_u+1))+d-sign(d)*sqrt(d*d+BI(i_u,i_u)*BI(i_u,i_u)*BI(i_u,i_u+1)*BI(i_u,i_u+1));
        else
            s=BI(i_l+1,i_l+1);
            [~,~,r0]=r_hp_mul(BI(i_l+2,i_l+2),BI(i_l+2,i_l+2));
            [~,~,r1]=r_hp_mul(BI(i_l+2,i_l+3),BI(i_l+2,i_l+3));
            [~,~,r2]=r_hp_add(r0,r1);
            [~,~,r3]=r_hp_mul(BI(i_l+1,i_l+1),BI(i_l+1,i_l+1));
            [~,~,r4]=r_hp_mul(BI(i_l+1,i_l+2),BI(i_l+1,i_l+2));
            [~,~,r5]=r_hp_add(r3,r4);
            [~,~,r6]=r_hp_sub(r2,r5);
            [~,~,d]=r_hp_div(r6,2);
%             d=((BI(i_l+2,i_l+2)^2+BI(i_l+2,i_l+3)^2)-(BI(i_l+1,i_l+1)^2+BI(i_l+1,i_l+2)^2))/2;
            [~,~,r7]=r_hp_mul(BI(i_l+1,i_l+1),BI(i_l+1,i_l+1));
            [~,~,r8]=r_hp_mul(BI(i_l+1,i_l+2),BI(i_l+1,i_l+2));
            [~,~,r9]=r_hp_add(r7,r8);
            [~,~,r10]=r_hp_add(r9,d);
            [~,~,d2]=r_hp_mul(d,d);
            [~,~,r11]=r_hp_mul(BI(i_l+2,i_l+2),BI(i_l+2,i_l+2));
            [~,~,r12]=r_hp_mul(BI(i_l+1,i_l+2),BI(i_l+1,i_l+2));
            [~,~,r13]=r_hp_mul(r11,r12);
            [~,~,r14]=r_hp_add(d2,r13);
            [~,r15]=r_hp_sqrt(r14);
            [~,~,r16]=r_hp_mul(sign(d),r15);
            [~,~,shift]=r_hp_sub(r10,r16);
%             shift=(BI(i_l+1,i_l+1)^2+BI(i_l+1,i_l+2)^2)+d-sign(d)*sqrt(d^2+BI(i_l+2,i_l+2)^2*BI(i_l+1,i_l+2)^2);
        end
        [~,~,shift2]=r_hp_mul(shift,shift);
        [~,~,s2]=r_hp_mul(s,s);
        [~,~,shiftdivs2]=r_hp_div(shift2,s2);
        if(shiftdivs2)<=epcl
            shift=0;
        end
    end
    if shift~=0
        if direction==1
            [BI,U,V]=QR_Wilkinson_shift_Iteration_once(BI,U,V,i_l,i_u,shift);
        else
            [BI,U]=QR_Wilkinson_shift_Iteration_once_upward(BI,U,i_l,i_u,shift);
        end
    else
        if direction==2
            [BI,U,V]=QR_zero_shift_once_iteration(BI,U,V,i_l,i_u);
        else
            [BI,U]=QR_zero_shift_once_iteration_upward(BI,U,i_l,i_u);
        end
    end
end
V=V';
end

function [H,beta,v,nx]=newht(x)
% 对向量进行Householder变换，返回Householder矩阵H，beta、v、nx为构建H的参数


% 计算向量范数
[n,~]=size(x);
[~,~,xtx]=i_hp_matrix_mul(complex(x'),complex(x));
[~,nx]=i_hp_sqrt(xtx);


y=zeros(n,1);
y(1)=nx;
% 计算v=x-y;
[~,~,v]=i_hp_matrix_sub(complex(x),complex(y));

% 计算beta=v'*x;
[~,~,beta]=i_hp_matrix_mul(complex(v'),complex(x));

% 计算beta=1/beta;
[~,~,beta]=i_hp_div(1,beta);
I=eye(n);

% 计算H=(I-beta*v*v');
[~,~,vvt]=i_hp_matrix_mul(complex(v),complex(v'));
[~,~,betavvt]=i_hp_matrix_mul_e(complex(vvt),complex(beta));
[~,~,H]=i_hp_matrix_sub(complex(I),complex(betavvt));

end

function [G,U1,U2,V1,V2]=mul_22_submatrix(G,U1,U2,V1,V2)
G=real(G);
g11=G(1,1);
g12=G(1,2);
g22=G(2,2);
[~,~,a]=r_hp_mul(g11,g11);
% a=g11^2;
[~,~,g122]=r_hp_mul(g12,g12);
[~,~,g222]=r_hp_mul(g22,g22);
[~,~,b]=r_hp_add(g122,g222);
% b=g12^2+g22^2;
[~,~,c]=r_hp_mul(g11,g12);
% c=g11*g12;
[~,~,c2]=r_hp_mul(2,c);
if abs(c2)>1e-3
    [~,~,bsa]=r_hp_sub(b,a);
    [~,~,theta]=r_hp_div(bsa,c2);
    % theta=(b-a)/(2*c);
    [~,~,theta2]=r_hp_mul(theta,theta);
    [~,~,result]=r_hp_add(1,theta2);
    [~,result1]=r_hp_sqrt(result);
    [~,~,fenmu]=r_hp_add(abs(theta),result1);
    [~,~,t]=r_hp_div(sign(theta),fenmu);
else
    t=1e-3;
end
% t=sign(theta)/(abs(theta)+sqrt(1+theta^2));
[~,~,t2]=r_hp_mul(t,t);
[~,~,result2]=r_hp_add(1,t2);
[~,fenmu]=r_hp_sqrt(result2);
[~,~,cs]=r_hp_div(1,fenmu);
% cs=1/sqrt(1+t^2);
[~,~,sn]=r_hp_mul(cs,t);
% sn=cs*t;
[~,~,G]=r_hp_matrix_mul(G,[cs,sn;-sn,cs]);
% G=G*[cs,sn;-sn,cs];
[V1,V2]=updatecsvv(cs,sn,V1,V2);
[~,~,G112]=r_hp_mul(G(1,1),G(1,1));
[~,~,G212]=r_hp_mul(G(2,1),G(2,1));
[~,~,result3]=r_hp_add(G112,G212);
[~,alpha]=r_hp_sqrt(result3);
% alpha=sqrt(G(1,1)^2+G(2,1)^2);
[~,~,G122]=r_hp_mul(G(1,2),G(1,2));
[~,~,G222]=r_hp_mul(G(2,2),G(2,2));
[~,~,result4]=r_hp_add(G122,G222);
[~,beta]=r_hp_sqrt(result4);
% beta=sqrt(G(1,2)^2+G(2,2)^2);
[~,~,c1]=r_hp_div(G(1,1),alpha);
% c1=G(1,1)/alpha;
[~,~,c2]=r_hp_div(G(2,2),beta);
% c2=G(2,2)/beta;
[~,~,s1]=r_hp_div(G(2,1),alpha);
% s1=G(2,1)/alpha;
[~,~,s2]=r_hp_div(G(1,2),beta);
% s2=G(1,2)/beta;
[~,~,G]=r_hp_matrix_mul([c1,s1;s2,c2],G);
% G=[c1,s1;s2,c2]*G;
[U1,U2]=updatecsvv_last(c1,s1,c2,s2,U1,U2);
end

function [tri,U,V]=QR_Wilkinson_shift_Iteration_once(tri,U,V,i_min,i_max,shift)
n=i_max+1;
[~,~,r0]=r_hp_mul(tri(i_min+1,i_min+1),tri(i_min+1,i_min+1));
[~,~,x]=r_hp_sub(r0,shift);
% x=tri(i_min+1,i_min+1)*tri(i_min+1,i_min+1)-shift;
[~,~,y]=r_hp_mul(tri(i_min+1,i_min+1),tri(i_min+1,i_min+2));
% y=tri(i_min+1,i_min+1)*tri(i_min+1,i_min+2);
for k=i_min:n-2
    [c,s,r]=csr(x,y);
    [~,~,A]=r_hp_matrix_mul([tri(k+1,k+1),tri(k+1,k+2);0,tri(k+2,k+2)],[c,s;-s,c]);
    x=A(1,1);tri(k+1,k+2)=A(1,2);y=A(2,1);tri(k+2,k+2)=A(2,2);
    [V(:,k+1),V(:,k+2)]=updatecsvv(c,s,V(:,k+1),V(:,k+2));
    if(k>0)
        tri(k,k+1)=r;
    end
    [c,s,r]=csr(x,y);
    tri(k+1,k+1)=r;
    [U(:,k+1),U(:,k+2)]=updatecsvv(c,s,U(:,k+1),U(:,k+2));
    if k~=n-2
        [~,~,A]=r_hp_matrix_mul([c,-s;s,c],[tri(k+1,k+2),0;tri(k+2,k+2),tri(k+2,k+3)]);
        x=A(1,1);y=A(1,2);tri(k+2,k+2)=A(2,1);tri(k+2,k+3)=A(2,2);
    else
        [~,~,A]=r_hp_matrix_mul([c,-s;s,c],[tri(n-1,n);tri(n,n)]);
        tri(n-1,n)=A(1);tri(n,n)=A(2);
    end
end
end

function [c,s,rr]=csr(x,y)
if(y==0)
    c=1;
    s=0;
    rr=x;
else
    if(abs(y)>abs(x))
        [~,~,tao]=r_hp_div(-x,y);
%         tao=-x/y;
%         s=sqrt(1+tao^2);
        [~,~,tao2]=r_hp_mul(tao,tao);
        [~,~,r0]=r_hp_add(1,tao2);
        [~,s]=r_hp_sqrt(r0);
        [~,~,rr]=r_hp_mul(-y,s);
%         rr=-y*s;
        [~,~,s]=r_hp_div(1,s);
%         s=1/s;
        [~,~,c]=r_hp_mul(s,tao);
%         c=s*tao;
    else
        [~,~,tao]=r_hp_div(-y,x);
%         tao=-y/x;
        [~,~,tao2]=r_hp_mul(tao,tao);
        [~,~,r0]=r_hp_add(1,tao2);
        [~,c]=r_hp_sqrt(r0);
%         c=sqrt(1+tao^2);
        [~,~,rr]=r_hp_mul(x,c);
%         rr=x*c;
        [~,~,c]=r_hp_div(1,c);
%         c=1/c;
        [~,~,s]=r_hp_mul(c,tao);
%         s=c*tao;
    end
end
end

function [tri,U]=QR_zero_shift_once_iteration_upward(tri,U,i_min,i_max)
n=i_max+1;
oldc=1;
x=tri(n,n);
y=tri(n-1,n);
for i=n-1:-1:i_min+1
    [c,s,r]=csr(x,y);
    [U(:,i+1),U(:,i)]=updatecsvv(c,s,U(:,i+1),U(:,i));
    if i~=n-1
        [~,~,tri(i+1,i+2)]=r_hp_mul(-olds,r);
%         tri(i+1,i+2)=-olds*r;
    end
    [~,~,x]=r_hp_mul(oldc,r);
%     x=oldc*r;
    [~,~,y]=r_hp_mul(-tri(i,i),s);
%     y=-tri(i,i)*s;
    [~,~,h]=r_hp_mul(tri(i,i),c);
%     h=tri(i,i)*c;
    [c,s,r]=csr(x,y);
    oldc=c;
    olds=s;
    tri(i+1,i+1)=r;
    x=h;
    if(i~=1)
        y=tri(i-1,i);
    end
end
[~,~,tri(i_min+1,i_min+2)]=r_hp_mul(-h,s);
% tri(i_min+1,i_min+2)=-h*s;
[~,~,tri(i_min+1,i_min+1)]=r_hp_mul(h,c);
% tri(i_min+1,i_min+1)=h*c;
end

function [tri,U]=QR_Wilkinson_shift_Iteration_once_upward(tri,U,i_min,i_max,shift)
n=i_max+1;
[~,~,ct]=r_hp_sub(i_min,shift);
% ct=i_min-shift;
[~,~,x]=r_hp_mul(tri(n,n),tri(n,n));
% x=tri(n,n)^2;
[~,~,y]=r_hp_mul(tri(n-1,n-1),tri(n-1,n));
% y=tri(n-1,n-1)*tri(n-1,n);
for k=n-1:-1:i_min+1
    [c,s,r]=csr(x,y);
    if k<n-1
        tri(k+1,k+2)=r;
    end
    [~,~,A]=r_hp_matrix_mul([c,-s;s,c],[tri(k+1,k+1),0;tri(k,k+1),tri(k,k)]);
    x=A(1,1);y=A(1,2);tri(k,k+1)=A(2,1);tri(k,k)=A(2,2);
    [U(:,k+1),U(:,k)]=updatecsvv(c,s,U(:,k+1),U(:,k));
    [c,s,r]=csr(x,y);
    tri(k+1,k+1)=r;
    if k~=ct+1
        [~,~,A]=r_hp_matrix_mul([tri(k,k+1),tri(k,k);0,tri(k-1,k)],[c,s;-s,c]);
        x=A(1,1);tri(k,k)=A(1,2);y=A(2,1);tri(k-1,k)=A(2,2);
    else
        [~,~,A]=r_hp_matrix_mul([tri(ct+1,ct+2),tri(ct+1,ct+1)],[c,s;-s,c]);
        tri(ct+1,ct+2)=A(1);tri(ct+1,ct+1)=A(2);
    end
end
end

function [tri,U,VtT]=QR_zero_shift_once_iteration(tri,U,Vt,i_min,i_max)
n=i_max+1;
oldc=1;
x=tri(i_min+1,i_min+1);
y=tri(i_min+1,i_min+2);
for i=i_min:n-1
    [c,s,r]=csr(x,y);
    VtT=Vt';
    [VtT(:,i+1),VtT(:,i+2)]=updatecsvv(c,s,VtT(:,i+1),VtT(:,i+2));
    if i~=i_min
        [~,~,tri(i,i+1)]=r_hp_mul(-olds,r);
%         tri(i,i+1)=-olds*r;
    end
    [~,~,x]=r_hp_mul(oldc,r);
%     x=oldc*r;
    [~,~,y]=r_hp_mul(-tri(i+2,i+2),s);
%     y=-tri(i+2,i+2)*s;
    [~,~,h]=r_hp_mul(tri(i+2,i+2),c);
%     h=tri(i+2,i+2)*c;
    [c,s,r]=csr(x,y);
    [U(:,i+1),U(:,i+2)]=updatecsvv(c,s,U(:,i+1),U(:,i+2));
    tri(i+1,i+1)=r;
    x=h;
    if i~=n-2
        y=tri(i+2,i+3);
    end
    oldc=c;
    olds=s;
end
[~,~,tri(n-1,n)]=r_hp_mul(-h,s);
% tri(n-1,n)=-h*s;
[~,~,tri(n,n)]=r_hp_mul(h,c);
tri(n,n)=h*c;
end

function [v1,v2]=updatecsvv(c,s,v1,v2)
[n,~]=size(v1);
for i=1:n
    t=v1(i);
    [~,~,ct]=i_hp_mul(c,t);
    [~,~,sv2i]=i_hp_mul(s,v2(i));
    [~,~,v1(i)]=i_hp_sub(ct,sv2i);
%     v1(i)=c*t-s*v2(i);
    [~,~,st]=i_hp_mul(s,t);
    [~,~,cv2i]=i_hp_mul(c,v2(i));
    [~,~,v2(i)]=i_hp_add(st,cv2i);
%     v2(i)=s*t+c*v2(i);
end
end

function [v1,v2]=updatecsvv_last(c1,s1,c2,s2,v1,v2)
[n,~]=size(v1);
for i=1:n
    t=v1(i);
    v1(i)=c1*t+s1*v2(i);
    v2(i)=s2*t+c2*v2(i);
end
end