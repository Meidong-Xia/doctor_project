function [U,A,V]=lansvd(mat)
[size_matm,size_matn]=size(mat);
if size_matm>size_matn
    islong=1;
else
    islong=0;
end
if islong
    mat=mat';
end
[U,A,V]=newbilan(mat);
T=U;A=A';U=V;V=T;
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
V=Vt';
%规范相位
for i=1:min(lenu,lenv)
    if(real(U(1,i))<0)
        U(:,i)=-U(:,i);
        V(:,i)=-V(:,i);
    end

end
end


function [U,B,V]=newbilan(mat)
% mat=mat';
[m,n]=size(mat);
if m>n
    mat=mat';
    isT=1;
else
    isT=0;
end
[m,n]=size(mat);
mini=min(m,n);
p0=zeros(m,1);
p0(1)=1;
beta=zeros(mini+1);
alpha=zeros(mini);
U=zeros(m,mini+1);
V=zeros(n,mini);
beta(1)=1;
U(:,1)=p0;
B=zeros(mini+1,mini);
r=zeros(n,mini+1);
P=zeros(m,mini);
for j=1:mini
    if j~=1
        r(:,j)=mat'*U(:,j)-beta(j)*V(:,j-1);
        for i=1:j-1
            r(:,j)=r(:,j)-(V(:,i)'*r(:,j))*V(:,i);
        end
        alpha(j)=norm(r(:,j));
        V(:,j)=r(:,j)/alpha(j);
        P(:,j)=mat*V(:,j)-alpha(j)*U(:,j);
        for i=1:j
            P(:,j)=P(:,j)-(U(:,i)'*P(:,j))*U(:,i);
        end
        beta(j+1)=norm(P(:,j));
        U(:,j+1)=P(:,j)/beta(j+1);
    else
        r(:,1)=mat'*U(:,1);
        alpha(j)=norm(r(:,1));
        V(:,j)=r(:,1)/alpha(j);
        P(:,j)=mat*V(:,j)-alpha(j)*U(:,j);
        P(:,j)=P(:,j)-(U(:,j)'*P(:,j))*U(:,j);
        beta(j+1)=norm(P(:,j));
        U(:,j+1)=P(:,j)/beta(j+1);
    end
end
for i=1:mini
    B(i,i)=alpha(i);
    B(i+1,i)=beta(i+1);
end
U=U(:,1:mini);
B=B(1:mini,:);
if isT==1
    T=U;
    B=B';
    U=V;
    V=T;
else
end
end

function [U,BI,V]=bi_diag_svd(U,BI,V)
underflow=1e-16;
tol=1e-10;
epcl=1e-12;
iter_num=0;
old_i_low=-1;
old_i_high=-1;
[m,n]=size(BI);
fudge=min(m,n);
maxit=n*n;
[low_bound_sigma,lad,~]=estimation(BI);
max_bound_sigma=max(abs(BI),[],'all');
thresh=max(tol*low_bound_sigma,maxit*underflow);
while true
    iter_num=iter_num+1;
    if iter_num>maxit
        break
    end
    i_u=n-1;
    while((i_u>=1)&&(abs(BI(i_u,i_u+1))<=thresh))
        i_u=i_u-1;
    end
    if i_u==0
        break
    end
    i_l=i_u-1;
    if i_l~=0
        while((abs(BI(i_l,i_l+1))>thresh)&&(i_u>=1))
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
            direction=1;
        else
            direction=2;
        end
    end
    old_i_low=i_l;
    old_i_high=i_u;
    direction=1;

    if direction==1
        if(abs(BI(i_u,i_u+1)/lad(i_u+1))<=tol)
            BI(i_u,i_u+1)=0;
        end
    end
    %compute shift
    if (fudge*tol*low_bound_sigma/max_bound_sigma<=epcl)
        shift=0;
    else
        if direction==1
            s=BI(i_u+1,i_u+1);
            d=((BI(i_u,i_u)*BI(i_u,i_u)+BI(i_u-1,i_u)*BI(i_u-1,i_u))-(BI(i_u+1,i_u+1)*BI(i_u+1,i_u+1)+BI(i_u,i_u+1)*BI(i_u,i_u+1)))/2;
            shift=(BI(i_u+1,i_u+1)*BI(i_u+1,i_u+1)+BI(i_u,i_u+1)*BI(i_u,i_u+1))+d-sign(d)*sqrt(d*d+BI(i_u,i_u)*BI(i_u,i_u)*BI(i_u,i_u+1)*BI(i_u,i_u+1));
        else
            s=BI(i_l+1,i_l+1);
            d=((BI(i_l+2,i_l+2)^2+BI(i_l+2,i_l+3)^2)-(BI(i_l+1,i_l+1)^2+BI(i_l+1,i_l+2)^2))/2;
            shift=(BI(i_l+1,i_l+1)^2+BI(i_l+1,i_l+2)^2)+d-sign(d)*sqrt(d^2+BI(i_l+2,i_l+2)^2*BI(i_l+1,i_l+2)^2);
        end
        if(shift*shift/(s*s))<=epcl
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
        if direction==1
            [BI,U,V]=QR_zero_shift_once_iteration(BI,U,V,i_l,i_u);
        else
            [BI,U]=QR_zero_shift_once_iteration_upward(BI,U,i_l,i_u);
        end
    end
end
V=V';
end

function [v1,v2]=updatecsvv(c,s,v1,v2)
[n,~]=size(v1);
for i=1:n
    t=v1(i);
    v1(i)=c*t-s*v2(i);
    v2(i)=s*t+c*v2(i);
end
end

function [sigma_low_bound,lamda,u]=estimation(B)
[m0,n0]=size(B);
n=min(m0,n0);
lamda=zeros(n,1);
lamda(n)=abs(B(n,n));
for j=n-2:-1:0
    lamda(j+1)=abs(B(j+1,j+1))*(lamda(j+2)/(lamda(j+2)+abs(B(j+1,j+2))));
end
u=zeros(n,1);
u(1)=abs(B(1,1));
for j=0:n-2
    u(j+2)=abs(B(j+2,j+2))*(u(j+1)/(u(j+1)+abs(B(j+1,j+2))));
end
B_Infinity=min(lamda);
B_1=min(u);
sigma_low_bound=min(B_Infinity,B_1);

end
function [G,U1,U2,V1,V2]=mul_22_submatrix(G,U1,U2,V1,V2)
G=real(G);
g11=G(1,1);
g12=G(1,2);
g22=G(2,2);
a=g11^2;
b=g12^2+g22^2;
c=g11*g12;
theta=(b-a)/(2*c);
t=sign(theta)/(abs(theta)+sqrt(1+theta^2));
cs=1/sqrt(1+t^2);
sn=cs*t;
G=G*[cs,sn;-sn,cs];
[V1,V2]=updatecsvv(cs,sn,V1,V2);
alpha=sqrt(G(1,1)^2+G(2,1)^2);
beta=sqrt(G(1,2)^2+G(2,2)^2);
c1=G(1,1)/alpha;
c2=G(2,2)/beta;
s1=G(2,1)/alpha;
s2=G(1,2)/beta;
G=[c1,s1;s2,c2]*G;
[U1,U2]=updatecsvv_last(c1,s1,c2,s2,U1,U2);
end

function [v1,v2]=updatecsvv_last(c1,s1,c2,s2,v1,v2)
[n,~]=size(v1);
for i=1:n
    t=v1(i);
    v1(i)=c1*t+s1*v2(i);
    v2(i)=s2*t+c2*v2(i);
end
end


function [c,s,rr]=csr(x,y)
if(y==0)
    c=1;
    s=0;
    rr=x;
else
    if(abs(y)>abs(x))
        tao=-x/y;
        s=sqrt(1+tao^2);
        rr=-y*s;
        s=1/s;
        c=s*tao;
    else
        tao=-y/x;
        c=sqrt(1+tao^2);
        rr=x*c;
        c=1/c;
        s=c*tao;
    end
end
end


function [tri,U,V]=QR_Wilkinson_shift_Iteration_once(tri,U,V,i_min,i_max,shift)
n=i_max+1;
x=tri(i_min+1,i_min+1)*tri(i_min+1,i_min+1)-shift;
y=tri(i_min+1,i_min+1)*tri(i_min+1,i_min+2);
for k=i_min:n-2
    [c,s,r]=csr(x,y);
    A=[tri(k+1,k+1),tri(k+1,k+2);0,tri(k+2,k+2)]*[c,s;-s,c];
    x=A(1,1);tri(k+1,k+2)=A(1,2);y=A(2,1);tri(k+2,k+2)=A(2,2);
    [V(:,k+1),V(:,k+2)]=updatecsvv(c,s,V(:,k+1),V(:,k+2));
    if(k>0)
        tri(k,k+1)=r;
    end
    [c,s,r]=csr(x,y);
    tri(k+1,k+1)=r;
    [U(:,k+1),U(:,k+2)]=updatecsvv(c,s,U(:,k+1),U(:,k+2));
    if k~=n-2
        A=[c,-s;s,c]*[tri(k+1,k+2),0;tri(k+2,k+2),tri(k+2,k+3)];
        x=A(1,1);y=A(1,2);tri(k+2,k+2)=A(2,1);tri(k+2,k+3)=A(2,2);
    else
        A=[c,-s;s,c]*[tri(n-1,n);tri(n,n)];
        tri(n-1,n)=A(1);tri(n,n)=A(2);
    end
end
end

function [tri,U]=QR_Wilkinson_shift_Iteration_once_upward(tri,U,i_min,i_max,shift)
n=i_max+1;
ct=i_min-shift;
x=tri(n,n)^2;
y=tri(n-1,n-1)*tri(n-1,n);
for k=n-1:-1:i_min
    [c,s,r]=csr(x,y);
    if k<n-1
        tri(k+1,k+2)=r;
    end
    A=[c,-s;s,c]*[tri(k+1,k+1),0;tri(k,k+1),tri(k,k)];
    x=A(1,1);y=A(1,2);tri(k,k+1)=A(2,1);tri(k,k)=A(2,2);
    [U(:,k+1),U(:,k)]=updatecsvv(c,s,U(:,k+1),U(:,k));
    [c,s,r]=csr(x,y);
    tri(k+1,k+1)=r;
    if k~=ct+1
        A=[tri(k,k+1),tri(k,k);0,tri(k-1,k)]*[c,s;-s,c];
        x=A(1,1);tri(k,k)=A(1,2);y=A(2,1);tri(k-1,k)=A(2,2);
    else
        A=[tri(ct+1,ct+2),tri(ct+1,ct+1)]*[c,s;-s,c];
        tri(ct+1,ct+2)=A(1);tri(ct+1,ct+1)=A(2);
    end
end
end


function [tri,U,VtT]=QR_zero_shift_once_iteration(tri,U,Vt,i_min,i_max)
n=i_max+1;
oldc=1;
x=tri(i_min+1,i_min+1);
y=tri(i_min+1,i_min+2);
for i=i_min:n-2
    [c,s,r]=csr(x,y);
    VtT=Vt';
    [VtT(:,i+1),VtT(:,i+2)]=updatecsvv(c,s,VtT(:,i+1),VtT(:,i+2));
    Vt=VtT';
    if i~=i_min
        tri(i,i+1)=-olds*r;
    end
    x=oldc*r;
    y=-tri(i+2,i+2)*s;
    h=tri(i+2,i+2)*c;
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
tri(n-1,n)=-h*s;
tri(n,n)=h*c;
end

function [tri,U]=QR_zero_shift_once_iteration_upward(tri,U,i_min,i_max)
n=i_max+1;
oldc=1;
x=tri(n,n);
y=tri(n-1,n);
for i=n-1:-1:i_min
    [c,s,r]=csr(x,y);
    [U(:,i+1),U(:,i)]=updatecsvv(c,s,U(:,i+1),U(:,i));
    if i~=n-1
        tri(i+1,i+2)=-olds*r;
    end
    x=oldc*r;
    y=-tri(i,i)*s;
    h=tri(i,i)*c;
    [c,s,r]=csr(x,y);
    oldc=c;
    olds=s;
    tri(i+1,i+1)=r;
    x=h;
    if(i~=1)
        y=tri(i-1,i);
    end
end
tri(i_min+1,i_min+2)=-h*s;
tri(i_min+1,i_min+1)=h*c;
end