function [U,Sigma,V,Plist,stop] =f_svd_precoding(H,P,Ns,flg_b,noise,flg_a,rho)
stop = false;
if flg_b == 1 % 64bit
    [U,Sigma,V]=svd(H);
    V = V(:,1:Ns);
    U = U(:,1:Ns);
    Sigma = Sigma(1:Ns,1:Ns);
    if flg_a == 1 % 平均功率分配
        Plist = sqrt(P/(Ns))*diag(ones(Ns,1));  
    elseif flg_a == 2 % 注水功率分配
        Plist = power_allocation(Sigma,P,noise,flg_b); 
    end
    V = V*Plist;
elseif flg_b == 2 % 32bit
    [U,Sigma,V]=cfp32_svd_household(H);
    acc = norm(H-U*Sigma*V','fro')/norm(H,'fro');
    if acc > 0.1 || isnan(acc)
        stop = true;
    end
    if nnz(Sigma) < Ns
        stop = true;
    end
    if stop
        Plist = [];
        return
    end
    V = V(:,1:Ns);
    U = U(:,1:Ns);
    Sigma = Sigma(1:Ns,1:Ns);
    if flg_a == 1 % 平均功率分配
        Plist = single(sqrt(single(P)/trace(single(V)*single(V'))))*single(diag(ones(Ns,1))); 
    elseif flg_a == 2 % 注水功率分配
        Plist = power_allocation(Sigma,P,noise,flg_b); 
    end
    V = single(V)*single(Plist);
elseif flg_b == 3 % 16bit
    [U,Sigma,V]=cfp16_mysvd_4(H);
    % [U,A,V,err]=cfp16_svd_lanczos(H);
    % err = norm(U*A*V'-H,"fro")/norm(H,"fro");
    acc = norm(H-U*Sigma*V','fro')/norm(H,'fro');
    if acc > 0.1 || isnan(acc)
        stop = true;
    end
    if nnz(Sigma) < Ns
        stop = true;
    end
    if stop
        Plist = [];
        return
    end
    V = V(:,1:Ns);
    U = U(:,1:Ns);
    Sigma = Sigma(1:Ns,1:Ns);
    if flg_a == 1 % 平均功率分配
        [~,~,V_head] = i_hp_matrix_mul(complex(V),complex(V'));
        power = 0;
        for k=1:size(V_head,1)
            [~,~,power] = i_hp_add(power,V_head(k,k));
        end
        [~,~,tmp] = i_hp_div(P,power);
        [~,p] = i_hp_sqrt(tmp);
        Plist = p*diag(ones(Ns,1)); 
    elseif flg_a == 2 % 注水功率分配
        Plist = power_allocation(Sigma,P,noise,flg_b); 
    elseif flg_a == 3 % 数值解法
        Plist = f_power_optimization(P,Sigma,noise,rho,1); 
    elseif flg_a == 4 % 二分查找
        Plist = f_power_optimization(P,Sigma,noise,rho,2); 
    end
    [~,~,V] = i_hp_matrix_mul(complex(V),complex(Plist));
end

end

function P_list=power_allocation(Sigma,P,noise,flg_b)
% Generate the water-filling power allocation matrix whose i th diagonal
% entry represent the power allocated to i th channel, given channel gains,
% total power, and noise power
%
% Based on:
%   https://zhuanlan.zhihu.com/p/502453127
%
% Inputs:
%   - Sigma: the channel gains matrix which usually generated by svd
%   - P: the total power
%   - noise: noise power
%   - flg_b: computation precision
%
% Outputs:
%   - P_list: the power allocation matrix
%
%
% log:
%   - initialized by Meidong Xia on 11/26/2023
%   - Using bisection algorithm by Meidong Xia on 01/21/2024
% 
%
% Rest of the code... 
size_0=size(Sigma,1);
lamda=zeros(1,size_0);
if flg_b == 1 % 64bit
    for i=1:size_0
        lamda(i)=noise/Sigma(i,i)^2;
    end
    w_pre = 0;
    w_next=200*P;
    while true
        w = (w_next+w_pre)/2;
        f = sum(max((w-lamda),0));
        if abs(f-P)<1e-5
            break
        end
        if f>P
            w_next = w; 
        else
            w_pre = w;       
        end
    end
    P_list=sqrt(diag(max((w-lamda),0)));
elseif flg_b==2 % 32bit
    for i=1:size_0
        lamda(i)=single(noise)/single(Sigma(i,i))^2;
    end
    w_pre = 0;
    w_next=single(200)*single(P);
    while true
        w = (single(w_next)+single(w_pre))/single(2);
        f = sum(max(single(w)-single(lamda),single(0)));
        if abs(single(f)-single(P))<1e-5
            break;
        end
        if f>P
            w_next = w; 
        else
            w_pre = w;       
        end
        
    end
    P_list=sqrt(diag(max(single(w)-single(lamda),single(0))));
elseif flg_b==3 % 16bit
    for i=1:size_0
        [~,~,sigma_tmp] = i_hp_mul(Sigma(i,i),Sigma(i,i));
        [~,~,lamda(i)]=i_hp_div(noise,sigma_tmp);
    end
    w_pre = 0;
    [~,~,w_next]=i_hp_mul(200,P);
    
    cnt = 1e2;
    while true
        [~,~,w_tmp] = i_hp_add(w_next,w_pre);
        [~,~,w] = i_hp_div(w_tmp,2);
        f = real(f_cal(w,lamda,0));
        [~,~,judge] = i_hp_sub(f,P);
        if real(abs(judge)) < 1e-5
            break;
        end
        if f>P
            w_next = w; 
        else
            w_pre = w;       
        end
        
        cnt = cnt-1;
        if cnt<0
            break
        end
    end
    P_list=diag(f_cal_tmp(w,lamda));
    for i=1:length(lamda)
        [~,P_list(i,i)] = i_hp_sqrt(P_list(i,i));
    end
end
end

function val = f_cal(w,lamda,P) % sum(max((w-lamda),0))-P
    tmp = f_cal_tmp(w,lamda);
    val = 0;
    for i=1:length(tmp)
        [~,~,val] = i_hp_add(val,tmp(i));
    end
    [~,~,val] = i_hp_sub(val,P);
end

function val = f_cal_tmp(w,lamda) % max((w-lamda),0)
    val = zeros(size(lamda));
    for i = 1:length(val)
        [~,~,tmp] = i_hp_sub(w,lamda(i));
        val(i) = max(real(tmp),0);
    end
end

function [U,A,V]=cfp32_svd_household(mat)
% 1.26更新 规范U、V相位
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
A=real(A);
[U,A,Vt]=bi_diag_svd(U,A,V);
[m0,n0]=size(A);
[lenu,~]=size(U);
[lenv,~]=size(Vt);
n=min(m0,n0);
for i=1:n
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

function [U,BI,V]=bi_diag_svd(U,BI,V)
% underflow=1e-4;
tol=1e-4;
epcl=1e-14;
iter_num=0;
old_i_low=-1;
old_i_high=-1;
[m,n]=size(BI);
fudge=min(m,n);
maxit=3*n*n;
[low_bound_sigma,lad,~]=estimation(BI);
max_bound_sigma=max(abs(BI),[],'all');
count=0;
count1=0;
% thresh=max(tol*low_bound_sigma,maxit*underflow);
while true
    iter_num=iter_num+1;
    if iter_num>maxit
        break
    end
    i_u=n-1;
    while((i_u>=1)&&(abs(BI(i_u,i_u+1))<=10^(-1.5)))
        i_u=i_u-1;
    end
    if i_u==0
        break
    end
    i_l=i_u-1;
    if i_l~=0
        while((abs(BI(i_l,i_l+1))>10^(-1.5))&&(i_u>=1))
            i_l=i_l-1;
            if i_l==0
                break
            end
        end
    end
    if i_u==i_l+1
        [BI(i_l+1:i_u+1,i_l+1:i_u+1),U(:,i_l+1),U(:,i_u+1),V(:,i_l+1),V(:,i_u+1)]=mul_22_submatrix(BI(i_l+1:i_u+1,i_l+1:i_u+1),U(:,i_l+1),U(:,i_u+1),V(:,i_l+1),V(:,i_u+1));
        count1=count1+1;
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
            s=single(BI(i_u+1,i_u+1));
            d=single(((BI(i_u,i_u)*BI(i_u,i_u)+BI(i_u-1,i_u)*BI(i_u-1,i_u))-(BI(i_u+1,i_u+1)*BI(i_u+1,i_u+1)+BI(i_u,i_u+1)*BI(i_u,i_u+1)))/2);
            shift=single((BI(i_u+1,i_u+1)*BI(i_u+1,i_u+1)+BI(i_u,i_u+1)*BI(i_u,i_u+1))+d-sign(d)*sqrt(d*d+BI(i_u,i_u)*BI(i_u,i_u)*BI(i_u,i_u+1)*BI(i_u,i_u+1)));
        else
            s=single(BI(i_l+1,i_l+1));
            d=single(((BI(i_l+2,i_l+2)^2+BI(i_l+2,i_l+3)^2)-(BI(i_l+1,i_l+1)^2+BI(i_l+1,i_l+2)^2))/2);
            shift=single((BI(i_l+1,i_l+1)^2+BI(i_l+1,i_l+2)^2)+d-sign(d)*sqrt(d^2+BI(i_l+2,i_l+2)^2*BI(i_l+1,i_l+2)^2));
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

function [H,beta,v,nx]=newht(x)
% nx=norm(x);
[n,~]=size(x);
xtx=single(x'*x);
nx=sqrt(single(xtx));
y=zeros(n,1);
y(1)=nx;
v=single(x-y);
beta=single(v'*x);
beta=single(1/beta);
I=eye(n);
H=(I-beta*(v*v'));
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
    UK(j:m,j:m)=single(UK(j:m,j:m)*UJH');
    U=single(U*UK);
    A(j:m,j:n)=single(UJH*A(j:m,j:n));
    if j<n
        [VJH,~,~,~]=newht(A(j,j+1:n)');
        VK(j+1:n,j+1:n)=single(VK(j+1:n,j+1:n)*VJH');
        V=single(V*VK);
        A(j:m,j+1:n)=single(A(j:m,j+1:n)*VJH');
    end
end
if ~is
    A=A';
    u=U;
    U=V;
    V=u;
end
end

function [v1,v2]=updatecsvv(c,s,v1,v2)
[n,~]=size(v1);
for i=1:n
    t=single(v1(i));
    v1(i)=single(c*t-s*v2(i));
    v2(i)=single(s*t+c*v2(i));
end
end

function [sigma_low_bound,lamda,u]=estimation(B)
[m0,n0]=size(B);
n=min(m0,n0);
lamda=zeros(n,1);
lamda(n)=abs(B(n,n));
for j=n-2:-1:0
    lamda(j+1)=single(abs(B(j+1,j+1))*(lamda(j+2)/(lamda(j+2)+abs(B(j+1,j+2)))));
end
u=zeros(n,1);
u(1)=single(abs(B(1,1)));
for j=0:n-2
    u(j+2)=single(abs(B(j+2,j+2))*(u(j+1)/(u(j+1)+abs(B(j+1,j+2)))));
end
B_Infinity=min(lamda);
B_1=min(u);
sigma_low_bound=min(B_Infinity,B_1);

end

function [G,U1,U2,V1,V2]=mul_22_submatrix(G,U1,U2,V1,V2)
G=single(real(G));
g11=G(1,1);
g12=G(1,2);
g22=G(2,2);
a=single(g11^2);
b=single(g12^2+g22^2);
c=single(g11*g12);
theta=single((b-a)/(2*c));
t=single(sign(theta)/(abs(theta)+sqrt(1+theta^2)));
cs=single(1/sqrt(1+t^2));
sn=single(cs*t);
G=single(G*[cs,sn;-sn,cs]);
[V1,V2]=updatecsvv(cs,sn,V1,V2);
alpha=single(sqrt(G(1,1)^2+G(2,1)^2));
beta=single(sqrt(G(1,2)^2+G(2,2)^2));
c1=single(G(1,1)/alpha);
c2=single(G(2,2)/beta);
s1=single(G(2,1)/alpha);
s2=single(G(1,2)/beta);
G=single([c1,s1;s2,c2]*G);
[U1,U2]=updatecsvv_last(c1,s1,c2,s2,U1,U2);
end

function [v1,v2]=updatecsvv_last(c1,s1,c2,s2,v1,v2)
[n,~]=size(v1);
for i=1:n
    t=v1(i);
    v1(i)=single(c1*t+s1*v2(i));
    v2(i)=single(s2*t+c2*v2(i));
end
end


function [c,s,rr]=csr(x,y)
if(y==0)
    c=1;
    s=0;
    rr=x;
else
    if(abs(y)>abs(x))
        tao=single(-x/y);
        s=single(sqrt(1+tao^2));
        rr=single(-y*s);
        s=single(1/s);
        c=single(s*tao);
    else
        tao=single(-y/x);
        c=single(sqrt(1+tao^2));
        rr=single(x*c);
        c=single(1/c);
        s=single(c*tao);
    end
end
end


function [tri,U,V]=QR_Wilkinson_shift_Iteration_once(tri,U,V,i_min,i_max,shift)
n=i_max+1;
x=single(tri(i_min+1,i_min+1)*tri(i_min+1,i_min+1)-shift);
y=single(tri(i_min+1,i_min+1)*tri(i_min+1,i_min+2));
for k=i_min:n-2
    [c,s,r]=csr(x,y);
    A=single([tri(k+1,k+1),tri(k+1,k+2);0,tri(k+2,k+2)]*[c,s;-s,c]);
    x=single(A(1,1));tri(k+1,k+2)=single(A(1,2));y=single(A(2,1));tri(k+2,k+2)=single(A(2,2));
    [V(:,k+1),V(:,k+2)]=updatecsvv(c,s,V(:,k+1),V(:,k+2));
    if(k>0)
        tri(k,k+1)=single(r);
    end
    [c,s,r]=csr(x,y);
    tri(k+1,k+1)=single(r);
    [U(:,k+1),U(:,k+2)]=updatecsvv(c,s,U(:,k+1),U(:,k+2));
    if k~=n-2
        A=single([c,-s;s,c]*[tri(k+1,k+2),0;tri(k+2,k+2),tri(k+2,k+3)]);
        x=single(A(1,1));y=single(A(1,2));tri(k+2,k+2)=single(A(2,1));tri(k+2,k+3)=single(A(2,2));
    else
        A=single([c,-s;s,c]*[tri(n-1,n);tri(n,n)]);
        tri(n-1,n)=single(A(1));tri(n,n)=single(A(2));
    end
end
end

function [tri,U]=QR_Wilkinson_shift_Iteration_once_upward(tri,U,i_min,i_max,shift)
n=i_max+1;
ct=single(i_min-shift);
x=single(tri(n,n)^2);
y=single(tri(n-1,n-1)*tri(n-1,n));
for k=n-1:-1:i_min
    [c,s,r]=csr(x,y);
    if k<n-1
        tri(k+1,k+2)=single(r);
    end
    A=single([c,-s;s,c]*[tri(k+1,k+1),0;tri(k,k+1),tri(k,k)]);
    x=single(A(1,1));y=single(A(1,2));tri(k,k+1)=single(A(2,1));tri(k,k)=single(A(2,2));
    [U(:,k+1),U(:,k)]=updatecsvv(c,s,U(:,k+1),U(:,k));
    [c,s,r]=csr(x,y);
    tri(k+1,k+1)=single(r);
    if k~=ct+1
        A=single([tri(k,k+1),tri(k,k);0,tri(k-1,k)]*[c,s;-s,c]);
        x=single(A(1,1));tri(k,k)=single(A(1,2));y=single(A(2,1));tri(k-1,k)=single(A(2,2));
    else
        A=single([tri(ct+1,ct+2),tri(ct+1,ct+1)]*[c,s;-s,c]);
        tri(ct+1,ct+2)=single(A(1));tri(ct+1,ct+1)=single(A(2));
    end
end
end


function [tri,U,VtT]=QR_zero_shift_once_iteration(tri,U,Vt,i_min,i_max)
n=i_max+1;
oldc=1;
x=single(tri(i_min+1,i_min+1));
y=single(tri(i_min+1,i_min+2));
for i=i_min:n-2
    [c,s,r]=csr(x,y);
    VtT=Vt';
    [VtT(:,i+1),VtT(:,i+2)]=updatecsvv(c,s,VtT(:,i+1),VtT(:,i+2));
    Vt=VtT';
    if i~=i_min
        tri(i,i+1)=single(-olds*r);
    end
    x=single(oldc*r);
    y=single(-tri(i+2,i+2)*s);
    h=single(tri(i+2,i+2)*c);
    [c,s,r]=csr(x,y);
    [U(:,i+1),U(:,i+2)]=updatecsvv(c,s,U(:,i+1),U(:,i+2));
    tri(i+1,i+1)=single(r);
    x=h;
    if i~=n-2
        y=single(tri(i+2,i+3));
    end
    oldc=single(c);
    olds=single(s);
end
tri(n-1,n)=single(-h*s);
tri(n,n)=single(h*c);
end

function [tri,U]=QR_zero_shift_once_iteration_upward(tri,U,i_min,i_max)
n=i_max+1;
oldc=1;
x=single(tri(n,n));
y=single(tri(n-1,n));
for i=n-1:-1:i_min
    [c,s,r]=csr(x,y);
    [U(:,i+1),U(:,i)]=updatecsvv(c,s,U(:,i+1),U(:,i));
    if i~=n-1
        tri(i+1,i+2)=single(-olds*r);
    end
    x=single(oldc*r);
    y=single(-tri(i,i)*s);
    h=single(tri(i,i)*c);
    [c,s,r]=single(csr(x,y));
    oldc=single(c);
    olds=single(s);
    tri(i+1,i+1)=single(r);
    x=single(h);
    if(i~=1)
        y=single(tri(i-1,i));
    end
end
tri(i_min+1,i_min+2)=single(-h*s);
tri(i_min+1,i_min+1)=single(h*c);
end

function [U,A,Vt]=cfp16_mysvd_4(mat)
    mat_half=hp_matrix_turn(mat);
    [U_single,A_single,Vt_single]=cfp32_svd_household(mat_half);
    U = hp_matrix_turn(double(U_single));
    A = hp_matrix_turn(double(A_single));
    Vt = hp_matrix_turn(double(Vt_single));
end


function Plist = f_power_optimization(P,Sigma,noise,rho,flg)
% Generate the power allocation matrix whose i th diagonal
% entry represent the power allocated to i th channel, given channel gains,
% total power, and noise power
%
% Based on:
%   
%
% Inputs:
%   - Sigma: the channel gains matrix which usually generated by svd
%   - P: the total power
%   - noise: noise power
%   - rho: the variance of computation error
%   - flg: 1: numerical method, 2: bisection method 
%
% Outputs:
%   - P_list: the power allocation matrix
%
%
% log:
%   - initialized by Meidong Xia on 01/21/2024
% 
%
% Rest of the code...

size_0=size(Sigma,1);
sigma = diag(Sigma); % 信道向量
eposilon = 1e-5;

if flg == 1 % numerical
    
    cvx_begin quiet
        variable p(size_0) nonnegative
        expression f
        for i = 1:size_0
            f = f+P*rho*inv_pos(p(i))+noise*inv_pos(p(i)*sigma(i)^2);
        end
        minimize(f)
        subject to
            sum(p) == P;
    cvx_end
    Plist = diag(sqrt(p));

elseif flg == 2 % bisection 
    mu_min = 0;
    mu_max = 100*P;
    while mu_max - mu_min > 1e-5
        mu = (mu_max+mu_min)/2;
        val = 0;
        for i = 1:size_0
            tmp = sqrt(P*rho*sigma(i)^2+noise);
            val = val + max(tmp/(sigma(i)*sqrt(mu)),eposilon);
        end
        if abs(val-P)<1e-5
            break
        end
        if val>P
            mu_min = mu;
        elseif val<P
            mu_max = mu;
        end

    end
    p = zeros(size_0,1);
    for i = 1:size_0
        tmp = sqrt(P*rho*sigma(i)^2+noise);
        p(i) = max(tmp/(sigma(i)*sqrt(mu)),eposilon);
    end
    Plist = diag(sqrt(p));
end
end