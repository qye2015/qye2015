%% joint opt with ABS+association, one drop, avg over fading (100)
% z is the off time
%  last update: 3/12/13

%% initialization
clear;clc;close all;
% number of macrocells to be counted
macrocount = 7;

densepico = [7];
densefemt = [32];

for ip=1:length(densepico)
    for ifm=1:length(densefemt)
            nploop = densepico(ip);
            nfloop = densefemt(ifm);


for loop=2:3
bias1 = 5;
bias2 = 10;
bias1b = 5;
bias2b = 10;

%coverage area
L0      = 500;
%number of BSs to be calculated 
RC      = 4;
%density of cellular
lambdaB = 1/(L0^2);
%density of users
lambdaU = 80/(L0^2);
%the radius of hexagonal model
rhex    = sqrt(2/3/sqrt(3)/lambdaB);
L       = L0*RC;

%# of Macros
nm=29;
%# of picos and femtos per macro 
nu=lambdaU*L*L;
np=round(lambdaB*nploop*L*L);
%np=160;
nf=round(lambdaB*nfloop*L*L);
nb=nm+np+nf;

bs1index=1:nm;
bs2index=nm+1:nm+np;
bs3index=nm+np+1:nb;

%% Deployment, BSs, users,
BS      = zeros(nm,2);
BS(2,:) = [0, sqrt(3)*rhex];
BS(3,:) = [3/2*rhex, sqrt(3)/2*rhex]; 
BS(4,:) = [1, -1].*BS(3,:);
BS(5,:) = [0, -sqrt(3)*rhex];
BS(6,:) = [-1, -1].*BS(3,:);
BS(7,:) = [-1, 1].*BS(3,:);
BS(8,:) = [0, 2*sqrt(3)*rhex];
BS(9,:)= [3/2*rhex, sqrt(3)*3/2*rhex]; 
BS(10,:)= [3*rhex, sqrt(3)*rhex]; 
BS(11,:)= [3*rhex, 0]; 
BS(12,:)= [1, -1].*BS(10,:);
BS(13,:)= [1, -1].*BS(9,:);
BS(14,:)= [-1,-1].*BS(8,:);
BS(15,:)= [-1,-1].*BS(9,:);
BS(16,:)= [-1,-1].*BS(10,:);
BS(17,:)= [-1, 1].*BS(11,:);
BS(18,:)= [-1, 1].*BS(10,:);
BS(19,:)= [-1, 1].*BS(9,:);
BS(20,:)= [3/2*rhex, sqrt(3)*5/2*rhex]; 
BS(21,:)= [3*rhex, 2*sqrt(3)*rhex];
BS(22,:)= [1, -1].*BS(20,:); 
BS(23,:)= [1, -1].*BS(21,:); 
BS(24,:)= [-1,-1].*BS(20,:); 
BS(25,:)= [-1,-1].*BS(21,:); 
BS(26,:)= [-1, 1].*BS(20,:); 
BS(27,:)= [-1, 1].*BS(21,:); 
BS(28,:)= [0, sqrt(3)*3*rhex]; 
BS(29,:)= [0, -sqrt(3)*3*rhex]; 

% deployment of users
Ur=rand(nu,2)*L-L/2;
BS(bs2index,:)=rand(length(bs2index),2)*L-L/2;
BS(bs3index,:)=rand(length(bs3index),2)*L-L/2;
% Ur=getCor(sqrt(nu))*L-L/2;
% BS(bs2index,:)=getCor(sqrt(np))*L-L/2;
% BS(bs3index,:)=getCor(sqrt(nf))*L-L/2;


% scatter(BS(bs1index,1),BS(bs1index,2),40,'red');
% hold on;
% scatter(BS(bs2index,1),BS(bs2index,2),30,'black');
% scatter(BS(bs3index,1),BS(bs3index,2),20,'green');
% scatter(Ur(:,1),Ur(:,2),10,'blue');

%% get the uesrs closest to BS 1
distcell= pdist2(BS(1:nm,:),Ur,'euclidean');

nearB   = zeros(nu,1);
for i=1:nu
    tempin = find(distcell(:,i)==min(distcell(:,i)));
    nearB(i)=tempin(1);
end


% only consider the users in the center cell
tempu = find (nearB>=1 & nearB <= macrocount);
%tempu = find (nearB <= 2 | nearB==4 | nearB ==6);
U = Ur(tempu,:);
ncu = length(tempu);

%% get the small cells closest to BS 1
distbs= pdist2(BS(1:nm,:),BS(nm+1:nb,:),'euclidean');

nearbs   = zeros(nb-nm,1);
for i=1:nb-nm
    tempinb = find(distbs(:,i)==min(distbs(:,i)));
    nearbs(i)=tempinb(1);
end

% only consider the users in the center cell
temppf = find (nearbs >= 1 & nearbs <=macrocount);
%temppf = find (nearbs <= 2 | nearbs==4 | nearbs ==6);
smallBS = BS(nm+temppf,:);
nsmall = length(temppf);

npico = length(find(temppf<=np));

%% paramenters to get SINR
%path loss factor
alpha=4; 
%with temperature 290K and BW=10MHz,the thermal noise
noise=4*1.38*10^(-23)*290*10*10^(6);
%power of each BSs.
index=zeros(1,nb);
index(bs1index)=1;index(bs2index)=2;index(bs3index)=3;
%power of BSs. Here assumes there are 3 kinds of BSs, in dBm
p=[40 1 0.1];
power=p(index);
%power of BSs. at ABSs
p2=[0 1 0.1];
%power of each BSs.
power2=p2(index);

%The distance
D = pdist2(BS,U,'euclidean');


%% run 100 times to average performance
avgfad = 5;
Rc = zeros(avgfad,ncu);
km = zeros(avgfad, nb);
kmb = km;
ko = km;
kob = km;
kb = km;
kbb = km;
Rm = Rc;
Rb = Rc;
Ral = Rc;
sumlog = zeros(1,avgfad);
zl = sumlog;
sumlogm=sumlog;
sumlogb=sumlog;
sumlogal = sumlog;
cvxbin = [1:macrocount temppf'+nm];
ncb=length(cvxbin);
%find according bias factor
indexbias = zeros(1,ncb);
indexbias(1:macrocount) = 1;
indexbias(macrocount+1:macrocount+npico) = 2;
indexbias(macrocount+npico+1:ncb) = 3;
bias=[1 bias1 bias2];
biasfac = bias(indexbias)'*ones(1,ncu);
biasb=[1 bias1 bias2];
biasfacb = biasb(indexbias)'*ones(1,ncu);
kal = zeros(avgfad, ncb);

for l=1:avgfad
    
    %channel gain
    hcell  = exprnd(1,nb,ncu);
   % hcell  = ones(nb,ncu);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %==========================================================================
    % at normal subframes
    % get the SINR
    [SINR, S, I]=getsinr(power', D, alpha, noise, hcell);
    cn=log2(1+SINR);
    c1 = cn(cvxbin,:);
    %==========================================================================
    % BS off
    %get the SINR
    [SINR2, S2, I2]=getsinr(power2', D, alpha, noise, hcell);
    cb=log2(1+SINR2);
    c2 = cb(cvxbin,:);

    % the total rate
    %c = z*c1 + (1-z)*c2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Solve the joint user association problem, which is convex
    %optimization by tool cvx
    A=ones(1,ncb);
    P=ones(ncu,1);
    %fprintf('minQ=%e\n',min(min(Q)));
    %c2(c2<10^(-200))=10^(-200);
    cvx_expert true
    cvx_clear
    cvx_begin
        variable x(ncb,ncu);
        variable y(ncb,ncu);
        variable z(1,1);
        maximize( sum(log(sum(c1.*x+c2.*y,1))));
        %minimize( -sum(entr(x*P))-sum(diag(x*Q')) );%mark
        %minimize( -sum(entr(x*P))-sum(sum(x.*Q)));
        %maximize(sum(log(A'*(x.*C))));
    subject to
        x*P<=z*A';
        y*P<=(1-z)*A';
        x>=0;
        x<=1;
        y>=0;
        y<=1;
        z>=0;
        z<=1;
    cvx_end;

    x2=x;
    x2(find(x<0.00001))=0;
    x2(find(x>0.99999))=1;

    y2=y;
    y2(find(y<0.00001))=0;
    y2(find(y>0.99999))=1;
    
    % find the load of each BSs
    for i = 1:ncb
        ko(l,i) = length(find(x2(i,:) >0));
    end
    % find the load of each BSs
    for i = 1:ncb
        kob(l,i) = length(find(y2(i,:) >0));
    end

    %the rate for each user
    Rc(l,:)=sum(c1.*x2+c2.*y2,1)';

    %fprintf('rounding sum rate = %f \n',sum(Rr));
    sumlog(l)=sum(log(Rc(l,:)));%mark   
    
    zl(l)=z;
    z
 
    
    %% on bs 
    % distributed algorithm
    ut=zeros(1,ncb);
    ut1=zeros(1,ncb);
    j=0;
    xl=zeros(ncb,ncu);
    K=ones(1,ncb);
    step=0.0001;
    Rl=zeros(1,ncu);
    sumlogl=10^100;
    metric=10^1000;
    best=-10^1000;
    bestxl=zeros(ncb,ncu);
    bestK=ones(1,ncb);
    bestut=zeros(1,ncb);
    cnt=0;
    R=c1;

    tic;
    while step > 10^-8
        ut=ut1;
        for i=1:ncu
            [maxj j]=max(log(R(:,i))'-ut);
            xl(:,i)=0;
            xl(j(1),i)=1;
        end
        for m=1:ncb
            K(m)=min(nu,exp(-1+ut(m)));
            ut1(m)=ut(m)+step*(sum(xl(m,:))-K(m));
        end
        kl=(sum(xl,2))';
        [nkl xln]=find(xl==1);
        nkl=nkl';
        Rl=(diag(R(nkl,:)))'./kl(nkl);
        sumlogl=sum(log10(Rl));
        cur=max(abs(K-sum(xl,2)'));
        if metric > cur
            metric = cur;
            step = step * 1.01;
        else
            step = step / 1.01;
        end
        if best < sumlogl
            best = sumlogl;
            bestxl = xl;
            bestK = K;
            bestut = ut;
        end
        cnt = cnt + 1;
    %fprintf('cnt=%d best=%f cur=%f sumlogl=%f step=%f\n', cnt, best, cur, sumlogl, step);
    end
    dual_time=toc;

    xl=bestxl;
    K=bestK;
    ut=bestut;
    kl=(sum(xl,2))';
    kal(l,:)=kl;
    [nkl xln]=find(xl==1);
    nkl=nkl';

    %fprintf('best=%f\n', metric);
    Rl=(diag(R(nkl,:)))'./kl(nkl);
    Ral(l,:)=Rl;
    %fprintf('sum rate of algorithm = %f\n',sum(Rl));
    sumlogal(l)=sum(log(Rl));%mark    
    fprintf('sumlog of algorithm= %f\n', sum(log(Rl)));
    %fprintf('dual_time=%f\n', dual_time);


    %% In max-received signal (Max-SINR)
    %find associated BSs, the closest BS
    % normal subframes ====================================================
    xm   = zeros(1,ncu);
    for i = 1:ncu
        tempin = find(c1(:,i) == max(c1(:,i)));
        xm(i) = tempin(1);
        %nearB(tempin(1),i)=1;
    end

    % find the load of each BSs
 
    for i = 1:nb
        km(l,i) = length(find(xm == i));
    end

    % ABSs =====================================================
    xmb   = zeros(1,ncu);
    for i = 1:ncu
        tempinb = find(c2(:,i) == max(c2(:,i)));
        xmb(i) = tempinb(1);
        %nearB(tempin(1),i)=1;
    end

    % find the load of each BSs

    for i = 1:nb
        kmb(l,i) = length(find(xmb == i));
    end

    % record values
    Cm1 = diag(c1(xm,:));
    Cm2 = diag(c2(xmb,:));
    Rm(l,:)=z*Cm1./km(l,xm)'+(1-z)*Cm2./kmb(l,xmb)';
    sumlogm(l)=sum(log(Rm(l,:)));
    
    %% In max biased received signal (Max-bias)
    %get sinr
    %==========================================================================
    % at normal subframes
    % get the biased SINR
    sbn=SINR(cvxbin,:).*biasfac;
    %==========================================================================
    % BS off
    %get the biased SINR
    sbb=SINR2(cvxbin,:).*biasfacb;
    
    %find associated BSs, the closest BS
    % normal subframes ====================================================
    xb   = zeros(1,ncu);
    for i = 1:ncu
        tempinb = find(sbn(:,i) == max(sbn(:,i)));
        xb(i) = tempinb(1);
        %nearB(tempin(1),i)=1;
    end

    % find the load of each BSs
 
    for i = 1:nb
        kb(l,i) = length(find(xb == i));
    end

    % ABSs =====================================================
    xbb   = zeros(1,ncu);
    for i = 1:ncu
        tempinbb = find(sbb(:,i) == max(sbb(:,i)));
        xbb(i) = tempinbb(1);
        %nearB(tempin(1),i)=1;
    end

    % find the load of each BSs

    for i = 1:nb
        kbb(l,i) = length(find(xbb == i));
    end

    % record values
    Cb1 = diag(c1(xb,:));
    Cb2 = diag(c2(xbb,:));
    Rb(l,:)=z*Cb1./kb(l,xb)'+(1-z)*Cb2./kbb(l,xbb)';
    sumlogb(l)=sum(log(Rb(l,:)));

    l
end
save(strcat(int2str(nploop),'_',int2str(nfloop),'_',int2str(loop),'.mat'),'Rm','Rc','Rb','Ral','km','kmb','ko','kob','kb','kbb','kal','zl');
end

    end
end

%save('data1.mat',Rm,Rc,Rb,km,kmb,ko,kob,kb,kbb,zd)
% 
% %% figure
% % association indicator at normal subframes: ain
% % association indicator at ABSs: aib
% ain = zeros(nb,nu);
% aib = ain;
% for i=1:nb
%     jn = find(x2(i,:)>0);
%     ain(i,jn) = 1;
%     jb = find(y2(i,:)>0);
%     aib(i,jb) = 1;    
% end
% 
% aib(1:nm,:) = 0;
% num_fracn = length(find(sum(ain,1)>1));
% num_fracb = length(find(sum(aib,1)>1));
% 
% num_samenb = length(find(sum(ain-aib,1)==0));
% 
% 
% scatter(BS(bs1index,1),BS(bs1index,2),40,'red');
% hold on;
% scatter(BS(bs2index,1),BS(bs2index,2),30,'black');
% scatter(BS(bs3index,1),BS(bs3index,2),20,'green');
% scatter(Ur(:,1),Ur(:,2),10,'blue');
% for i=1:nu
%     nbni = find(ain(:,i)==1);
%     for j=1:length(nbni)
%         plot([BS(nbni(j),1),Ur(i,1)],[BS(nbni(j),2),Ur(i,2)],'r');
%     end
% end
% 
% 
% for i=1:nu
%     nbbi = find(aib(:,i)==1);
%     for j=1:length(nbbi)
%         plot([BS(nbbi(j),1),Ur(i,1)],[BS(nbbi(j),2),Ur(i,2)],'r--');
%     end
% end
% 
% 
% sameuser = find(sum(ain-aib,1)==0);
% for i=1:num_samenb
%     nbbi = find(aib(:,sameuser(i))==1);
%     for j=1:length(nbbi)
%         plot([BS(nbbi(j),1),Ur(sameuser(i),1)],[BS(nbbi(j),2),Ur(sameuser(i),2)],'g--');
%     end
% end

