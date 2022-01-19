%% Variational Bayesian for Gaussian Mixture Model
% concentration parameter is very important!
K = 6; % identified No. cluster by DP
inter = 13;
nfft=3000;  % FFT length
load modalproperty

[m,n] = size(F); [r,~,~] = size(M);
F = reshape(F,[m*n,1]);D = reshape(D,[m*n,1]);
std_F = reshape(std_F,[m*n,1]);std_D = reshape(std_D,[m*n,1]);

Mo = permute(M,[1 3 2]); std_M = permute(std_M,[1 3 2]);
Moo = reshape(Mo,[r,m*n]); std_M = reshape(std_M,[r,m*n]); 

X = [F D]'; 
a = find(X(1,:)>0 & X(1,:)<inter); % remove 0 elements
X = X(:,a); 
[dim N] = size(X);
d = dim;
STD = [std_F std_D]';
STD = STD(:,a);
Moo = Moo(:,a); std_M = std_M(:,a);

% VB fitting
[y1, model, L] = mixGaussVb(X,K); % m is mean; v is No. in each cluster


%% find weight w
model.weight = zeros(1, size(model.R, 2));
for i=1:size(model.R, 2)
    model.weight(1, i) = sum(model.R(:,i) / size(model.R, 1)); % average the model. R
end
rule = find(model.weight(1, :)>(1/3*m/N));
w = model.weight(1, rule); % screen out cluster with less than 1/3

%% find mean 
FD = model.m(:,rule); % mean of frequency and damping 

uni_labels = unique(y1);
No = size(uni_labels, 2);
S_FD = zeros(2,No); % standard derivation
for i=1:No
    a = find(y1==i);
    S_FD(1,i) = mean(STD(1,a),2);
    S_FD(2,i) = mean(STD(2,a),2);
end
S_FD = S_FD(:,rule);

%% find mean of mode shape
Mm = zeros(r,No,2);
std_Mm = zeros(r,No,2);
for j=1:2
    if j==2
        for i=1:No
             a = find(y1==i); 
             Mm(:,i,j)=mean(abs(imag(Moo(:,a))),2);
             std_Mm(:,i,j)=mean(abs(imag(std_M(:,a))),2);
         end
     else
        for i=1:No
            a = find(y1==i); 
             Mm(:,i,j)=mean(abs(real(Moo(:,a))),2);
             std_Mm(:,i,j)=mean(abs(real(std_M(:,a))),2);
        end 
     end
end
% assign +/-
Ms=zeros(r,No,2);
std_Ms=zeros(r,No,2);
for i=1:No
    a = find(y1==i); 
    Ms(:,i,1)=sign(real(Moo(:,a(end))));
    Ms(:,i,2)=sign(imag(Moo(:,a(end))));
    std_Ms(:,i,1)=sign(real(std_M(:,a(end))));
    std_Ms(:,i,2)=sign(imag(std_M(:,a(end))));
end
Mm=Mm.*Ms;
Mm=complex(Mm(:,:,1),Mm(:,:,2));
Mm=Mm(:,rule);
std_Mm=std_Mm.*std_Ms;
std_Mm=complex(std_Mm(:,:,1),std_Mm(:,:,2));
std_Mm=std_Mm(:,rule);

%% Mac value
Mmac=zeros(1,No);

for i=1:No
    a = find(y1==i);
    pe = Moo(:,a);
    aa = length(a);
    mac = zeros(aa-1,1);
    for j=1:aa-1;
        mac(j)=abs((pe(:,j)')*pe(:,j+1))^2/(((pe(:,j)')*pe(:,j))*((pe(:,j+1)')*pe(:,j+1)));
    end
     Mmac(1,i)=mean(mac);
end
Mmac = Mmac(:,rule);

%% find covariance
for i=1:No
    constant(i) = (model.kappa(i)+1)./(model.kappa(i)*(model.v(i)-d+1));
    model.T(:,:,i) = constant(i)*(model.U(:,:,i)'*model.U(:,:,i)); % covariance 
end
model.T = model.T(:,:,rule);
%% SUMMARY
[ave,h] = sort(FD(1,:));
cov = model.T(:,:,h);
S_FD = S_FD(:,h);
result = [FD(:,h); S_FD; Mmac(:,h); w(:,h)] % frequency,damping,std_freq,std_damp,MAC,weighting
Mm=Mm(:,h);
std_Mm=std_Mm(:,h);
save('modeshape','Mm');

%% plot error ellipse
figure(6);
plotClass(X,y1);
ylim([-0.5,3]);
xlim([0 inter]);
% set(gca,'ytick',[-1:0.5:2]);
hold on
for i = 1:K
%     MyEllipse(model.T(:,:,i), FD(:,i),'conf', 0.9,'style','r','intensity',w(i), 'facefill',0.8)
    plot(FD(1,i), FD(2,i),'xk','linewidth',2,'markersize',10); 
end

%% 3D Gaussian distribution
figure(7)
[X1,X2] = meshgrid(linspace(-1,inter,500)',linspace(-1,inter,500)');
XX = [X1(:) X2(:)];
for i = 1:K
    p = mvnpdf(XX,FD(:,i)',model.T(:,:,i));
    Z = reshape(p,500,500);
    mesh(X1,X2,Z); % or use surf
    colormap ( gray ); %  colormap( gray ) 
    grid on; % 
    hold on
    ylim([-1,3]);
    xlim([0,inter]);
    
%     ylabel('Damping ratio (%)');
%     xlabel('Frequency');
%     zlabel('Probability');
end
tx=xlabel('Frequency (Hz)');
ty=ylabel('Damping ratio (%)');
zlabel('Probability density');
set(tx,'Rotation',15)
set(ty,'Rotation',-45)
view(gca,[-30 50]);
% view(-20,50)
rotate3d on
set(gca,'fontsize',20);
%% lpot final
        figure(8)
        SWW = zeros(length(sw1),1,3); 
        SWW(:,:,1)=sw1; SWW(:,:,2)=sw2; SWW(:,:,3)=sw3; 
        color = 'krgmcy';  % cm = jet(size(uni_labels, 2));  set(gca,'colororder',cm)
        m = length(color);
        for i=1:3
            semilogy((0:length(sw1)-1)/nfft*fs,SWW(:,:,i),color(mod(i-1,m)+1),'Linewidth',1.5);
            hold on
        end 
        hold on
        xlim([0 13]);
        ylim([10^(-6) 10^(-4)]);
        range = get(gca,'ylim');
%         for i = 1:size(result,2)
%             plot([result(1,i),result(1,i)],range,'color', [0.7,0.7,0.7],'Linewidth',1.5)
%         end
        fac=0.1;
        for i=1:K
            plot([result(1,i)-fac result(1,i)-fac],range,'color', [0.8,0.8,0.8],'Linewidth',0.5)
            plot([result(1,i)+fac result(1,i)+fac],range,'color', [0.8,0.8,0.8],'Linewidth',0.5)  
            gray = [0,0.,0.];
            h1 = fill([result(1,i)-fac result(1,i)-fac result(1,i)+fac result(1,i)+fac], [range fliplr(range)],gray,...
                'FaceColor',[0.3,0.3,0.3],...
                'EdgeColor','none','FaceAlpha',0.5); % FaceAlpha control transparency 
            hold on
        end
        hold on
  
        xlabel('Frequency (Hz)');
        ylabel('Singular Values');
        xlim([0 13]);xticks([0:3:13]);
        set(gca,'FontSize',22);

