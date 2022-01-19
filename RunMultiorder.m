
% SSI_cov
% multi order computation
ref = 1:12;    % referrence location of  sensors
q = 60;       % SSI number of block row; 
n_m = 130;       % max model order; 
n_mim = 40;     % min model order
n_b = 100;    % number of data block
fs = 100;    % sampling rate
step=2;
nfft=3000;  % FFT length

load shearframe_data
data = tdata;
fs=fs;

[samples channel]=size(data);
Yref=zeros(samples,length(ref));
for i=1:length(ref);
    Yref(:,i)= data(:,ref(1,i));
end
% Overplot of the first two singular values of the PSD matrix:
Gyy = zeros(channel, length(ref), nfft/2+1);
for i=1:channel % response
    for j=1:length(ref) % reference
        Gyy(i,j,:) = cpsd(data(:,i), Yref(:,j), nfft , nfft/2, nfft, fs);
    end
end

sw1 = zeros(size(Gyy,3),1);
sw2 = zeros(size(Gyy,3),1);
sw3 = zeros(size(Gyy,3),1);

for i=1:size(Gyy,3)
    [~, sg, ~] = svd(Gyy(:,:,i));
    sw1(i) = sqrt(sg(1,1));
    sw2(i) = sqrt(sg(2,2));
    sw3(i) = sqrt(sg(3,3));
 
end
SW=[sw1 sw2 sw3 ]; % aggiungi sw4

figure(1)
semilogy((0:length(sw1)-1)/nfft*fs,sw1,'k','Linewidth',1.5);
hold on
semilogy((0:length(sw1)-1)/nfft*fs,sw2,'r','Linewidth',1.5);
hold on
semilogy((0:length(sw1)-1)/nfft*fs,sw3,'g','Linewidth',1.5);
xlabel('Frequency (Hz)');
ylabel('Singular Values');
xlim([0 13]);
ylim([1*10^(-6) 10^(-4)]);
aa=axis;


tic
[modal,phi,std_modal,std_phi]...
    = SSI_cov_uncertainty(n_mim,n_m,q,n_b,data,ref,fs);
toc

F=zeros(length(n_mim:step:n_m),n_m/2);   % frequency poles matrix
D=zeros(length(n_mim:step:n_m),n_m/2); % damping poles matrix
M=zeros(channel,n_m/2,length(n_mim:step:n_m)); % mode shape pole matrix
std_F=zeros(length(n_mim:step:n_m),n_m/2);    % frequency std pole matrix
std_D=zeros(length(n_mim:step:n_m),n_m/2);     % damping std pole matrix
std_M=zeros(channel,n_m/2,length(n_mim:step:n_m)); % mode shape std pole matrix

k=0;
for i=n_mim:step:n_m
     k=k+1;
[m,n]=find(modal==i);
F(k,1:length(m))=modal(m,2);
D(k,1:length(m))=100*modal(m,3);
M(1:channel,1:length(m),k)=phi(:,m);
std_F(k,1:length(m))=std_modal(m,2);
std_D(k,1:length(m))=std_modal(m,3);
std_M(1:channel,1:length(m),k)=std_phi(:,m);
end
%%% Hard criteria remove poles with damping<0 and damping >10
a=D<=0;
b=find(a);
D(b)=0;
F(b)=0;
std_F(b)=0;
std_D(b)=0;
a= D>10;
b=find(a);
D(b)=0;
F(b)=0;
std_F(b)=0;
std_D(b)=0;
a=F>0;

MM=zeros(size(M));
for i=1:channel;
    MM(i,:,:)=a';
end
M=MM.*M;
std_M=MM.*std_M;


% all uncertainty in stabilization diagram after damping criteria
w=n_mim-step;
b=0;
for i=1:(n_m-n_mim)/step+1
    w=w+2;
    b=b+1;
    NN_(1:length(nonzeros(F(i,:))),b)=w;
end
NN_=nonzeros(NN_(:));          % model order column vector
FF=reshape(nonzeros(F'),[],1);      % frequency column vector
std_FF=reshape(nonzeros(std_F'),[],1);  % frequency std column vector
frequency = [FF,NN_];
std_f = [std_FF,NN_];
% frequency = [modal(:,2),modal(:,1)];
% std_f = [std_modal(:,2),modal(:,1)];

frequency1 = frequency(:,1)+std_f(:,1);  % upper bound
frequency2 = frequency(:,1)-std_f(:,1);   % lower bound

figure(2)
[haxes hline1 hline2]=plotyy(F,n_mim:step:n_m,(0:length(sw1)-1)/nfft*fs,SW,'plot','semilogy');
set(hline1,'LineStyle','none');
set(hline1,'Marker','o','markersize',7,'Linewidth',1);
set(hline1,'Color','r');
set(hline1,'MarkerSize',7);
set(hline2,'Linewidth',1);

hold on 
for i = 1:size(frequency1,1)
    plot([frequency1(i,1) frequency2(i,1)], [frequency(i,2) frequency(i,2)],'color','b')
    plot([frequency1(i,1) frequency1(i,1)], [frequency(i,2)+0.7 frequency(i,2)-0.7],'color','b')
    plot([frequency2(i,1) frequency2(i,1)], [frequency(i,2)+0.7 frequency(i,2)-0.7],'color','b')
end
% xlim(haxes(1),[0 fsamp/2]);%fsamp/2
% xlim(haxes(2),[0 fsamp/2]);%fsamp/2
xlim(haxes(1),[0.1 13]);%fsamp/2
xlim(haxes(2),[0.1 13]);%fsamp/2
ylim(haxes(1),[n_mim n_m]);
ylim(haxes(2),[aa(1,3) aa(1,4)]);
xlabel('Frequency (Hz)');
ylabel(haxes(1),'Model order');
ylabel(haxes(2),'SV');
set(haxes(1),'YTick',(0:step*5:n_m));
set(gca,'FontSize',20);

% Open the mode validation criteria tool:
mode_criteria




