
% Uncertainty criteria
% set Coefficient of variation COV = std/frequency < 0.015
n=0.015;
COV = std_f(:,1)./frequency(:,1); 
% COV = std_d(:,1); 
I = find(COV<n);
[m,~] = size(frequency);
AB = (1:m)';
Im = setdiff(AB,I);  % false mode 

% True mode
COV_true = COV(I,:);
NNN_true=NNN(I,1);
frequency_true = frequency(I,:);
damping_true = damping(I,:);
modeshape_true = MMM(:,I);

std_f_true=std_FFF(I,1);
std_d_true=std_DDD(I,1);
std_m_true=std_MMM(:,I);

frequency1_true = frequency1(I,:);
frequency2_true = frequency2(I,:);

% pole matrix after hard criteria
F=zeros(length(n_mim:step:n_m),n_m/2);   % frequency poles matrix after hard criteria
D=zeros(length(n_mim:step:n_m),n_m/2); % damping poles matrix after hard criteria
M=zeros(channel,n_m/2,length(n_mim:step:n_m)); % mode shape pole matrix
std_F=zeros(length(n_mim:step:n_m),n_m/2);    % frequency std pole matrix after hard criteria
std_D=zeros(length(n_mim:step:n_m),n_m/2);     % damping std pole matrix after hard criteria
std_M=zeros(channel,n_m/2,length(n_mim:step:n_m)); % mode shape std pole matrix

k=0;
for i=n_mim:step:n_m
     k=k+1;
[m,n]=find(NNN_true==i);
F(k,1:length(m))=frequency_true(m,1);
D(k,1:length(m))=damping_true(m,1);
M(1:channel,1:length(m),k)=modeshape_true(:,m);
std_F(k,1:length(m))=std_f_true(m,1);
std_D(k,1:length(m))=std_d_true(m,1);
std_M(1:channel,1:length(m),k)=std_m_true(:,m);
end
% save('modalproperty','F','D','M','std_F','std_D','std_M','sw1','sw2','sw3','fs');

figure(6)
[haxes hline1 hline2]=plotyy(F,n_mim:step:n_m,(0:length(sw1)-1)/nfft*fs,SW,'plot','semilogy');
set(hline1,'LineStyle','none');
set(hline1,'Marker','o','markersize',7,'Linewidth',1);
set(hline1,'Color','r');
set(hline2,'Linewidth',1);

hold on 
for i = 1:size(frequency1_true,1)
    plot([frequency1_true(i,1) frequency2_true(i,1)], [frequency_true(i,2) frequency_true(i,2)], 'color','b')
    plot([frequency1_true(i,1) frequency1_true(i,1)], [frequency_true(i,2)+0.7 frequency_true(i,2)-0.7],  'color','b')
    plot([frequency2_true(i,1) frequency2_true(i,1)], [frequency_true(i,2)+0.7 frequency_true(i,2)-0.7],  'color','b')
end

xlim(haxes(1),[0.1 13]);%fsamp/2
xlim(haxes(2),[0.1 13]);%fsamp/2
ylim(haxes(1),[n_mim n_m]);
ylim(haxes(2),[aa(1,3) aa(1,4)]);
xlabel('Frequency (Hz)');
ylabel(haxes(1),'Model order');
ylabel(haxes(2),'SV');
set(haxes(1),'YTick',(0:step*5:n_m));
set(gca,'FontSize',20);

 
