
[r c]=size(F);
[rm cm km]=size(M);

mpc=zeros(size(F));
for i=1:r;
    for j=1:c
        if F(i,j)>0;
            Crr=(real(M(:,j,i))'*real(M(:,j,i))); % variance of real part
            Cii=(imag(M(:,j,i))'*imag(M(:,j,i))); % variance of imaginary part
            Cri=(real(M(:,j,i))'*imag(M(:,j,i))); % covariance
            mu=(Cii-Crr)/(2*Cri);
            b=mu+(sign(mu))*sqrt(mu^2+1);
            t=atan(b);
            g1=Crr+Cri*(2*(mu^2+1)*(sin(t))^2-1)/mu;
            g2=Cii-Cri*(2*(mu^2+1)*(sin(t))^2-1)/mu;
            mpc(i,j)=(2*(g1/(g1+g2)-0.5))^2; % modal phase collinearity
        end
    end
end

            
            
mpd=zeros(size(F));

for i=1:r;
    for j=1:c;
        if F(i,j)>0;
            [u s v]=svd([real(M(:,j,i)) imag(M(:,j,i))]);% mean phase deviation
            a=0;
            b=0;
            for k=1:rm;
                a=a+norm(M(k,j,i))*acos(abs((((real(M(k,j,i))))*v(2,2)-(imag(M(k,j,i)))*v(1,2))/(sqrt(v(1,2)^2+v(2,2)^2)*norm(M(k,j,i)))));
                b=b+norm(M(k,j,i));
            end
            mpd(i,j)=radtodeg((a/b));
        end
    end
end

 mq=(1-mpc+(mpd)/45)/2;
 
for i=1:r;
    for j=1:c;
        if mq(i,j)>0.7
%             ((mpd(i,j)/90)>0.3 & mpc(i,j)<0.7); % This is up to the user should be chosen between 0 and 1, see my paper published in 2017. 
            F(i,j)=0;
            D(i,j)=0;
            M(:,j,i)=0;
            std_F(i,j)=0;
            std_D(i,j)=0;
            std_M(:,j,i)=0;
        end
    end
end

% Converting column vector 
w=n_mim-step;
b=0;
for i=1:(n_m-n_mim)/step+1
    w=w+2;
    b=b+1;
    NNN(1:length(nonzeros(F(i,:))),b)=w;
end
NNN=nonzeros(NNN(:));          % model order column vector
FFF=reshape(nonzeros(F'),[],1);      % frequency column vector
DDD=reshape(nonzeros(D'),[],1);      % damping column vector
MMM=reshape(nonzeros(M),channel,[]);  % mode shape matrix
% 
std_FFF=reshape(nonzeros(std_F'),[],1);  % frequency std column vector
std_DDD=reshape(nonzeros(std_D'),[],1);  % damping std column vector
std_MMM=std_M;
std_MMM(:,all(std_MMM==0,1))= []; % remove 0 column vector % mode shape std matrix
 
% 
frequency = [FFF,NNN];damping = [DDD,NNN];
std_f = [std_FFF,NNN];std_d=[std_DDD,NNN];

frequency1 = frequency(:,1)+std_f(:,1);  % upper bound
frequency2 = frequency(:,1)-std_f(:,1);   % lower bound

figure(3)
[haxes hline1 hline2]=plotyy(F,n_mim:step:n_m,(0:length(sw1)-1)/nfft*fs,SW,'plot','semilogy');
set(hline1,'LineStyle','none');
set(hline1,'Marker','o','markersize',7,'Linewidth',1);
set(hline1,'Color','r');
set(hline2,'Linewidth',1);

hold on 
for i = 1:size(frequency1,1)
    plot([frequency1(i,1) frequency2(i,1)], [frequency(i,2) frequency(i,2)],'color','b')
    plot([frequency1(i,1) frequency1(i,1)], [frequency(i,2)+0.7 frequency(i,2)-0.7],'color','b')
    plot([frequency2(i,1) frequency2(i,1)], [frequency(i,2)+0.7 frequency(i,2)-0.7],'color','b')
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




                
