

%% Automatically find No. Clusters
%% Collapse Gibbs sampling for Dirichelt process gaussian mixture model
% load modal_property
inter = 13; % frequency of interest
[m,n] = size(F);
freq = reshape(F,[m*n,1]);
damp = reshape(D,[m*n,1]);
dpX = [freq damp]';
a = find(dpX(1,:)>0 & dpX(1,:)<inter);
dpX = dpX(:,a);

[dim N] = size(dpX);
d = dim;

[y,model,w,llh,cluster] = mixGaussGb(dpX); % llh log likelihood


figure(5);
plot(cluster,'ro-','linewidth',1.5);
ylim([0 10]);
set(gca, 'YTick', 0:1:20);
set(gca,'fontsize',23);
xlabel('Iterations');
ylabel('Number of clusters');
% nnn = size(cluster,2);
avercluster = mean(cluster(:,end),2)

