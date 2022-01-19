function [Modmean] = mean_modal_vector(M1,z,ifr,channel)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Mm=zeros(channel,2);
Ms=zeros(channel,2);
Mm(:,1)=mean(abs(real(M1(:,z,ifr))),3);
Mm(:,2)=mean(abs(imag(M1(:,z,ifr))),3);
Ms(:,1)=sign(real(M1(:,z,ifr(end))));
Ms(:,2)=sign(imag(M1(:,z,ifr(end))));
Mm=Mm.*Ms;
Modmean=complex(Mm(:,1),Mm(:,2));
end
