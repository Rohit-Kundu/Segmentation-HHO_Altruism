% Main paper:
% Segmentation of Brain MRI using an Altruistic Harris Hawks' Optimization algorithm
% Rajarshi Bandopadhyay, Rohit Kundu, Diego Oliva, Ram Sarkar
% _____________________________________________________

function imgOut=imageGRAY(img,Rvec)
% imgOut=img;
limites=[0 Rvec 255];
tamanho=size(img);
imgOut(:,:)=img*0;
cores=[ 0   0   0;
        255 0   0;
        0   255 0;
        0   0   255;
        255 255 0;
        0   255 255;
        255 0   255;
        255 255 255];
        
%cores=colormap(lines)*255;
%close figure(1);
%close all;
%tic
k=1;
    for i= 1:tamanho(1,1)
        for j=1:tamanho(1,2)
            while(k<size(limites,2))
                if(img(i,j)>=limites(1,k) && img(i,j)<=limites(1,k+1))
                    imgOut(i,j,1)=limites(1,k);
%                     imgOut(i,j,2)=cores(k,2);
%                     imgOut(i,j,3)=cores(k,3);
                end
                k=k+1;
            end
            k=1;
        end
    end
