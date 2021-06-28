% Main paper:
% Segmentation of Brain MRI using an Altruistic Harris Hawks' Optimization algorithm
% Rajarshi Bandopadhyay, Rohit Kundu, Diego Oliva, Ram Sarkar
% _____________________________________________________

function PSNRV = PSNR(origImg, distImg)

origImg = double(origImg);
distImg = double(distImg);

[M N] = size(origImg);
error = origImg - distImg;
MSE = sum(sum(error .* error)) / (M * N);

if(MSE > 0)
    PSNRV = 10*log(255*255/MSE) / log(10);
else
    PSNRV = 99;
end
