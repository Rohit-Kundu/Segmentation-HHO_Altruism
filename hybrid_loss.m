% Main paper:
% Segmentation of Brain MRI using an Altruistic Harris Hawks' Optimization algorithm
% Rajarshi Bandopadhyay, Rohit Kundu, Diego Oliva, Ram Sarkar
% _____________________________________________________

function [fit] = hybrid_loss(x,h)
    fit = 0.5*cross_entropy(x,h)+0.5*New_loss(x,h);
end
