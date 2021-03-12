function [fit] = hybrid_loss(x,h)
    fit = 0.5*cross_entropy(x,h)+0.5*New_loss(x,h);
end
