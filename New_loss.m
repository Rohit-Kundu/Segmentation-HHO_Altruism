% Main paper:
% Segmentation of Brain MRI using an Altruistic Harris Hawks' Optimization algorithm
% Rajarshi Bandopadhyay, Rohit Kundu, Diego Oliva, Ram Sarkar
% _____________________________________________________

function [fit] = New_loss(x,h)
    
    st = size(x,2);
    for i = 1:st+1
        
        if i == 1  % The 1st Part
            ti = x(i);
	        ti_1 = 1;
            nu(i) = -mEin(ti_1,ti,h)*(1-mEin(ti_1,ti,h)/mZero(ti_1,ti,h)) - log(mEin(ti_1,ti,h))*log((1-mEin(ti_1,ti,h))/mZero(ti_1,ti,h));
            
        
        elseif i> st % Last Part
            ti=256;
	        ti_1=x(i-1);
            nu(i) = -mEin(ti_1,ti,h)*(1-mEin(ti_1,ti,h)/mZero(ti_1,ti,h)) - log(mEin(ti_1,ti,h))*log((1-mEin(ti_1,ti,h))/mZero(ti_1,ti,h));
        
        else % Original
            ti=x(i);
	        ti_1=x(i-1);
            nu(i) = -mEin(ti_1,ti,h)*(1-mEin(ti_1,ti,h)/mZero(ti_1,ti,h)) - log(mEin(ti_1,ti,h))*log((1-mEin(ti_1,ti,h))/mZero(ti_1,ti,h));
        end
    end

    sumNU=sum(nu);

    if (isnan(sumNU)==1)
        fit=-1*sumNU;  
    else
        fit=-1*sumNU;    
    end

end

function [u] = mEin(aa,bb,ih)
    bm_1=bb-1;
    dif_bm_1=abs(aa-bb);
    inten= linspace(aa,bm_1,dif_bm_1);
    inten=inten';
    H=ih(aa:bm_1);
    
%     if size(inten,1)~=size(H,1)
%         size(inten,1)
%         size(H,1)
%         size(inten,1)
%     end
%     if size(inten,2)~=size(H,2)
%         size(inten,2)
%         size(H,2)
%     end
%     
    uu=inten.*H;
    u=sum(uu);
end

function [u]=mZero(aa,bb,ih)
    bm_1=bb-1;
    h=ih(aa:bm_1);
    u=  sum(h);
end
