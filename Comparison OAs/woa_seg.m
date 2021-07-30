
% The Whale Optimization Algorithm
function [Leader_score,Convergence_curve,Rabbit_Location, metrics, time, Iout]=woa_seg(SearchAgents_no,dim,Max_iter,lb,ub,Im)

 format long;
    format compact; 
    rand('seed', sum(100 * clock));

[ih, ~] = imhist(Im(:,:,1));% histogram  Check Normalizaftion
max_nfes = SearchAgents_no * Max_iter;
[sz1,sz2] = size(Im(:,:,1));% im size
fa=0; % flag



% initialize position vector and score for the leader
Rabbit_Location=zeros(1,dim);
Leader_score=inf; %change this to -inf for maximization problems


%Initialize the positions of search agents
Positions=initializationWOA(SearchAgents_no,dim,ub,lb);

Convergence_curve=zeros(1,Max_iter);

t=0;% Loop counter

% Main loop
while t<Max_iter
    for i=1:size(Positions,1)
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        % Calculate objective function for each search agent
        indivi = sort(fix(Positions(i,:)));
        fitness(i) = CE(indivi,ih);
        fitness = fitness';
%         indivi = sort(fix(Positions(i,:)));
%         fitness = CE(indivi,ih);
%         fitness = fitness';
%         fitness=fobj(Positions(i,:));
        
        % Update the leader
        if fitness(i) < Leader_score % Change this to > for maximization problem
            Leader_score=fitness(i); % Update alpha
            Rabbit_Location=Positions(i,:);
        end
        
    end
    
    a=2-t*((2)/Max_iter); % a decreases linearly fron 2 to 0 in Eq. (2.3)
    
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1+t*((-1)/Max_iter);
    
    % Update the Position of search agents 
    for i=1:size(Positions,1)
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper
        
        
        b=1;               %  parameters in Eq. (2.5)
        l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
        
        p = rand();        % p in Eq. (2.6)
        
        for j=1:size(Positions,2)
            
            if p<0.5   
                if abs(A)>=1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand=abs(C*X_rand(j)-Positions(i,j)); % Eq. (2.7)
                    Positions(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)
                    
                elseif abs(A)<1
                    D_Leader=abs(C*Rabbit_Location(j)-Positions(i,j)); % Eq. (2.1)
                    Positions(i,j)=Rabbit_Location(j)-A*D_Leader;      % Eq. (2.2)
                end
                
            elseif p>=0.5
              
                distance2Leader=abs(Rabbit_Location(j)-Positions(i,j));
                % Eq. (2.5)
                Positions(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Rabbit_Location(j);
                
            end
            
        end
    end
    t=t+1;
    Convergence_curve(t)=Leader_score;
%     [t Leader_score];
end

bsf_solution=Rabbit_Location;
bsf_solution = fix(bsf_solution);
BThresholds = sort(bsf_solution);

gBestR = BThresholds;   
Iout = imageGRAY(Im,gBestR); %Segmented Image

%% Metric calculation

window = 3;
Imd = double(Im);
Ioutd = double(Iout);
psnr = PSNR(Im, Iout); %1
SSIM = ssim(Im, Iout); %2
FSIM = FeatureSIM(Im, Iout); %3
UIQI = img_qi(Im, Iout, window); %4
QILV = qilv(Im, Iout, [window, window]); %5
HPSI = HaarPSI(Imd, Ioutd); %6
metrics = [psnr, SSIM, FSIM, UIQI, QILV, HPSI];

time = toc;
%converH = converH';
bsf_solution = sort(BThresholds);
Rabbit_Location=bsf_solution;

% Leader_score
% 
% bsf_solution=Rabbit_Location;
% bsf_solution = fix(bsf_solution);
%     BThresholds = sort(bsf_solution);
% 
%     gBestR = BThresholds;   
%     Iout = imageGRAY(Im,gBestR); %Segmented Image
%     
%     %% Metric calculation
%     
%     window = 3;
%     Imd = double(Im);
%     Ioutd = double(Iout);
%     psnr = PSNR(Im, Iout); %1
%     SSIM = ssim(Im, Iout); %2
%     FSIM = FeatureSIM(Im, Iout); %3
%     UIQI = img_qi(Im, Iout, window); %4
%     QILV = qilv(Im, Iout, [window, window]); %5
%     HPSI = HaarPSI(Imd, Ioutd); %6
%     metrics = [psnr, SSIM, FSIM, UIQI, QILV, HPSI];
% 
%     time = toc;
% %     converH = converH';
%     bsf_solution = sort(BThresholds);
% bsf_solution = fix(bsf_solution);
% BThresholds = sort(bsf_solution);
% 
% gBestR = BThresholds;   
% Iout = imageGRAY(Im,gBestR); %Segmented Image
% 
% %% Metric calculation
% 
% window = 3;
% Imd = double(Im);
% Ioutd = double(Iout);
% psnr = PSNR(Im, Iout); %1
% SSIM = ssim(Im, Iout); %2
% FSIM = FeatureSIM(Im, Iout); %3
% UIQI = img_qi(Im, Iout, window); %4
% QILV = qilv(Im, Iout, [window, window]); %5
% HPSI = HaarPSI(Imd, Ioutd); %6
% metrics = [psnr, SSIM, FSIM, UIQI, QILV, HPSI];
% 
% time = toc;
% %converH = converH';
% bsf_solution = sort(BThresholds);
% Rabbit_Location=bsf_solution;

end



function Positions=initializationWOA(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end

end

