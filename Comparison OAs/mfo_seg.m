
function [Best_flame_score,Convergence_curve,Best_flame_pos,metrics, time, Iout]=mfo_seg(N,dim,Max_iteration,lb,ub,Im)
format long;
format compact; 
rand('seed', sum(100 * clock));

[ih, ~] = imhist(Im(:,:,1));% histogram  Check Normalizaftion
% max_nfes = SearchAgents_no * Max_iter;
[sz1,sz2] = size(Im(:,:,1));% im size
fa=0; % flag

% display('MFO is optimizing your problem');

%Initialize the positions of moths
Moth_pos=initializationMFO(N,dim,ub,lb);

Convergence_curve=zeros(1,Max_iteration);

Iteration=1;

% Main loop
while Iteration<Max_iteration+1
    
    % Number of flames Eq. (3.14) in the paper
    Flame_no=round(N-Iteration*((N-1)/Max_iteration));
    
    for i=1:size(Moth_pos,1)
        
        % Check if moths go out of the search spaceand bring it back
        Flag4ub=Moth_pos(i,:)>ub;
        Flag4lb=Moth_pos(i,:)<lb;
        Moth_pos(i,:)=(Moth_pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
        
        % Calculate the fitness of moths
        indivi = sort(fix(Moth_pos(i,:)));
        Moth_fitness(1,i) = CE(indivi,ih);
%         Moth_fitness = Moth_fitness';
%         Moth_fitness(1,i)=fobj(Moth_pos(i,:));  
        
    end
       
    if Iteration==1
        % Sort the first population of moths
        [fitness_sorted I]=sort(Moth_fitness);
        sorted_population=Moth_pos(I,:);
        
        % Update the flames
        best_flames=sorted_population;
        best_flame_fitness=fitness_sorted;
    else
        
        % Sort the moths
        double_population=[previous_population;best_flames];
        double_fitness=[previous_fitness best_flame_fitness];
        
        [double_fitness_sorted I]=sort(double_fitness);
        double_sorted_population=double_population(I,:);
        
        fitness_sorted=double_fitness_sorted(1:N);
        sorted_population=double_sorted_population(1:N,:);
        
        % Update the flames
        best_flames=sorted_population;
        best_flame_fitness=fitness_sorted;
    end
    
    % Update the position best flame obtained so far
    Best_flame_score=fitness_sorted(1);
    Best_flame_pos=sorted_population(1,:);
      
    previous_population=Moth_pos;
    previous_fitness=Moth_fitness;
    
    % a linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a=-1+Iteration*((-1)/Max_iteration);
    
    for i=1:size(Moth_pos,1)
        
        for j=1:size(Moth_pos,2)
            if i<=Flame_no % Update the position of the moth with respect to its corresponsing flame
                
                % D in Eq. (3.13)
                distance_to_flame=abs(sorted_population(i,j)-Moth_pos(i,j));
                b=1;
                t=(a-1)*rand+1;
                
                % Eq. (3.12)
                Moth_pos(i,j)=distance_to_flame*exp(b.*t).*cos(t.*2*pi)+sorted_population(i,j);
            end
            
            if i>Flame_no % Upaate the position of the moth with respct to one flame
                
                % Eq. (3.13)
                distance_to_flame=abs(sorted_population(i,j)-Moth_pos(i,j));
                b=1;
                t=(a-1)*rand+1;
                
                % Eq. (3.12)
                Moth_pos(i,j)=distance_to_flame*exp(b.*t).*cos(t.*2*pi)+sorted_population(Flame_no,j);
            end
            
        end
        
    end
    
    Convergence_curve(Iteration)=Best_flame_score;
    
    % Display the iteration and best optimum obtained so far
%     if mod(Iteration,50)==0
%         display(['At iteration ', num2str(Iteration), ' the best fitness is ', num2str(Best_flame_score)]);
%     end
    Iteration=Iteration+1; 
end

bsf_solution=Best_flame_pos;
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
Best_flame_pos=bsf_solution;






function X=initializationMFO(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    X=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        X(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end