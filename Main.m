% Authors: Rajarshi Banerjee, Rohit Kundu, Diego Oliva, Ram Sarkar.

% Developed in MATLAB R2018b

clear;
warning('off','all');
close all;
clc;

fprintf('\nProcess HHO_Altruism starts...\n')
Alg_name = 'HHO+altruism'; %Algorithm name
num_runs= 2; %Number of runs
pop=50; %Population 


%% File for results
date = datestr(datetime('now','Format','dd-MM-yyyy'));
namexc = strcat(Alg_name,'_Mth_',date(1:11),'.xlsx');
titles = {'BstFit','Bstsol','PSNR','SSIM','FSIM','UIQI','QILV','HPSI','Time'};
%xlswrite(namexc,titles,1,'A1')
fprintf("%s", namexc);
%show figs  
flagx = 1;

%% Gray-level values
MinValue = 1;
MaxValue = 256;

% Thresholds levels
num_thr = [2,3,4,5];
ThNum = numel(num_thr);

%% Magnetic Resonance image matrix 
load('Harvard_MRI.mat')
names = {'h22','h42','h62','h82'};
%names={'h22'};

MR_Images = numel(names); 

for I = 1:MR_Images
    nut=cell2str(names(1,I));
    nut = nut(3:5);
    fprintf(' \n\n----- Image = %s -----\n', nut)
    fil = 1;
    dataBest = cell(1*ThNum, 17);
    Im = eval(cell2mat(names(1, I)));
    
    for thr = num_thr
        T=num2str(thr);
        fprintf(' \nThreshold = %s', T)
        levels = thr;
        itermax = 1000;
        TempVSolu = 0;

        %% Result Vector
        BestFitStore = []; % Best Fitness vector
        BestSolStore = []; % Best Solution vector
        PSNRvect = [];
        FSIMvect = [];
        SSIMvect = [];
        UIQIvect = [];
        QILVvect = [];
        HPSIvect = [];
        TimeVect = [];
    
        for k=1:num_runs
            %=HHO(pop, itermax, lb, ub, levels+1, Im)
            %[BestFitt, temmmmp, BestSolu, metrics, time, Iout] = lshade(pop, levels+1, itermax, MinValue, MaxValue, Im);
            %return
            %pop+2
            [BestFitt, BestSolu, metrics, time, Iout] = HHO(pop, itermax, MinValue, MaxValue,levels+1, Im);
            %function [bsf_fit_var, converH, bsf_solution, metrics, time, Iout] = lshade(SearchAgents_no, problem_size=dims, Max_iter, lb, ub, Im)
            %          Rabbit_Energy,NULL, Rabbit_Locn, metrics, time, Io 
            
            if isinf(BestFitt)
                continue        
            end

            BestFitStore = [BestFitStore, BestFitt];%concatenation of vector by column
            BestSolu = BestSolu(1:levels);%strip data for CE detail
            BestSolStore = [BestSolStore, mat2str(BestSolu)];

            PSNRvect = [PSNRvect, metrics(1)];  
            SSIMvect = [SSIMvect, metrics(2)];
            FSIMvect = [FSIMvect, metrics(3)];
            UIQIvect = [UIQIvect, metrics(4)];  
            QILVvect = [QILVvect, metrics(5)];
            HPSIvect = [HPSIvect, metrics(6)];

            TimeVect = [TimeVect, time];
            TempVSolu = TempVSolu + BestSolu;
            
        end  % runs   
        
        %% Show Visual Results
        if flagx==1
            L = strcat('Image=',nut,' Thr=',T);
            Axes = axes(figure);
            %Image = imshow(Iout, 'Parent', Axes);
            title(Axes, L);
        end
        
        %% Average solution
        TempVSolu=TempVSolu/num_runs;
        [BestFitnexx, posx]=min(BestFitStore);

        dataBest{fil, 1} = BestFitnexx;
        dataBest{fil, 2} = Strxtract(BestSolStore,posx);
        dataBest{fil, 3} = mean(PSNRvect);
        mean(PSNRvect)
        dataBest{fil, 4} = mean(SSIMvect);
        mean(SSIMvect)
        dataBest{fil, 5} = mean(FSIMvect);
        dataBest{fil, 6} = mean(UIQIvect);    
        dataBest{fil, 7} = mean(QILVvect);   
        dataBest{fil, 8} = mean(HPSIvect);
        dataBest{fil, 9} = mean(TimeVect);

    %---------------------------------------------------*

        dataBest{fil, 10} = BestFitStore;
        dataBest{fil, 11} = PSNRvect;
        dataBest{fil, 12} = SSIMvect;   
        dataBest{fil, 13} = FSIMvect;
        dataBest{fil, 14} = UIQIvect;    
        dataBest{fil, 15} = QILVvect;  
        dataBest{fil, 16} = HPSIvect;  
        dataBest{fil, 17} = TimeVect; 

        fil = fil+1; 

    end% th
    
    %% Save data in matrix
    aname2=strcat(Alg_name,'_',nut,'_MTH',date(1:11),'.mat');
    save(aname2,'dataBest');
    
    %% Save data in file
    sheet = 1; 
    xlRange = strcat('A',num2str(I*(ThNum+2)));
    %dlmwrite(namexc,dataBest,sheet,xlRange)
    %dlmwrite(namexc,dataBest,sheet,xlRange)
    %xlswrite(dataBest, namexc,'Sheet',sheet,'Range',xlRange);
end % images
%writecell
fprintf('\n\nProcess HHO ends...\n')
