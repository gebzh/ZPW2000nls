%% simulation of Fig 5.13
addpath(strcat(pwd,'/methods/'));
addpath(strcat(pwd,'/supplementaryFunctions/'));

global Timer;
global f0_std;
global f1_std;
TIMER_ALL = 0;

% loading data
if isempty(f0_std) || isempty(f1_std) 
    f0_std = [1701.4; 1698.7; 2001.4; 1998.7; 2301.4; 2298.7; 2601.4; 2598.7];
    f1_std = (10.3:1.1:29)';
end

fs = 8000;
Td = 0.7;
Aadjs = (0.4:0.1:0.9);

errorcount = zeros(6,2);
errorcount_nosa = zeros(6,2);
notok = zeros(8,18,8,18,6);

i = 0;
for p1 = 1:8    
    for q1 = 1:18
        for p2 = 1:8
            for q2 = 1:18
                i = 0;
                for Aadj = Aadjs
                    i = i + 1;
                    x = Generate2000Signal(f0_std(p1),f1_std(q1),fs,Td) +...
                        Aadj*Generate2000Signal(f0_std(p2),f1_std(q2),fs,Td);
                    [f0_hat,f1_hat] = NLSM_Based_Algorithm_plus(x,fs);
                    if (abs(f0_hat-f0_std(p1))>0.5 || abs(f1_hat-f1_std(q1))>0.5)
                        fprintf('p1 = %d;q1 = %d;p2 = %d;q2 = %d;Aadj = %f;Ê±NLSM³öÏÖÒëÂë´íÎó£¡\n',...
                            p1,q1,p2,q2,Aadj);
                        errorcount(i,1) = errorcount(i,1) + 1;
                        if floor((p1-1)/2)~=floor((p2-1)/2)
                            errorcount_nosa(i,1) = errorcount_nosa(i,1) + 1;
                        end
                        notok(p1,q1,p2,q2,i) = 1;
%                         breakpointForSectionImplementation( x )
                    end
                    TIMER_ALL = TIMER_ALL + Timer;
                    [f0_hat,f1_hat] = YangFan2010(x,fs);
                    if (abs(f0_hat-f0_std(p1))>0.5 || abs(f1_hat-f1_std(q1))>0.5)
                        fprintf('p1 = %d;q1 = %d;p2 = %d;q2 = %d;Aadj = %f;Ê±YF³öÏÖÒëÂë´íÎó£¡\n',...
                            p1,q1,p2,q2,Aadj);
                        errorcount(i,2) = errorcount(i,2) + 1;
                        if floor((p1-1)/2)~=floor((p2-1)/2)
                            errorcount_nosa(i,2) = errorcount_nosa(i,2) + 1;
                        end
                    end
                    TIMER_ALL = TIMER_ALL + Timer;
                end   
            end
        end
    end
end
disp(TIMER_ALL);
msgbox('Operation Completed');
save(strcat(pwd,'/data_files/Fig 5_13 simulation results.mat'),...
'Aadjs','notok','errorcount','errorcount_nosa');
% costs about  seconds in total
%% simulation of Fig 5.13
addpath(strcat(pwd,'/methods/'));
addpath(strcat(pwd,'/supplementaryFunctions/'));

global Timer;
global f0_std;
global f1_std;
TIMER_ALL = 0;

% loading data
if isempty(f0_std) || isempty(f1_std) 
    f0_std = [1701.4; 1698.7; 2001.4; 1998.7; 2301.4; 2298.7; 2601.4; 2598.7];
    f1_std = (10.3:1.1:29)';
end

fs = 8000;
Td = 0.7;
Aadjs = (0.4:0.1:0.9);

errorcount = zeros(6,2);
errorcount_nosa = zeros(6,2);
notok = zeros(8,18,8,18,6);

i = 0;
for p1 = 1:8    
    for q1 = 1:18
        for p2 = 1:8
            for q2 = 1:18
                i = 0;
                for Aadj = Aadjs
                    i = i + 1;
                    x = Generate2000Signal(f0_std(p1),f1_std(q1),fs,Td) +...
                        Aadj*Generate2000Signal(f0_std(p2),f1_std(q2),fs,Td);
                    [f0_hat,f1_hat] = Elias_Aboutanios2004(x,fs);
                    if (abs(f0_hat-f0_std(p1))>0.5 || abs(f1_hat-f1_std(q1))>0.5)
                        fprintf('p1 = %d;q1 = %d;p2 = %d;q2 = %d;Aadj = %f;Ê±NLSM³öÏÖÒëÂë´íÎó£¡\n',...
                            p1,q1,p2,q2,Aadj);
                        errorcount(i,1) = errorcount(i,1) + 1;
                        if floor((p1-1)/2)~=floor((p2-1)/2)
                            errorcount_nosa(i,1) = errorcount_nosa(i,1) + 1;
                        end
                        notok(p1,q1,p2,q2,i) = 1;
%                         breakpointForSectionImplementation( x )
                    end
                    TIMER_ALL = TIMER_ALL + Timer;
%                     [f0_hat,f1_hat] = Rapid_RELAX_real_Dichotomous_Search(x,fs);
%                     if (abs(f0_hat-f0_std(p1))>0.5 || abs(f1_hat-f1_std(q1))>0.5)
%                         fprintf('p1 = %d;q1 = %d;p2 = %d;q2 = %d;Aadj = %f;Ê±YF³öÏÖÒëÂë´íÎó£¡\n',...
%                             p1,q1,p2,q2,Aadj);
%                         errorcount(i,2) = errorcount(i,2) + 1;
%                         if floor((p1-1)/2)~=floor((p2-1)/2)
%                             errorcount_nosa(i,2) = errorcount_nosa(i,2) + 1;
%                         end
%                     end
%                     TIMER_ALL = TIMER_ALL + Timer;
                end   
            end
        end
    end
end
disp(TIMER_ALL);
msgbox('Operation Completed');
save(strcat(pwd,'/data_files/Fig 5_13 simulation results part2.mat'),...
'Aadjs','notok','errorcount','errorcount_nosa');
% costs about  seconds in total