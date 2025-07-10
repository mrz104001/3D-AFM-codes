% %trying to read the plotter saving files from the Zurich Instrument
% %lock-in. it autodetects if you've calculated the data. If so, it will skip
% %the calculation part and do post-processing.
% 
% %the procedure of the program:
% %1. reads the data from .m file saved with the matlab.
% %2. ! need to know the meaning of auxin1 and auxin2
% %3. cut the matrix into forward and backward scanning for x scan.
% %4. use only the forward scan to further get the Z ramp data. cut.
% %5. caculate everything based on the cut data from 4.
% 
% s=getenv('COMPUTERNAME');
% if strcmp(s,'DESKTOP-LSI2BC0')
%     cd 'C:\work\postdoc\WCM\projects\2021.11.1 AqpZ protein and water\matlab code\force_curve';
% else
%     if (strcmp(s,'DESKTOP-5EBBHSB'))
%         cd 'D:\work\matlab coding\force_curve';
%     else
%         cd 'C:\Runze\matlab coding\force_curve';
%     end
% end
% 
% clc;
% close all;
% clearvars -except path;
% 
% %trying to read the reading path from userPath.mat, the var is readingpath.
% %the use of readingpath record the the working folder.
% load('userPath.mat','readingpath');
% 
% if (exist('path','var')==0)||~ischar(path)
%     if (exist('readingpath','var')==1)
%         path=readingpath;
%     else
%         path = 'C:\work\postdoc\WCM\projects\2021.11.1 AqpZ protein and water\*.mat';
%     end
%     clear readingpath;
% end
% 
% 
% %parameter settings.
% f1=1;   %if 1, save the raw data for debug, 0 for no.
% f2=1;   %if 1, average phase and 2D plot it.
% f3=0;   %if 1, force the program to reload and recalculate despite the same file name.
% f4=0;   %defines the scan direction, 0 forward, 1 backward.
% missing_line_flag=0;    % record the status of missing line detection.
% 
% smooth=2;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %define the coefficient for transferring data
% v2x=12.6*2/2.5;     %voltage to nm for Z piezo, in nm/V, *2 because the Z piezo driver calculates (A-B)*2.
% v2A=20;         %voltage to nm for deflection, in nm/V
% xv2x=7.51*5;    %voltage to nm for x piezo, in nm/V, *5 because the x piezo driver calculates *5.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %trying to load the workspace saved for the 3D scanning.
% if (exist('path','var')==0)||~ischar(path)
% path='C:\work\postdoc\WCM\projects\2021.11.1 AqpZ protein and water\12.1 A-phase method\*.mat';
% end
% 
% [name0,path0]=uigetfile(path);
% if name0~=0
% path=[path0,name0];
% end
% 
% disp('Loading data, please wait.');
% 
% load(path);
% 
% readingpath=path;
% save('userPath.mat','readingpath','-append');
% clear readingpath;
% 
% %after loading from the path, need to define the data in auxin1 and auxin2,
% %one is Z ramp, the other is x scan.
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %trying to cut the matrix loaded from the file.
% %based on my current understanding, the forward and backward here means the
% %scan direction in x.
% 
% %the cutting here is already wrong, fix it.
% %just used a simple switch for the two aux channels temporarily.
% %need an elegant way to do it right.
% 
% %use fft to find the higher frequency ramp and assign the low frequency as
% %x ramp and high one as z ramp.
% 
% %this function is not working at all, will try redo the function later.
% %[auxin1,auxin2]=x_z_ramp(auxin1,auxin2,5000,5);
% % 
% % auxin10=auxin1;
% % auxin1=auxin2;
% % auxin2=auxin10;
% 
% disp('Loading data finished, please continue to the next section...');
% fprintf('------------------------------------------------------------------\n');
% fprintf('\n');
% 
% %need to determine the range of the data before cutting it.
% %% try to cut the matrix into useful data before cutting the period.
% 
% % optional, cut the data if strange thing happens.
% if (0)
%     N=37967;
%     auxin1=auxin1(N:end);
%     auxin2=auxin2(N:end);
%     r=r(N:end);
%     phase=phase(N:end);
% end
% 
% % ask if determine the pixel size automatically from the data file.
% answer = questdlg('Automatically determine pixel size?', ...
%         'Auto?', ...
%         'Yes','No','Yes');
% 
% if (isfield(para,'pixel')&&~isempty(para.pixel))&&strcmp(answer,'Yes')
%     % check the para for the pixel value.
%     N=str2double(para.pixel);
% else
%     % ask for input of pixel number.
%     disp("Manually input of pixel size, wait for input");
%     Ans=inputdlg("pixel number", ...
%     "Input",[1 40],"50");
%     N=str2double(Ans);
% end
% 
% % try to locate the period for x scan to locate the scanning data.
% % if x scan is not good, won't affect the z ramp calculation much.
% disp("Finding local max for x scan signal...");
% [~,max1]=period_divide(auxin2,N,0);
% disp("Finding local min for x scan signal...");
% [~,min1]=period_divide(-auxin2(max1(1):max1(end)),N,0);
% min1=min1+max1(1)-1;
% 
% % check if the period is right for the x scan signal.
% if (length(max1)<N||length(min1)<N-1)
%     answer = questdlg('possible missing lines detected, proceed?', ...
% 	'Proceed?', ...
% 	'Yes','No','Yes');
%     if strcmp(answer,'No')
%         return;
%     end
%     %get the start and end point for the x scan cutting.
%     period_max=mean(max1(2:end)-max1(1:end-1));
%     auxin11=auxin1(max(1,ceil((max1(1)-period_max))):min(length(auxin1),ceil(max1(end)+period_max)));
%     missing_line_flag=1;
% else
%     %get the start and end point for the x scan cutting.
%     period_max=mean(max1(2:end)-max1(1:end-1));
%     auxin11=auxin1(max(1,ceil((max1(1)-period_max))):min(length(auxin1),ceil(max1(end)+period_max)));
% end
% 
% % optional, cut the meaningful part of the auxin11.
% % try to plot the auxin1 for the choosing of the cutting point, then ask
% % for input of the useful part.
% if missing_line_flag==1         % check if line missing.
%     answer = questdlg('Possible missing lines, plot auxin1 to choose useful part of the data?', ...
%         'Proceed?', ...
%         'Yes','No','Yes');
%     if strcmp(answer,'Yes')
%         f=figure;
%         N_downsample=100;
%         plot(downsample(auxin1,N_downsample));
%         title('auxin1 for measuring region, the code will resume after the figure is closed.');
%         uiwait(f);
%         prompt = {'Enter region start point:','Enter region end point:'};
%         dlgtitle = 'Useful region';
%         fieldsize = [1 45; 1 45];
%         definput = {'0',num2str(length(auxin1)/N_downsample)};
%         answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
%         N_data_start=str2double(answer{1})*N_downsample+1;
%         N_data_end=str2double(answer{2})*N_downsample;
%         auxin11=auxin11(N_data_start:N_data_end);
% 
%         % need to deal with max1 and min1 to remove the outlier for the
%         % data end.
%         max1(max1>N_data_end)=[];
%         min1(min1>N_data_end)=[];
%     end
% end
% 
% L=length(auxin11);
% auxin111=smoothdata(auxin11,'movmean',ceil(L/N^2/2));
% auxin11=auxin11-auxin111;
% 
% disp("Finding local max for z ramp signal...");
% [min0,max0]=period_divide(auxin11,N^2,1);
% 
% if ceil((max1(1)-period_max))>=1
%     min00=min0+ceil((max1(1)-period_max))-1;
%     max00=max0+ceil((max1(1)-period_max))-1;
% else
%     min00=min0;
%     max00=max0;
% end
% 
% min0=min00;
% max0=max00;
% 
% % as for 4/26, these two flags are just for denote purpose, no actual use.
% missing_data_flag=0;
% missing_line_direction=0;   % 0 means first few lines are missing, 1 means last few lines are missing.
% 
% if length(min0)<2*N^2
%     answer = questdlg('Missing data detected, proceed?', ...
% 	'Proceed?', ...
% 	'Yes','No','Yes');
%     if strcmp(answer,'No')
%         return;
%     end
%     missing_data_flag=1;
%     missing_line_flag=1;
% end
% 
% % done here, every single pixel is determined.
% 
% fprintf('\n');
% disp('Period max and min determined, please move on to the next section.');
% 
% if missing_line_flag==1
%     answer = questdlg('Missing lines, choose the missing direction:', ...
%         'Ask for direction?', ...
%         'First','Last','First');
%     switch answer
%         case 'First'
%             missing_line_direction=0;
%         case 'Last'
%             missing_line_direction=1;
%     end
%     fprintf('\n');
%     disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
%     if missing_line_direction == 0
%         disp('Detected missing lines, data will be aligned to the end!');
%     else
%         disp('Detected missing lines, data will be aligned to the front!');
%     end
%     disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
% end
% 
% fprintf('------------------------------------------------------------------\n');
% fprintf('\n');


%% start part for the retrace calculation, should input the before and after avg mat file from trace calculation.

%% temperorily isolate this part for coding, add the auto fix function.
%use square divide function to cut the matrix into square. need to check
%before next step.
[N_forward,N_backward,missing_line_flag]=square_divide(auxin2,max0,min0,max1,min1,missing_line_flag);

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% N_forward, column 1, start of the forward, column 2, end of the forward,
% column 3, cycle number in the forward.0
%N=size(N_forward,1);

fprintf('\n');
disp("Pre-cutting for matrix done. Please check N_backward for the squareness manually.");

% check the oddness in the N_forward matrix and give a warning.
AA_warning=find(N_backward(:,3)~=N);
%fprintf("The following index might have problem: %s.\n",mat2str(AA_warning'));

% try to write auto fix for obvious division errors.
AA_warning1_fixed=[];
AA_warning2=[];
for i=1:length(AA_warning)
    % go through the problematic points and try to fix it.
    j=AA_warning(i);
    % for the first one, directly cut the size down to N.
    if j==1
        if N_backward(1,3)>N
            N_backward(1,1)=N_backward(1,2)-N+1;
            N_backward(1,3)=N;
            AA_warning1_fixed(end+1)=1;
        else
            AA_warning2(end+1)=1;
        end
        continue;
    end

    % for every other cycle, if the length is N-1, check the previous
    % N_backward, if N+1, fix the cut number.
    if j>1&&(N_backward(j,3)==N-1)&&(N_backward(j-1,3)==N+1)
        N_backward(j,1)=N_backward(j,1)-1;
        N_backward(j,3)=N_backward(j,3)+1;
        N_backward(j-1,2)=N_backward(j-1,2)-1;
        N_backward(j-1,3)=N_backward(j-1,3)-1;
        AA_warning1_fixed(end+1)=j;
        continue;
    end
    AA_warning2(end+1)=j;
end
if (~isempty(AA_warning))
    disp(' ');
    disp('Problem with cutting:');
end
if (~isempty(AA_warning1_fixed))
    fprintf("The following index is auto corrected: %s. Please check again.\n",mat2str(AA_warning1_fixed));
end
if (~isempty(AA_warning2))
    fprintf("The following index needs manual correction: %s.\n",mat2str(AA_warning2));
end
fprintf('------------------------------------------------------------------\n');
fprintf('\n');


%% IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% test region before next step.
% advise to open the N_backward matrix and manually adjust it to square.
close all
% plot function to see cutting situation.
i=47;        % the backward cycle to be checked.
dN=5;       % larger range to show the misalignment of the cutting.
N_backward(i,3)=N_backward(i,2)-N_backward(i,1)+1;
N1=min0(N_backward(i,1)-dN);
N2=max0(min([N_backward(i,2)+dN,length(max0)]));
figure;
plot(auxin1(N1:N2));
hold on;
yyaxis right;
plot(auxin2(N1:N2));
%xline(min0(min0_ind));

% try to plot the number of lines for addition of missing cycles.
for j=N_backward(i,1):N_backward(i,2)
    label=sprintf('%d',j-N_backward(i,1)+1);
    xline(max0(j)-N1,'-',label);
end
hold off;
legend("z ramp","x scan");

%% USE ACCORDINGLY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% this is the function to add certain min and max position into the matrix.
i_add=N_backward(i,1)+49;      % the number here is the line number before the missing line in the plot.
m_base=N1;
min_add=267340+m_base;
max_add=267769+m_base;

% adding the element into the min0 and max0, cautious!!!!!!!!!!!!
min0=[min0(1:i_add-1);min_add;min0(i_add:end)];
max0=[max0(1:i_add-1);max_add;max0(i_add:end)];

% add 1 to everything in N_backward after the added line.
N_backward(i,2)=N_backward(i,2)+1;
N_backward(i,3)=N_backward(i,3)+1;
for j=i+1:length(N_backward(:,1))
    N_backward(j,1)=N_backward(j,1)+1;
    N_backward(j,2)=N_backward(j,2)+1;
end


%% Dealing with the loss of data in the last retrace.
N_backward(end,2)=length(max0);

%% IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Run before next step.
fprintf('\n');
AA_warning=find(N_backward(:,3)~=N);
if (~isempty(AA_warning))
    w_message=sprintf("The following index might have problem: %s. Refer to AA_warning for later use.",mat2str(AA_warning'));
    disp(w_message);
else
    disp("Everything seems OK, please proceed...");
end
fprintf('------------------------------------------------------------------\n');
fprintf('\n');

%% cut the matrix based on the revised cutting information.
% now we have 5 channels to cut. just in case we don't have full pixels,
% give the undetected as null values.

% a pre-fix for the N_backward, align the lines to the end.
lines=size(N_backward,1);
if lines<N_backward(1,3)
    N_backward=[nan(N_backward(1,3)-lines,3);N_backward];
end

fprintf('\n');
auxin2_backward=flip(matrix_cut_step(min0,max0,auxin2,N_backward),2);
auxin1_backward=flip(matrix_cut_step(min0,max0,auxin1,N_backward),2);
r_backward=flip(matrix_cut_step(min0,max0,r,N_backward),2);
phase_backward=flip(matrix_cut_step(min0,max0,phase,N_backward),2);
%z_backward=matrix_cut_step(min0,max0,z,N_backward);

disp('Cutting for x scan, z ramp, r and phase done...');
fprintf('------------------------------------------------------------------\n');
fprintf('\n');



%% averaging section.
% try to figure out the averaging function I programed before and do a bit
% wise averaging for the data.

% test the average function.
% instructions on step_3D_avg.
% three input variables. first one is the input matrix, nxn cell matrix.
% second one is the method option, depends on the method, the third input
% and return value is different.
%%% for method 1, simple cut, define the step resolution and cut based on
% this, the third value is the resolution, and the return value is a struct
% variable, containing cut point, step average, pixel start point and
% total pixel.
%%% for method 2, planning to divide the matrix based solely on the steps,
% see what I do with the single point spectrum for manually saved files.
%%% for method 3, cut the input matrix based on the given matrix. the third
% input would be the position matrix from method 1 or 2, the output would
% be only the cutted matrix.
fprintf('\n');
% ask for input of two variables, the z resolution time, and bin size.
Ans=inputdlg(["z resolution time (1-10)","bin size (1-10)"], ...
    "Input",[1 40],["1" "3"]);
z_resolution_time=round(str2double(Ans{1}));
bin_size=round(str2double(Ans{2}));

%dummy protection. forbid overlimit values.
if (z_resolution_time<1||z_resolution_time>10)
    z_resolution_time=1;
end
if (bin_size<1||bin_size>10)
    bin_size=1;
end

z_resolution=2/4096*z_resolution_time;


% the tested good results are 2-4 bin size for resolution time 1.



disp('Averaging the force curve based on height resolution...');
auxin1_avg=step_3D_average(auxin1_backward,2,z_resolution);

r_avg_bin1=step_3D_average(r_backward,3,auxin1_avg.position);
phase_avg_bin1=step_3D_average(phase_backward,3,auxin1_avg.position);
z_avg_bin1=step_3D_average(auxin1_backward,3,auxin1_avg.position);

% try to bin from the contact point.
r_avg=bin_avg(r_avg_bin1,bin_size);
phase_avg=bin_avg(phase_avg_bin1,bin_size);
z_avg=bin_avg(z_avg_bin1,bin_size);

disp('Averaging done, need to check the averaging quality before next step.');
disp('Try using single_cell_variable_readout app to explore the cut matrix.');
%the averaging seems OK, but still need to program another way of cutting
%the matrix. also, need good data to test the display of the matrix.
fprintf('------------------------------------------------------------------\n');
fprintf('\n');

%% get rid of the z drift or not, currently no.
% try to get all the ideal minimum point for each pixel based on amplitude.
remove_LP_flag=0;

if remove_LP_flag==1
% first, get the min point of amplitude for each pixel.
r_cell=r_avg;   %using r_cell for easy replacement of the processing cell variable.
z_cell=z_avg;

[m,n]=size(r_cell);
r_min=nan(m,n);
for i=1:m
    for j=1:n
        r_min(i,j)=min(r_cell{i,j});
    end
end

r_min_L=mean(r_min,'all');        % get the largest min point of all pixels.

% try to find the pixel point for each pixel that equals the largest min
% point.

% the method is to get the first one that's larger and linear fit to get
% the more accurate index out.
ind_min=nan(m,n);
z_min=nan(m,n);
ind_min0=nan(m,n);
for i=1:m
    for j=1:n
        % ind_min00=find(r_cell{i,j}<=r_min_L,1);
        % if (isempty(ind_min00))
        %     ind_min0(i,j)=length(r_cell{i,j});
        % else
        %     ind_min0(i,j)=ind_min00;
        % end
        % if (ind_min0(i,j)==1)
        %     continue;
        %     % this is causing nan problem.
        % end
        % % do a linear fit of ind_min0 and ind_min0-1 to get the more
        % % accurate index.
        % y=polyfit(r_cell{i,j}(ind_min0(i,j)-1:ind_min0(i,j)),[ind_min0(i,j)-1,ind_min0(i,j)],1);
        % ind_min(i,j)=polyval(y,r_min_L)-length(r_cell{i,j});
        % % get the accurate z_min point for each pixel.
        % z_min(i,j)=interp1(z_cell{i,j},polyval(y,r_min_L));

        % new way of searching the index.
        fprintf('Searching for amplitude align index, %0.1f%% in progress.\n',(j+(i-1)*n)/(m*n)*100);
        try
            ind_min(i,j)=index_search(r_cell{i,j},r_min_L);
        catch
            fprintf('Something is wrong with index i=%d, j=%d.\n',i,j);
        end
    end
end

% this part was used to get the topograph out for the same setpoint.
% poly_num=1;
% 
% % try to do polynomial fit and remove the background.
% z_min_corr=nan(m,n);
% for i=1:m
%     % try to remove the nan value from the matrix and restore it
%     % afterwards.
%     z_line=z_min(i,:);
%     z_line1=z_line(~isnan(z_line));
% 
%     % poly fit to remove the background drift.
%     y=polyfit(1:length(z_line1),z_line1,poly_num);
%     z_line(~isnan(z_line))=polyval(y,1:length(z_line1));
%     z_line_corr1=z_line;
% 
%     % next fit sine wave.
%     z_line=z_min(i,:)-z_line;
%     % solve the nan problem.
%     nan_ind0=find(isnan(z_line));
%     for j=1:length(nan_ind0)
%         nan_ind=nan_ind0(j);
%         if nan_ind-1>0
%             z_line(nan_ind)=z_line(nan_ind-1);
%         else
%             z_line(nan_ind)=z_line(nan_ind+1);
%         end
%     end
%     z_line_corr2=zeros(1,length(z_line));
%     for j=1:2
%         y=sineFit(1:length(z_line),z_line,0);
%         ft=fittype('offs+amp*sin(2*pi*f*x+phi)','independent',{'x'},'coefficients',{'offs','amp','f','phi'});
%         z_line_corr2=z_line_corr2+feval(ft,y(1),y(2),y(3),y(4),1:length(z_line));
%         z_line=z_line-feval(ft,y(1),y(2),y(3),y(4),1:length(z_line));
%     end
% 
%     % try to use sine wave to remove the background drift again.
%     z_min_corr(i,:)=z_line_corr1+z_line_corr2;
% end
% 
% z_min1=z_min-z_min_corr;
% 
% z_min2=z_min1/z_real_resolution;


    ind_shift=ind_min;%-z_min2;
    ind_shift=ind_shift-min(min(ind_shift));


    try
        MIJ.closeAllWindows();
        MIJ.exit();
    catch
    end

    javaaddpath 'C:\Program Files\MATLAB\R2023b\java\mij.jar'
    javaaddpath 'C:\imaging software\ImageJ\ij.jar'

    MIJ.start();
    MIJ.createImage(ind_min);


    disp('The z drift has been calculated, applying to the raw data.');

    % adding nan in the very end of the cell variable for r, phase, z.
    % this is modifying the avg matrix directly.


    for i=1:m
        for j=1:n
            if (~isnan(ind_shift(i,j)))
                r_avg_bin1{i,j}=[r_avg_bin1{i,j},nan(1,round(ind_shift(i,j)*bin_size))];
                phase_avg_bin1{i,j}=[phase_avg_bin1{i,j},nan(1,round(ind_shift(i,j)*bin_size))];
                z_avg_bin1{i,j}=[z_avg_bin1{i,j},nan(1,round(ind_shift(i,j)*bin_size))];
            end
        end
    end
end

r_avg=bin_avg(r_avg_bin1,bin_size);
phase_avg=bin_avg(phase_avg_bin1,bin_size);
z_avg=bin_avg(z_avg_bin1,bin_size);

disp('Recalculation for bin avg done.');
fprintf('------------------------------------------------------------------\n');
fprintf('\n');



%% try to convert the 3d cell variable to a 3d matrix, easier to plot in softwares.
% the format of the axis is y,x, z.

% the flip here is to flip the z direction so that the z is smaller for
% smaller indices.

r_avg_m=flip(cell_2_3dmatrx(r_avg),3);
phase_avg_m=flip(cell_2_3dmatrx(phase_avg),3);
z_avg_m=flip(cell_2_3dmatrx(z_avg),3);

r_avg_bin1_m=flip(cell_2_3dmatrx(r_avg_bin1),3);
phase_avg_bin1_m=flip(cell_2_3dmatrx(phase_avg_bin1),3);
z_avg_bin1_m=flip(cell_2_3dmatrx(z_avg_bin1),3);

disp('Converting the r and phase into 3d matrix, ready to be processed by display_3D_curve_all.mlapp.');
fprintf('------------------------------------------------------------------\n');
fprintf('\n');

%% saving part, temporarily save the processed data to a file.
% save into two parts, before averaging and after averaging.
fprintf('\n');
disp('Saving all the calculated information...');
% before averaging part. all the backward, para and path.
path1=[path(1:end-4),'_processed_before_avg_backward.mat'];
save(path1,"auxin1_backward","auxin2_backward","r_backward","phase_backward","para","path","N_backward","N_backward","missing_line_flag","min0","max0","min1","max1",'-mat');

% after averaging part, all the avg and m data.
path1=[path(1:end-4),'_processed_after_avg_backward.mat'];
save(path1,"auxin1_avg","auxin1_avg","r_avg","r_avg_bin1","phase_avg","phase_avg_bin1","r_avg_m","r_avg_bin1_m","phase_avg_m","phase_avg_bin1_m","z_avg","z_avg_m","z_avg_bin1_m","para","path","z_resolution","bin_size",'-mat');

disp('All data saved for later reference...');

fprintf('------------------------------------------------------------------\n');
fprintf('\n');



%% SECTION II: to get the matrix ready for the IgorPro software.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% after processing, the data will not be raw after this.

%% get rid of the z error effect in the data.
lines=sum(isnan(N_backward),'all')/3;
N=size(N_backward,1);

r_m=r_avg_m(lines+1:N,:,:);
phase_m=phase_avg_m(lines+1:N,:,:);
z_m=z_avg_m(lines+1:N,:,:);

[r_avg_m_corr,phase_avg_m_corr,z_avg_m_corr]=z_corr(r_m,phase_m,z_m);

% save the corrected data into file.
path1=[path(1:end-4),'_processed_after_avg_backward.mat'];
save(path1,"phase_avg_m_corr","r_avg_m_corr","z_avg_m_corr",'-mat','-append');

disp('<a href="matlab:opentoline(''force_curve_3D_5_channels_stepZ_for_missing_lines_retrace.m'',670)">Next code</a>');
beep;

%% matrix processing, assuming the name of the data matrix is given.
% the processed matrix names are phase_avg_m_corr, r_avg_m_corr,
% z_avg_m_corr.


list0={'phase_avg_m_corr','r_avg_m_corr','z_avg_m_corr'};

tic;
for ii=1:length(list0)
    s=eval(list0{ii});

    % determine if the name contains phase.
    tf=strfind(list0{ii},'phase');

    % do a reverse for the phase data.
    if ~isempty(tf)
        s=-s;
    end


    % to record the start and end point indices in the data matrix.
    start_index=[];
    end_index=[];
    % try to find start and end points for the matrix.
    for i=1:size(s,1)
        for j=1:size(s,2)
            ans0=find(~isnan(s(i,j,:)));
            if ~isempty(ans0)
                start_index(i,j)=min(ans0);
                end_index(i,j)=max(ans0);
            else
                start_index(i,j)=size(s,3)/2;
                end_index(i,j)=size(s,3)/2;
            end
        end
    end

    N=50;
    z_start_bound=prctile(start_index(:),N);
    z_end_bound=prctile(end_index(:),100-N);

    % average to get rid of the phase difference in each pixel.
    if ~isempty(tf)
        z_lower_bound=roundn(z_start_bound+(z_end_bound-z_start_bound)/2,0);
        z_upper_bound=z_end_bound;
        s1=[];

        for j=1:size(s,1)
            for i=1:size(s,2)
                s1(j,i,:)=s(j,i,:)-mean(s(j,i,z_lower_bound:z_upper_bound),'all','omitmissing');
            end
        end
        s=s1;
    end

    % cut the matrix s in z direction.

    z_lower_bound=z_start_bound;
    z_upper_bound=z_end_bound;

    s1=s(:,:,z_lower_bound:z_upper_bound);
    s=s1;

    
    % try to code another way of dealing with nan values in the matrix, meshgrid.
    s1=s;   % temporarily stores the s value to be used in this section.

    for i=1:size(s,1)
        s1(i,:,:)=fillmissing2(squeeze(s1(i,:,:)),'v4');
        tt=toc;
        fprintf("Removing NAN from %s:%d/%d,%d/%d, remaning time: %0.1f s.\n",list0{ii},i,size(s,1),ii,length(list0),tt/((ii-1)*size(s,1)+i)*((length(list0)+1-ii)*size(s,1)-i));
    end

    if ii==3
        s3=s;
        s=s1;
    end


    % save the s matrix into the mat file with corresponding name.
    var_name=[list0{ii},'_filled'];
    assignin("base",var_name,s);

    path1=[path(1:end-4),'_processed_after_avg_backward.mat'];
    save(path1,var_name,'-mat','-append');

    if ii==3
        var_name=[list0{ii},'_aligned'];
        assignin("base",var_name,s3);
    end
    % save the data into file.
    path1=[path(1:end-4),'_processed_after_avg_backward.mat'];
    save(path1,var_name,'-mat','-append');

end

disp('<a href="matlab:opentoline(''force_curve_3D_5_channels_stepZ_for_missing_lines_retrace.m'',760)">Next code</a>');
beep;

%% after everything, calculate the force and energy.

 
answer = questdlg('Extend the calculation for E and F?', ...
        'Extended calculation?', ...
        'Yes','No','Yes');

v2nm_A=10;      % nm/V
v2nm_z=11*2;    % nm/V
Q=10;
k=3;            % N/m

F=[];
E=[];
F_extended=[];
E_extended=[];
error=[];
trust_values=[];

tic;
for i=1:size(r_avg_m_corr,1)
    for j=1:size(r_avg_m_corr,1)
        r_line=flip(squeeze(r_avg_m_corr(i,j,:)));
        phase_line=-flip(squeeze(phase_avg_m_corr(i,j,:)));
        z_line=-flip(squeeze(z_avg_m_corr(i,j,:)));

        % add a nan removing part. try to do it clean.

        
        ind=(i-1)*size(r_avg_m_corr,1)+j;
        ind_total=size(r_avg_m_corr,1)*size(r_avg_m_corr,2);

        fprintf('Single force curve calculation: %d/%d.',ind,ind_total);
        try
            [F0,E0,~]=force_cal(r_line.*1e-9*v2nm_A,phase_line,z_line.*v2nm_z*1e-9,Q,k);
            if strcmp(answer,'Yes')
                [F0_extended,E0_extended,trust_value]=force_cal_extended(r_line.*1e-9*v2nm_A,phase_line,z_line.*v2nm_z*1e-9,Q,k);
            end
            if (sum(imag(E0))~=0)
                temp.pos=[i,j];
                temp.error_code='Complex number detected.';
                error{end+1}=temp;
            end
        catch
        end

        t1=toc;
        fprintf('Remaining time: %0.2f min.\n',t1/ind*(ind_total-ind)/60);

        F{i,j}=F0';
        E{i,j}=E0';
        if strcmp(answer,'Yes')
            F_extended{i,j}=F0_extended';
            E_extended{i,j}=E0_extended';
            trust_values(i,j)=trust_value;
        end
    end
end


E=cell_2_3dmatrx(E);
F=cell_2_3dmatrx(F);
E=flip(E,3);
F=flip(F,3);

if strcmp(answer,'Yes')
    E_extended=cell_2_3dmatrx(E_extended);
    F_extended=cell_2_3dmatrx(F_extended);
    E_extended=flip(E_extended,3);
    F_extended=flip(F_extended,3);
end

% save the data into file.
path1=[path(1:end-4),'_processed_after_avg_backward.mat'];
save(path1,"E","F","error",'-mat','-append');

if strcmp(answer,'Yes')
    path1=[path(1:end-4),'_processed_after_avg_backward.mat'];
    save(path1,"E_extended","F_extended","trust_values",'-mat','-append');
end


% after saving all the variables, save another file containing the key
% matrices.
path1=[path(1:end-4),'_r_phase_z+E_F_backward.mat'];
save(path1,"r_avg_m_corr","phase_avg_m_corr","z_avg_m_corr","E","F","error",'-mat');

if strcmp(answer,'Yes')
    path1=[path(1:end-4),'_r_phase_z+E_F_backward.mat'];
    save(path1,"E_extended","F_extended","trust_values",'-mat','-append');
end

% for the target of z_avg_m_corr_filled data matrix.
s=z_avg_m_corr_filled;

% select the bottom plane as the base for the calculation.
s0=squeeze(s(:,:,1));


% do a line subtraction of polynomial 2.
s1=s0.*0;
x=1:size(s0,2);

for i=1:size(s0,1)
    s1(i,:)=s0(i,:)-polyval(polyfit(x,s0(i,:),2),x);
end

lpm=s1;     % give a lowest point mapping.

% give the subtracted z matrix.
s2=s.*0;
for i=1:size(s,3)
    s2(:,:,i)=squeeze(s(:,:,i))-s0+s1;
end

z_avg_m_corr_filled_sub=s2;


% save the data into file.
path1=[path(1:end-4),'_processed_after_avg_backward.mat'];
save(path1,"lpm","z_avg_m_corr_filled_sub",'-mat','-append');

% after saving all the variables, save another file containing the key
% matrices.
path1=[path(1:end-4),'_r_phase_z+E_F_backward.mat'];
save(path1,"lpm","z_avg_m_corr_filled_sub",'-mat','-append');

disp('force_curve_3D_5_channels_stepZ_for_missing_lines_retrace.m is complete, move on to next data.');
beep;

%% test function.


v2nm_A=10;      % nm/V
v2nm_z=11*2;    % nm/V
Q=10;
k=3;            % N/m

F=[];
E=[];

for i=1:size(r_avg_m,1)
    for j=1:size(r_avg_m,2)
        r_line=flip(squeeze(r_avg_m(i,j,:)));
        phase_line=flip(squeeze(phase_avg_m(i,j,:)));
        z_line=flip(squeeze(z_avg_m(i,j,:)));
        try
            [F0,E0,~]=force_cal(r_line.*1e-9*v2nm_A,phase_line,z_line.*v2nm_z*1e-9,Q,k);
        catch
        end
        F{i,j}=F0;
        E{i,j}=E0;
    end
end

%% test function region 1, try to deal with the drift. Abolished.
% suppose the constant drift in x direction, given the time duration in the
% para settings, do correction for bidirection scanning. also require the
% pixel information in para.
N=50;

flip_x=0;

if flip_x==1
    phase_avg_m=flip(phase_avg_m,3);
end
% size input.
drift_velocity=3;  % in nm/min.
scan_size=3;    % in nm.
scan_time=str2double(para.frame_time);  % get scan time from para, in s.
scan_pixel=N;   % get scan pixels from previous settings.
bidirection_flag=1; % determine the scanning type.

drift_pixel_per_line=(scan_time*drift_velocity/60)/scan_pixel/(scan_size/scan_pixel);
drift_pixel_total=round(drift_pixel_per_line*scan_pixel);

% for now, only works for bidirection scanning, only deal with phase matrix
% for now. possibly works for monodirection scanning.
phase_avg_m_drift_corrected=nan(size(phase_avg_m)+[0,0,drift_pixel_total]);
if bidirection_flag==1
    % for every scan line, do the drift correction.
    for i=1:scan_pixel
        phase0=phase_avg_m(i,:,:);
        pixel0=1+round((i-1)*drift_pixel_per_line);
        phase_avg_m_drift_corrected(i,:,pixel0:pixel0+scan_pixel-1)=phase0;
    end
end

if flip_x==1
    phase_avg_m=flip(phase_avg_m,3);
    phase_avg_m_drift_corrected=flip(phase_avg_m_drift_corrected,3);
end




%% calculating the force part.
%the force calculation meets some problems, first how to attach them into a
%bigger matrix, second the force calculation tends to give f1 as in 1x1
%matrix, ending in error, why?
clear f_m f_error f_cell;

f_m=nan(size(r_avg_m));

for i=1:size(r_avg_m,1)
    for j=1:size(r_avg_m,2)
        r_c=flip(squeeze(r_avg_m(i,j,:)));
        phase_c=flip(squeeze(phase_avg_m(i,j,:)));
        z_c=flip(squeeze(z_avg_m(i,j,:)));
        
        % processing for the phase. need to customize for each data. or set a
        % standard.
        phase_c=-phase_c;
        
        % try to find the first data in phase that is not nan.
        phase_c_tm=phase_c(~isnan(phase_c));
        phase_c=phase_c-phase_c_tm(1)+90;
        
        % parameter part.
        Q=10;
        k_stiff=3;
        
        r_c1=r_c(~isnan(r_c));
        phase_c1=phase_c(~isnan(phase_c));
        z_c1=z_c(~isnan(z_c));
        
        try
            [f,~,~]=force_cal(r_c1.*1e-9*20,phase_c1,z_c1.*v2x*1e-9,Q,k_stiff);
        catch
            f=[];
        end
        
        f=real(f);
        % how to align the calculated force, now I see that the reduced
        % number of data point is about the same overall. 22 or 23 for most
        % cases, this means probably the reduced number is associated with
        % the amplitude in the data. now need to figure out whether align it to
        % the approach point or the retract point.
        f_error0=length(r_c1)-length(f);
        f_cell{i,j}=f;
        f_tmp=nan(size(r_c));
        f_tmp1=[nan(length(r_c1)-length(f),1);f];
        f_tmp2=[f;nan(length(r_c1)-length(f),1)];
        
        f_tmp(~isnan(r_c))=f_tmp1;
        f_m(i,j,:)=f_tmp;
        f_error(i,j)=f_error0;
    end
end
f_m=flip(f_m,3);

% dock the calculated force data into the amplitude corrected matrix.
f_m_corr=nan(size(r_avg_m_corr));

for i=1:size(r_avg_m_corr,1)
    for j=1:size(r_avg_m_corr,2)
        % get the length of nan before the real data.
        n0=find(~isnan(squeeze(r_avg_m_corr(i,j,:)))>0,1);
        f_m0=squeeze(f_m(i,j,:));
        f_m0=f_m0(~isnan(f_m0));
        n1=length(f_m0);
        f_m_corr(i,j,n0:n0+n1-1)=f_m0;
    end
end

%% for the calculated force, align the force in the far end.



%%
f=real(f);

if m>n
    phase_app_avg=phase_app_avg(:,1:n);
    f=f(:,1:n);
    r_app_avg=r_app_avg(:,1:n);
    aux1_app_avg=aux1_app_avg(:,1:n);
    x_app=x_app(1:n);
end
if m<n
    phase_app_avg=[phase_app_avg,nan(length(phase_app_avg(:,1),n-m))];
    f=[f,nan(length(f(:,1),n-m))];
    r_app_avg=[r_app_avg,nan(length(r_app_avg(:,1),n-m))];
    aux1_app_avg=[aux1_app_avg,nan(length(aux1_app_avg(:,1),n-m))];
    x_app=[x_app,nan(n-m)];
end
phase_app_avg_m(k,:,:)=phase_app_avg;

%f should be enlarged to fill in the f_app_m matrix.
[~,m0,~]=size(f_app_m);

f_app_m(k,:,:)=[f;nan(m0-length(f(:,1)),n)];

r_app_m(k,:,:)=r_app_avg;
zmod_app_m(k,:,:)=aux1_app_avg;
x_app_m(k,:)=x_app;
disp(k);


%clear all the data that are used before.
if f1==0
    clear aux1 aux10 aux2 aux20 phase0 r0 i j k m n;
end

%trying to write the calculated results to file.
path1=path(1:end-4);
path1=[path1,'_results_f_phase.mat'];

save(path1,'f_app_m','r_app_m','phase_app_avg_m','zmod_app_m','x_app_m');

if exist('para','var')
    save(path1,'para','-append');
end

% removing all the parameters, this sentence should be at the very last.
if f1==0
    clear path0 path1 name0 f1 f2 f3 f4 v2x v2A name0 path0 N1 N2 xv2x range smooth;
end

%should be at the very end, save path to reading path in file.
readingpath=path;
save('userPath.mat','readingpath','-append');
clear readingpath;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test function part.
%% try to get the raw data.
s01=[];
s0=r_forward;
for i=1:size(s0,1)
    for j=1:size(s0,2)
        s01{i,j}=s0{i,j}';
    end
end

s=cell_2_3dmatrx(s01);
s=flip(s,3);


% try to do a filter instead of doing average.
s1=smoothdata(s,3,'sgolay',20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% SECTION III: local function region.
% period divide based on auxin2, give the min and max point to cut all
% other matrix.
%% period divide function.
function [min2,max2]=period_divide(auxin2,N,find_min)
L=length(auxin2);
period2=find_period_max_fun_3D_cut(auxin2(),L/N,0);

%try to plot and evaluate the accuracy of the function.
% disp('Plotting.');
% plot(auxin2());
% hold on;
% %xline(period1.loc,'r');
% xline(period2.loc,'b');
% hold off;

max2=period2.loc;

%divide all the data based on max2, not including the first one.
auxin2_table=[];
for i=1:length(max2)-1
    auxin2_table{i}=auxin2(max2(i):max2(i+1));
end

%try to use the above data to find the min2.
min2=[];
if find_min==1  %if checked, try to find the minimum position.
    for i=1:length(auxin2_table)
        fprintf('Finding local min for z ramp signal, %0.1f%% in progress.\n',i/length(auxin2_table)*100);
        try
            min2(i)=find_period_max_fun_3D_cut(auxin2_table{i},length(auxin2_table{i})/2,1);
        catch
            assignin('base',"i",i);
            assignin('base','auxin_table',auxin2_table);
            assignin('base','period2',period2);
        end
    end

    min2=max2(1:end-1)+min2';

    %then we have all the max and min. except the first min.
    %using the max(1) to deduce the min(1).
    min2=[(max2(1)-max2(2)+min2(1));min2];
end
end

%% matrix cut function based on the input min and max position.
function auxin11_forward=matrix_cut_step(min0,max0,auxin1,N_forward)
% need to input min, max, matrix, pixel(only support square).
%pixel_x=100;

% use auxin1 as example, output to auxin11_forward, auxin11_backward.
pixel_x=N_forward(end,3);

auxin11_forward=cell(pixel_x);



%try to assign the data into each pixel.
%for the cell format, {j,i}, i is x, and j is y.
for j=1:pixel_x
    for i=1:pixel_x
        if isnan(N_forward(j,1))
            auxin11_forward{j,i}=[];
        else
            auxin11_forward{j,i}=auxin1(min0(N_forward(j,1)+i-1):max0(N_forward(j,1)+i-1));
        end
    end
end
end

%% averaging function to average.
%step wise 3D averaging.
function auxin1_cut=step_3D_average(auxin1_forward,cut_method,cut_matrix)
%might require input, a matrix, method, second matrix(cut position).
%method note: 1. simple cut to get the point
[m,n]=size(auxin1_forward);

%first function: simple cut to get the cut point.
if cut_method==1
    for i=1:m
        for j=1:n
            s1=auxin1_forward{i,j};
            % 1. simple way.
            N=cut_matrix;   %determines how many points I want to keep. later change to height size.

            %now just a simple way of averaging.
            min1=min(s1);
            max1=max(s1);

            

            s_th=min1:(max1-min1)/N:max1;
            s_p=[]; %the calculated s_p is the position of the threshold cut.
            for ii=1:length(s_th)-1
                s_p(ii)=find(s1>=s_th(ii),1);
            end
            s_p(1)=1;
            s_p=[s_p,length(s1)];


            %average according to the above threshold point.
            for ii=1:length(s_p)-1
                s1_avg(ii)=mean(s1(s_p(ii):s_p(ii+1)));
            end

            auxin1_cut.position{i,j}=s_p;
            auxin1_cut.avg{i,j}=s1_avg;
        end
    end
end

%second function: cut the matrix based on height resolution.
if cut_method==2
    for i=1:m
        for j=1:n
            s1=auxin1_forward{i,j};
            if isempty(s1)
                continue;
            end
            s_p0=1;
            
            % smoothing the s1 before the cutting, get right position.
            for ii=1:2
                % 1. simple way.
                % need to verify that the retraction part is removed correctly.

                %now just a simple way of averaging.
                min1=min(s1);
                max1=max(s1);

                N=round((max1-min1)/cut_matrix);    %assuming the cut_matrix is the resolution, get a integer step number.
                try 
                    s1=smoothdata(s1,'movmean',length(s1)/N);
                catch
                    disp(i);
                    disp(j);
                end

                if (s1(1)>min1+cut_matrix)&&ii>1
                    s_p0=find(s1<min1+cut_matrix,1);
                end
            end

            s1=s1(s_p0:end);
    
            assignin('base',"s_p0",s_p0);
            assignin('base',"i",i);
            assignin('base',"j",j);

            s_th=min1:(max1-min1)/N:max1;

            % go through all the data to find the index.
            s_p=[]; % record all the index that meets the standard.
            for ii=1:length(s_th)-1
                try 
                    s_p(ii)=find(s1>=s_th(ii),1);
                catch
                    % assignin('base','i',i);
                    % assignin('base','j',j);
                    % assignin('base','ii',ii);
                    % try to display the values where something is wrong.
                    fprintf('The index that went wrong is i=%d, j=%d.\n\n',i,j);
                end
            end

            s_p(1)=1;
            s_p=[s_p,length(s1)];

            s_p=s_p+s_p0-1;


            %average according to the above threshold point.
            s1_avg=[];
            s1=auxin1_forward{i,j};
            for ii=1:length(s_p)-1
                s1_avg(ii)=mean(s1(s_p(ii):s_p(ii+1)));
            end

            auxin1_cut.position{i,j}=s_p;
            auxin1_cut.avg{i,j}=s1_avg;
        end
    end
end
%third function: cut the matrix based on the cut_matrix.
%the cut matrix should be a cell variable that is given in the cut.position
%in the previous two functions.
if cut_method==3
    auxin1_cut=[];
    for i=1:m
        for j=1:n
            s1=auxin1_forward{i,j};
            s1_avg=[];
            for ii=1:length(cut_matrix{i,j})-1
                s1_avg(ii)=mean(s1(cut_matrix{i,j}(ii):cut_matrix{i,j}(ii+1)));
            end
            auxin1_cut{i,j}=s1_avg;
        end
    end
end
end

%% square distribution function, try to cut the matrix into square pixels.
% this is the function to determine the square distribution of the matrix
function [N_forward,N_backward,missing_line_flag]=square_divide(auxin2,max00,min00,max1,min1,missing_line_flag)
% probably need a special variable to record everything. do that later.

% need the input of x scan signal, assuming auxin2.
% need the input of max and min for x scan, z ramp. 1 for x scan, 0 for z
% ramp.
max0=max00;
min0=min00;

max11=[];
max0_ind=[];
X_mean=[];
% for the calculation of the max point.
for i=1:length(max1)
    % for every given position of the max point for x scan, find the first
    % point bigger than the max value.
    N_max=find((max0>max1(i)),1);   %locate the peak position in z ramp peak index.
    for j=-2:1:2    %calculate the 5 average numbers around the peak point.
        if(N_max+j>0&&N_max+j<=length(min0))
            X_mean(j+3)=mean(auxin2(min0(N_max+j):max0(N_max+j)));
        else
            X_mean(j+3)=nan;
        end
    end

    % try to find the top two largest value for X_mean, the min for the
    % second one is the exact peak.
    [~,N]=maxk(X_mean,2);
    max11(i)=min0(N_max-3+max(N));  %record the index directly for the x scan and z ramp.
    max0_ind(i)=N_max-4+max(N);     %record the max index for z ramp, this max point still belongs to the forward motion part.
end
max11=max11';
max0_ind=max0_ind';

%i=1;
% for the calculation of the min point.
min11=[];
min0_ind=[];
X_mean=[];
for i=1:length(min1)
    % for every given position of the min point for x scan, find the first
    % point smaller than the min value.
    N_min=find((min0>min1(i)),1);   %locate the deep position in z ramp deep index.
    for j=-2:1:2
        if (N_min+j>0&&N_min+j<=length(max0))
            X_mean(j+3)=mean(auxin2(min0(N_min+j):max0(N_min+j)));
        else
            X_mean(j+3)=nan;
        end
    end

    % try to find the top two smallest value for X_mean, the min for the
    % second one is the exact deep.
    [~,N]=mink(X_mean,2);
    min11(i)=min0(N_min-3+max(N));  %record the index directly for the x scan and z ramp.
    min0_ind(i)=N_min-4+max(N);     %record the min index for z ramp, this max point belongs to the next forward motion part.
end
min11=min11';
min0_ind=min0_ind';


%evaluate the peak and deep position for the z ramp.
% the period evaluation for the max and min alone.
max_period=max0_ind(2:end)-max0_ind(1:end-1);
min_period=min0_ind(2:end)-min0_ind(1:end-1);
x_min=[];
for i=1:length(min0_ind)
    X_min(i)=mean(auxin2(min0(min0_ind(i)):max0(min0_ind(i))));
end
X_max=[];
for i=1:length(max0_ind)
    X_max(i)=mean(auxin2(max0(max0_ind(i)):max0(max0_ind(i))));
end

% based on the above data, the estimated period would be period0.
period0=ceil(mean([max_period;min_period]));
X_min0=mean(X_min);
X_max0=mean(X_max);

% further data to be added.

% the min and max located above is pretty accurate, now need to add the
% first and last min position.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEALING WITH THE FIRST MIN POSITION.
% for the first min position. use the current first max, jump half period
% backward.
N_first=max0_ind(1)-ceil(period0/2);

if missing_line_flag~=0&&N_first<1
    N_first1=[];
    missing_line_flag=2;    % 2 means extra retrace situation for the data front.
else
    % get the 15 neighboring heights to determine the real first step.
    X_mean=[];
    for j=-4:1:10
        X_mean(j+5)=mean(auxin2(min0(N_first+j):max0(N_first+j)));
    end
    X_step=X_mean(2:end)-X_mean(1:end-1);
    X_step1=isoutlier(X_step);
    N_first1=find(X_step1==0,1)-1+N_first-4;
end

%min0_ind=[N_first1;min0_ind];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEALING WITH THE LAST MIN POSITION.
% for the last min position, use the current last max, jump half period
% forward.
N_last=min(length(min0),max0_ind(end)+ceil(period0/2));

%get the 14 neighboring heights to determine the real first step.
X_mean=[];
for j=-10:1:4
    if (N_last+j>length(max0))
        break;
    end
    X_mean(j+11)=mean(auxin2(min0(N_last+j):max0(N_last+j)));
end
X_step=X_mean(2:end)-X_mean(1:end-1);
X_step1=isoutlier(X_step);
N_last1=find(X_step1==1,1)-1+N_last-10;
if (isempty(N_last1))
    N_last1=length(X_step1)+N_last-10;
end

min0_ind=[N_first1;min0_ind;N_last1];
min11=[min0(N_first1);min11;min0(N_last1)];

if missing_line_flag==2
    N_forward=[];
    for i=1:length(max0_ind)-1
        N_forward(i,1)=min0_ind(i)+1;
        N_forward(i,2)=max0_ind(i+1);
        N_forward(i,3)=N_forward(i,2)-N_forward(i,1)+1;
    end

    N_backward=[];
    for i=1:length(max0_ind)
        N_backward(i,1)=max0_ind(i)+1;
        N_backward(i,2)=min0_ind(i);
        N_backward(i,3)=N_backward(i,2)-N_backward(i,1)+1;
    end
else
    N_forward=[];
    for i=1:length(max0_ind)
        N_forward(i,1)=min0_ind(i)+1;
        N_forward(i,2)=max0_ind(i);
        N_forward(i,3)=N_forward(i,2)-N_forward(i,1)+1;
    end

    N_backward=[];
    for i=1:length(max0_ind)
        N_backward(i,1)=max0_ind(i)+1;
        N_backward(i,2)=min0_ind(i+1);
        N_backward(i,3)=N_backward(i,2)-N_backward(i,1)+1;
    end
end

end

%% z axis alignment coding.
function r_avg_m=cell_2_3dmatrx(r_avg)
% using r_avg cell variable as a reference. assuming we have the r_avg.
% note that the format for the 3d matrix is y,z,x. the format for the cell
% variable is y,x.

% go through every element of the cell variable, get the length of the z
% vector.
[m,n]=size(r_avg);
cell_length=zeros(m,n);
for j=1:m
    for i=1:n
        cell_length(j,i)=length(r_avg{j,i});
    end
end

k=ceil(mean(cell_length,'all')*1.2);
%k=ceil(max(cell_length,[],'all'));

% for now, trim every vector to the length of k.
r_avg_m=nan(m,n,k);
for j=1:m
    for i=1:n
        s0=r_avg{j,i};
        if length(s0)>=k
            s0=s0(end-k+1:end);
        else
            s0=[nan(1,k-length(s0)),s0];
        end
        r_avg_m(j,i,:)=s0;
    end
end
end

%% bin from the contact point.
function r_avg = bin_avg(r_avg_bin1,bin_size)
% assuming that I already have the r_avg_bin1, try to bin it with the
% bin_size.
[m,n]=size(r_avg_bin1);
r_avg=[];
for i=1:m
    for j=1:n
        s1=r_avg_bin1{i,j};
        s1=flip(s1);
        s2=[];
        for ii=1:bin_size:length(s1)
            if (ii+bin_size-1<=length(s1))
                s11=s1(ii:ii+bin_size-1);
                s11=s11(~isnan(s11));
                if isempty(s11)
                    s11=nan;
                end
                s2(end+1)=mean(s11);
            else
                s11=s1(ii:end);
                s11=s11(~isnan(s11));
                if isempty(s11)
                    s11=nan;
                end
                s2(end+1)=mean(s11);
            end
        end
        r_avg{i,j}=flip(s2);
    end
end
end

%% index searching by exponential fitting.
function ind=index_search(r,r0)
warning('off','all');
r=r(:);
r00=r-max(r);
y=fit((1:length(r00))',r00,'exp1');
ft=fittype('a+b*exp(c*x)');
y1=fit((1:length(r))',r,ft,'StartPoint',[max(r),y.a,y.b]);
y1.a=y1.a-r0;
ind = fzero(y1, length(r));
y1.a=y1.a+r0;
warning('on','all');
end

%% z error fixing from amplitude data and applied to amplitde and phase.
function [r,phase,z]=z_corr(r0,phase0,z0)

lpf=0;
s1=r0;
s11=phase0;
s2=ones(size(s1));
s22=ones(size(s11));
N=10;
k=1;
tic;

for i=1:size(s1,1)
    for j=1:size(s1,2)
        s3=squeeze(s1(i,j,:));
        s33=squeeze(s11(i,j,:));
        s4=s3(~isnan(s3));
        s4=lowpass(s4,1,N);
        s3(~isnan(s3))=s4;
        s2(i,j,:)=s3;
        if lpf==1
            s33(~isnan(s33))=lowpass(s33(~isnan(s33)),1,N);
            s22(i,j,:)=s33;
        end
        tt=toc;
        fprintf('Low pass smoothing amplitude: %0.2f %%. Remaining time: %0.1f s.\n',round(100*k/(size(s1,1)*size(s1,2)),2),(size(s1,1)*size(s1,2)-k)*tt/k);
        k=k+1;
    end
end

%%

s0=s2;
phase_avg_bin1_m_lpf=s22;   % filtered phase data, not used for now.

s=s0;
s=permute(s,[3 2 1]);
s1=reshape(s,[size(s,1),size(s,2)*size(s,3)]);

s1=permute(s1,[2 1]);
%%
% use anything larger than this z index to calculate an average amplitude
% to normalize the amplitude data. the ratio of the data would be in
% s1.

lower_bound=0.5;

% try to get the largest few points in each spectrum.
s1_max=mean(s1(:,(round(end*lower_bound):end)),2,'omitmissing');
s1=s1./s1_max;

%%
% find the min data for each pixel.
s1_min=min(s1,[],2);
s1_min_all=max(s1_min);

offset=zeros(size(s1_min));
for i=1:size(s1_min,1)
    s1_x=s1(i,:);
    try
        offset(i)=find(s1_x>s1_min_all,1);
    catch
        if isnan(s1_x)
            offset(i)=0;
        else
            offset(i)=offset(i-1);
        end
    end
end

%%
% the noted part is the additional alignment using cross-correlation, not
% good.
% 
% offset_max=max(offset);
% s2=zeros(size(s1)+[0,offset_max]);
% for i=1:size(s1_min,1)
%     s2(i,:)=[nan(1,offset_max-offset(i)),s1(i,:),nan(1,offset(i))];
% end
% 
% % need to get the s2 for later computation.
% for i=0:size(s,2)-1
%     s21=s2(1+i*size(s,3):size(s,2)+i*size(s,3),:);
%     s21_ref=s21(1,offset_max:end);
%     s3=s21(:,offset_max:end);
%     offset_corr=[];
%     % the offset_max should be the meeting point of everyone.
%     for j=1:size(s21,1)
%         s_x=s3(j,:);
%         s_x=s_x(~isnan(s_x));
%         if isempty(s_x)
%             offset_corr(j)=0;
%             continue;
%         end
%         s2_ref0=s21_ref(1:length(s_x));
%         [r,lags]=xcorr(s2_ref0,s_x);
%         offset_corr(j)=lags(r==max(r));
%     end
%     offset(1+i*size(s,3):size(s,2)+i*size(s,3))=offset(1+i*size(s,3):size(s,2)+i*size(s,3));%+offset_corr';
% end
% 
% offset=offset-min(offset);

% correct for r
s=r0;
s=permute(s,[3 2 1]);
s1=reshape(s,[size(s,1),size(s,2)*size(s,3)]);
s1=permute(s1,[2 1]);

offset_max=max(offset);
s2=zeros(size(s1)+[0,offset_max]);
for i=1:size(offset,1)
    s2(i,:)=[nan(1,offset_max-offset(i)),s1(i,:),nan(1,offset(i))];
end

N=offset_max;
s2=permute(s2,[2 1]);
s1=reshape(s2,size(s)+[N,0,0]);
r=permute(s1,[3 2 1]);


% correct for phase
s=phase0;
s=permute(s,[3 2 1]);
s1=reshape(s,[size(s,1),size(s,2)*size(s,3)]);
s1=permute(s1,[2 1]);

offset_max=max(offset);
s2=zeros(size(s1)+[0,offset_max]);
for i=1:size(offset,1)
    s2(i,:)=[nan(1,offset_max-offset(i)),s1(i,:),nan(1,offset(i))];
end

N=offset_max;
s2=permute(s2,[2 1]);
s1=reshape(s2,size(s)+[N,0,0]);
phase=permute(s1,[3 2 1]);


% correct for z
s=z0;
s=permute(s,[3 2 1]);
s1=reshape(s,[size(s,1),size(s,2)*size(s,3)]);
s1=permute(s1,[2 1]);

offset_max=max(offset);
s2=zeros(size(s1)+[0,offset_max]);
for i=1:size(offset,1)
    s2(i,:)=[nan(1,offset_max-offset(i)),s1(i,:),nan(1,offset(i))];
end

N=offset_max;
s2=permute(s2,[2 1]);
s1=reshape(s2,size(s)+[N,0,0]);
z=permute(s1,[3 2 1]);
end