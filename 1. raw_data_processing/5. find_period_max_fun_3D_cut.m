%a function designed for force_curve_3D_5_channels function.
%used to find the period in the set function. 8/7/23.
function period = find_period_max_fun_3D_cut(aux1,wind_size,single_flag)
%% 

%parameter description:
%1. aux1, the original signal;
%2. sampling_rate, used to determine the smooth window size.
%3. wind_size, ask for input a estimation of the period.

pks_win=0.5;    %period seperation threshold, if the two peaks are too close to each other. remove the second one?


if single_flag==1
    N=0.01;
    s11=smoothdata(medfilt1(aux1,ceil(wind_size/10)),'movmean',ceil(wind_size/10));
    s12=diff(smoothdata(s11,'movmean',ceil(wind_size/10)));
    s12(1:ceil(N*end))=nan;
    s12(ceil(end-N*end):end)=nan;
    s13=diff(s12);
    period=find(s13==max(s13),1);
    return;
end

%try to use the same function to find the period first,then use the period
%to get the filter window for the function.

%smoothing part of the function.

%if method is -1, then the window size of the filter is not given, try
%initialize the function with starting parameters.

aux12=smoothdata(medfilt1(aux1,ceil(wind_size/10)),'movmean',ceil(wind_size/10));

%use the peak width and prominence to determine if the peak is normal.
[pks,y,w,p]=findpeaks(aux12);
y_n=find(pks<max(pks)/2);
%y_n=find(isoutlier(w).*isoutlier(p));

%need to run the data set again to find the too close peaks.
y(y_n)=[];
w(y_n)=[];
p(y_n)=[];
y_period=y(2:end)-y(1:end-1);
y_period0=mean(y_period);
y_s=y_period<0.5*y_period0;


y(y_s)=[];
w(y_s)=[];
p(y_s)=[];

%this gives the outlier of the peaks determined, but another problem is if
%I need to add points before and after the first and last point. the idea
%would be to add 1* period to the first or last point.

y_period=y(2:end)-y(1:end-1);
y_period0=mean(y_period);
% s0=['The period is ',num2str(y_period0),'.'];
% disp(s0);

%try to find the real peaks for the period.
y1=y;
for i=1:length(y)
    pks_win_l=floor(y(i)-w(i)*pks_win);
    if pks_win_l<1
        pks_win_l=1;
    end
    pks_win_r=ceil(y(i)+w(i)*pks_win);
    if pks_win_r>length(aux1)
        pks_win_r=length(aux1);
    end
    y1(i)=find(aux1(pks_win_l:pks_win_r)==max(aux1(pks_win_l:pks_win_r)),1,'last')+pks_win_l-1;
end

period={};
period.loc=y1;
period.loc_des='the location of the period maximum.';
period.hw=w;
period.hw_des='the half width of the period';
period.period_avg=y_period0;
period.period_avg_des='the average period';

%I think I should add a part to test if the data is as I required, if
%something went wrong, let the user know. Not yet done.
% figure;
% plot(aux1);
% hold on;
% plot(y1,aux1(y1),'o');
% hold off;
% 
% figure;
% plot(aux12);
% hold on;
% plot(y,aux12(y),'o');
% hold off;
end