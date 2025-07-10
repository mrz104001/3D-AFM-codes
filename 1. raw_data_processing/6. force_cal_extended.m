%trying to build a function to calculate force via A and phi.
%at least the r and phi data should be inputted, the other five can be set
%as 1 or left blank. h can be inputted with r and phase.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT !!!!
% phase should be in degrees, being 90 in the far side.
% r should be in nm and everything should be given in real value to get
% real force values in nN, except w0, w could be left blank.
% all the data should be column vectors.
% the input data should be full for now.

% 1.  r and phase could be 2d or 1d matrix. be sure to align the z direction in
% the first dimension, Nx1 vector.

% 2.  the direction of the approach. as the row number grows, the tip should be
% closer to the surface.

% 3.  the phase should be 90 degrees, and decreasing with approaching.

% 4.  requirement on z? the zero point would be which side?

%define the input parameters.



function [f,f2,trust_value] = force_cal_extended(r,phase,h,Q,k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testing parameters.
% v2nm_A=10;      % nm/V
% v2nm_z=11*2;    % nm/V
% 
% 
% n1=4;
% n2=1;
% r=flip(squeeze(r_avg_m(n1,n2,:))*v2nm_A*1e-9);
% phase=flip(squeeze(phase_avg_m(n1,n2,:)));
% phase=90+mean(phase(1:20),'omitmissing')-phase;
% 
% h=-flip(squeeze(z_avg_m(n1,n2,:))*v2nm_z*1e-9);
% h=h-min(h);
% 
% k=3;
% Q=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% flip the direction of the height change.
h=h-min(h);

h00=h;
sfun0=0;
f=[];
f2=[];
f20=[];

if sfun0==0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %using Holscher's method to reconstruct the force.
    %also assuming the approaching measurements, with increasing row number,
    %the tip is closer to the surface.
    %requiring the r and z be in the same unit so that subtraction could be
    %done.
    
    %first, calculate the integral.
    %take the first line as example.
    f=[];
    for i=1:length(r(1,:))
        r1=flip(r(:,i));
        phase1=flip(phase(:,i));
        phase0=max(phase1);
        phase1=phase1+90-phase0;
        h1=flip(h00(:,i));

        % double expand on r2 and phase2.
        r2=double_expand(r1);
        phase2=double_expand(phase1);
        % x1=1:length(r1);
        % x2= 1:length(r2);
        % x1(isnan(h1))=[];
        % h1(isnan(h1))=[];
        % h2=interp1(x1,h1,x2','linear','extrap');
        % % h2(1:length(h1)) = h1;

        h2=double_expand(h1);




        kappa1=(k*max(r1)/Q.*sqrt(r1./2).*cos(phase1/180*pi));
        kappa=(k*max(r2)/Q.*sqrt(r2./2).*cos(phase2/180*pi));
        f1=[];
        
        % rewrite the intergral part for the conservative energy.

        for i1=1:length(r1)
            % check if can get full range of z - z+2A.
            if (h2(i1)+2*r2(i1))>max(h2)
                break;
            end
            if isnan(r2(i1))%||isnan(h2(i1))
                f1(i1)=nan;
                continue;
            end

            i2=find(h2>(h2(i1)+2*r2(i1)),1)-1;
            %disp(i1);

            f1(i1)=0;
            
            for i3=i1+1:i2
                f1(i1)=f1(i1)+kappa(i3)/sqrt(h2(i3)-h2(i1))*(h2(i3+1)-h2(i3));
            end
            
            % i1=i1;
            % % stop.
        end

        % give a trust value.
        for i1=1:length(r1)
            if (h1(i1)+2*r1(i1))>max(h1)
                break;
            end
        end

        trust_value=i1;

        % assignin("base",'f1',r);


        f1=sgolayfilt(f1,1,3);
        for i1=1:length(f1)-1
            %find out whether the minus sign is needed.
            f(i1,i)=-(f1(i1+1)-f1(i1))/(h2(i1+1)-h2(i1));
        end
    end
end
f=flip(f,1);
f2=flip(f1',1);



%% currently not using.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sfun0==1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %using Payam's method to reconstruct the force.

    %since it is done via two integrals, for the loop, the formula is divided
    %into several parts for the calculation.
    [m,n]=size(r);

    A0=max(max(r));
    phase0=max(max(phase));
    phase=phase+90-phase0;



    phase=smoothdata(phase,1,'movmean',20);

    %gives the option to use default values for h.
    if h==1
        h=1:m;
    else
        h0=min(h);
        h1=max(h);
        ha=(h1-h0)/(m-1);
        h=0:ha:(h1-h0);
    end

    f2=zeros(m,n);
    f20=zeros(m,n);
    f1=zeros(m,n);
    f3=zeros(m,n);
    A=zeros(1,m);


    %use the first two terms of the formula to calculate force?
    for i=1:n
        for j=1:m
            A(j)=A0*cos(phase(j,i)/180*3.14)/(2*Q*r(j,i));
            f20(j,i)=A(j)*(h(2)-h(1))*2*k;
        end
        f20(:,i)=f20(:,i)-mean(f20(1:20,i));
        for j=1:m
            for j1=1:j
                f2(j,i)=f2(j,i)+f20(j1,i);
            end
        end
        for j=m:-1:1
            if j>1
                d=h(j-1);
            end
            if j==1
                f1(j,i)=NaN;
                continue;
            end
            for j1=m:-1:j
                B=r(j1,i)^1.5/sqrt(2*(h(j1)-d+r(j1,i)));
                f1(j,i)=f1(j,i)+B*A(j1)*(h(2)-h(1));
            end
        end
        for j=m:-1:1
            if j==m
                f3(j,i)=(f1(j,i)-f1(j-1,i))/(h(2)-h(1))*2*k;
            else
                f3(j,i)=(f1(j+1,i)-f1(j,i))/(h(2)-h(1))*2*k;
            end
        end
        f=f2-f3;
    end

    %fitting the f to a function to subtract background.
    f=f(2:end,:);

    phase_fun=@(p,x) p(1).*exp(x.*p(2))+p(3);
    startvalue=[1 1 0];
    options=optimset('TolFun',1e-15,'TolX',1e-15,'Display','none');
    [p,~,~,~,~,~,~]=lsqcurvefit(phase_fun,startvalue,(2:m)'./m,mean(f,2)*1e10,[1e-50 0 0],[],options);
    f0=phase_fun(p,(2:m)./m)/1e10;

    for i=1:n
        f_s(:,i)=f(:,i)-f0';
    end

    figure;
    plot(f_s);

    %try fitting phase to get rid of the phase shifts.
    phase0=smoothdata(phase,1,"movmean",50);
    phase=phase-phase0;

    f000=smoothdata(f,1,"movmean",20);
    f=f-f000;

end


