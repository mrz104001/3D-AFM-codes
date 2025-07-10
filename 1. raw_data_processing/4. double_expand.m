%% for r1, how to add further data points to fill up for the calculation of energy.

function r2=double_expand(r1)
% try to fit exponential.

UpLim = [90 10000 90 10000 90];
LoLim = [-90 0 -90 0 0 ];
StPoint = [0.5 100 0.5 100 0];

foption = fitoptions('method','NonlinearLeastSquares','Upper',UpLim,...
        'Lower', LoLim, 'StartPoint', StPoint, 'normalize', 'off');

fittype1 = fittype('a1 * exp(-x / b1) + a2 * exp(-x/b2) + c ',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a1' 'b1' 'a2' 'b2' 'c'});

x=(1:length(r1))';
y=r1;



x(isnan(y))=[];
y(isnan(y))=[];

% do a normalization.
[yy,nc,ns]=normalize(y);


exps0 = fit(x, yy,fittype1,foption);

x_add=(x-min(x)+1)+max(x);
x_add=[x;x_add];

y1=feval(exps0,x_add).*ns+nc;

% figure;
% plot(x,y,'o');
% hold on;
% plot(x_add,y1);
% hold off;

r2=nan(max(x_add),1);
r2(x_add)=y1;
r2(x)=y;
end