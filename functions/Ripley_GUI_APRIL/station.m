function station=station(data, doplot)
% station=station(data, doplot)
% turns a distribution in space
% Input:
% data - a matrix of format [X Y Z] representing the distribution
% doplot - whether to plot a visualization of how data was turned (0 or 1)
% Output:
% station - a matrix of format [X Y Z] with the turned plot
%
% License: RipleyGUI is distributed free under the conditions that
% (1) it shall not be incorporated in software that is subsequently sold; 
% (2) the authorship of the software shall be acknowledged in any publication that uses results generated by the software; 
% (3) this notice shall remain in place in each source file. 
if (nargin < 2) || isempty(doplot)
    doplot = 1;
end
data=sortrows(data);
siz=size(data);
dim=siz(1);
s=(max(data)-min(data));
if siz(2)>2
    data2d=zeros(dim,2);
    for i=1:3
        if not(s(i)==min(s))
            data2d(:,i)=data(:,i);
        end
    end
end
if s(2)>s(1)
    temp=data(:,1);
    data(:,1)=data(:,2);
    data(:,2)=temp;
    clear temp
end
s=(max(data)-min(data));
xmin=min(data(:,1));xmax=max(data(:,1));
ymin=min(data(:,2));ymax=max(data(:,2));
% zmin=min(data(:,3));zmax=max(data(:,3));
l=sqrt((xmin-xmax)^2+(ymin-ymax)^2);
f=asin(s(2)/l);
datarot=zeros(dim,3);
if f>0.08*pi
    if and(ismember(xmin,data(1:round(length(data)/2),1)),ismember(ymin,data(1:round(length(data)/2),2)))
        for i=1:length(data)
            datarot(i,1)=data(i,1)*cos(f)+data(i,2)*sin(f);
            datarot(i,2)=data(i,1)*sin(f)-data(i,2)*cos(f);
        end
    else
        for i=1:length(data)
            datarot(i,1)=data(i,1)*cos(f)-data(i,2)*sin(f);
            datarot(i,2)=data(i,1)*sin(f)+data(i,2)*cos(f);
        end
    end
    datarot(:,3)=data(:,3);
else
    datarot=data;
end
% a1=s(1)*s(2);
% sr=(max(datarot)-min(datarot));
% a2=sr(1)*sr(2);
% station.prop=a2/a1;
if doplot
    subplot(1,2,1)
    scatter(data(:,1),data(:,2)),axis equal,hold on
    plot([min(data(:,1)):1:max(data(:,1))],min(data(:,2)),'r')
    plot([min(data(:,1)):1:max(data(:,1))],max(data(:,2)),'r')
    subplot(1,2,2)
    scatter(datarot(:,1),datarot(:,2)),axis equal,hold on
    plot([min(datarot(:,1)):1:max(datarot(:,1))],min(datarot(:,2)),'r')
    plot([min(datarot(:,1)):1:max(datarot(:,1))],max(datarot(:,2)),'r')
end
station=datarot;