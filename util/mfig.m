
function h=mfig(varargin)
if nargin>0
    n=varargin{1};
    h=figure(n);
    set(h,'menubar','none','toolbar','none');
else 
    figure('menubar','none','toolbar','none')
end
set(gcf,'PaperPositionMode','auto');
set(gcf,'color','w')