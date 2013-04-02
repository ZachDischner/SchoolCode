function plotZ(x,y)
idx = [1:105:length(x)];
plot(x(idx),y(idx,1),'r.',x(idx),y(idx,2),'b--',x(idx),y(idx,3),'g-','MarkerSize',2)

end