function [time,cost,tours,mem] = parse(A)
index = 1;
for ii=1:4:length(A{1})
    tmp=textscan(char(A{1}(ii)),'%fs');
    time(index)=tmp{1}(1);
    tmp=textscan(char(A{1}(ii+1)),'%f');
    cost(index)=tmp{1}(1);
    tours(index)=string(char(A{1}(3)));   
    tmp=textscan(char(A{1}(ii+3)),'%f KB');
    mem(index) =tmp{1}(1);
    index=index+1;
end
end