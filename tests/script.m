clear;
x = dir();
index = 1;
for ii=1:length(x)
    [~,~,ext] = fileparts(x(ii).name);
    if strcmp(ext,'.1')
        x(ii).name
        fileID = fopen(x(ii).name,'r');
        formatSpec = '%s';
        A = textscan(fileID,formatSpec,'Delimiter','\n');
        [tmp(index).time tmp(index).cost tmp(index).tours tmp(index).mem] = parse(A);
        index=index+1;
        fclose(fileID);
    end
end

for ii=1:length(tmp)
    Means(1).time(ii) = mean(tmp(ii).time);
    Means(1).mem(ii) = mean(tmp(ii).mem);
end


% all(Info(1).cost == Info(1).cost(1)) %confirmation that all costs are the same

clear Info;
clear A;
clear index;
clear ans;
clear ext;
clear fileID;
clear ii;
clear x;
clear formatSpec;
clear tmp;

x = dir();
index = 1;
for ii=1:length(x)
    [~,~,ext] = fileparts(x(ii).name);
    if strcmp(ext,'.2')
        fileID = fopen(x(ii).name,'r');
        formatSpec = '%s';
        A = textscan(fileID,formatSpec,'Delimiter','\n');
        [tmp(index).time tmp(index).cost tmp(index).tours tmp(index).mem] = parse(A);
        index=index+1;
        fclose(fileID);
    end
end

for ii=1:length(tmp)
    Means(2).time(ii) = mean(tmp(ii).time);
    Means(2).mem(ii) = mean(tmp(ii).mem);
end

clear Info;
clear A;
clear index;
clear ans;
clear ext;
clear fileID;
clear ii;
clear x;
clear formatSpec;
clear tmp;

x = dir();
index = 1;
for ii=1:length(x)
    [~,~,ext] = fileparts(x(ii).name);
    if strcmp(ext,'.4')
        fileID = fopen(x(ii).name,'r');
        formatSpec = '%s';
        A = textscan(fileID,formatSpec,'Delimiter','\n');
        [tmp(index).time tmp(index).cost tmp(index).tours tmp(index).mem] = parse(A);
        index=index+1;
        fclose(fileID);
    end
end

for ii=1:length(tmp)
    Means(3).time(ii) = mean(tmp(ii).time);
    Means(3).mem(ii) = mean(tmp(ii).mem);
end

clear Info;
clear A;
clear index;
clear ans;
clear ext;
clear fileID;
clear ii;
clear x;
clear formatSpec;
clear tmp;

x = dir();
index = 1;
for ii=1:length(x)
    [~,~,ext] = fileparts(x(ii).name);
    if strcmp(ext,'.8')
        fileID = fopen(x(ii).name,'r');
        formatSpec = '%s';
        A = textscan(fileID,formatSpec,'Delimiter','\n');
        [tmp(index).time tmp(index).cost tmp(index).tours tmp(index).mem] = parse(A);
        index=index+1;
        fclose(fileID);
    end
end

for ii=1:length(tmp)
    Means(4).time(ii) = mean(tmp(ii).time);
    Means(4).mem(ii) = mean(tmp(ii).mem);
end

clear Info;
clear A;
clear index;
clear ans;
clear ext;
clear fileID;
clear ii;
clear x;
clear formatSpec;
clear tmp;

x = dir();
index = 1;
for ii=1:length(x)
    [~,~,ext] = fileparts(x(ii).name);
    if strcmp(ext,'.serial')
        fileID = fopen(x(ii).name,'r');
        formatSpec = '%s';
        A = textscan(fileID,formatSpec,'Delimiter','\n');
        [tmp(index).time tmp(index).cost tmp(index).tours tmp(index).mem] = parse(A);
        index=index+1;
        fclose(fileID);
    end
end

for ii=1:length(tmp)
    Means(5).time(ii) = mean(tmp(ii).time);
    Means(5).mem(ii) = mean(tmp(ii).mem);
end

clear Info;
clear A;
clear index;
clear ans;
clear ext;
clear fileID;
clear ii;
clear x;
clear formatSpec;
clear tmp;

for i=1:6
    for ii=1:4
        a(i).time(ii) = Means(ii).time(i);
        a(i).mem(ii) = Means(ii).mem(i);
        
        a(i).speedup(ii)  = Means(5).time(i) / Means(ii).time(i);
        
    end
end



for ii=1:4
    speedup_mean(ii)= 0;
    speedup_max(ii) = 0;
    speedup_min(ii) = 10000;
    
    for i=1:6
        speedup_mean(ii) = speedup_mean(ii)+ a(i).speedup(ii);
        if(speedup_max(ii) <  a(i).speedup(ii))
            speedup_max(ii)  = a(i).speedup(ii);
        end
        if(speedup_min(ii) >  a(i).speedup(ii))
            speedup_min(ii)  = a(i).speedup(ii);
        end
    end
    speedup_mean(ii) = speedup_mean(ii)/7;
end

% Create a figure and plot the data
figure;
colors = {'b-', 'g-', 'r-', 'c-', 'm-', 'y-', 'k-', 'b--', 'g--', 'r--'};

hold on;
for i=1:6
    plot([1, 2, 3, 4], a(i).speedup, colors{i});
end


title('OMP Results');
xlabel('Number of Threads');
ylabel('Speedup referenced to 1 thread run');


legend('gen19', 'gen20', 'gen22', 'gen24','gen26', 'gen30');
figure;
hold on;
for i=1:6
    plot([1, 2, 3, 4], a(i).mem, colors{i});
end
legend('gen19', 'gen20', 'gen22', 'gen24','gen26', 'gen30');
title('OMP Memory Results');
xlabel('Number of Threads');
ylabel('Memory KB');