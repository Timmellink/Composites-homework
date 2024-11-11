function array = zArray(h,n,z)
% input for z is vector of angles 
%% create an array ranging from h = -h/2 to h = h/2, divided in n sublengths
k = h/n;
zn = -h/2:k:h/2;
txt = ["length of zn = "+length(zn)];
disp(txt)
array = cell(2,length(n));
i = 1;
for x = zn
    if i<=length(z)
        array{i,1} = zn;
        cell_x = {x z(i)};
        array{i,2} = cell_x;
    else
        array{i,1} = zn;
        ind = length(z)+5-i
        cell_x = {x, z(ind)};
        array{i,2} = cell_x;
    end
    i = i+1
end
