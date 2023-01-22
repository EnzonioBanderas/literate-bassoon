function fl = func_dirl2fl(dirl)
%FUNC_DIRL2FL Summary of this function goes here
%   Detailed explanation goes here

n = length(dirl);
fl = cell(n, 1);
for i = 1:n
    fl{i} = fullfile(dirl(i).folder, dirl(i).name);
end

end

