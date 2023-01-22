
a = 500
b = 1000
ctr=1
i=0:0.05:2
area=zeros(1,numel(i)) % preallocate
for i=i
    elements = a * i 
    area(ctr) = b + elements 
    ctr=ctr+1
end

