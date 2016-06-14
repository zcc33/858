canc=[];
[n,m] = size(Data.GE);
for i = 1:m
    a=Data.sample(i);
    b=a{1};
    if b(end-1)=='0'
        canc(i) = 1;
    else
        canc(i) = 0;
    end
end