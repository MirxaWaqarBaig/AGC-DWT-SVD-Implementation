function [I]=GHE_gray(I)
[R, C] = size(I);

cnt = zeros(1, 256);
for i = 1 : R
    for j = 1 : C
        cnt(1, I(i, j) + 1) = cnt(1, I(i, j) + 1) + 1;
    end
end

f = zeros(1, 256);
f = double(f); cnt = double(cnt);


for i = 1 : 256
    f(1, i) = cnt(1, i) / (R * C);
end


for i = 2 : 256
    f(1, i) = f(1, i - 1) + f(1, i);
end


for i = 1 : 256
    f(1, i) = f(1, i) * 255;
end


I = double(I);
for i = 1 : R
    for j = 1 : C
        I(i, j) = f(1, I(i, j) + 1);
    end
end