function [] = printMatrix(A,width,decimals)

for i = 1:size(A,1)
    for j = 1:size(A,2)
        fprintf(['%' num2str(width) '.' num2str(decimals) 'f\t'],A(i,j))
    end
    fprintf('\n')
end