function z = applypriori(z,prioris)
%applypriori apply some a-priori information
%   此处显示详细说明
% 'sparsity',40,'real-valuedness','nonnegativity','support'};
d = find(strcmp(prioris,'support'));
if d
    z = zerosadd(zerosadd(z,2,2),1,2);
end
a=find(strcmp(prioris,'real-valuedness'));
if a
   z = real(z); 
end
b = find(strcmp(prioris,'nonnegativity'));
if b
   z = max(z,0);
end
c = find(strcmp(prioris,'sparsity'));
if c
    tempz = sum(abs(z).^2,2);
    [~,inx] = sort(tempz,'descend');
    y = zeros(size(z));
    y(inx(1:prioris{c+1}),:) = z(inx(1:prioris{c+1}),:);
    z = y;
end
end

