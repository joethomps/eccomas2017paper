function [J U] = lanczos(La,x1)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[n w] = size(La);
a = zeros(1,n);
z = zeros(n,n-1);
b = zeros(1,n-1);
X = [x1,zeros(n,n-1)];

a(1) = X(:,1).'*La*X(:,1);
z(:,1) = a(1).*X(:,1) - La*X(:,1);
b(1) = norm(z(:,1));
X(:,2) = z(:,1)./b(1);
for i = 2:n-1
    a(i) = X(:,i).'*La*X(:,i);
    z(:,i) = a(i).*X(:,i) - b(i-1).*X(:,i-1) - La*X(:,i);
    b(i) = norm(z(:,i));
    X(:,i+1) = z(:,i)./b(i);
end
a(n) = X(:,n).'*La*X(:,n);

J = diag(a,0)-diag(b,-1)-diag(b,1);
U = X.';