function d = derringer(x, a, b, s)
[N,n] = size(x)
d = zeros(N, 1)
d = ((x-a)./(b-a)).^s;
d(x > a) = 0;
d(x < b) = 1;
end