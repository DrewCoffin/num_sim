function fib = fibonacci(n)
fib = zeros(n:1);
fib(1)=1;
fib(2)=1;
k=3;
for k = 3:n
    fib(k) = fib(k-2) + fib(k-1);
end
