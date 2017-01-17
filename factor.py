def factor(n):
    factors = []
    d = 2
    while n>1:
        while n%d == 0:
            factors.append(d)
            n /= d
        else:
            break
    d = 3
    while n>1:
        while n%d == 0:
            factors.append(d)
            n /= d
        d = d + 2
        if d*d > n:
            if n > 1: factors.append(n)
            break
    return factors    

num = input('Enter number to factor: ')
print factor(num)
