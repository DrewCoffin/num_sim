function[] = randplot(num, scale, color)
hold off
arr = (rand([1,num])-0.5)*scale;
bar = rand([1,num])*(scale/10);
x = [1:num];
errorbar(x,arr,bar);
hold on
avg = ones([1,num])*mean(arr);
plot(x,avg)
sig = ones([1,num])*std(arr);
yup = avg + sig;
ydown = yup - 2*sig;
plot(x,yup)
plot(x,ydown)
x2 = [x, fliplr(x)];
fill(x2,[avg,yup], color, 'FaceAlpha', .05)
fill(x2,[avg,ydown], color, 'FaceAlpha', .05)
end