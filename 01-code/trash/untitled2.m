
x = chebfun('x');

F = @(x) diff(besselj(2,x),x);
f = chebfun(F,[-100,100]);

figure
r = roots(f,'complex'); 
hold off, plot(r,'.')



    temp1=roots(diff(chebfun(@(t) besselj(2,t),[0,600]))+0.00001,'complex');
