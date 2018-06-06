clear all
syms w(r) k M y(t)
k = 10^9 
M = 500
cond = w(6) == 1
ode2 = w + (r * diff(w, r))/2 == -(k * (w(r) - 1))/(M * r^4) 
sol1 = simplify(dsolve(ode2, cond))
ezplot(sol1,[6,40])


