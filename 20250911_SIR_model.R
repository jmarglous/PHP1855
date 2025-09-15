beta <- 1/300
r <- .2
t_max <- 30 

S <- numeric(t_max)
I <- numeric(t_max)
R <- numeric(t_max)
N <- numeric(t_max)


S[1] <- 299
I[1] <- 1
R[1] <- 0
N[1] <- S[1] + I[1] + R[1] 

for (t in 1:(t_max - 1)) {
  S_current <- S[t]
  I_current <- I[t]
  R_current <- R[t]
  
  S_new <- S_current - beta * S_current * I_current
  I_new <- I_current + beta * S_current * I_current - r * I_current
  R_new <- R_current + r * I_current  
  N_new <- S_new + I_new + R_new
  
  S[t+1] <- S_new
  I[t+1] <- I_new
  R[t+1] <- R_new
  N[t+1] <- N_new
}

print(S)
print(I)
print(R)

data.frame(S, I, R, N)
plot(S, type = 'l', lwd = 2, col = 'green')
lines(I, lwd = 2, col = 'red')
lines(R, lwd = 2, col = 'blue')
lines(N, lwd = 2, col = 'black')


legend("right", legend = c("Susceptible (S)", "Infected (I)", "Recovered (R)", "total (N)"),
       col = c("green", "red",'blue', 'black'), pch = 1)

max(I)
