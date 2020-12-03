#MODEL BUILDING (up to line 353)

population = 3200000
I0 = 2734
S0 = population - I0
R0 = 2226

setwd("~/R/my_r_progs/M3550_Project")
library(R0)

df1 = read.csv('rt.csv')
df1 = df1[which(df1$region == "UT"),]
df1$date = as.Date(df1$date)
df1[is.na.data.frame(df1)]<-0
df1 = df1[which(df1$date >= "2020-3-6" & df1$date <= "2020-10-31"),]

df2 = read.csv('Overview_Cumulative COVID-19 Cases with Estimated Recoveries_2020-11-01.csv')
df2$Date = as.Date(df2$Date)
df2 = df2[which(df2$Date >= "2020-3-6" & df2$Date <= "2020-10-31"),]

covid = data.frame(df1$date,df2$Total,df1$new_deaths,df2$Died,
                   df2$Estimated.Recovered..,df2$Estimated.Active)
names(covid) = c('date','total.infected','new.deaths','total.deaths',
                 'total.recovered','active.infected')
plot(covid$active.infected)
new.recovered = covid$total.recovered
i = length(new.recovered)
while(i>1){
  new.recovered[i] = new.recovered[i] - new.recovered[i-1]
  i = i-1
}
new.removed = new.recovered + covid$new.deaths
total.removed = new.removed
i=1
while(i<length(total.removed)){
  total.removed[i+1] = total.removed[i+1] + total.removed[i]
  i=i+1
}
covid$total.removed = total.removed
covid$susceptible = population - (covid$active.infected + covid$total.removed)

init.df = covid[which(covid$date <= "2020-5-1"),]

GT = generation.time("gamma", c(17,17/2))
Rt = est.R0.TD(init.df$active.infected, GT, n.t0 = init.df$active.infected[1],
               t = init.df$date)
Rt$R[length(Rt$R)] = (Rt$R[length(Rt$R)-1]^2)/Rt$R[length(Rt$R)-2]
Rt.Rat = Rt$R  #Ratio of subsequent Rt to current Rt
t = 1
while(t < length(Rt.Rat)){
  Rt.Rat[t] = Rt.Rat[t+1]/Rt.Rat[t]
  t = t+1
}
Rt.Rat[t] = Rt.Rat[t-1]
plot(Rt.Rat)

#Our constants for this model are population, a, and lag
#Our initial conditions are I0, S0, R0, and Rt0

model.df = covid[which(covid$date >= "2020-5-1"),]

a = 1/25
lag = 23
I = model.df$active.infected
S = model.df$susceptible
R = model.df$total.removed
R.rat.t = 1:length(model.df$date)
R.t.col = 1:length(model.df$date)
I[1]
S[1]
R[1]
Rt0 = Rt$R[length(Rt$R)]
t = 1
R.t = Rt0
R.t
R.rat = .99
t.end = t + 14
while(t < t.end){
  R.rat.t[t] = R.rat
  R.t = R.t*R.rat
  R.t.col[t] = R.t
  r = a*R.t/S[t]
  S[t+1] = S[t] - r*I[t]*S[t]
  I[t+1] = I[t] + r*I[t]*S[t] - a*I[t]
  R[t+1] = R[t] + a*I[t]
  t = t+1
}
plot(model.df$susceptible~model.df$date)
points(S~model.df$date,col = 'red')
plot(model.df$active.infected~model.df$date)
points(I~model.df$date,col = 'red')
plot(model.df$total.removed~model.df$date)
points(R~model.df$date,col = 'red')
t.end
S[t.end]
I[t.end]
R[t.end]
model.df$active.infected[t.end]
R.t
model.df$date[t.end]

R.rat = 1.02
t.end = t + lag
while(t < t.end){
  R.rat.t[t] = R.rat
  R.t = R.t*R.rat
  R.t.col[t] = R.t
  r = a*R.t/S[t]
  S[t+1] = S[t] - r*I[t]*S[t]
  I[t+1] = I[t] + r*I[t]*S[t] - a*I[t]
  R[t+1] = R[t] + a*I[t]
  t = t+1
}
plot(model.df$susceptible~model.df$date)
points(S~model.df$date,col = 'red')
plot(model.df$active.infected~model.df$date)
points(I~model.df$date,col = 'red')
plot(model.df$active.infected,xlim = c(1,t.end),
     ylim = c(model.df$active.infected[1],model.df$active.infected[t.end+1]))
points(I,col = 'red')
plot(model.df$total.removed~model.df$date)
points(R~model.df$date,col = 'red')
t.end
S[t.end]
I[t.end]
R[t.end]
model.df$active.infected[t.end]
R.t
model.df$date[t.end]

#Then we handle a time period with a constant R.t
t.end = t+29
r = a*R.t/S[t]
while(t < t.end){
  R.rat.t[t] = 1
  R.t.col[t] = R.t
  S[t+1] = S[t] - r*I[t]*S[t]
  I[t+1] = I[t] + r*I[t]*S[t] - a*I[t]
  R[t+1] = R[t] + a*I[t]
  t = t+1
}
plot(model.df$susceptible~model.df$date)
points(S~model.df$date,col = 'red')
plot(model.df$active.infected~model.df$date)
points(I~model.df$date,col = 'red')
plot(model.df$active.infected,xlim = c(1,t.end),
     ylim = c(model.df$active.infected[1],model.df$active.infected[t.end+1]))
points(I,col = 'red')
plot(model.df$total.removed~model.df$date)
points(R~model.df$date,col = 'red')
t.end
S[t.end]
I[t.end]
R[t.end]
model.df$active.infected[t.end]
R.t
model.df$date[t.end]

R.rat = .956
t.end = t + lag
while(t < t.end){
  R.rat.t[t] = R.rat
  R.t = R.t*R.rat
  R.t.col[t] = R.t
  r = a*R.t/S[t]
  S[t+1] = S[t] - r*I[t]*S[t]
  I[t+1] = I[t] + r*I[t]*S[t] - a*I[t]
  R[t+1] = R[t] + a*I[t]
  t = t+1
}
plot(model.df$susceptible~model.df$date)
points(S~model.df$date,col = 'red')
plot(model.df$active.infected~model.df$date)
points(I~model.df$date,col = 'red')
plot(model.df$active.infected,xlim = c(1,t.end),
     ylim = c(model.df$active.infected[1],model.df$active.infected[t.end+1]+1500))
points(I,col = 'red')
plot(model.df$total.removed~model.df$date)
points(R~model.df$date,col = 'red')
t.end
S[t.end]
I[t.end]
R[t.end]
model.df$active.infected[t.end]
R.t
model.df$date[t.end]

a = 1/17

#Then we handle a time period with a constant R.t
t.end = t+19
r = a*R.t/S[t]
while(t < t.end){
  R.rat.t[t] = 1
  R.t.col[t] = R.t
  S[t+1] = S[t] - r*I[t]*S[t]
  I[t+1] = I[t] + r*I[t]*S[t] - a*I[t]
  R[t+1] = R[t] + a*I[t]
  t = t+1
}
plot(model.df$susceptible~model.df$date)
points(S~model.df$date,col = 'red')
plot(model.df$active.infected~model.df$date)
points(I~model.df$date,col = 'red')
plot(model.df$active.infected,xlim = c(1,t.end),
     ylim = c(model.df$active.infected[1],model.df$active.infected[t.end+1]+5000))
points(I,col = 'red')
plot(model.df$total.removed~model.df$date)
points(R~model.df$date,col = 'red')
t.end
S[t.end]
I[t.end]
R[t.end]
model.df$active.infected[t.end]
R.t
model.df$date[t.end]

a = 1/24

#Most schools opened within ~12 days from the end of the last period
R.rat = 1.037
t.end = t + 13 + lag
while(t < t.end){
  R.rat.t[t] = R.rat
  R.t = R.t*R.rat
  R.t.col[t] = R.t
  r = a*R.t/S[t]
  S[t+1] = S[t] - r*I[t]*S[t]
  I[t+1] = I[t] + r*I[t]*S[t] - a*I[t]
  R[t+1] = R[t] + a*I[t]
  t = t+1
}
plot(model.df$susceptible~model.df$date)
points(S~model.df$date,col = 'red')
plot(model.df$active.infected~model.df$date)
points(I~model.df$date,col = 'red')
plot(model.df$active.infected,xlim = c(1,t.end),
     ylim = c(model.df$active.infected[1],model.df$active.infected[t.end+1]+2500))
points(I,col = 'red')
plot(model.df$total.removed~model.df$date)
points(R~model.df$date,col = 'red')
t.end
S[t.end]
I[t.end]
R[t.end]
model.df$active.infected[t.end]
R.t
model.df$date[t.end]

#Varying R.t
R.rat = .966
t.end = t +  lag
while(t < t.end){
  R.rat.t[t] = R.rat
  R.t = R.t*R.rat
  R.t.col[t] = R.t
  r = a*R.t/S[t]
  S[t+1] = S[t] - r*I[t]*S[t]
  I[t+1] = I[t] + r*I[t]*S[t] - a*I[t]
  R[t+1] = R[t] + a*I[t]
  t = t+1
}
plot(model.df$susceptible~model.df$date)
points(S~model.df$date,col = 'red')
plot(model.df$active.infected~model.df$date)
points(I~model.df$date,col = 'red')
plot(model.df$active.infected,xlim = c(1,t.end),
     ylim = c(model.df$active.infected[1],model.df$active.infected[t.end+1]+2500))
points(I,col = 'red')
plot(model.df$total.removed~model.df$date)
points(R~model.df$date,col = 'red')
t.end
S[t.end]
I[t.end]
R[t.end]
model.df$active.infected[t.end]
R.t
model.df$date[t.end]

#Varying R.t
R.rat = 1.032
t.end = t +  16
while(t < t.end){
  R.rat.t[t] = R.rat
  R.t = R.t*R.rat
  R.t.col[t] = R.t
  r = a*R.t/S[t]
  S[t+1] = S[t] - r*I[t]*S[t]
  I[t+1] = I[t] + r*I[t]*S[t] - a*I[t]
  R[t+1] = R[t] + a*I[t]
  t = t+1
}
plot(model.df$susceptible~model.df$date)
points(S~model.df$date,col = 'red')
plot(model.df$active.infected~model.df$date)
points(I~model.df$date,col = 'red')
plot(model.df$active.infected,xlim = c(1,t.end),
     ylim = c(model.df$active.infected[1],model.df$active.infected[t.end]))
points(I,col = 'red')
plot(model.df$total.removed~model.df$date)
points(R~model.df$date,col = 'red')
t.end
S[t.end]
I[t.end]
R[t.end]
model.df$active.infected[t.end]
R.t
model.df$date[t.end]

R.rat.t[t] = R.rat
R.t.col[t] = R.t*R.rat

plot(R.rat.t~model.df$date)
lines(R.rat.t~model.df$date)
hist(R.rat.t)
boxplot(R.rat.t)
summary(R.rat.t)

plot(R.t.col~model.df$date)
lines(R.t.col~model.df$date)
plot(log(R.t.col)~model.df$date)
lines(log(R.t.col)~model.df$date)
hist(R.t.col)
hist(log(R.t.col))
boxplot(R.t.col)
boxplot(log(R.t.col))
summary(R.t.col)
summary(log(R.t.col))

plot(model.df$susceptible~model.df$date,xlab = "Date",
     ylab = "Susceptible",main = "Data & Model Comparison (Susceptible)",
     cex=.75,lty=1)
points(S~model.df$date,col = 'red',cex=.75)
plot(model.df$active.infected~model.df$date,xlab = "Date",
     ylab = "Infected",main = "Data & Model Comparison (Infected)",
     cex=.75,lty=1)
points(I~model.df$date,col = 'red',cex=.75)
plot(model.df$total.removed~model.df$date,xlab = "Date",
     ylab = "Removed",main = "Data & Model Comparison (Removed)",
     cex=.75,lty=1)
points(R~model.df$date,col = 'red',cex=.75)

#Dates:
#5-1 to 5-15 varying (decreasing, R.rat = .99, partial lag = 14)
#5-15 to 6-7 varying (increasing, R.rat = 1.02, full lag)
#6-7 to 7-6 constant (>1)
#7-6 to 7-29 varying (decreasing, R.rat = .956, full lag)
    #a increases from 1/25 to 1/17
#7-29 to 8-17 constant (<1)
    #a decreases from 1/17 to 1/24
#8-17 to 9-22 varying (increasing, R.rat = 1.037, full lag + 13)
#9-22 to 10-15 varying (decreasing, R.rat = .966, full lag)
#10-15 to 10-31 varying (increasing, R.rat = 1.032, partial lag = 16)

#Investigating patterns in R.rat and R.t
#shapiro.test(R.rat.t)
#shapiro.test(R.t.col)
#shapiro.test(log(R.t.col))
#Nothing interesting here

#PROJECTIONS
#First projection out to line 890
p.t.end = 50
x = 1:p.t.end

#Central Projection
R.t.1 = 1:p.t.end
R.t.1a = 1:p.t.end
R.t.1b = 1:p.t.end
R.t.1c = 1:p.t.end
R.rat = 1
R.t.1[1] = R.rat*R.t.col[t]
R.t.1a[1] = R.rat*R.t.col[t]
R.t.1b[1] = R.rat*R.t.col[t]
R.t.1c[1] = R.rat*R.t.col[t]
i=1
while(i<8){
  R.t.1[i+1] = R.rat*R.t.1[i]
  R.t.1a[i+1] = R.rat*R.t.1a[i]
  R.t.1b[i+1] = R.rat*R.t.1b[i]
  R.t.1c[i+1] = R.rat*R.t.1c[i]
  i = i+1
}
i1 = i
i.lag = i+lag+13
R.rat.a = .99
R.rat.b = .975
R.rat.c = .96
while(i1 < p.t.end){
  R.t.1[i1+1] = R.t.1[i1]
  i1 = i1+1
}
while(i < i.lag){
  R.t.1a[i+1] = R.rat.a*R.t.1a[i]
  R.t.1b[i+1] = R.rat.b*R.t.1b[i]
  R.t.1c[i+1] = R.rat.c*R.t.1c[i]
  i = i+1
}
while(i < p.t.end){
  R.t.1a[i+1] = R.t.1a[i]
  R.t.1b[i+1] = R.t.1b[i]
  R.t.1c[i+1] = R.t.1c[i]
  i = i+1
}
#plot(R.t.1)
#plot(R.t.1a)
#plot(R.t.1b)
#plot(R.t.1c)
S1 = 1:p.t.end
S1[1] = S[t]
S1a = 1:p.t.end
S1a[1] = S[t]
S1b = 1:p.t.end
S1b[1] = S[t]
S1c = 1:p.t.end
S1c[1] = S[t]
I1 = 1:p.t.end
I1[1] = I[t]
I1a = 1:p.t.end
I1a[1] = I[t]
I1b = 1:p.t.end
I1b[1] = I[t]
I1c = 1:p.t.end
I1c[1] = I[t]
R1 = 1:p.t.end
R1[1] = R[t]
R1a = 1:p.t.end
R1a[1] = R[t]
R1b = 1:p.t.end
R1b[1] = R[t]
R1c = 1:p.t.end
R1c[1] = R[t]
p.t = 1
S1.ref = S1[1]
S1a.ref = S1a[1]
S1b.ref = S1b[1]
S1c.ref = S1c[1]
while(p.t < p.t.end){
  if(R.t.1[p.t+1] == R.t.1[p.t]){
    r1 = a*R.t.1[p.t]/S1.ref
  } else {
    S1.ref = S1[p.t]
    r1 = a*R.t.1[p.t]/S1.ref
  }
  if(R.t.1a[p.t+1] == R.t.1a[p.t]){
    r1a = a*R.t.1a[p.t]/S1a.ref
  } else {
    S1a.ref = S1a[p.t]
    r1a = a*R.t.1a[p.t]/S1a.ref
  }
  if(R.t.1b[p.t+1] == R.t.1b[p.t]){
    r1b = a*R.t.1b[p.t]/S1b.ref
  } else {
    S1b.ref = S1b[p.t]
    r1b = a*R.t.1b[p.t]/S1b.ref
  }
  if(R.t.1c[p.t+1] == R.t.1c[p.t]){
    r1c = a*R.t.1c[p.t]/S1c.ref
  } else {
    S1c.ref = S1c[p.t]
    r1c = a*R.t.1c[p.t]/S1c.ref
  }
  S1[p.t+1] = S1[p.t] - r1*I1[p.t]*S1[p.t]
  S1a[p.t+1] = S1a[p.t] - r1a*I1a[p.t]*S1a[p.t]
  S1b[p.t+1] = S1b[p.t] - r1b*I1b[p.t]*S1b[p.t]
  S1c[p.t+1] = S1c[p.t] - r1c*I1c[p.t]*S1c[p.t]
  I1[p.t+1] = I1[p.t] + r1*I1[p.t]*S1[p.t] - a*I1[p.t]
  I1a[p.t+1] = I1a[p.t] + r1a*I1a[p.t]*S1a[p.t] - a*I1a[p.t]
  I1b[p.t+1] = I1b[p.t] + r1b*I1b[p.t]*S1b[p.t] - a*I1b[p.t]
  I1c[p.t+1] = I1c[p.t] + r1c*I1c[p.t]*S1c[p.t] - a*I1c[p.t]
  R1[p.t+1] = R1[p.t] + a*I1[p.t]
  R1a[p.t+1] = R1a[p.t] + a*I1a[p.t]
  R1b[p.t+1] = R1b[p.t] + a*I1b[p.t]
  R1c[p.t+1] = R1c[p.t] + a*I1c[p.t]
  p.t = p.t+1
}
plot(x, S1, type="o", col="red", cex=.75, lty=1, ylim=c(2825829,3081985),
     xlab = "Days since Oct. 31", ylab = "Susceptible",
     main = "Central Projections (Susceptible)")
max(S1)
min(S1)
S1[p.t.end]
points(x, S1a, col="purple",cex=.75)
lines(x, S1a, col="purple",lty=1)
max(S1a)
min(S1a)
S1a[p.t.end]
points(x, S1b, col="blue",cex=.75)
lines(x, S1b, col="blue",lty=1)
max(S1b)
min(S1b)
S1b[p.t.end]
points(x, S1c, col="green",cex=.75)
lines(x, S1c, col="green",lty=1)
max(S1c)
min(S1c)
S1c[p.t.end]
plot(x, I1, type="o", col="red", cex=.75, lty=1, ylim=c(30463.87,139218.1),
     xlab = "Days since Oct. 31", ylab = "Infected",
     main = "Central Projections (Infected)")
max(I1)
min(I1)
I1[p.t.end]
points(x, I1a, col="purple",cex=.75)
lines(x, I1a, col="purple",lty=1)
max(I1a)
min(I1a)
I1a[p.t.end]
points(x, I1b, col="blue",cex=.75)
lines(x, I1b, col="blue",lty=1)
max(I1b)
min(I1b)
I1b[p.t.end]
points(x, I1c, col="green",cex=.75)
lines(x, I1c, col="green",lty=1)
max(I1c)
min(I1c)
I1c[p.t.end]
plot(x, R1, type="o", col="red", cex=.75, lty=1, ylim=c(87551.33,234952.6),
     xlab = "Days since Oct. 31", ylab = "Removed",
     main = "Central Projections (Removed)")
max(R1)
min(R1)
R1[p.t.end]
points(x, R1a, col="purple",cex=.75)
lines(x, R1a, col="purple",lty=1)
max(R1a)
min(R1a)
R1a[p.t.end]
points(x, R1b, col="blue",cex=.75)
lines(x, R1b, col="blue",lty=1)
max(R1b)
min(R1b)
R1b[p.t.end]
points(x, R1c, col="green",cex=.75)
lines(x, R1c, col="green",lty=1)
max(R1c)
min(R1c)
R1c[p.t.end]

#Upper Projection
R.t.2 = 1:p.t.end
R.t.2a = 1:p.t.end
R.t.2b = 1:p.t.end
R.t.2c = 1:p.t.end
R.rat = 1.02
R.t.2[1] = R.rat*R.t.col[t]
R.t.2a[1] = R.rat*R.t.col[t]
R.t.2b[1] = R.rat*R.t.col[t]
R.t.2c[1] = R.rat*R.t.col[t]
i=1
while(i<8){
  R.t.2[i+1] = R.rat*R.t.2[i]
  R.t.2a[i+1] = R.rat*R.t.2a[i]
  R.t.2b[i+1] = R.rat*R.t.2b[i]
  R.t.2c[i+1] = R.rat*R.t.2c[i]
  i = i+1
}
i2 = i
i.lag = i+lag+13
R.rat.a = .99
R.rat.b = .975
R.rat.c = .96
while(i2 < p.t.end){
  R.t.2[i2+1] = R.t.2[i2]
  i2 = i2+1
}
while(i < i.lag){
  R.t.2a[i+1] = R.rat.a*R.t.2a[i]
  R.t.2b[i+1] = R.rat.b*R.t.2b[i]
  R.t.2c[i+1] = R.rat.c*R.t.2c[i]
  i = i+1
}
while(i < p.t.end){
  R.t.2a[i+1] = R.t.2a[i]
  R.t.2b[i+1] = R.t.2b[i]
  R.t.2c[i+1] = R.t.2c[i]
  i = i+1
}
#plot(R.t.2)
#plot(R.t.2a)
#plot(R.t.2b)
#plot(R.t.2c)
S2 = 1:p.t.end
S2[1] = S[t]
S2a = 1:p.t.end
S2a[1] = S[t]
S2b = 1:p.t.end
S2b[1] = S[t]
S2c = 1:p.t.end
S2c[1] = S[t]
I2 = 1:p.t.end
I2[1] = I[t]
I2a = 1:p.t.end
I2a[1] = I[t]
I2b = 1:p.t.end
I2b[1] = I[t]
I2c = 1:p.t.end
I2c[1] = I[t]
R2 = 1:p.t.end
R2[1] = R[t]
R2a = 1:p.t.end
R2a[1] = R[t]
R2b = 1:p.t.end
R2b[1] = R[t]
R2c = 1:p.t.end
R2c[1] = R[t]
p.t = 1
S2.ref = S2[1]
S2a.ref = S2a[1]
S2b.ref = S2b[1]
S2c.ref = S2c[1]
while(p.t < p.t.end){
  if(R.t.2[p.t+1] == R.t.2[p.t]){
    r2 = a*R.t.2[p.t]/S2.ref
  } else {
    S2.ref = S2[p.t]
    r2 = a*R.t.2[p.t]/S2.ref
  }
  if(R.t.2a[p.t+1] == R.t.2a[p.t]){
    r2a = a*R.t.2a[p.t]/S2a.ref
  } else {
    S2a.ref = S2a[p.t]
    r2a = a*R.t.2a[p.t]/S2a.ref
  }
  if(R.t.2b[p.t+1] == R.t.2b[p.t]){
    r2b = a*R.t.2b[p.t]/S2b.ref
  } else {
    S2b.ref = S2b[p.t]
    r2b = a*R.t.2b[p.t]/S2b.ref
  }
  if(R.t.2c[p.t+1] == R.t.2c[p.t]){
    r2c = a*R.t.2c[p.t]/S2c.ref
  } else {
    S2c.ref = S2c[p.t]
    r2c = a*R.t.2c[p.t]/S2c.ref
  }
  S2[p.t+1] = S2[p.t] - r1*I2[p.t]*S2[p.t]
  S2a[p.t+1] = S2a[p.t] - r2a*I2a[p.t]*S2a[p.t]
  S2b[p.t+1] = S2b[p.t] - r2b*I2b[p.t]*S2b[p.t]
  S2c[p.t+1] = S2c[p.t] - r2c*I2c[p.t]*S2c[p.t]
  I2[p.t+1] = I2[p.t] + r2*I2[p.t]*S2[p.t] - a*I2[p.t]
  I2a[p.t+1] = I2a[p.t] + r2a*I2a[p.t]*S2a[p.t] - a*I2a[p.t]
  I2b[p.t+1] = I2b[p.t] + r2b*I2b[p.t]*S2b[p.t] - a*I2b[p.t]
  I2c[p.t+1] = I2c[p.t] + r2c*I2c[p.t]*S2c[p.t] - a*I2c[p.t]
  R2[p.t+1] = R2[p.t] + a*I2[p.t]
  R2a[p.t+1] = R2a[p.t] + a*I2a[p.t]
  R2b[p.t+1] = R2b[p.t] + a*I2b[p.t]
  R2c[p.t+1] = R2c[p.t] + a*I2c[p.t]
  p.t = p.t+1
}
plot(x, S2, type="o", col="red", cex=.75, lty=1, ylim=c(2729201,3081985),
     xlab = "Days since Oct. 31", ylab = "Susceptible",
     main = "Upper Projections (Susceptible)")
max(S2)
min(S2)
S2[p.t.end]
points(x, S2a, col="purple",cex=.75)
lines(x, S2a, col="purple",lty=1)
max(S2a)
min(S2a)
S2a[p.t.end]
points(x, S2b, col="blue",cex=.75)
lines(x, S2b, col="blue",lty=1)
max(S2b)
min(S2b)
S2b[p.t.end]
points(x, S2c, col="green",cex=.75)
lines(x, S2c, col="green",lty=1)
max(S2c)
min(S2c)
S2c[p.t.end]
plot(x, I2, type="o", col="red", cex=.75, lty=1, ylim=c(30463.87,237878.7),
     xlab = "Days since Oct. 31", ylab = "Infected",
    main = "Upper Projections (Infected)")
max(I2)
min(I2)
I2[p.t.end]
points(x, I2a, col="purple",cex=.75)
lines(x, I2a, col="purple",lty=1)
max(I2a)
min(I2a)
I2a[p.t.end]
points(x, I2b, col="blue",cex=.75)
lines(x, I2b, col="blue",lty=1)
max(I2b)
min(I2b)
I2b[p.t.end]
points(x, I2c, col="green",cex=.75)
lines(x, I2c, col="green",lty=1)
max(I2c)
min(I2c)
I2c[p.t.end]
plot(x, R2, type="o", col="red", cex=.75, lty=1, ylim=c(87551.33,293930.9),
     xlab = "Days since Oct. 31", ylab = "Removed",
     main = "Upper Projections (Removed)")
max(R2)
min(R2)
R2[p.t.end]
points(x, R2a, col="purple",cex=.75)
lines(x, R2a, col="purple",lty=1)
max(R2a)
min(R2a)
R2a[p.t.end]
points(x, R2b, col="blue",cex=.75)
lines(x, R2b, col="blue",lty=1)
max(R2b)
min(R2b)
R2b[p.t.end]
points(x, R2c, col="green",cex=.75)
lines(x, R2c, col="green",lty=1)
max(R2c)
min(R2c)
R2c[p.t.end]

#Lower Projection
R.t.3 = 1:p.t.end
R.t.3a = 1:p.t.end
R.t.3b = 1:p.t.end
R.t.3c = 1:p.t.end
R.rat = .98
R.t.3[1] = R.rat*R.t.col[t]
R.t.3a[1] = R.rat*R.t.col[t]
R.t.3b[1] = R.rat*R.t.col[t]
R.t.3c[1] = R.rat*R.t.col[t]
i=1
while(i<8){
  R.t.3[i+1] = R.rat*R.t.3[i]
  R.t.3a[i+1] = R.rat*R.t.3a[i]
  R.t.3b[i+1] = R.rat*R.t.3b[i]
  R.t.3c[i+1] = R.rat*R.t.3c[i]
  i = i+1
}
i3 = i
i.lag = i+lag+13
R.rat.a = .99
R.rat.b = .975
R.rat.c = .96
while(i3 < p.t.end){
  R.t.3[i3+1] = R.t.3[i3]
  i3 = i3+1
}
while(i < i.lag){
  R.t.3a[i+1] = R.rat.a*R.t.3a[i]
  R.t.3b[i+1] = R.rat.b*R.t.3b[i]
  R.t.3c[i+1] = R.rat.c*R.t.3c[i]
  i = i+1
}
while(i < p.t.end){
  R.t.3a[i+1] = R.t.3a[i]
  R.t.3b[i+1] = R.t.3b[i]
  R.t.3c[i+1] = R.t.3c[i]
  i = i+1
}
#plot(R.t.3)
#plot(R.t.3a)
#plot(R.t.3b)
#plot(R.t.3c)
S3 = 1:p.t.end
S3[1] = S[t]
S3a = 1:p.t.end
S3a[1] = S[t]
S3b = 1:p.t.end
S3b[1] = S[t]
S3c = 1:p.t.end
S3c[1] = S[t]
I3 = 1:p.t.end
I3[1] = I[t]
I3a = 1:p.t.end
I3a[1] = I[t]
I3b = 1:p.t.end
I3b[1] = I[t]
I3c = 1:p.t.end
I3c[1] = I[t]
R3 = 1:p.t.end
R3[1] = R[t]
R3a = 1:p.t.end
R3a[1] = R[t]
R3b = 1:p.t.end
R3b[1] = R[t]
R3c = 1:p.t.end
R3c[1] = R[t]
p.t = 1
S3.ref = S1[1]
S3a.ref = S1a[1]
S3b.ref = S1b[1]
S3c.ref = S1c[1]
while(p.t < p.t.end){
  if(R.t.3[p.t+1] == R.t.3[p.t]){
    r3 = a*R.t.3[p.t]/S3.ref
  } else {
    S3.ref = S3[p.t]
    r3 = a*R.t.3[p.t]/S3.ref
  }
  if(R.t.3a[p.t+1] == R.t.3a[p.t]){
    r3a = a*R.t.3a[p.t]/S3a.ref
  } else {
    S3a.ref = S3a[p.t]
    r3a = a*R.t.3a[p.t]/S3a.ref
  }
  if(R.t.3b[p.t+1] == R.t.3b[p.t]){
    r3b = a*R.t.3b[p.t]/S3b.ref
  } else {
    S3b.ref = S3b[p.t]
    r3b = a*R.t.3b[p.t]/S3b.ref
  }
  if(R.t.3c[p.t+1] == R.t.3c[p.t]){
    r3c = a*R.t.3c[p.t]/S3c.ref
  } else {
    S3c.ref = S3c[p.t]
    r3c = a*R.t.3c[p.t]/S3c.ref
  }
  S3[p.t+1] = S3[p.t] - r1*I3[p.t]*S3[p.t]
  S3a[p.t+1] = S3a[p.t] - r3a*I3a[p.t]*S3a[p.t]
  S3b[p.t+1] = S3b[p.t] - r3b*I3b[p.t]*S3b[p.t]
  S3c[p.t+1] = S3c[p.t] - r3c*I3c[p.t]*S3c[p.t]
  I3[p.t+1] = I3[p.t] + r3*I3[p.t]*S3[p.t] - a*I3[p.t]
  I3a[p.t+1] = I3a[p.t] + r3a*I3a[p.t]*S3a[p.t] - a*I3a[p.t]
  I3b[p.t+1] = I3b[p.t] + r3b*I3b[p.t]*S3b[p.t] - a*I3b[p.t]
  I3c[p.t+1] = I3c[p.t] + r3c*I3c[p.t]*S3c[p.t] - a*I3c[p.t]
  R3[p.t+1] = R3[p.t] + a*I3[p.t]
  R3a[p.t+1] = R3a[p.t] + a*I3a[p.t]
  R3b[p.t+1] = R3b[p.t] + a*I3b[p.t]
  R3c[p.t+1] = R3c[p.t] + a*I3c[p.t]
  p.t = p.t+1
}
plot(x, S3, type="o", col="red", cex=.75, lty=1, ylim=c(2883509,3081985),
     xlab = "Days since Oct. 31", ylab = "Susceptible",
     main = "Lower Projections (Susceptible)")
max(S3)
min(S3)
S3[p.t.end]
points(x, S3a, col="purple",cex=.75)
lines(x, S3a, col="purple",lty=1)
max(S3a)
min(S3a)
S3a[p.t.end]
points(x, S3b, col="blue",cex=.75)
lines(x, S3b, col="blue",lty=1)
max(S3b)
min(S3b)
S3b[p.t.end]
points(x, S3c, col="green",cex=.75)
lines(x, S3c, col="green",lty=1)
max(S3c)
min(S3c)
S3c[p.t.end]
plot(x, I3, type="o", col="red", cex=.75, lty=1, ylim=c(24012.24,88228.36),
     xlab = "Days since Oct. 31", ylab = "Infected",
     main = "Lower Projections (Infected)")
max(I3)
min(I3)
I3[p.t.end]
points(x, I3a, col="purple",cex=.75)
lines(x, I3a, col="purple",lty=1)
max(I3a)
min(I3a)
I3a[p.t.end]
points(x, I3b, col="blue",cex=.75)
lines(x, I3b, col="blue",lty=1)
max(I3b)
min(I3b)
I3b[p.t.end]
points(x, I3c, col="green",cex=.75)
lines(x, I3c, col="green",lty=1)
max(I3c)
min(I3c)
I3c[p.t.end]
plot(x, R3, type="o", col="red", cex=.75, lty=1, ylim=c(87551.33,200662.9),
     xlab = "Days since Oct. 31", ylab = "Removed",
     main = "Lower Projections (Removed)")
max(R3)
min(R3)
R3[p.t.end]
points(x, R3a, col="purple",cex=.75)
lines(x, R3a, col="purple",lty=1)
max(R3a)
min(R3a)
R3a[p.t.end]
points(x, R3b, col="blue",cex=.75)
lines(x, R3b, col="blue",lty=1)
max(R3b)
min(R3b)
R3b[p.t.end]
points(x, R3c, col="green",cex=.75)
lines(x, R3c, col="green",lty=1)
max(R3c)
min(R3c)
R3c[p.t.end]




p.t.end = 100
x = 1:p.t.end

#Central Projection
R.t.1 = 1:p.t.end
R.t.1a = 1:p.t.end
R.t.1b = 1:p.t.end
R.t.1c = 1:p.t.end
R.rat = 1
R.t.1[1] = R.rat*R.t.col[t]
R.t.1a[1] = R.rat*R.t.col[t]
R.t.1b[1] = R.rat*R.t.col[t]
R.t.1c[1] = R.rat*R.t.col[t]
i=1
while(i<8){
  R.t.1[i+1] = R.rat*R.t.1[i]
  R.t.1a[i+1] = R.rat*R.t.1a[i]
  R.t.1b[i+1] = R.rat*R.t.1b[i]
  R.t.1c[i+1] = R.rat*R.t.1c[i]
  i = i+1
}
i1 = i
i.lag = i+lag+13
R.rat.a = .99
R.rat.b = .975
R.rat.c = .96
while(i1 < p.t.end){
  R.t.1[i1+1] = R.t.1[i1]
  i1 = i1+1
}
while(i < i.lag){
  R.t.1a[i+1] = R.rat.a*R.t.1a[i]
  R.t.1b[i+1] = R.rat.b*R.t.1b[i]
  R.t.1c[i+1] = R.rat.c*R.t.1c[i]
  i = i+1
}
while(i < p.t.end){
  R.t.1a[i+1] = R.t.1a[i]
  R.t.1b[i+1] = R.t.1b[i]
  R.t.1c[i+1] = R.t.1c[i]
  i = i+1
}
#plot(R.t.1)
#plot(R.t.1a)
#plot(R.t.1b)
#plot(R.t.1c)
S1 = 1:p.t.end
S1[1] = S[t]
S1a = 1:p.t.end
S1a[1] = S[t]
S1b = 1:p.t.end
S1b[1] = S[t]
S1c = 1:p.t.end
S1c[1] = S[t]
I1 = 1:p.t.end
I1[1] = I[t]
I1a = 1:p.t.end
I1a[1] = I[t]
I1b = 1:p.t.end
I1b[1] = I[t]
I1c = 1:p.t.end
I1c[1] = I[t]
R1 = 1:p.t.end
R1[1] = R[t]
R1a = 1:p.t.end
R1a[1] = R[t]
R1b = 1:p.t.end
R1b[1] = R[t]
R1c = 1:p.t.end
R1c[1] = R[t]
p.t = 1
S1.ref = S1[1]
S1a.ref = S1a[1]
S1b.ref = S1b[1]
S1c.ref = S1c[1]
while(p.t < p.t.end){
  if(R.t.1[p.t+1] == R.t.1[p.t]){
    r1 = a*R.t.1[p.t]/S1.ref
  } else {
    S1.ref = S1[p.t]
    r1 = a*R.t.1[p.t]/S1.ref
  }
  if(R.t.1a[p.t+1] == R.t.1a[p.t]){
    r1a = a*R.t.1a[p.t]/S1a.ref
  } else {
    S1a.ref = S1a[p.t]
    r1a = a*R.t.1a[p.t]/S1a.ref
  }
  if(R.t.1b[p.t+1] == R.t.1b[p.t]){
    r1b = a*R.t.1b[p.t]/S1b.ref
  } else {
    S1b.ref = S1b[p.t]
    r1b = a*R.t.1b[p.t]/S1b.ref
  }
  if(R.t.1c[p.t+1] == R.t.1c[p.t]){
    r1c = a*R.t.1c[p.t]/S1c.ref
  } else {
    S1c.ref = S1c[p.t]
    r1c = a*R.t.1c[p.t]/S1c.ref
  }
  S1[p.t+1] = S1[p.t] - r1*I1[p.t]*S1[p.t]
  S1a[p.t+1] = S1a[p.t] - r1a*I1a[p.t]*S1a[p.t]
  S1b[p.t+1] = S1b[p.t] - r1b*I1b[p.t]*S1b[p.t]
  S1c[p.t+1] = S1c[p.t] - r1c*I1c[p.t]*S1c[p.t]
  I1[p.t+1] = I1[p.t] + r1*I1[p.t]*S1[p.t] - a*I1[p.t]
  I1a[p.t+1] = I1a[p.t] + r1a*I1a[p.t]*S1a[p.t] - a*I1a[p.t]
  I1b[p.t+1] = I1b[p.t] + r1b*I1b[p.t]*S1b[p.t] - a*I1b[p.t]
  I1c[p.t+1] = I1c[p.t] + r1c*I1c[p.t]*S1c[p.t] - a*I1c[p.t]
  R1[p.t+1] = R1[p.t] + a*I1[p.t]
  R1a[p.t+1] = R1a[p.t] + a*I1a[p.t]
  R1b[p.t+1] = R1b[p.t] + a*I1b[p.t]
  R1c[p.t+1] = R1c[p.t] + a*I1c[p.t]
  p.t = p.t+1
}
plot(x, S1, type="o", col="red", cex=.75, lty=1, ylim=c(2075558,3081985),
     xlab = "Days since Oct. 31", ylab = "Susceptible",
     main = "Central Projections (Susceptible)")
max(S1)
min(S1)
S1[p.t.end]
points(x, S1a, col="purple",cex=.75)
lines(x, S1a, col="purple",lty=1)
max(S1a)
min(S1a)
S1a[p.t.end]
points(x, S1b, col="blue",cex=.75)
lines(x, S1b, col="blue",lty=1)
max(S1b)
min(S1b)
S1b[p.t.end]
points(x, S1c, col="green",cex=.75)
lines(x, S1c, col="green",lty=1)
max(S1c)
min(S1c)
S1c[p.t.end]
plot(x, I1, type="o", col="red", cex=.75, lty=1, ylim=c(8652.792,366498.8),
     xlab = "Days since Oct. 31", ylab = "Infected",
     main = "Central Projections (Infected)")
max(I1)
min(I1)
I1[p.t.end]
points(x, I1a, col="purple",cex=.75)
lines(x, I1a, col="purple",lty=1)
max(I1a)
min(I1a)
I1a[p.t.end]
points(x, I1b, col="blue",cex=.75)
lines(x, I1b, col="blue",lty=1)
max(I1b)
min(I1b)
I1b[p.t.end]
points(x, I1c, col="green",cex=.75)
lines(x, I1c, col="green",lty=1)
max(I1c)
min(I1c)
I1c[p.t.end]
plot(x, R1, type="o", col="red", cex=.75, lty=1, ylim=c(87551.33,757943),
     xlab = "Days since Oct. 31", ylab = "Removed",
     main = "Central Projections (Removed)")
max(R1)
min(R1)
R1[p.t.end]
points(x, R1a, col="purple",cex=.75)
lines(x, R1a, col="purple",lty=1)
max(R1a)
min(R1a)
R1a[p.t.end]
points(x, R1b, col="blue",cex=.75)
lines(x, R1b, col="blue",lty=1)
max(R1b)
min(R1b)
R1b[p.t.end]
points(x, R1c, col="green",cex=.75)
lines(x, R1c, col="green",lty=1)
max(R1c)
min(R1c)
R1c[p.t.end]

#Upper Projection
R.t.2 = 1:p.t.end
R.t.2a = 1:p.t.end
R.t.2b = 1:p.t.end
R.t.2c = 1:p.t.end
R.rat = 1.02
R.t.2[1] = R.rat*R.t.col[t]
R.t.2a[1] = R.rat*R.t.col[t]
R.t.2b[1] = R.rat*R.t.col[t]
R.t.2c[1] = R.rat*R.t.col[t]
i=1
while(i<8){
  R.t.2[i+1] = R.rat*R.t.2[i]
  R.t.2a[i+1] = R.rat*R.t.2a[i]
  R.t.2b[i+1] = R.rat*R.t.2b[i]
  R.t.2c[i+1] = R.rat*R.t.2c[i]
  i = i+1
}
i2 = i
i.lag = i+lag+13
R.rat.a = .99
R.rat.b = .975
R.rat.c = .96
while(i2 < p.t.end){
  R.t.2[i2+1] = R.t.2[i2]
  i2 = i2+1
}
while(i < i.lag){
  R.t.2a[i+1] = R.rat.a*R.t.2a[i]
  R.t.2b[i+1] = R.rat.b*R.t.2b[i]
  R.t.2c[i+1] = R.rat.c*R.t.2c[i]
  i = i+1
}
while(i < p.t.end){
  R.t.2a[i+1] = R.t.2a[i]
  R.t.2b[i+1] = R.t.2b[i]
  R.t.2c[i+1] = R.t.2c[i]
  i = i+1
}
#plot(R.t.2)
#plot(R.t.2a)
#plot(R.t.2b)
#plot(R.t.2c)
S2 = 1:p.t.end
S2[1] = S[t]
S2a = 1:p.t.end
S2a[1] = S[t]
S2b = 1:p.t.end
S2b[1] = S[t]
S2c = 1:p.t.end
S2c[1] = S[t]
I2 = 1:p.t.end
I2[1] = I[t]
I2a = 1:p.t.end
I2a[1] = I[t]
I2b = 1:p.t.end
I2b[1] = I[t]
I2c = 1:p.t.end
I2c[1] = I[t]
R2 = 1:p.t.end
R2[1] = R[t]
R2a = 1:p.t.end
R2a[1] = R[t]
R2b = 1:p.t.end
R2b[1] = R[t]
R2c = 1:p.t.end
R2c[1] = R[t]
p.t = 1
S2.ref = S2[1]
S2a.ref = S2a[1]
S2b.ref = S2b[1]
S2c.ref = S2c[1]
while(p.t < p.t.end){
  if(R.t.2[p.t+1] == R.t.2[p.t]){
    r2 = a*R.t.2[p.t]/S2.ref
  } else {
    S2.ref = S2[p.t]
    r2 = a*R.t.2[p.t]/S2.ref
  }
  if(R.t.2a[p.t+1] == R.t.2a[p.t]){
    r2a = a*R.t.2a[p.t]/S2a.ref
  } else {
    S2a.ref = S2a[p.t]
    r2a = a*R.t.2a[p.t]/S2a.ref
  }
  if(R.t.2b[p.t+1] == R.t.2b[p.t]){
    r2b = a*R.t.2b[p.t]/S2b.ref
  } else {
    S2b.ref = S2b[p.t]
    r2b = a*R.t.2b[p.t]/S2b.ref
  }
  if(R.t.2c[p.t+1] == R.t.2c[p.t]){
    r2c = a*R.t.2c[p.t]/S2c.ref
  } else {
    S2c.ref = S2c[p.t]
    r2c = a*R.t.2c[p.t]/S2c.ref
  }
  S2[p.t+1] = S2[p.t] - r1*I2[p.t]*S2[p.t]
  S2a[p.t+1] = S2a[p.t] - r2a*I2a[p.t]*S2a[p.t]
  S2b[p.t+1] = S2b[p.t] - r2b*I2b[p.t]*S2b[p.t]
  S2c[p.t+1] = S2c[p.t] - r2c*I2c[p.t]*S2c[p.t]
  I2[p.t+1] = I2[p.t] + r2*I2[p.t]*S2[p.t] - a*I2[p.t]
  I2a[p.t+1] = I2a[p.t] + r2a*I2a[p.t]*S2a[p.t] - a*I2a[p.t]
  I2b[p.t+1] = I2b[p.t] + r2b*I2b[p.t]*S2b[p.t] - a*I2b[p.t]
  I2c[p.t+1] = I2c[p.t] + r2c*I2c[p.t]*S2c[p.t] - a*I2c[p.t]
  R2[p.t+1] = R2[p.t] + a*I2[p.t]
  R2a[p.t+1] = R2a[p.t] + a*I2a[p.t]
  R2b[p.t+1] = R2b[p.t] + a*I2b[p.t]
  R2c[p.t+1] = R2c[p.t] + a*I2c[p.t]
  p.t = p.t+1
}
plot(x, S2, type="o", col="red", cex=.75, lty=1, ylim=c(1507013,3081985),
     xlab = "Days since Oct. 31", ylab = "Susceptible",
     main = "Upper Projections (Susceptible)")
max(S2)
min(S2)
S2[p.t.end]
points(x, S2a, col="purple",cex=.75)
lines(x, S2a, col="purple",lty=1)
max(S2a)
min(S2a)
S2a[p.t.end]
points(x, S2b, col="blue",cex=.75)
lines(x, S2b, col="blue",lty=1)
max(S2b)
min(S2b)
S2b[p.t.end]
points(x, S2c, col="green",cex=.75)
lines(x, S2c, col="green",lty=1)
max(S2c)
min(S2c)
S2c[p.t.end]
plot(x, I2, type="o", col="red", cex=.75, lty=1, ylim=c(14487.73,673512.7),
     xlab = "Days since Oct. 31", ylab = "Infected",
     main = "Upper Projections (Infected)")
max(I2)
min(I2)
I2[p.t.end]
points(x, I2a, col="purple",cex=.75)
lines(x, I2a, col="purple",lty=1)
max(I2a)
min(I2a)
I2a[p.t.end]
points(x, I2b, col="blue",cex=.75)
lines(x, I2b, col="blue",lty=1)
max(I2b)
min(I2b)
I2b[p.t.end]
points(x, I2c, col="green",cex=.75)
lines(x, I2c, col="green",lty=1)
max(I2c)
min(I2c)
I2c[p.t.end]
plot(x, R2, type="o", col="red", cex=.75, lty=1, ylim=c(87551.33,1297350),
     xlab = "Days since Oct. 31", ylab = "Removed",
     main = "Upper Projections (Removed)")
max(R2)
min(R2)
R2[p.t.end]
points(x, R2a, col="purple",cex=.75)
lines(x, R2a, col="purple",lty=1)
max(R2a)
min(R2a)
R2a[p.t.end]
points(x, R2b, col="blue",cex=.75)
lines(x, R2b, col="blue",lty=1)
max(R2b)
min(R2b)
R2b[p.t.end]
points(x, R2c, col="green",cex=.75)
lines(x, R2c, col="green",lty=1)
max(R2c)
min(R2c)
R2c[p.t.end]

#Lower Projection
R.t.3 = 1:p.t.end
R.t.3a = 1:p.t.end
R.t.3b = 1:p.t.end
R.t.3c = 1:p.t.end
R.rat = .98
R.t.3[1] = R.rat*R.t.col[t]
R.t.3a[1] = R.rat*R.t.col[t]
R.t.3b[1] = R.rat*R.t.col[t]
R.t.3c[1] = R.rat*R.t.col[t]
i=1
while(i<8){
  R.t.3[i+1] = R.rat*R.t.3[i]
  R.t.3a[i+1] = R.rat*R.t.3a[i]
  R.t.3b[i+1] = R.rat*R.t.3b[i]
  R.t.3c[i+1] = R.rat*R.t.3c[i]
  i = i+1
}
i3 = i
i.lag = i+lag+13
R.rat.a = .99
R.rat.b = .975
R.rat.c = .96
while(i3 < p.t.end){
  R.t.3[i3+1] = R.t.3[i3]
  i3 = i3+1
}
while(i < i.lag){
  R.t.3a[i+1] = R.rat.a*R.t.3a[i]
  R.t.3b[i+1] = R.rat.b*R.t.3b[i]
  R.t.3c[i+1] = R.rat.c*R.t.3c[i]
  i = i+1
}
while(i < p.t.end){
  R.t.3a[i+1] = R.t.3a[i]
  R.t.3b[i+1] = R.t.3b[i]
  R.t.3c[i+1] = R.t.3c[i]
  i = i+1
}
#plot(R.t.3)
#plot(R.t.3a)
#plot(R.t.3b)
#plot(R.t.3c)
S3 = 1:p.t.end
S3[1] = S[t]
S3a = 1:p.t.end
S3a[1] = S[t]
S3b = 1:p.t.end
S3b[1] = S[t]
S3c = 1:p.t.end
S3c[1] = S[t]
I3 = 1:p.t.end
I3[1] = I[t]
I3a = 1:p.t.end
I3a[1] = I[t]
I3b = 1:p.t.end
I3b[1] = I[t]
I3c = 1:p.t.end
I3c[1] = I[t]
R3 = 1:p.t.end
R3[1] = R[t]
R3a = 1:p.t.end
R3a[1] = R[t]
R3b = 1:p.t.end
R3b[1] = R[t]
R3c = 1:p.t.end
R3c[1] = R[t]
p.t = 1
S3.ref = S1[1]
S3a.ref = S1a[1]
S3b.ref = S1b[1]
S3c.ref = S1c[1]
while(p.t < p.t.end){
  if(R.t.3[p.t+1] == R.t.3[p.t]){
    r3 = a*R.t.3[p.t]/S3.ref
  } else {
    S3.ref = S3[p.t]
    r3 = a*R.t.3[p.t]/S3.ref
  }
  if(R.t.3a[p.t+1] == R.t.3a[p.t]){
    r3a = a*R.t.3a[p.t]/S3a.ref
  } else {
    S3a.ref = S3a[p.t]
    r3a = a*R.t.3a[p.t]/S3a.ref
  }
  if(R.t.3b[p.t+1] == R.t.3b[p.t]){
    r3b = a*R.t.3b[p.t]/S3b.ref
  } else {
    S3b.ref = S3b[p.t]
    r3b = a*R.t.3b[p.t]/S3b.ref
  }
  if(R.t.3c[p.t+1] == R.t.3c[p.t]){
    r3c = a*R.t.3c[p.t]/S3c.ref
  } else {
    S3c.ref = S3c[p.t]
    r3c = a*R.t.3c[p.t]/S3c.ref
  }
  S3[p.t+1] = S3[p.t] - r1*I3[p.t]*S3[p.t]
  S3a[p.t+1] = S3a[p.t] - r3a*I3a[p.t]*S3a[p.t]
  S3b[p.t+1] = S3b[p.t] - r3b*I3b[p.t]*S3b[p.t]
  S3c[p.t+1] = S3c[p.t] - r3c*I3c[p.t]*S3c[p.t]
  I3[p.t+1] = I3[p.t] + r3*I3[p.t]*S3[p.t] - a*I3[p.t]
  I3a[p.t+1] = I3a[p.t] + r3a*I3a[p.t]*S3a[p.t] - a*I3a[p.t]
  I3b[p.t+1] = I3b[p.t] + r3b*I3b[p.t]*S3b[p.t] - a*I3b[p.t]
  I3c[p.t+1] = I3c[p.t] + r3c*I3c[p.t]*S3c[p.t] - a*I3c[p.t]
  R3[p.t+1] = R3[p.t] + a*I3[p.t]
  R3a[p.t+1] = R3a[p.t] + a*I3a[p.t]
  R3b[p.t+1] = R3b[p.t] + a*I3b[p.t]
  R3c[p.t+1] = R3c[p.t] + a*I3c[p.t]
  p.t = p.t+1
}
plot(x, S3, type="o", col="red", cex=.75, lty=1, ylim=c(2446718,3081985),
     xlab = "Days since Oct. 31", ylab = "Susceptible",
     main = "Lower Projections (Susceptible)")
max(S3)
min(S3)
S3[p.t.end]
points(x, S3a, col="purple",cex=.75)
lines(x, S3a, col="purple",lty=1)
max(S3a)
min(S3a)
S3a[p.t.end]
points(x, S3b, col="blue",cex=.75)
lines(x, S3b, col="blue",lty=1)
max(S3b)
min(S3b)
S3b[p.t.end]
points(x, S3c, col="green",cex=.75)
lines(x, S3c, col="green",lty=1)
max(S3c)
min(S3c)
S3c[p.t.end]
plot(x, I3, type="o", col="red", cex=.75, lty=1, ylim=c(6132.136,182788.3),
     xlab = "Days since Oct. 31", ylab = "Infected",
     main = "Lower Projections (Infected)")
max(I3)
min(I3)
I3[p.t.end]
points(x, I3a, col="purple",cex=.75)
lines(x, I3a, col="purple",lty=1)
max(I3a)
min(I3a)
I3a[p.t.end]
points(x, I3b, col="blue",cex=.75)
lines(x, I3b, col="blue",lty=1)
max(I3b)
min(I3b)
I3b[p.t.end]
points(x, I3c, col="green",cex=.75)
lines(x, I3c, col="green",lty=1)
max(I3c)
min(I3c)
I3c[p.t.end]
plot(x, R3, type="o", col="red", cex=.75, lty=1, ylim=c(87551.33,479510.1),
     xlab = "Days since Oct. 31", ylab = "Removed",
     main = "Lower Projections (Removed)")
max(R3)
min(R3)
R3[p.t.end]
points(x, R3a, col="purple",cex=.75)
lines(x, R3a, col="purple",lty=1)
max(R3a)
min(R3a)
R3a[p.t.end]
points(x, R3b, col="blue",cex=.75)
lines(x, R3b, col="blue",lty=1)
max(R3b)
min(R3b)
R3b[p.t.end]
points(x, R3c, col="green",cex=.75)
lines(x, R3c, col="green",lty=1)
max(R3c)
min(R3c)
R3c[p.t.end]
