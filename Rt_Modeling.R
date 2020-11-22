#MODEL BUILDING

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

plot(log(model.df$susceptible)~model.df$date)
points(log(S)~model.df$date,col = 'red')
plot(log(model.df$active.infected)~model.df$date)
points(log(I)~model.df$date,col = 'red')
plot(log(model.df$total.removed)~model.df$date)
points(log(R)~model.df$date,col = 'red')

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
p.t.end = 100
x = 1:p.t.end

#Define p.R.t (projected R.t)
p.R.t = 1:p.t.end
R.rat = 1.032
p.R.t[1] = R.rat*R.t.col[t]
i=1
while(i<7){
  p.R.t[i+1] = R.rat*p.R.t[i]
  i = i+1
}
while(i<p.t.end){
  p.R.t[i+1] = p.R.t[i]
  i = i+1
}
#plot(p.R.t)
#1st Simulation: R.t changes through remainder of lag, then stays constant
p.S = 1:p.t.end
p.S[1] = S[t]
p.I = 1:p.t.end
p.I[1] = I[t]
p.R = 1:p.t.end
p.R[1] = R[t]
p.t = 1
S.ref = p.S[1]
while(p.t < p.t.end){
  if(p.R.t[p.t+1] == p.R.t[p.t]){
    r = a*p.R.t[p.t]/S.ref
  } else {
    S.ref = p.S[p.t]
    r = a*p.R.t[p.t]/S.ref
  }
  p.S[p.t+1] = p.S[p.t] - r*p.I[p.t]*p.S[p.t]
  p.I[p.t+1] = p.I[p.t] + r*p.I[p.t]*p.S[p.t] - a*p.I[p.t]
  p.R[p.t+1] = p.R[p.t] + a*p.I[p.t]
  p.t = p.t+1
}
plot(x, p.S, type="o", col="red", cex=.75, lty=1, ylim=c(1148446,3081985))
max(p.S)
min(p.S)
p.S[p.t.end]
plot(x,p.I,type="o", col="red", cex=.75, lty=1, ylim=c(30079.51,648454.2))
max(p.I)
min(p.I)
p.I[p.t.end]
plot(x,p.R,type="o", col="red", cex=.75, lty=1, ylim=c(87551.33,1418950))
max(p.R)
min(p.R)
p.R[p.t.end]

#2nd Simulation: Extreme decline following announcement
p.R.t = 1:p.t.end
R.rat = 1.032
p.R.t[1] = R.rat*R.t.col[t]
i=1
while(i<7){
  p.R.t[i+1] = R.rat*p.R.t[i]
  i = i+1
}
while(i<9){
  p.R.t[i+1] = p.R.t[i]
  i = i+1
}
R.rat = .956
i.end = i+lag
while(i<i.end){
  p.R.t[i+1] = R.rat*p.R.t[i]
  i = i+1
}
while(i<p.t.end){
  p.R.t[i+1] = p.R.t[i]
  i = i+1
}
plot(p.R.t)

p.S = 1:p.t.end
p.S[1] = S[t]
p.I = 1:p.t.end
p.I[1] = I[t]
p.R = 1:p.t.end
p.R[1] = R[t]
p.t = 1
S.ref = p.S[1]
while(p.t < p.t.end){
  if(p.R.t[p.t+1] == p.R.t[p.t]){
    r = a*p.R.t[p.t]/S.ref
  } else {
    S.ref = p.S[p.t]
    r = a*p.R.t[p.t]/S.ref
  }
  p.S[p.t+1] = p.S[p.t] - r*p.I[p.t]*p.S[p.t]
  p.I[p.t+1] = p.I[p.t] + r*p.I[p.t]*p.S[p.t] - a*p.I[p.t]
  p.R[p.t+1] = p.R[p.t] + a*p.I[p.t]
  p.t = p.t+1
}
plot(p.S)
plot(p.I)
plot(p.R)

#3rd Simulation: moderate decline following announcement
p.R.t = 1:p.t.end
R.rat = 1.032
p.R.t[1] = R.rat*R.t.col[t]
i=1
while(i<7){
  p.R.t[i+1] = R.rat*p.R.t[i]
  i = i+1
}
while(i<9){
  p.R.t[i+1] = p.R.t[i]
  i = i+1
}
R.rat = .966
i.end = i+lag
while(i<i.end){
  p.R.t[i+1] = R.rat*p.R.t[i]
  i = i+1
}
while(i<p.t.end){
  p.R.t[i+1] = p.R.t[i]
  i = i+1
}
plot(p.R.t)

p.S = 1:p.t.end
p.S[1] = S[t]
p.I = 1:p.t.end
p.I[1] = I[t]
p.R = 1:p.t.end
p.R[1] = R[t]
p.t = 1
S.ref = p.S[1]
while(p.t < p.t.end){
  if(p.R.t[p.t+1] == p.R.t[p.t]){
    r = a*p.R.t[p.t]/S.ref
  } else {
    S.ref = p.S[p.t]
    r = a*p.R.t[p.t]/S.ref
  }
  p.S[p.t+1] = p.S[p.t] - r*p.I[p.t]*p.S[p.t]
  p.I[p.t+1] = p.I[p.t] + r*p.I[p.t]*p.S[p.t] - a*p.I[p.t]
  p.R[p.t+1] = p.R[p.t] + a*p.I[p.t]
  p.t = p.t+1
}
plot(p.S)
plot(p.I)
plot(p.R)

#4th Simulation: minimal decline following announcement
p.R.t = 1:p.t.end
R.rat = 1.032
p.R.t[1] = R.rat*R.t.col[t]
i=1
while(i<7){
  p.R.t[i+1] = R.rat*p.R.t[i]
  i = i+1
}
while(i<9){
  p.R.t[i+1] = p.R.t[i]
  i = i+1
}
R.rat = .99
i.end = i+lag
while(i<i.end){
  p.R.t[i+1] = R.rat*p.R.t[i]
  i = i+1
}
while(i<p.t.end){
  p.R.t[i+1] = p.R.t[i]
  i = i+1
}
plot(p.R.t)

p.S = 1:p.t.end
p.S[1] = S[t]
p.I = 1:p.t.end
p.I[1] = I[t]
p.R = 1:p.t.end
p.R[1] = R[t]
p.t = 1
S.ref = p.S[1]
while(p.t < p.t.end){
  if(p.R.t[p.t+1] == p.R.t[p.t]){
    r = a*p.R.t[p.t]/S.ref
  } else {
    S.ref = p.S[p.t]
    r = a*p.R.t[p.t]/S.ref
  }
  p.S[p.t+1] = p.S[p.t] - r*p.I[p.t]*p.S[p.t]
  p.I[p.t+1] = p.I[p.t] + r*p.I[p.t]*p.S[p.t] - a*p.I[p.t]
  p.R[p.t+1] = p.R[p.t] + a*p.I[p.t]
  p.t = p.t+1
}
plot(p.S)
plot(p.I)
plot(p.R)

#5th Simulation: Extreme decline following announcement (lag + 13)
p.R.t = 1:p.t.end
R.rat = 1.032
p.R.t[1] = R.rat*R.t.col[t]
i=1
while(i<7){
  p.R.t[i+1] = R.rat*p.R.t[i]
  i = i+1
}
while(i<9){
  p.R.t[i+1] = p.R.t[i]
  i = i+1
}
R.rat = .956
i.end = i+lag+13
while(i<i.end){
  p.R.t[i+1] = R.rat*p.R.t[i]
  i = i+1
}
while(i<p.t.end){
  p.R.t[i+1] = p.R.t[i]
  i = i+1
}
#plot(p.R.t)
p.S = 1:p.t.end
p.S[1] = S[t]
p.I = 1:p.t.end
p.I[1] = I[t]
p.R = 1:p.t.end
p.R[1] = R[t]
p.t = 1
S.ref = p.S[1]
while(p.t < p.t.end){
  if(p.R.t[p.t+1] == p.R.t[p.t]){
    r = a*p.R.t[p.t]/S.ref
  } else {
    S.ref = p.S[p.t]
    r = a*p.R.t[p.t]/S.ref
  }
  p.S[p.t+1] = p.S[p.t] - r*p.I[p.t]*p.S[p.t]
  p.I[p.t+1] = p.I[p.t] + r*p.I[p.t]*p.S[p.t] - a*p.I[p.t]
  p.R[p.t+1] = p.R[p.t] + a*p.I[p.t]
  p.t = p.t+1
}
points(x, p.S, col="blue")
lines(x, p.S, col="blue",lty=1)
max(p.S)
min(p.S)
p.S[p.t.end]
points(x, p.I, col="blue")
lines(x, p.I, col="blue",lty=1)
max(p.I)
min(p.I)
p.I[p.t.end]
points(x, p.R, col="blue")
lines(x, p.R, col="blue",lty=1)
max(p.R)
min(p.R)
p.R[p.t.end]

#6th Simulation: moderate decline following announcement (lag + 13)
p.R.t = 1:p.t.end
R.rat = 1.032
p.R.t[1] = R.rat*R.t.col[t]
i=1
while(i<7){
  p.R.t[i+1] = R.rat*p.R.t[i]
  i = i+1
}
while(i<9){
  p.R.t[i+1] = p.R.t[i]
  i = i+1
}
R.rat = .966
i.end = i+lag+13
while(i<i.end){
  p.R.t[i+1] = R.rat*p.R.t[i]
  i = i+1
}
while(i<p.t.end){
  p.R.t[i+1] = p.R.t[i]
  i = i+1
}
#plot(p.R.t)
p.S = 1:p.t.end
p.S[1] = S[t]
p.I = 1:p.t.end
p.I[1] = I[t]
p.R = 1:p.t.end
p.R[1] = R[t]
p.t = 1
S.ref = p.S[1]
while(p.t < p.t.end){
  if(p.R.t[p.t+1] == p.R.t[p.t]){
    r = a*p.R.t[p.t]/S.ref
  } else {
    S.ref = p.S[p.t]
    r = a*p.R.t[p.t]/S.ref
  }
  p.S[p.t+1] = p.S[p.t] - r*p.I[p.t]*p.S[p.t]
  p.I[p.t+1] = p.I[p.t] + r*p.I[p.t]*p.S[p.t] - a*p.I[p.t]
  p.R[p.t+1] = p.R[p.t] + a*p.I[p.t]
  p.t = p.t+1
}
points(x, p.S, col="purple")
lines(x, p.S, col="purple",lty=1)
max(p.S)
min(p.S)
p.S[p.t.end]
points(x, p.I, col="purple")
lines(x, p.I, col="purple",lty=1)
max(p.I)
min(p.I)
p.I[p.t.end]
points(x, p.R, col="purple")
lines(x, p.R, col="purple",lty=1)
max(p.R)
min(p.R)
p.R[p.t.end]

#7th Simulation: minimal decline following announcement (lag + 13)
p.R.t = 1:p.t.end
R.rat = 1.032
p.R.t[1] = R.rat*R.t.col[t]
i=1
while(i<7){
  p.R.t[i+1] = R.rat*p.R.t[i]
  i = i+1
}
while(i<9){
  p.R.t[i+1] = p.R.t[i]
  i = i+1
}
R.rat = .99
i.end = i+lag+13
while(i<i.end){
  p.R.t[i+1] = R.rat*p.R.t[i]
  i = i+1
}
while(i<p.t.end){
  p.R.t[i+1] = p.R.t[i]
  i = i+1
}
#plot(p.R.t)
p.S = 1:p.t.end
p.S[1] = S[t]
p.I = 1:p.t.end
p.I[1] = I[t]
p.R = 1:p.t.end
p.R[1] = R[t]
p.t = 1
S.ref = p.S[1]
while(p.t < p.t.end){
  if(p.R.t[p.t+1] == p.R.t[p.t]){
    r = a*p.R.t[p.t]/S.ref
  } else {
    S.ref = p.S[p.t]
    r = a*p.R.t[p.t]/S.ref
  }
  p.S[p.t+1] = p.S[p.t] - r*p.I[p.t]*p.S[p.t]
  p.I[p.t+1] = p.I[p.t] + r*p.I[p.t]*p.S[p.t] - a*p.I[p.t]
  p.R[p.t+1] = p.R[p.t] + a*p.I[p.t]
  p.t = p.t+1
}
points(x, p.S, col="green")
lines(x, p.S, col="green",lty=1)
max(p.S)
min(p.S)
p.S[p.t.end]
points(x, p.I, col="green")
lines(x, p.I, col="green",lty=1)
max(p.I)
min(p.I)
p.I[p.t.end]
points(x, p.R, col="green")
lines(x, p.R, col="green",lty=1)
max(p.R)
min(p.R)
p.R[p.t.end]
