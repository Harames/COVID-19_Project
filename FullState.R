setwd("~/R/my_r_progs/M3550_Project")
library(R0)

#Create data.frames for each data set.
#Then work those data frames to have the same length.
df1 = read.csv('rt.csv')
head(df1)
names(df1)
df1 = df1[which(df1$region == "UT"),]
df1$date = as.Date(df1$date)
df1[is.na.data.frame(df1)]<-0
df1 = df1[which(df1$date >= "2020-3-6" & df1$date <= "2020-10-31"),]
tail(df1)

df2 = read.csv('Overview_Cumulative COVID-19 Cases with Estimated Recoveries_2020-11-01.csv')
head(df2)
names(df2)
df2$Date = as.Date(df2$Date)
df2 = df2[which(df2$Date >= "2020-3-6" & df2$Date <= "2020-10-31"),]
tail(df2)

plot(df2$Estimated.Active~df2$Date)
plot(df2$Total~df2$Date)
plot(df2$Died~df2$Date)
plot(df2$Estimated.Recovered~df2$Date)

#Select the needed columns
covid = data.frame(df1$date,df2$Total,df1$new_deaths,df2$Died,
                   df2$Estimated.Recovered..,df2$Estimated.Active)
names(covid) = c('date','total.infected','new.deaths','total.deaths',
                 'total.recovered','active.infected')

#Generate additional columns
tot_pop = 3210000
current.susceptible = tot_pop - covid$total.infected
plot(current.susceptible)

new.infected = covid$total.infected
i = length(new.infected)
while(i>1){
  new.infected[i] = new.infected[i] - new.infected[i-1]
  i = i-1
}
plot(new.infected)

new.recovered = covid$total.recovered
i = length(new.recovered)
while(i>1){
  new.recovered[i] = new.recovered[i] - new.recovered[i-1]
  i = i-1
}
plot(new.recovered)

new.removed = new.recovered + covid$new.deaths
plot(new.removed)

total.removed = new.removed
i=1
while(i<length(total.removed)){
  total.removed[i+1] = total.removed[i+1] + total.removed[i]
  i=i+1
}
plot(total.removed)

#Add new columns to data.frame
covid = data.frame(covid,current.susceptible,new.infected,
                   new.recovered,new.removed,total.removed)

#Now we must generate approximations for Rt
attach(covid)
#cGT = generation.time("gamma",c(5.2,2.8))  #Standard
k=5.2/2.8
cGT = generation.time("gamma",c(17,17/k))
cGT[is.na(cGT)] <- 0
cGT
plot(cGT)
Rt = est.R0.TD(active.infected,cGT,n.t0 = active.infected[1],t=date)
Rt
Rt$R
plot(Rt)
weekly.Rt = smooth.Rt(Rt,7)
weekly.Rt
plot(weekly.Rt)
month.Rt = smooth.Rt(Rt,28)
month.Rt
plot(month.Rt)

Rt.col = covid$total.infected
i = 1
while (i<=length(Rt$R)) {
  Rt.col[i] = Rt$R[i]
  i = i+1
}
while(i<=length(date)){
  Rt.col[i] = Rt$R[length(Rt$R)]
  i = i+1
}
Rt.col[length(Rt.col)]=Rt.col[length(Rt.col)-1]
plot(Rt.col)
covid$Rt.col = Rt.col

smooth.Rt.col = covid$total.infected
i=1
while(i<=length(weekly.Rt$R)){
  index=1
  while(index<=7){
    smooth.Rt.col[7*(i-1)+index]=weekly.Rt$R[i]
    index=index+1
  }
  i=i+1
}
i=7*(i-1)+1
while(i<=length(covid$date)){
  smooth.Rt.col[i] = weekly.Rt$R[length(weekly.Rt$R)]
  i = i+1
}
plot(smooth.Rt.col)
covid$smooth.Rt.col = smooth.Rt.col

long.Rt.col = covid$total.infected
i=1
while(i<=length(month.Rt$R)){
  index=1
  while(index<=28){
    long.Rt.col[28*(i-1)+index]=month.Rt$R[i]
    index=index+1
  }
  i=i+1
}
i=28*(i-1)+1
while(i<=length(covid$date)){
  long.Rt.col[i] = month.Rt$R[length(month.Rt$R)]
  i = i+1
}
plot(long.Rt.col)
covid$long.Rt.col = long.Rt.col

covid = data.frame(covid[which(covid$date >= "2020-5-1"),])
head(covid)
tail(covid)

#Can't properly attach covid df, so we have to use covid$ from her on
a = 1/20
r = covid$new.deaths
current.Inf = covid$active.infected
current.Sus = covid$current.susceptible
current.Rem = covid$total.removed
current.Inf[1]
current.Sus[1]
current.Rem[1]
t = 1
while(t<length(current.Inf)){
  r[t] = a*covid$Rt.col[t]/current.Sus[t]
  current.Sus[t+1] = current.Sus[t] - r[t]*current.Inf[t]*current.Sus[t]
  current.Inf[t+1] = current.Inf[t] + r[t]*current.Inf[t]*current.Sus[t] - a*current.Inf[t]
  current.Rem[t+1] = current.Rem[t] + a*current.Inf[t]
  t = t+1
}
r[length(r)]=r[length(r)-1]
#plot(r)
plot(covid$active.infected~covid$date)
points(current.Inf~covid$date, col='red')

t = 1
while(t<length(current.Inf)){
  r[t] = a*covid$smooth.Rt.col[t]/current.Sus[t]
  #TODO: update each model on each step through the loop
  current.Sus[t+1] = current.Sus[t] - r[t]*current.Inf[t]*current.Sus[t]
  current.Inf[t+1] = current.Inf[t] + r[t]*current.Inf[t]*current.Sus[t] - a*current.Inf[t]
  current.Rem[t+1] = current.Rem[t] + a*current.Inf[t]
  t = t+1
}
r[length(r)]=r[length(r)-1]
#plot(r)
plot(covid$active.infected~covid$date)
points(current.Inf~covid$date, col='red')

t = 1
while(t<length(current.Inf)){
  r[t] = a*covid$long.Rt.col[t]/current.Sus[t]
  #TODO: update each model on each step through the loop
  current.Sus[t+1] = current.Sus[t] - r[t]*current.Inf[t]*current.Sus[t]
  current.Inf[t+1] = current.Inf[t] + r[t]*current.Inf[t]*current.Sus[t] - a*current.Inf[t]
  current.Rem[t+1] = current.Rem[t] + a*current.Inf[t]
  t = t+1
}
r[length(r)]=r[length(r)-1]
#plot(r)
plot(covid$active.infected~covid$date)
points(current.Inf~covid$date, col='red')
