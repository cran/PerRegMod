
#""""""""""""""""""""""""Detecting periodicity
#install.packages("expm", repos = "https://cran.r-project.org")
#install.packages("sn", repos = "https://cran.r-project.org")
#install.packages("readxl", repos = "https://cran.r-project.org")

#install.packages("expm")
#install.packages("sn")
#install.packages("readxl")


#library(PerRegMod)
requireNamespace("PerRegMod")
requireNamespace("expm")
requireNamespace("readxl")
requireNamespace("sn")

library(PerRegMod)
library(expm)
library(readxl)
library(sn)# for skew normal density



n=200
s=2
x1=rnorm(n,0,1)
x2=rnorm(n,0,2)
x3=rnorm(n,0,3)
x4=rnorm(n,0,2.7)
y=rnorm(n,0,2.5)
x=list(x1,x2,x3,x4)
model=lm(y~x1+x2+x3+x4)
z=model$residuals
check_periodicity(x,y,s)

####""""""""""""""""""""""""" LSE method
set.seed(4)
s=2
n=200
m=n/s
p=3
mu=c(2,6)
beta1=c(3,5)
beta2=c(1,2.5)
beta3=c(-2,1)
x1=runif(n,0,5)
x2=runif(n,0,10)
x3=runif(n,0,15)
y=rep(0,n)
for (i in 1:s) {
  q=seq(i,n,s)
  y[q]=mu[i] + beta1[i] * x1[q] + beta2[i] * x2[q] + beta3[i] * x3[q] + rsn(m,alpha=10)
}
x=list(x1,x2,x3)
lm_per(x,y,s)


##""""""""""""""""""""""" AE
set.seed(2)
s=2
n=200
m=n/s
p=3
mu=c(2,6)
beta1=c(3,5)
beta2=c(1,2.5)
beta3=c(-2,1)
x1=runif(n,0,5)
x2=runif(n,0,10)
x3=runif(n,0,15)
y=rep(0,n)
for (i in 1:s) {
  q=seq(i,n,s)
  y[q]=mu[i] + beta1[i] * x1[q] + beta2[i] * x2[q] + beta3[i] * x3[q] + rnorm(m,0,1)
}
x=list(x1,x2,x3)
lm_per_AE(x,y,s)

######"""""""""""""""""""""""""" real data examples
#####"
#Data <- system.file("data", "weather_data_casablanca-settat.xlsx", package = "PerRegMod")
#load("~/PerRegMod_4.4.1/PerRegMod/data/weather_data_casablanca_settat.RData")
#Data=weather_data_casablanca_settat

#file_path <- system.file("data", "weather_data_casablanca_settat.xlsx", package = "PerRegMod")
#Data <- readxl::read_excel(file_path)
#file_path <- system.file("inst/extdata", "weather_data_casablanca_settat.xlsx", package = "PerRegMod")
#Data <- readxl::read_excel(file_path)
#Data=read_xlsx("inst/extdata/weather_data_casablanca_settat.xlsx")
file_path <- system.file("extdata", "weather_data_casablanca_settat.xlsx", package = "PerRegMod")

# Check if the file path is found, then read the file
if (file_path != "") {
  Data <- readxl::read_excel(file_path)
} else {
  stop("File not found. Please check that the file is correctly placed in 'inst/extdata' within the package directory.")
}
months <- c('January','February','March','April','May','June',
            'July','August','September','October','November','December')
years <- c('2010','2011','2012','2013','2014','2015','2016',
           '2017','2018','2019','2020')

#########""""""""""""for temperature
y=rep(0,11 * 12)
ly=length(years)
lm=length(months)
for (i in 1:ly ) {
  for (j in 1:lm) {
    y[(i - 1) * lm + j] = mean(Data$Temperature..C[Data$Year == years[i] &
                                                     Data$Month == months[j] ]  )
  }

}

#########""""""""""""for Humidity
x1=rep(0,11 * 12)
ly=length(years)
lm=length(months)
for (i in 1:ly ) {
  for (j in 1:lm) {
    x1[(i - 1) * lm + j]=mean(Data$Humidity..[Data$Year==years[i] &
                                                Data$Month==months[j] ]  )
  }

}

#########""""""""""""for dew point
x2=rep(0,11 * 12)
ly=length(years)
lm=length(months)
for (i in 1:ly ) {
  for (j in 1:lm) {
    x2[(i - 1) * lm + j]=mean(Data$Dew.Point..C[Data$Year==years[i] &
                                                  Data$Month==months[j] ]  )
  }

}

#########""""""""""""for wind
x3=rep(0,11 * 12)
ly=length(years)
lm=length(months)
for (i in 1:ly ) {
  for (j in 1:lm) {
    x3[(i - 1) * lm + j]=mean(Data$Wind.Speed.km.h[Data$Year==years[i] &
                                                     Data$Month==months[j] ]  )
  }

}

###""""""""""""""""""""""""example 1

#""""""""""check periodicity
x=list(x2,x3)# dew point and wind
n=length(y)
s=12
m=n/s
check_periodicity(x,y,s)

#""""""""estimating parameters using LSE
lm_per(x,y,s)
#""""""""estimating parameters using AE
lm_per_AE(x,y,s)

####"""""""""""""""""""""Example 2
##""""""check periodicity
x=list(x1,x3)# humidity and wind
n=length(y)
s=12
m=n/s
check_periodicity(x,y,s)

#""""""""estimating parameters using LSE
lm_per(x,y,s)
#""""""""estimating parameters using AE
lm_per_AE(x,y,s)

