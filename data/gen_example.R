library(dplyr)
library(lubridate)

# one week, hourly measurements
start_time <- ymd_hms("2024-01-01 00:00:00")
end_time <- start_time + days(7)
time_series <- seq(start_time, end_time, by = "hour")

# number of data points
n <- length(time_series)

# data generation
set.seed(42)  

df <- tibble(
  time = time_series,
  V = runif(n, 10,15),  # Wind between 10 and 15 m/s
  Ta = runif(n, 260, 280), # Air temperature between 260 K and 280 K (-13.15 C and 10.85 C)
  RH = runif(n, 10, 20),  # Relative humidity between 10 and 20 percent
  p = runif(n, 650, 660), # Atmospheric pressure between 650 and 660 mbars
  Zv = sample(seq(2, 3, by = 0.5), n, replace = TRUE),  # Height of wind measurements between 2 and 2.5 m
  Zt = sample(seq(1, 1.5, by = 0.5), n, replace = TRUE),   # Height of temperature measurements between 1 and 1.5 m
  Zrh = sample(seq(1, 1.5, by = 0.5), n, replace = TRUE),  # Heigth if relative humidity measurements between 1 and 1.5 m
  Ts = runif(n, 260.15, 273.15), # Surface temperature between 260.15 K and 273.15 K (-10 C and 0 C) snow surface temperature is maximum 0 C (melting)
  z0 = runif(n, 0.001, 0.001),  # dynamic roughness length (m)
  z0t = runif(n, 0.0001, 0.0001), # thermal (m)
  z0q = runif(n, 0.0001, 0.0001)  # humidity (m)
)


z=(1:nrow(df))
result <- apply(t(z),2, function(x) MObulk(df$V[x], df$Ta[x], df$RH[x], df$p[x], 
                                               df$Zv[x],df$Zt[x],df$Zrh[x], df$Ts[x], 
                                               df$z0[x], df$z0t[x], df$z0q[x]))
df$H=result[1,]
df$LE=result[2,]
df$dH=result[7,]
df$dLE=result[8,]
