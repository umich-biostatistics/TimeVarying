# TimeVarying

Package description:

# Install

```
library(devtools)
install_github("umich-biostatistics/TimeVarying")

library(TimeVarying)
```

# Learn

Open R/Rstudio and run the following:

```
?simul

# generate the simuluation data
data <- simul(N = 2000, N_Strata = 5, p=5)

?SGD

result <- SGD(data$delta, data$z, data$time, facility = data$facility,
                knot = 10, M_stop = 10000, tol = 10^(-6), rate = 0.001)

?NR_new

result <- NR_new(data$delta, data$z, data$time, knot = 10 ,
                  facility = data$facility, M_stop = 10000, tol = 10^(-6),
                  rate = 0.001)
```
