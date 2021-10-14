## Basic usage of simulate_infections
# Simulating infection incidence corresponding to a drop in Re value from 2.3
# to 0.5, at the half of the time period, then recovering the Re values using the
# estimate_Re function

Re_evolution <- c(rep(2.3, 100), rep(0.5, 100))
simulated_incidence_1 <- simulate_infections(
  Rt = Re_evolution
)
Re_recovered_1 <- estimate_Re(simulated_incidence_1)


## Advanced usage of simulate_infections
# Simulating infection incidence using the same Re progression as above, but a
# assuming a constant import of 100 cases per day

imported_infections <- rep(100, length(Re_evolution))
simulated_incidence_2 <- simulate_infections(
  Rt = Re_evolution,
  imported_infections = imported_infections
)
Re_recovered_2 <- estimate_Re(simulated_incidence_2)

