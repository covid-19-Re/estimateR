## Obtaining the discretized probabiility vector for different distributions

gamma_distribution <- list(name= "gamma", shape = 3.2, scale = 1.3)
delay_distribution_vector_1 <- build_delay_distribution(gamma_distribution)

normal_distribution <- list(name = "norm", mean = 5, sd = 2)
delay_distribution_vector_2 <- build_delay_distribution(normal_distribution)

uniform_distribution <- list(name = "unif", min = 0.5, max = 4)
delay_distribution_vector_3 <- build_delay_distribution(uniform_distribution)
