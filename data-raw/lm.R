set.seed(2017)

sim_lm <- function(X, Beta){
	X %*% Beta + rnorm(nrow(X))
}

#### Two group ####

# Two group with and without effect
d1 <- data.frame(Group=c(rep("Ctrl", 5), rep("Trt", 5)))
d1$Expression <- sim_lm(X= model.matrix(~Group, data=d1), Beta=c(5, 2))

d2 <- data.frame(Group=c(rep("Ctrl", 5), rep("Trt", 5)))
d2$Expression <- sim_lm(X= model.matrix(~Group, data=d2), Beta=c(5, 0))

#### Interactions ####

# Interaction effect
d3 <- data.frame(expand.grid(Group=c(rep("Ctrl", 5), rep("Trt", 5)), Time=c("0h", "1h")))
d3$Expression <- sim_lm(X= model.matrix(~Group*Time, data=d3), Beta=c(5, 2, 3, 0))

# Interaction effect
d4 <- data.frame(expand.grid(Group=c(rep("Ctrl", 5), rep("Trt", 5)), Time=c("0h", "1h")))
d4$Expression <- sim_lm(X= model.matrix(~Group*Time, data=d4), Beta=c(5, 2, 3, 4))

#### Batch ####

d5 <- data.frame(expand.grid(Group=c(rep("Ctrl", 5), rep("Trt", 5)), Batch=c("A", "B")))
d5$Expression <- sim_lm(X=model.matrix(~Group+Batch, data=d5), Beta=c(5, 2, 8))

### Three Group ####

# Two group with and without effect
d6 <- data.frame(Group=c(rep("A", 5), rep("B", 5), rep("C", 5)))
d6$Expression <- sim_lm(X=model.matrix(~Group, data=d6), Beta=c(5, 1, 4))

#### Format output ####

# As list
lmExamples <- list(twoGroup1=d1,
									 twoGroup2=d2,
									 interaction1=d3,
									 interaction2=d4,
									 batchEffect=d5,
									 threeGroup=d6)

# Save
devtools::use_data(lmExamples, overwrite = TRUE)
