data <- read.csv("temp.csv", header =FALSE)
D<-data[2]
#head(D)
block_length<- ceiling(nrow(D)^(1/3))
blocks <- D[1:block_length, 1]
for (i in 2:(length(D[, 1]) - (block_length-1))) {
blocks <- rbind(blocks, D[i:(i + (block_length-1)), 1])
}
mbb_size<-floor(nrow(D)/block_length)
# MOVING BLOCK BOOTSTRAP
xbar <- NULL
for (i in 1:1000) {
take.blocks <- sample(1:(length(D[, 1]) - (block_length-1)), mbb_size,
replace = TRUE)
newdat = c(t(blocks[take.blocks, ]))
xbar[i] = mean(newdat)
}
###hist(xbar)
#cat("the mean estimate is", mean(xbar), "the sample
#standard deviation is ",
#sd(xbar))

cat(mean(xbar), " ", sd(xbar), " ")
z_value = mean(xbar)/sd(xbar)
cat(2*pnorm(-abs(z_value)), " ")

D<-data[3]
#head(D)
block_length<- ceiling(nrow(D)^(1/3))
blocks <- D[1:block_length, 1]
for (i in 2:(length(D[, 1]) - (block_length-1))) {
blocks <- rbind(blocks, D[i:(i + (block_length-1)), 1])
}
mbb_size<-floor(nrow(D)/block_length)
# MOVING BLOCK BOOTSTRAP
xbar <- NULL
for (i in 1:1000) {
take.blocks <- sample(1:(length(D[, 1]) - (block_length-1)), mbb_size,
replace = TRUE)
newdat = c(t(blocks[take.blocks, ]))
xbar[i] = mean(newdat)
}
###hist(xbar)
#cat("the mean estimate is", mean(xbar), "the sample
#standard deviation is ",
#sd(xbar))

cat(mean(xbar), " ", sd(xbar), " ")
z_value = mean(xbar)/sd(xbar)
cat(2*pnorm(-abs(z_value)))
