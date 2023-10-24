# trig-reg example

R <- 136 # tank radius, ft
H <- 66.4 # shell height, ft
Y <- 36e3 # yield strength, psi
E <- 30e6 # young's modulus, psi

# load in data and process ----
df <- read.csv("XYZ Example 1.csv")

df$x <- (df$X - mean(df$X)) * 3.28084 # re-center + meter-ft conversion
df$y <- (df$Y - mean(df$Y)) * 3.28084
df$z <- df$Z * 3.28084
theta <- atan2(df$y, df$x) # determine angle about center, radians
df$azimuth <- (theta + (theta < 0) * 2 * pi) # only positive angles
df <- df[order(df$azimuth),] # reorder by angle

# cosine curve fit ----

lm <- lm(z ~ cos(azimuth) + sin(azimuth), data = df)
df$zhat <- fitted(lm) # cosine curve fit
df$z1 <- residuals(lm) # out of plane settlement

# settlement data vs cosine curve plot
plot(df$azimuth, df$z*12, pch=16, cex = 0.5, col = "gray", xlab="angle (rad)", ylab="elevation (in)")
lines(df$azimuth, df$zhat*12, lwd = 2, col = "black")
legend("topright", legend = c("settlement data", "cosine curve fit"), pch = c(16, NA), lwd = c(NA, 2), col = c("gray", "black"), bty = "n")
title(main = "settlement data and cosine curve")

# out-of-plane settlement ----

# forward variable selection

# create data frame of basis cos and sin functions, B
# (easier to define them as columns in a matrix than defining them individually)
maxfreq = floor(2 * pi * R / (2 * 20)) # max frequency = circum / (2 * 20 ft)
B <- as.data.frame(matrix(ncol = 1 + 2 * maxfreq, nrow = nrow(df)))
B[,1] <- df$z1
colnames(B)[1] <- "z1"
for(freq in 1:maxfreq){
  B[,2 * (freq - 1) + 2] <- cos(freq * df$azimuth)
  B[,2 * (freq - 1) + 3] <- sin(freq * df$azimuth)
  colnames(B)[2 * (freq - 1) + 2] <- paste0("c", freq)
  colnames(B)[2 * (freq - 1) + 3] <- paste0("s", freq)
  # B$c1 is cos(azimuth), B$s2 is sin(2*azimuth), and so on...
}

# try adding each frequency and determine when to stop using more cos/sin terms
adjR2 <- vector(mode="numeric", length = maxfreq)
for(freq in seq(from = 1, to = maxfreq)){
  vars <- colnames(B)[seq(from = 2, to = 2 * (freq - 1) + 3)]
  frm <- as.formula(paste("z1 ~ -1 + ", paste(vars, collapse = " + ")))
  
  lmx <- lm(frm, data = B)
  adjR2[freq] <- summary(lmx)$adj.r.squared
}

# find the maximum frequency that increases adj. R^2
endfreq <- max(which(diff(adjR2) > 0)) + 1
vars <- colnames(B)[seq(from = 2, to = 2 * (endfreq - 1) + 3)]
frm <- as.formula(paste("z1 ~ -1 + ", paste(vars, collapse = " + ")))
lmx <- lm(frm, data = B)

df$z1hat <- fitted(lmx)

# out-of-plane-settlement plot
plot(df$azimuth, df$z1*12, pch=16, cex = 0.5, col = "gray", xlab="angle (rad)", ylab="elevation (in)")
lines(df$azimuth, df$z1hat*12, lwd = 2, col = "black")
legend("topright", legend = c("out-of-plane settlement", paste0(endfreq, "-mode trig-reg fit")), pch = c(16, NA), lwd = c(NA, 2), col = c("gray", "black"), bty = "n")
title(main = "out-of-plane settlement and trig-reg fit")

# allowable out-of-plane settlement ----

# calculate second derivative of trig-reg fit
coeffs <- lmx$coefficients
d2 <- vector(mode = "numeric", length = nrow(df))
for(freq in 1:endfreq){
  d2 <- d2 + freq^2 * coeffs[2 * (freq - 1) + 1] * B[,2 * (freq - 1) + 2] * -1 / R^2
  d2 <- d2 + freq^2 * coeffs[2 * (freq - 1) + 2] * B[,2 * (freq - 1) + 3] * -1 / R^2
}

# allowable second derivative based on modified Marr
d2max <- 11 * Y / (12 * H * E) # in/in^2

# second derivative and allowable plot
plot(df$azimuth, -d2 * 12 / 144, type = "l", lwd = 2, col = "red", xlab="angle (rad)", ylab="elevation (in)",
     ylim = c(min(-d2 * 12 / 144, -d2max), max(-d2 * 12 / 144, d2max)))
lines(c(0,2*pi), c(d2max, d2max), lty = 2, lwd = 2, col = "red")
lines(c(0,2*pi), c(-d2max, -d2max), lty = 2, lwd = 2, col = "red")
legend("topleft", legend = c("-d2 of trig-reg fit", "allowable d2"), lwd = c(2, 2), lty = c(1, 2), col = c("red", "red"))
title(main = "trig-reg fit allowable second derivative")

# overlay out-of-plane settlement with trig-reg second derivative
scalefactor <- signif(diff(range(df$z1 * 12)) / diff(range(d2 * 12 / 144)), 1)/2

plot(df$azimuth, df$z1*12, pch=16, cex = 0.5, col = "gray", xlab="angle (rad)", ylab="elevation (in)")
lines(df$azimuth, df$z1hat*12, lwd = 2, col = "black")
lines(df$azimuth, -d2 * 12 / 144 * scalefactor, lwd = 2, col = "red")
lines(c(0,2*pi), c(d2max, d2max) * scalefactor, lty = 2, lwd = 2, col = "red")
lines(c(0,2*pi), c(-d2max, -d2max) * scalefactor, lty = 2, lwd = 2, col = "red")
legend("topright", legend = c("out-of-plane settlement", paste0(endfreq, "-mode trig-reg fit"), "-d2 of trig-reg fit", "allowable d2"),
       pch = c(16, NA, NA, NA), lwd = c(NA, 2, 2, 2), lty = c(NA, 1, 1, 2), col = c("gray", "black", "red", "red"), bty = "n")
title(main = "out-of-plane settlement: trig-reg fit and second derivative")

