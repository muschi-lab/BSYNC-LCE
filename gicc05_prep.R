###################################################
# Split the input in equally spaced intervals, and calculate their duration and age
interval <- (tail(input[,1], 1) - head(input[,1], 1)) / bsync_settings$N
nodes <- seq(head(input[,1], 1), tail(input[,1], 1), by = interval)

# Estimate GICC05 relative counting error at each node
rce <- as.data.frame(approx(GICC05_MCE[,1], GICC05_MCE[,2], nodes))
rce <- round(diff(rce[,2]), digits = 1)
# Allow faster compression/expansion of the GICC05 timescale
rce <- bsync_settings$f * rce 

# MCE at each node
mce <- approx(GICC05_MCE[,1], GICC05_MCE[,2], nodes)$y
# MCE of first node
mce1 <- mce[1]
# Explore alignments outside the GICC05 bounds 
mce <- mce * bsync_settings$flex 
###################################################
