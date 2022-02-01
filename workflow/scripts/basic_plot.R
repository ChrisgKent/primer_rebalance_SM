data = read.csv(snakemake@input[[1]])


fit = lm(y~x,
         data = data)

png(filename=snakemake@output[[1]])
par(mfcol = c(1,2))
plot(x=data$x,y=data$y)
abline(fit)
plot(fit$residuals)
abline(h=0)
dev.off()
