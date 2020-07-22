library(deldir)
x <- 1000*runif(10)
y <- 1000*runif(10)
vt <- deldir(x, y)

png("tessellation")
par(mar=rep(1, 4))
plot(x, y, axes=FALSE, ann=FALSE, pch=16)
plot(vt, wlines="tess", lty="solid", add=TRUE)
box()
dev.off()

w <- runif(10, 1, 100)
library("gpclib")

unitSquare <- as(list(x = c(0, 0, 1000, 1000, 0), y = c(0, 1000, 1000, 0, 0)),
                 "gpc.poly")
awvt <- awv(list(x=x, y=y), w, unitSquare)


png("awv.png")
par(mar=rep(1, 4))
plot(x, y, axes=FALSE, ann=FALSE, pch=16)
lapply(awvt, plot, add=TRUE)
text(x, y, round(w), pos=3)
box()
dev.off()



temp <- seq(.1, .9, length=10)
target <- temp/sum(temp)


treemap <- allocate(letters[1:10], 
                      list(x=x, y=y), 
                      w, unitSquare, target)
drawRegions(treemap, label=TRUE)







