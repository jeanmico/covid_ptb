# zph issue debug
frq_sex = births %>% tabyl(inf_sex, preterm)
cox_simple <- coxph(Surv(gest_weeks, preterm) ~ inf_sex, data = births)
simple_sr <- cox.zph(cox_simple)

xx <- simple_sr$x
yy = simple_sr$y
df = 4
nsmo = 5
d<- nrow(yy)
nvar <- ncol(yy)
pred.x <- seq(from = min(xx), to = max(xx), length = nsmo)
temp <- c(pred.x, xx)
lmat <- splines::ns(temp, df = 5, intercept = TRUE) # df = 1 is too small, df = 4 gives NA/NaN/Inf error



pmat <- lmat[1:nsmo, ]
xmat <- lmat[-(1:nsmo), ]
qmat <- qr(xmat)


dfggcoxzph(simple_sr)

tmpsr = data.frame(simple_sr$x, simple_sr$y)



tstplt <- ggplot(tmpsr, aes(x = simple_sr.x, y = inf_sex)) + geom_line()
tstplt
