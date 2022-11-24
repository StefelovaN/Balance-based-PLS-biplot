##################################################################################################################################
### Construction of balance-based PLS biplot to examine the association between 24-hour movement behaviours and health outcome ###
##################################################################################################################################

library(compositions)
library(ggplot2)
library(pls)


# loading data
MBdata = read.csv("MBdata.csv", sep=";")
head(MBdata)
n = nrow(MBdata)

# compositional explanatory variables (MB composition)
MB = MBdata[, 1:5] # Choose the appropriate columns.
head(MB)
sum(MB==0) # Are there any 0? If so, they need to be imputed.
D = ncol(MB)
cn = colnames(MB)

# additional real-valued covariates
Real = MBdata[, 6:7] # Choose the appropriate columns.
head(Real)
# Real = log(Real) # ?log-transformation
p1 = ncol(Real)
pn = colnames(Real)

# Response variable
Res = MBdata[, 8] # Choose the appropriate column.
# Res = log(Res) # ?log-transformation





#####################################################
# Matrices of logcontrast coefficients for balances #
#####################################################

GG = list()
nB = character()

codes = matrix(0, D-1, D)
for(i in 1:(D-1)){
  codes[i,] = c(rep(0,i-1),1,rep(-1,D-i))
}   
G = t(gsi.buildilrBase(t(codes)))
cnB = character()
for (i in 1:(D-1)) {
  num = paste(cn[codes[i,]==1], collapse='.')
  denum = paste(cn[codes[i,]==-1], collapse='.')
  cnB[i] = paste(num, denum, sep='_')
  if (i==1) cnB2 = paste(denum, num, sep='_')
}
GG[[1]] = G
nB = c(nB, cnB[1])

codes = matrix(c(1,1,-1,-1,-1,
                 1,-1,0,0,0,
                 0,0,1,-1,-1,
                 0,0,0,1,-1), D-1, D, byrow=T)
G = t(gsi.buildilrBase(t(codes)))
cnB = character()
for (i in 1:(D-1)) {
  num = paste(cn[codes[i,]==1], collapse='.')
  denum = paste(cn[codes[i,]==-1], collapse='.')
  cnB[i] = paste(num, denum, sep='_')
  if (i==1) cnB2 = paste(denum, num, sep='_')
}
GG[[2]] = G
nB = c(nB, cnB[1])

codes = matrix(c(1,1,1,-1,-1,
                 1,-1,-1,0,0,
                 0,1,-1,0,0,
                 0,0,0,1,-1), D-1, D, byrow=T)
G = t(gsi.buildilrBase(t(codes)))
cnB = character()
for (i in 1:(D-1)) {
  num = paste(cn[codes[i,]==1], collapse='.')
  denum = paste(cn[codes[i,]==-1], collapse='.')
  cnB[i] = paste(num, denum, sep='_')
  if (i==1) cnB2 = paste(denum, num, sep='_')
}
GG[[3]] = G
nB = c(nB, cnB[1])

codes = matrix(c(1,1,1,1,-1,
                 1,-1,-1,-1,0,
                 0,1,-1,-1,0,
                 0,0,1,-1,0), D-1, D, byrow=T)
G = t(gsi.buildilrBase(t(codes)))
cnB = character()
for (i in 1:(D-1)) {
  num = paste(cn[codes[i,]==1], collapse='.')
  denum = paste(cn[codes[i,]==-1], collapse='.')
  cnB[i] = paste(num, denum, sep='_')
  if (i==1) cnB2 = paste(denum, num, sep='_')
}
GG[[4]] = G
nB = c(nB, cnB[1])

codes = matrix(c(1,-1,-1,-1,1,
                 1,0,0,0,-1,
                 0,1,-1,-1,0,
                 0,0,1,-1,0), D-1, D, byrow=T)
G = t(gsi.buildilrBase(t(codes)))
cnB = character()
for (i in 1:(D-1)) {
  num = paste(cn[codes[i,]==1], collapse='.')
  denum = paste(cn[codes[i,]==-1], collapse='.')
  cnB[i] = paste(num, denum, sep='_')
  if (i==1) cnB2 = paste(denum, num, sep='_')
}
GG[[5]] = G
nB = c(nB, cnB[1])

codes = matrix(c(1,-1,-1,1,1,
                 0,1,-1,0,0,
                 1,0,0,-1,-1,
                 0,0,0,1,-1), D-1, D, byrow=T)
G = t(gsi.buildilrBase(t(codes)))
cnB = character()
for (i in 1:(D-1)) {
  num = paste(cn[codes[i,]==1], collapse='.')
  denum = paste(cn[codes[i,]==-1], collapse='.')
  cnB[i] = paste(num, denum, sep='_')
  if (i==1) cnB2 = paste(denum, num, sep='_')
}
GG[[6]] = G
nB = c(nB, cnB[1])

codes = matrix(c(1,-1,1,1,1,
                 1,0,-1,-1,-1,
                 0,0,1,-1,-1,
                 0,0,0,1,-1), D-1, D, byrow=T)
G = t(gsi.buildilrBase(t(codes)))
cnB = character()
for (i in 1:(D-1)) {
  num = paste(cn[codes[i,]==1], collapse='.')
  denum = paste(cn[codes[i,]==-1], collapse='.')
  cnB[i] = paste(num, denum, sep='_')
  if (i==1) cnB2 = paste(denum, num, sep='_')
}
GG[[7]] = G
nB = c(nB, cnB[1])

lB = length(nB)



H = matrix(0, lB, D)
for (i in 1:lB) {
  H[i, ] = GG[[i]][1, ]
}





##################
# PLS regression #
##################

#centering the variables
y = as.vector(scale(Res, scale = F))
X = as.matrix(as.data.frame(acomp(MB)-mean(acomp(MB))))
Xr = as.matrix(scale(Real, scale = F))

C = clr(X)

set.seed(1) # for reproducibility

misclass = mvr(y~cbind(C, Xr), ncomp = 5, validation = "CV")

# CV RMSEP and R^2
pls::RMSEP(misclass, estimate = "CV", intercept = FALSE)
pls::R2(misclass, estimate = "CV", intercept = FALSE)

# optimal number of PLS components based on the randomization test approach
comp = selectNcomp(misclass)

# % of explained variability
summary(mvr(y~cbind(C, Xr), ncomp = comp))

set.seed(11)
# bootstrap
kk = 1000
BE = matrix(rep(NA, kk*(lB+p1)), ncol = (lB+p1))
for(k in 1:kk){
  ind = sample(1:n, n, replace = T)  
  testy1 = Res[ind]
  testX1 = MB[ind, ]
  testX1r = Real[ind, ]
  testy = as.vector(scale(testy1, scale = F)) 
  testX = as.matrix(as.data.frame(acomp(testX1)-mean(acomp(testX1))))
  testXr = as.matrix(scale(testX1r, scale = F))
  testC = clr(testX)
  resst = mvr(testy~cbind(testC, testXr), ncomp = comp)
  tbeta0 = as.vector(coef(resst))
  tbeta_clr = tbeta0[1:D]
  tbeta_real = tbeta0[(D+1):(D+p1)]
  tbeta_bal = H %*% tbeta_clr
  tbeta = c(tbeta_bal, tbeta_real)
  BE[k,] = tbeta
}

# significance of bootstrap standardized regression coefficients
alpha0 = 0.05

beta = colMeans(BE)
names(beta) = c(nB, pn)
sd = apply(BE, 2, sd)
names(sd) = c(nB, pn)
me = apply(BE, 2, function(x) mean(x/sd(x)))
names(me) = c(nB, pn)
lw = apply(BE, 2, function(x) quantile(x/sd(x), alpha0/2))
names(lw) = c(nB, pn)
up = apply(BE, 2, function(x) quantile(x/sd(x), 1-alpha0/2))
names(up) = c(nB, pn)

signary = rep(0, lB+p1)
names(signary) = names(me)
for(i in 1:(lB+p1)){
  if((lw[i] > 0) & (up[i] > 0)){signary[i] = 1}
  if((lw[i] < 0) & (up[i] < 0)){signary[i] = -1}
}

round(me, 2)
round(lw, 2)
round(up, 2)
signary





##############
# PLS biplot #
##############

#centering the variables
y = as.vector(scale(Res, scale = F))
X = as.matrix(as.data.frame(acomp(MB)-mean(acomp(MB))))
Xr = as.matrix(scale(Real, scale = F))

C = clr(X)
resst = mvr(y~cbind(C, Xr), ncomp = comp) 

# scores and loadings
T = matrix(as.numeric(resst$scores), n, comp)[, 1:2]
P0 = matrix(as.numeric(resst$loadings), D+p1, comp)[, 1:2]
P_clr = P0[1:D, ]
P_real = P0[(D+1):(D+p1), ]
P_bal = H %*% P_clr
P = rbind(P_bal, P_real)

# scaling constant
scon = min(max(T[,1])/max(P[,1]), min(T[,1])/min(P[,1]), max(T[,2])/max(P[,2]), min(T[,2])/min(P[,2]))
P = scon*P

colnames(T) = c("Comp1", "Comp2")
colnames(P) = c("Comp1", "Comp2")

T = as.data.frame(T)
P = as.data.frame(P)
P$Variable = c(nB, pn)
P$Angle = ((180/pi)*atan(P$Comp2/P$Comp1))
P$hAdj = (1-1.1*sign(P$Comp1))/2

pdf("PLSbiplot.pdf", height = 9, width = 9)
g = ggplot(T, aes(x = Comp1, y = Comp2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  geom_point(aes(colour = Res), size = 2) +
  geom_segment(data = P, aes(x = 0, y = 0, xend = Comp1, yend = Comp2),
               arrow = arrow(length = unit(1, "picas")), 
               colour = "grey15", size = 1) +
  geom_text(data = P, aes(label = Variable, x = Comp1, y = Comp2, angle = Angle, hjust = hAdj), 
            size = 5, colour = c("darkblue","grey40","darkred")[factor(signary)], fontface = "bold") +
  scale_color_gradient(name = "Name", low = "grey10", high = "yellow") +
  scale_x_continuous("PLS comp. 1 (scores)",  sec.axis = sec_axis(~ . / scon, name = "PLS comp. 1 (loadings)")) +
  scale_y_continuous("PLS comp. 2 (scores)",  sec.axis = sec_axis(~ . / scon, name = "PLS comp. 2 (loadings)")) +
  theme(panel.background=element_blank(),
        axis.line = element_line(color="black"),
        axis.title.x.bottom = element_text(size = 15, face = "bold", margin = margin(t = 10)),
        axis.title.y.left = element_text(size = 15, face = "bold", margin = margin(r = 10)),
        axis.title.x.top = element_text(size = 15, face = "bold", margin = margin(b = 10)),
        axis.title.y.right = element_text(size = 15, face = "bold", margin = margin(l = 10)),
        legend.position = c(0.935, 0.125),
        legend.title = element_text(size = 15, face = "bold"), 
        legend.text = element_text(size = 13),
        legend.key = element_blank())
g
dev.off()
