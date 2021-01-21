#####################################################################################################################################
### Construction of balance-based PLS biplot to to examine the association between 24-hour movement behaviours and health outcome ###
#####################################################################################################################################

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



############
# Balances #
############

BAL = list()
CODES = list()
VV = list()
nB = character()
nB2 = character()

codes = matrix(0, D-1, D)
for(i in 1:(D-1)){
  codes[i,] = c(rep(0,i-1),1,rep(-1,D-i))
}   
V = t(gsi.buildilrBase(t(codes)))
cnB = character()
for (i in 1:(D-1)) {
  num = paste(cn[codes[i,]==1], collapse='.')
  denum = paste(cn[codes[i,]==-1], collapse='.')
  cnB[i] = paste(num, denum, sep='_')
  if (i==1) cnB2 = paste(denum, num, sep='_')
}
bal = log(as.matrix(MB))%*%t(V)
colnames(bal) = cnB
BAL[[1]] = bal
CODES[[1]] = codes
VV[[1]] = V
nB = c(nB, cnB[1])
nB2 = c(nB2, cnB2)

codes = matrix(c(1,1,-1,-1,-1,
                 1,-1,0,0,0,
                 0,0,1,-1,-1,
                 0,0,0,1,-1), D-1, D, byrow=T)
V = t(gsi.buildilrBase(t(codes)))
cnB = character()
for (i in 1:(D-1)) {
  num = paste(cn[codes[i,]==1], collapse='.')
  denum = paste(cn[codes[i,]==-1], collapse='.')
  cnB[i] = paste(num, denum, sep='_')
  if (i==1) cnB2 = paste(denum, num, sep='_')
}
bal = log(as.matrix(MB))%*%t(V)
colnames(bal) = cnB
BAL[[2]] = bal
CODES[[2]] = codes
VV[[2]] = V
nB = c(nB, cnB[1])
nB2 = c(nB2, cnB2)

codes = matrix(c(1,1,1,-1,-1,
                 1,-1,-1,0,0,
                 0,1,-1,0,0,
                 0,0,0,1,-1), D-1, D, byrow=T)
V = t(gsi.buildilrBase(t(codes)))
cnB = character()
for (i in 1:(D-1)) {
  num = paste(cn[codes[i,]==1], collapse='.')
  denum = paste(cn[codes[i,]==-1], collapse='.')
  cnB[i] = paste(num, denum, sep='_')
  if (i==1) cnB2 = paste(denum, num, sep='_')
}
bal = log(as.matrix(MB))%*%t(V)
colnames(bal) = cnB
BAL[[3]] = bal
CODES[[3]] = codes
VV[[3]] = V
nB = c(nB, cnB[1])
nB2 = c(nB2, cnB2)

codes = matrix(c(1,1,1,1,-1,
                 1,-1,-1,-1,0,
                 0,1,-1,-1,0,
                 0,0,1,-1,0), D-1, D, byrow=T)
V = t(gsi.buildilrBase(t(codes)))
cnB = character()
for (i in 1:(D-1)) {
  num = paste(cn[codes[i,]==1], collapse='.')
  denum = paste(cn[codes[i,]==-1], collapse='.')
  cnB[i] = paste(num, denum, sep='_')
  if (i==1) cnB2 = paste(denum, num, sep='_')
}
bal = log(as.matrix(MB))%*%t(V)
colnames(bal) = cnB
BAL[[4]] = bal
CODES[[4]] = codes
VV[[4]] = V
nB = c(nB, cnB[1])
nB2 = c(nB2, cnB2)

codes = matrix(c(1,-1,-1,-1,1,
                 1,0,0,0,-1,
                 0,1,-1,-1,0,
                 0,0,1,-1,0), D-1, D, byrow=T)
V = t(gsi.buildilrBase(t(codes)))
cnB = character()
for (i in 1:(D-1)) {
  num = paste(cn[codes[i,]==1], collapse='.')
  denum = paste(cn[codes[i,]==-1], collapse='.')
  cnB[i] = paste(num, denum, sep='_')
  if (i==1) cnB2 = paste(denum, num, sep='_')
}
bal = log(as.matrix(MB))%*%t(V)
colnames(bal) = cnB
BAL[[5]] = bal
CODES[[5]] = codes
VV[[5]] = V
nB = c(nB, cnB[1])
nB2 = c(nB2, cnB2)

codes = matrix(c(1,-1,-1,1,1,
                 0,1,-1,0,0,
                 1,0,0,-1,-1,
                 0,0,0,1,-1), D-1, D, byrow=T)
V = t(gsi.buildilrBase(t(codes)))
cnB = character()
for (i in 1:(D-1)) {
  num = paste(cn[codes[i,]==1], collapse='.')
  denum = paste(cn[codes[i,]==-1], collapse='.')
  cnB[i] = paste(num, denum, sep='_')
  if (i==1) cnB2 = paste(denum, num, sep='_')
}
bal = log(as.matrix(MB))%*%t(V)
colnames(bal) = cnB
BAL[[6]] = bal
CODES[[6]] = codes
VV[[6]] = V
nB = c(nB, cnB[1])
nB2 = c(nB2, cnB2)

codes = matrix(c(1,-1,1,1,1,
                 1,0,-1,-1,-1,
                 0,0,1,-1,-1,
                 0,0,0,1,-1), D-1, D, byrow=T)
V = t(gsi.buildilrBase(t(codes)))
cnB = character()
for (i in 1:(D-1)) {
  num = paste(cn[codes[i,]==1], collapse='.')
  denum = paste(cn[codes[i,]==-1], collapse='.')
  cnB[i] = paste(num, denum, sep='_')
  if (i==1) cnB2 = paste(denum, num, sep='_')
}
bal = log(as.matrix(MB))%*%t(V)
colnames(bal) = cnB
BAL[[7]] = bal
CODES[[7]] = codes
VV[[7]] = V
nB = c(nB, cnB[1])
nB2 = c(nB2, cnB2)

lB = length(BAL)



##################
# PLS regression #
##################

#centering the variables
y = as.vector(scale(Res, scale = F))
X = as.data.frame(scale(acomp(MB), scale = F))
Xr = as.matrix(scale(Real, scale = F))

V = VV[[1]]
Z = log(as.matrix(X))%*%t(V)

set.seed(1) # for reproducibility

misclass = mvr(y~cbind(Z, Xr), ncomp = 5, method = "kernel", validation = "CV")

# CV RMSEP and R^2
pls::RMSEP(misclass, estimate = "CV", intercept = FALSE)
pls::R2(misclass, estimate = "CV", intercept = FALSE)

# optimal number of PLS components based on the randomization test approach
comp = selectNcomp(misclass)


# bootstrap
kk = 1000
BE = matrix(rep(NA, kk*(lB+p1)), ncol = (lB+p1))
for(k in 1:kk){
  ind = sample(1:n, n, replace = T)  
  testy1 = Res[ind]
  testX1 = MB[ind, ]
  testX1r = Real[ind, ]
  testy = as.vector(scale(testy1, scale = F)) 
  testX = as.data.frame(scale(acomp(testX1), scale = F))
  testXr = as.matrix(scale(testX1r, scale = F))
  tbeta = rep(NA, lB+p1)
  for (ll in 1:lB) {
    selV = VV[[ll]]
    testZ = log(as.matrix(testX))%*%t(selV)
    resst = mvr(testy~cbind(testZ, testXr), ncomp = comp, method = "kernel")
    tbst = as.vector(coef(resst))
    tbeta[ll] = tbst[1]
  }
  tbeta[(lB+1):(lB+p1)] = tbst[D:(D-1+p1)]
  BE[k,] = tbeta
}

# significance of bootstrap standardized regression coefficients
beta = apply(BE, 2, mean)
sdbeta = apply(BE, 2, sd)
alpha1 = -qnorm(0.025)
a0 = beta/sdbeta
a = c(a0, -a0[1:lB])
names(a) = c(nB, pn, nB2)
stdRCo = (a[order(-abs(a))])
stdRCo
stdRCos = stdRCo[abs(stdRCo)>alpha1]
stdRCos

signary0 = rep(0, lB+p1)
for(i in 1:(lB+p1)){
  if((a[i]) > alpha1){signary0[i] = 1}
  if((a[i]) < -alpha1){signary0[i] = -1}
}
signary = c(signary0, -signary0[1:lB])



##############
# PLS biplot #
##############

selV = VV[[1]]
Z = log(as.matrix(X))%*%t(selV)
resst = mvr(y~cbind(Z, Xr), ncomp = comp, method = "kernel") 

# scores and loadings
G = matrix(as.numeric(resst$scores), n, comp)[, 1:2]
H = matrix(0, lB+p1, 2)
h = matrix(as.numeric(resst$loadings), D-1+p1, comp)[, 1:2]
H[1,] = h[1,]
H[(lB+1):(lB+p1), ] = h[D:(D-1+p1), ]

for (i in 2:lB){
  selV = VV[[i]]
  Z = log(as.matrix(X))%*%t(selV)
  resst = mvr(y~cbind(Z, Xr), ncomp = comp, method = "kernel") 
  h = as.numeric(resst$loadings)[c(1, D+p1)]
  H[i,] = h
}
H = rbind(H, -H[1:lB,])
H = 1*H # Scale loadings so that the arrows (points, respectively) are more visible in the biplot?

colnames(G) = c("Comp1", "Comp2")
colnames(H) = c("Comp1", "Comp2")

G = as.data.frame(G)
H = as.data.frame(H)
H$Variable = c(nB, pn, nB2)
H$Angle = ((180/pi)*atan(H$Comp2/H$Comp1))
H$Adj = (1-1.125*sign(H$Comp1))/2


pdf("PLSbiplot.pdf", height = 9, width = 9)
g = ggplot(G, aes(x = Comp1, y = Comp2)) +
  geom_point(aes(colour = Res), size = 2) +
  geom_segment(data = H, aes(x = 0, y = 0, xend = Comp1, yend = Comp2),
               arrow = arrow(length = unit(1/2, "picas")), 
               colour = c("darkblue","grey","darkred")[factor(signary)]) +
  geom_text(data = H, aes(label = Variable, x = Comp1, y = Comp2, angle = Angle, hjust = Adj), 
            size = 5, colour = c("darkblue","grey","darkred")[factor(signary)]) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_gradient(name = "Response", low = "black", high = "gold") +
  xlab("PLS comp. 1") + ylab("PLS comp. 2") +
  theme(panel.background=element_blank(),
        axis.line = element_line(color="black"),
        axis.text = element_blank(), 
        axis.ticks=element_blank(),
        axis.title = element_text(size = 15, vjust = 2, face = "bold"),
        legend.position = c(0.95, 0.125),
        legend.title = element_text(size = 15, face = "bold"), 
        legend.text = element_text(size = 13),
        legend.key = element_blank())
g + coord_fixed(ratio = 1) # Change ratio between y-axis and x-axis scale?
dev.off()
