load("~/RemaisWork/Schisto/Stanford/Prawn_fisheries_Public_health/Sensitivity_Analysis/sens_sims.RData")
library(sensitivity)
library(ggplot2)

# Assess PRCC of prawn model outputs ###############
# profit sensitivity  
profit.pcc = pcc(X = as.data.frame(lhcpars[,which(colnames(lhcpars) %in% c(names(par.aqua), 'c', 'p', 'eta'))]),
                 y = as.data.frame(lhcpars.aqua.fill$profit),
                 rank = TRUE)

profit.pcc.df = profit.pcc$PRCC
profit.pcc.df$cor = cor(lhcpars[,which(colnames(lhcpars) %in% c(names(par.aqua), 'c', 'p', 'eta'))], lhcpars.aqua.fill$profit)
profit.pcc.df$var = c(names(par.aqua), 'c', 'p', 'eta')

ggplot(profit.pcc.df, aes(x = cor, y = original, label = var)) + 
  geom_point(pch = 16) +
  geom_text(aes(label = var), hjust = -0.5, vjust = -0.5)

profit.lhsprcc = ggplot(profit.pcc.df, aes(x = var, y = original)) +
  theme_bw()+
  theme(axis.text = element_text(size = 12),  #increase axis label size
        axis.title = element_text(size = 15), #increase axis title size
        axis.title.x = element_blank()) + 
  scale_y_continuous(limits = c(-0.5,0.5), breaks = c(-0.5,0,0.5))+
  scale_x_discrete(labels = c("a.p" = expression(italic(a[P])),
                              "b.p" = expression(italic(b[P])),
                              "c" = expression(italic(c)),
                              "d" = expression(italic(d)),
                              "eta" = expression(italic(zeta)),
                              "gam" = expression(italic(gamma)),
                              "k" = expression(italic(kappa)),
                              "linf" = expression(italic(L[infinity])),
                              "muP" = expression(italic(mu[P])),
                              "om" = expression(italic(omega)),
                              "p" = expression(italic(p)))) +
  geom_bar(fill = 'darkgreen', stat = 'identity', width = 0.25)+
  geom_hline(yintercept = 0, lty = 2)+
  annotate('text', x = 1.2, y = 0.5, label = 'A)', size = 8) +
  labs(x = 'Parameter', y = 'Profit PRCC')

profit.lhsprcc

# roi sensitivity  
roi.pcc = pcc(X = as.data.frame(lhcpars[,which(colnames(lhcpars) %in% c(names(par.aqua), 'c', 'p', 'eta'))]),
              y = as.data.frame(lhcpars.aqua.fill$roi),
              rank = TRUE)

roi.pcc.df = roi.pcc$PRCC
roi.pcc.df$var = c(names(par.aqua), 'c', 'p', 'eta')

roi.lhsprcc = ggplot(roi.pcc.df, aes(x = var, y = original)) +
  theme_bw()+
  theme(axis.text = element_text(size = 12),  #increase axis label size
        axis.title = element_text(size = 15)) + #increase axis title size
  scale_y_continuous(limits = c(-1,1), breaks = c(-1,-0.5,0,0.5,1))+
  scale_x_discrete(labels = c("a.p" = expression(italic(a[P])),
                              "b.p" = expression(italic(b[P])),
                              "c" = expression(italic(c)),
                              "d" = expression(italic(d)),
                              "eta" = expression(italic(zeta)),
                              "gam" = expression(italic(gamma)),
                              "k" = expression(italic(kappa)),
                              "linf" = expression(italic(L[infinity])),
                              "muP" = expression(italic(mu[P])),
                              "om" = expression(italic(omega)),
                              "p" = expression(italic(p)))) +
  geom_bar(fill = 'darkgreen', stat = 'identity', width = 0.25)+
  geom_hline(yintercept = 0, lty = 2)+
  labs(title = 'roi sensitivity', x = 'Parameter', y = 'PRCC')

roi.lhsprcc

# harvest time sensitivity  
ht.pcc = pcc(X = as.data.frame(lhcpars[,which(colnames(lhcpars) %in% c(names(par.aqua), 'c', 'p', 'eta'))]),
             y = as.data.frame(lhcpars.aqua.fill$h.t),
             rank = TRUE)

ht.pcc.df = ht.pcc$PRCC
ht.pcc.df$var = c(names(par.aqua), 'c', 'p', 'eta')

ht.lhsprcc = ggplot(ht.pcc.df, aes(x = var, y = original)) +
  theme_bw()+
  theme(axis.text = element_text(size = 12),  #increase axis label size
        axis.title = element_text(size = 15)) + #increase axis title size
  scale_y_continuous(limits = c(-1,1), breaks = c(-1,-0.5,0,0.5,1))+
  scale_x_discrete(labels = c("a.p" = expression(italic(a[P])),
                              "b.p" = expression(italic(b[P])),
                              "c" = expression(italic(c)),
                              "d" = expression(italic(d)),
                              "eta" = expression(italic(zeta)),
                              "gam" = expression(italic(gamma)),
                              "k" = expression(italic(kappa)),
                              "linf" = expression(italic(L[infinity])),
                              "muP" = expression(italic(mu[P])),
                              "om" = expression(italic(omega)),
                              "p" = expression(italic(p)))) +
  geom_bar(fill = 'darkgreen', stat = 'identity', width = 0.25)+
  #geom_errorbar(x = var, )
  geom_hline(yintercept = 0, lty = 2)+
  labs(title = 'ht sensitivity', x = 'Parameter', y = 'PRCC')

ht.lhsprcc

# Assess PRCC of epi model outputs ###########
#Equilibrium mean worm burden  
W.pcc = pcc(X = as.data.frame(lhcpars[,which(colnames(lhcpars) %in% 
                                               names(par.snails)[-which(names(par.snails) %in% 
                                                                          c('A', 'H', 'psi1', 'psi2', 'psi3', 'phi'))])]),
            y = as.data.frame(lhcpars.epi.fill$W),
            rank = TRUE)

W.pcc.df = W.pcc$PRCC
W.pcc.df$cor = cor(lhcpars[,which(colnames(lhcpars) %in% 
                                    names(par.snails)[-which(names(par.snails) %in% 
                                                               c('A', 'H', 'psi1', 'psi2', 'psi3', 'phi'))])], 
                   lhcpars.epi.fill$W)
W.pcc.df$var = names(par.snails)[-which(names(par.snails) %in% 
                                          c('A', 'H', 'psi1', 'psi2', 'psi3', 'phi'))]

ggplot(W.pcc.df, aes(x = cor, y = original, label = var)) + 
  geom_point(pch = 16) +
  geom_text(aes(label = var), hjust = -0.5, vjust = -0.5)

W.lhsprcc = ggplot(W.pcc.df, aes(x = var, y = original)) +
  theme_bw()+
  theme(axis.text = element_text(size = 12),  #increase axis label size
        axis.title = element_text(size = 15)) + #increase axis title size
  scale_y_continuous(limits = c(-0.5,0.5), breaks = c(-0.5,0,0.5))+
  scale_x_discrete(labels = c("beta" = expression(italic(beta)),
                              "f" = expression(italic(f)),
                              "g1" = expression(g[1]),
                              "g2" = expression(italic(g[2])),
                              "Kn" = expression(italic(K)),
                              "lambda" = expression(italic(lambda)),
                              "m" = expression(italic(m)),
                              "muH" = expression(italic(mu[H])),
                              "muI" = expression(italic(mu[I])),
                              "muN1" = expression(italic(mu[N1])),
                              "muN2" = expression(italic(mu[N2])),
                              "muN3" = expression(italic(mu[N3])),
                              "muW" = expression(italic(mu[W])),
                              "sigma" = expression(italic(sigma)),
                              "theta" = expression(italic(theta)),
                              "z" = expression(italic(z)))) +
  geom_bar(fill = 'purple', stat = 'identity', width = 0.25)+
  #geom_errorbar(x = var, )
  geom_hline(yintercept = 0, lty = 2)+
  labs(x = 'Parameter', y = 'Equilibrium W PRCC')

W.lhsprcc

#Equilibrium total snail population 
Nt.pcc = pcc(X = as.data.frame(lhcpars[,which(colnames(lhcpars) %in% 
                                                names(par.snails)[-which(names(par.snails) %in% 
                                                                           c('A', 'H', 'psi1', 'psi2', 'psi3', 'phi'))])]),
             y = as.data.frame(lhcpars.epi.fill$N.t/area),
             rank = TRUE)

Nt.pcc.df = Nt.pcc$PRCC
Nt.pcc.df$cor = cor(lhcpars[,which(colnames(lhcpars) %in% 
                                    names(par.snails)[-which(names(par.snails) %in% 
                                                               c('A', 'H', 'psi1', 'psi2', 'psi3', 'phi'))])], 
                    lhcpars.epi.fill$N.t/area)
Nt.pcc.df$var = names(par.snails)[-which(names(par.snails) %in% 
                                           c('A', 'H', 'psi1', 'psi2', 'psi3', 'phi'))]

ggplot(Nt.pcc.df, aes(x = cor, y = original, label = var)) + 
  geom_point(pch = 16) +
  geom_text(aes(label = var), hjust = -0.5, vjust = -0.5)

Nt.lhsprcc = ggplot(Nt.pcc.df, aes(x = var, y = original)) +
  theme_bw()+
  theme(axis.text = element_text(size = 12),  #increase axis label size
        axis.title = element_text(size = 15)) + #increase axis title size
  scale_y_continuous(limits = c(-0.5,0.5), breaks = c(-0.5,0,0.5))+
  scale_x_discrete(labels = c("beta" = expression(italic(beta)),
                              "f" = expression(italic(f)),
                              "g1" = expression(g[1]),
                              "g2" = expression(italic(g[2])),
                              "Kn" = expression(italic(K)),
                              "lambda" = expression(italic(lambda)),
                              "m" = expression(italic(m)),
                              "muH" = expression(italic(mu[H])),
                              "muI" = expression(italic(mu[I])),
                              "muN1" = expression(italic(mu[N1])),
                              "muN2" = expression(italic(mu[N2])),
                              "muN3" = expression(italic(mu[N3])),
                              "muW" = expression(italic(mu[W])),
                              "sigma" = expression(italic(sigma)),
                              "theta" = expression(italic(theta)),
                              "z" = expression(italic(z)))) +
  geom_bar(fill = 'black', stat = 'identity', width = 0.25)+
  #geom_errorbar(x = var, )
  geom_hline(yintercept = 0, lty = 2)+
  labs(title = 'Nt sensitivity', x = 'Parameter', y = 'PRCC')

Nt.lhsprcc

#Equilibrium infected snail population 
It.pcc = pcc(X = as.data.frame(lhcpars[,which(colnames(lhcpars) %in% 
                                                names(par.snails)[-which(names(par.snails) %in% 
                                                                           c('A', 'H', 'psi1', 'psi2', 'psi3', 'phi'))])]),
             y = as.data.frame(lhcpars.epi.fill$I.t/area),
             rank = TRUE)

It.pcc.df = It.pcc$PRCC
It.pcc.df$cor = cor(lhcpars[,which(colnames(lhcpars) %in% 
                                     names(par.snails)[-which(names(par.snails) %in% 
                                                                c('A', 'H', 'psi1', 'psi2', 'psi3', 'phi'))])], 
                    lhcpars.epi.fill$I.t/area)
It.pcc.df$var = names(par.snails)[-which(names(par.snails) %in% 
                                           c('A', 'H', 'psi1', 'psi2', 'psi3', 'phi'))]

ggplot(It.pcc.df, aes(x = cor, y = original, label = var)) + 
  geom_point(pch = 16) +
  geom_text(aes(label = var), hjust = -0.5, vjust = -0.5)

ggplot(cbind(lhcpars.epi.fill, lhcpars.epi), aes(x = muW, y = I.t/area)) + 
  theme_bw() +
  geom_point(pch = 16, col = "red")+
  stat_smooth(method = "lm")

ggplot(cbind(lhcpars.epi.fill, lhcpars.epi), aes(x = muW, y = N.t/area)) + 
  theme_bw() +
  geom_point(pch = 16)+
  stat_smooth(method = "lm")


It.lhsprcc = ggplot(It.pcc.df, aes(x = var, y = original)) +
  theme_bw()+
  theme(axis.text = element_text(size = 12),  #increase axis label size
        axis.title = element_text(size = 15), #increase axis title size
        axis.title.x = element_blank()) + 
  scale_y_continuous(limits = c(-0.5,0.5), breaks = c(-0.5,0,0.5))+
  scale_x_discrete(labels = c("beta" = expression(italic(beta)),
                              "f" = expression(italic(f)),
                              "g1" = expression(g[1]),
                              "g2" = expression(italic(g[2])),
                              "Kn" = expression(italic(K)),
                              "lambda" = expression(italic(lambda)),
                              "m" = expression(italic(m)),
                              "muH" = expression(italic(mu[H])),
                              "muI" = expression(italic(mu[I])),
                              "muN1" = expression(italic(mu[N1])),
                              "muN2" = expression(italic(mu[N2])),
                              "muN3" = expression(italic(mu[N3])),
                              "muW" = expression(italic(mu[W])),
                              "sigma" = expression(italic(sigma)),
                              "theta" = expression(italic(theta)),
                              "z" = expression(italic(z)))) +
  geom_bar(fill = 'red', stat = 'identity', width = 0.25)+
  geom_hline(yintercept = 0, lty = 2)+
  annotate('text', x = 1.2, y = 0.5, label = 'B)', size = 8) +
  labs(x = 'Parameter', y = expression(paste('Infected snail density (',italic(N[I]) ,') PRCC')))

It.lhsprcc

# Assess PRCC of combined model (with intervention regime) outputs ###########
# Ending mean worm burden  
Wall.pcc = pcc(X = as.data.frame(lhcfin[,which(colnames(lhcfin) %in% 
                                                 names(par.all)[-which(names(par.all) %in% 
                                                                         c('A', 'H', 'psi1', 'psi2', 'psi3', 'phi'))])]),
               y = as.data.frame(lhcfin$W.all),
               rank = TRUE)

Wall.pcc.df = Wall.pcc$PRCC
Wall.pcc.df$cor = cor(lhcfin[,which(colnames(lhcfin) %in% 
                                      names(par.all)[-which(names(par.all) %in% 
                                                              c('A', 'H', 'psi1', 'psi2', 'psi3', 'phi'))])], 
                      lhcfin$W.all)
Wall.pcc.df$var = names(par.all)[-which(names(par.all) %in% 
                                          c('A', 'H', 'psi1', 'psi2', 'psi3', 'phi'))]

ggplot(Wall.pcc.df, aes(x = cor, y = original, label = var)) + 
  geom_point(pch = 16) +
  geom_text(aes(label = var), hjust = -0.5, vjust = -0.5)

Wall.lhsprcc = ggplot(Wall.pcc.df, aes(x = var, y = original)) +
  theme_bw()+
  theme(axis.text = element_text(size = 12),  #increase axis label size
        axis.title = element_text(size = 15)) + #increase axis title size
  scale_y_continuous(limits = c(-0.75,0.75), breaks = c(-0.75,-0.25,-0.5,0,0.25,0.5,0.75))+
  geom_bar(fill = 'purple', stat = 'identity', width = 0.25)+
  scale_x_discrete(labels = c("a.p" = expression(italic(a[P])),
                              "a.s" = expression(italic(a[N])),
                              "ar.slp" = expression(italic(alpha[m])),
                              "b.p" = expression(italic(b[P])),
                              "b.s" = expression(italic(b[N])),
                              "beta" = expression(italic(beta)),
                              "c" = expression(italic(c)),
                              "d" = expression(italic(d)),
                              "eta" = expression(italic(zeta)),
                              "f" = expression(italic(f)),
                              "g1" = expression(g[1]),
                              "g2" = expression(italic(g[2])),
                              "gam" = expression(italic(gamma)),
                              "k" = expression(italic(kappa)),
                              "Kn" = expression(italic(K)),
                              "lambda" = expression(italic(lambda)),
                              "linf" = expression(italic(L[infinity])),
                              "m" = expression(italic(m)),
                              "muH" = expression(italic(mu[H])),
                              "muI" = expression(italic(mu[I])),
                              "muN1" = expression(italic(mu[N1])),
                              "muN2" = expression(italic(mu[N2])),
                              "muN3" = expression(italic(mu[N3])),
                              "muP" = expression(italic(mu[P])),
                              "muW" = expression(italic(mu[W])),
                              "om" = expression(italic(omega)),
                              "p" = expression(italic(p)),
                              "sigma" = expression(italic(sigma)),
                              "th" = expression(italic(Th[m])),
                              "theta" = expression(italic(theta)),
                              "z" = expression(italic(z)))) +
  geom_hline(yintercept = 0, lty = 2)+
  annotate('text', x = 1.2, y = 0.75, label = 'C)', size = 8) +
  labs(x = 'Parameter', y = expression(paste('Final mean worm burden (', italic(W), ') PRCC')))

Wall.lhsprcc

#Ending infected snail density  
I_tall.pcc = pcc(X = as.data.frame(lhcfin[,which(colnames(lhcfin) %in% 
                                                   names(par.all)[-which(names(par.all) %in% 
                                                                           c('A', 'H', 'psi1', 'psi2', 'psi3'))])]),
                 y = as.data.frame(lhcfin$I.t.all),
                 rank = TRUE)

I_tall.pcc.df = I_tall.pcc$PRCC
I_tall.pcc.df$var = names(par.all)[-which(names(par.all) %in% 
                                            c('A', 'H', 'psi1', 'psi2', 'psi3'))]

I_tall.lhsprcc = ggplot(I_tall.pcc.df, aes(x = var, y = original)) +
  theme_bw()+
  scale_y_continuous(limits = c(-1,1), breaks = c(-1,-0.5,0,0.5,1))+
  geom_bar(fill = 'black', stat = 'identity', width = 0.25)+
  #geom_errorbar(x = var, )
  geom_hline(yintercept = 0, lty = 2)+
  labs(title = 'Infected snail density sensitivity', x = 'Parameter', y = 'PRCC')

I_tall.lhsprcc

#Ending total snail density  
N_tall.pcc = pcc(X = as.data.frame(lhcfin[,which(colnames(lhcfin) %in% 
                                                   names(par.all)[-which(names(par.all) %in% 
                                                                           c('A', 'H', 'psi1', 'psi2', 'psi3'))])]),
                 y = as.data.frame(colSums(lhcpars.all.sims[3651, 2:9, ])),
                 rank = TRUE)

N_tall.pcc.df = N_tall.pcc$PRCC
N_tall.pcc.df$var = names(par.all)[-which(names(par.all) %in% 
                                            c('A', 'H', 'psi1', 'psi2', 'psi3'))]

N_tall.lhsprcc = ggplot(N_tall.pcc.df, aes(x = var, y = original)) +
  theme_bw()+
  scale_y_continuous(limits = c(-0.75,0.75), breaks = c(-0.75,-0.25,-0.5,0,0.25,0.5,0.75))+
  geom_bar(fill = 'black', stat = 'identity', width = 0.25)+
  #geom_errorbar(x = var, )
  geom_hline(yintercept = 0, lty = 2)+
  labs(x = 'Parameter', y = 'Final snail density PRCC')

N_tall.lhsprcc

#Something weird going on here, PRCC is near 0 for all variables
plot(lhcfin$N.t, lhcfin$N.t.all, pch = 16, cex = 0.7)
plot(lhcfin$N.t.all, colSums(lhcpars.all.sims[3651, 2:9, ]), pch = 16, cex = 0.7,
     xlab = 'fin snail pop', ylab = 'median last year snail pop')

lhcfin$N.t.all.10yr = colSums(lhcpars.all.sims[3651, 2:9, ])
pairs(lhcfin[,c(1:8, 49)], pch = 18, cex = 0.7)
pairs(lhcfin[,c(1:8, 50)], pch = 18, cex = 0.7)

pairs(lhcfin[,c(11:17, 50)], pch = 18, cex = 0.7)
pairs(lhcfin[,c(21:34, 50)], pch = 18, cex = 0.7)

# Combine plots of key model outcomes ##################
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

fig1.layout = matrix(c(1,2,
                       3,3), ncol = 2, byrow = T)

windows(width = 350, height = 200)
multiplot(profit.lhsprcc, It.lhsprcc, 
          Wall.lhsprcc, 
          layout = fig1.layout)

