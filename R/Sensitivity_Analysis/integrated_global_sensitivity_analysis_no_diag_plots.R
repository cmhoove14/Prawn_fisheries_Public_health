load("Sensitivity_Analysis/sens_sims_no_diag.RData")
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
  scale_y_continuous(limits = c(-0.85,0.85), breaks = seq(-0.8, 0.8, 0.2))+
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
  annotate('text', x = 1.2, y = 0.8, label = 'A)', size = 8) +
  labs(x = 'Parameter', 
       y = expression(paste('Profit (', Pi,'(T,',P[0], ')) PRCC')))

profit.lhsprcc


#Equilibrium total snail population ###########
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
        axis.title = element_text(size = 15), #increase axis title size
        axis.title.x = element_blank()) + 
  scale_y_continuous(limits = c(-0.85,0.85), breaks = seq(-0.8, 0.8, 0.2))+
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
                              "theta1" = expression(italic(theta[1])),
                              "theta2" = expression(italic(theta[2])),
                              "z" = expression(italic(z)))) +
  geom_bar(fill = 'black', stat = 'identity', width = 0.25)+
  annotate('text', x = 1.2, y = 0.8, label = "B)", size = 8)  +
  geom_hline(yintercept = 0, lty = 2)+
  labs(x = 'Parameter', y = 'Total Snail Density (N) PRCC')

Nt.lhsprcc


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
  scale_y_continuous(limits = c(-0.85,0.85), breaks = seq(-0.8, 0.8, 0.2))+
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
                              "theta1" = expression(italic(theta[1])),
                              "theta2" = expression(italic(theta[2])),
                              "z" = expression(italic(z)))) +
  geom_hline(yintercept = 0, lty = 2)+
  annotate('text', x = 2, y = 0.8, label = 'C)', size = 8) +
  labs(x = 'Parameter', y = expression(paste('Final mean worm burden (', italic(W), ') PRCC')))

Wall.lhsprcc


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
multiplot(profit.lhsprcc, Nt.lhsprcc, 
          Wall.lhsprcc, 
          layout = fig1.layout)
