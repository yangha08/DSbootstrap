

##### PP1

set.seed(1)

sim.pp1.1 <- rpoispp(4000, win=disc(radius=0.05, centre=c(0.79, 0.79)))
sim.pp1.2 <- rpoispp(4000, win=disc(radius=0.05, centre=c(0.69, 0.69)))
sim.pp1.3 <- rpoispp(4000, win=disc(radius=0.05, centre=c(0.81, 0.65)))
sim.pp1.4 <- rpoispp(4000, win=disc(radius=0.15, centre=c(0.25, 0.25)))
sim.pp1.0 <- rpoispp(100)

sim.pp1 <- superimpose(sim.pp1.0, sim.pp1.1, sim.pp1.2, sim.pp1.3, sim.pp1.4)
sim.pp1$n   # 461
sim.pp1.df <- data.frame(x=sim.pp1$x, y=sim.pp1$y)

plot.pp1 = ggplot() +
  geom_point(mapping = aes(x = x, y = y), data = sim.pp1.df) +
  ggtitle("pp1(461 pts, cluster int = 4000, noise int = 100)") +
  theme(plot.title = element_text(hjust = 0.5))

plot.pp1

dens.pp1.alp1 = density.ppp(sim.pp1, dimyx= c(256, 256), 
                            varcov = Hpi.dir(sim.pp1.df, alpha = 1))
dens.pp1.alp2 = density.ppp(sim.pp1, dimyx= c(256, 256), 
                            varcov = Hpi.dir(sim.pp1.df, alpha = 2))
dens.pp1.alp3 = density.ppp(sim.pp1, dimyx= c(256, 256), 
                            varcov = Hpi.dir(sim.pp1.df, alpha = 3))
dens.pp1.alp5 = density.ppp(sim.pp1, dimyx= c(256, 256), 
                            varcov = Hpi.dir(sim.pp1.df, alpha = 5))
dens.pp1.alp10 = density.ppp(sim.pp1, dimyx= c(256, 256), 
                             varcov = Hpi.dir(sim.pp1.df, alpha = 10))
dens.pp1.alp20 = density.ppp(sim.pp1, dimyx= c(256, 256), 
                             varcov = Hpi.dir(sim.pp1.df, alpha = 20))
dens.pp1.alp50 = density.ppp(sim.pp1, dimyx= c(256, 256), 
                             varcov = Hpi.dir(sim.pp1.df, alpha = 50))
dens.pp1.alp100 = density.ppp(sim.pp1, dimyx= c(256, 256), 
                              varcov = Hpi.dir(sim.pp1.df, alpha = 100))

cont.pp1.1 <- plot_ly(z = ~dens.pp1.alp1$v, type = "contour") %>%
  layout(annotations = list(list(x = 0.37 , y = 1.1, text = "alpha=1", showarrow = F, xref='paper', yref='paper'))) 
cont.pp1.1 <- hide_colorbar(cont.pp1.1)
cont.pp1.2 <- plot_ly(z = ~dens.pp1.alp2$v, type = "contour") %>%
  layout(annotations = list(list(x = 0.5 , y = 1.1, text = "alpha=2", showarrow = F, xref='paper', yref='paper')))
cont.pp1.2 <- hide_colorbar(cont.pp1.2)
cont.pp1.3 <- plot_ly(z = ~dens.pp1.alp3$v, type = "contour") %>%
  layout(annotations = list(list(x = 0.5 , y = 1.1, text = "alpha=3", showarrow = F, xref='paper', yref='paper')))
cont.pp1.3 <- hide_colorbar(cont.pp1.3)
cont.pp1.5 <- plot_ly(z = ~dens.pp1.alp5$v, type = "contour") %>%
  layout(annotations = list(list(x = 0.62 , y = 1.1, text = "alpha=5", showarrow = F, xref='paper', yref='paper')))
cont.pp1.5 <- hide_colorbar(cont.pp1.5)
cont.pp1.10 <- plot_ly(z = ~dens.pp1.alp10$v, type = "contour") %>%
  layout(annotations = list(list(x = 0.37 , y = 1.05, text = "alpha=10", showarrow = F, xref='paper', yref='paper')))
cont.pp1.10 <- hide_colorbar(cont.pp1.10)
cont.pp1.20 <- plot_ly(z = ~dens.pp1.alp20$v, type = "contour") %>%
  layout(annotations = list(list(x = 0.5 , y = 1.05, text = "alpha=20", showarrow = F, xref='paper', yref='paper')))
cont.pp1.20 <- hide_colorbar(cont.pp1.20)
cont.pp1.50 <- plot_ly(z = ~dens.pp1.alp50$v, type = "contour") %>%
  layout(annotations = list(list(x = 0.5 , y = 1.05, text = "alpha=50", showarrow = F, xref='paper', yref='paper')))
cont.pp1.50 <- hide_colorbar(cont.pp1.50)
cont.pp1.100 <- plot_ly(z = ~dens.pp1.alp100$v, type = "contour") %>%
  layout(annotations = list(list(x = 0.67 , y = 1.05, text = "alpha=100", showarrow = F, xref='paper', yref='paper'))) %>%
  colorbar(title = "Intensity")


cont.pp1 <- subplot(cont.pp1.1, cont.pp1.2, cont.pp1.3, cont.pp1.5,
                    cont.pp1.10, cont.pp1.20, cont.pp1.50, cont.pp1.100, nrows = 2, widths = c(0.235, 0.265, 0.265, 0.235), 
                    shareX = T, shareY = T, margin = 0.03) %>% 
  layout(title = 'pp1, cluster int = 4000, noise int = 100', margin = 0.01)

cont.pp1




##### PP2

set.seed(1)

sim.pp2.1 <- rpoispp(4000, win=disc(radius=0.05, centre=c(0.79, 0.79)))
sim.pp2.2 <- rpoispp(4000, win=disc(radius=0.05, centre=c(0.69, 0.69)))
sim.pp2.3 <- rpoispp(4000, win=disc(radius=0.05, centre=c(0.81, 0.65)))
sim.pp2.4 <- rpoispp(4000, win=disc(radius=0.15, centre=c(0.25, 0.25)))
sim.pp2.0 <- rpoispp(100)

sim.pp2 <- superimpose(sim.pp2.0, sim.pp2.1, sim.pp2.2, sim.pp2.3, sim.pp2.4)
sim.pp2$n   # 748
sim.pp2.df <- data.frame(x=sim.pp2$x, y=sim.pp2$y)

plot.pp2 = ggplot() +
  geom_point(mapping = aes(x = x, y = y), data = sim.pp2.df) +
  ggtitle("pp2(748 pts, cluster int = 4000, noise int = 100)") +
  theme(plot.title = element_text(hjust = 0.5))

plot.pp2

dens.pp2.alp1 = density.ppp(sim.pp2, dimyx= c(256, 256), 
                            varcov = Hpi.dir(sim.pp2.df, alpha = 1))
dens.pp2.alp2 = density.ppp(sim.pp2, dimyx= c(256, 256), 
                            varcov = Hpi.dir(sim.pp2.df, alpha = 2))
dens.pp2.alp3 = density.ppp(sim.pp2, dimyx= c(256, 256), 
                            varcov = Hpi.dir(sim.pp2.df, alpha = 3))
dens.pp2.alp5 = density.ppp(sim.pp2, dimyx= c(256, 256), 
                            varcov = Hpi.dir(sim.pp2.df, alpha = 5))
dens.pp2.alp10 = density.ppp(sim.pp2, dimyx= c(256, 256), 
                             varcov = Hpi.dir(sim.pp2.df, alpha = 10))
dens.pp2.alp20 = density.ppp(sim.pp2, dimyx= c(256, 256), 
                             varcov = Hpi.dir(sim.pp2.df, alpha = 20))
dens.pp2.alp50 = density.ppp(sim.pp2, dimyx= c(256, 256), 
                             varcov = Hpi.dir(sim.pp2.df, alpha = 50))
dens.pp2.alp100 = density.ppp(sim.pp2, dimyx= c(256, 256), 
                              varcov = Hpi.dir(sim.pp2.df, alpha = 100))


cont.pp2.1 <- plot_ly(z = ~dens.pp2.alp1$v, type = "contour") %>%
  layout(annotations = list(list(x = 0.37 , y = 1.1, text = "alpha=1", showarrow = F, xref='paper', yref='paper'))) 
cont.pp2.1 <- hide_colorbar(cont.pp2.1)
cont.pp2.2 <- plot_ly(z = ~dens.pp2.alp2$v, type = "contour") %>%
  layout(annotations = list(list(x = 0.5 , y = 1.1, text = "alpha=2", showarrow = F, xref='paper', yref='paper')))
cont.pp2.2 <- hide_colorbar(cont.pp2.2)
cont.pp2.3 <- plot_ly(z = ~dens.pp2.alp3$v, type = "contour") %>%
  layout(annotations = list(list(x = 0.5 , y = 1.1, text = "alpha=3", showarrow = F, xref='paper', yref='paper')))
cont.pp2.3 <- hide_colorbar(cont.pp2.3)
cont.pp2.5 <- plot_ly(z = ~dens.pp2.alp5$v, type = "contour") %>%
  layout(annotations = list(list(x = 0.62 , y = 1.1, text = "alpha=5", showarrow = F, xref='paper', yref='paper')))
cont.pp2.5 <- hide_colorbar(cont.pp2.5)
cont.pp2.10 <- plot_ly(z = ~dens.pp2.alp10$v, type = "contour") %>%
  layout(annotations = list(list(x = 0.37 , y = 1.05, text = "alpha=10", showarrow = F, xref='paper', yref='paper')))
cont.pp2.10 <- hide_colorbar(cont.pp2.10)
cont.pp2.20 <- plot_ly(z = ~dens.pp2.alp20$v, type = "contour") %>%
  layout(annotations = list(list(x = 0.5 , y = 1.05, text = "alpha=20", showarrow = F, xref='paper', yref='paper')))
cont.pp2.20 <- hide_colorbar(cont.pp2.20)
cont.pp2.50 <- plot_ly(z = ~dens.pp2.alp50$v, type = "contour") %>%
  layout(annotations = list(list(x = 0.5 , y = 1.05, text = "alpha=50", showarrow = F, xref='paper', yref='paper')))
cont.pp2.50 <- hide_colorbar(cont.pp2.50)
cont.pp2.100 <- plot_ly(z = ~dens.pp2.alp100$v, type = "contour") %>%
  layout(annotations = list(list(x = 0.67 , y = 1.05, text = "alpha=100", showarrow = F, xref='paper', yref='paper'))) %>%
  colorbar(title = "Intensity")


cont.pp2 <- subplot(cont.pp2.1, cont.pp2.2, cont.pp2.3, cont.pp2.5,
                    cont.pp2.10, cont.pp2.20, cont.pp2.50, cont.pp2.100, nrows = 2, widths = c(0.235, 0.265, 0.265, 0.235), 
                    shareX = T, shareY = T, margin = 0.03) %>% 
  layout(title = 'pp2, cluster int = 4000, noise int = 400', margin = 0.01)

cont.pp2


