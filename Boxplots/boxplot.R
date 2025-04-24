require(ggplot2)
require(magrittr)
require(gridExtra)
load("Result1.RData")
load("Result1th.RData")
load("Result2.RData")
load("Result2th.RData")
load("Result3.RData")
load("Result4.RData")
load("Result5.RData")
load("Result6.RData")

class <- rep(c("Diggle", "Scott", "Flscv", "lscv",  "lcv",
               "F, SCV, Ori", "F, PI, Ori", 
               "Oracle", "Empirical",
               "F, SCV, F-S", "F, PI, F-S",
               "D, SCV, F-S", "D, PI, F-S",
               "F, SCV, 1", "F, PI, 1",
               "F, SCV, 5", "F, PI, 5"), each = 100)

####### Lambda 1, 2 #######


ISE1 <- c(ISE.1b.diggle,
          ISE.1b.scott,
          ISE.1b.mat.lscv,
          ISE.1b.diag.lscv,
          ISE.1b.lcv,
          ISE.mat.1b.scv,
          ISE.mat.1b.pi ,
          ISE.1b.oracle,
          ISE.1b.empboot,
          ISE.mat.1b.scv.FS,
          ISE.diag.1b.scv.FS,
          ISE.mat.1b.pi.FS,
          ISE.diag.1b.pi.FS,
          ISE.mat.1b.SB.scv.1 ,
          ISE.mat.1b.SB.pi.1 ,
          ISE.mat.1b.SB.scv.5 ,
          ISE.mat.1b.SB.pi.5)


MISE1 <- c(mean(ISE.1b.diggle),
           mean(ISE.1b.scott),
           mean(ISE.1b.mat.lscv),
           mean(ISE.1b.diag.lscv),
           mean(ISE.1b.lcv),
           mean(ISE.mat.1b.scv),
           mean(ISE.mat.1b.pi),
           mean(ISE.1b.oracle),
           mean(ISE.1b.empboot),
           mean(ISE.mat.1b.scv.FS),
           mean(ISE.diag.1b.scv.FS),
           mean(ISE.mat.1b.pi.FS),
           mean(ISE.diag.1b.pi.FS),
           mean(ISE.mat.1b.SB.scv.1),
           mean(ISE.mat.1b.SB.pi.1),
           mean(ISE.mat.1b.SB.scv.5),
           mean(ISE.mat.1b.SB.pi.5))




ISE1df = data.frame(ISE1, class)

ISE1df$class <- factor(ISE1df$class, levels = c("Diggle", "Scott", "Flscv", "lscv",  "lcv",
                                                "F, SCV, Ori", "F, PI, Ori", 
                                                "Oracle", "Empirical",
                                                "F, SCV, F-S", "F, PI, F-S",
                                                "D, SCV, F-S", "D, PI, F-S",
                                                "F, SCV, 1", "F, PI, 1",
                                                "F, SCV, 5", "F, PI, 5"))

which(MISE1[-8] == min(MISE1[-8]))


f1 = ISE1df %>% 
  
  # Build the boxplot. In the 'fill' argument, give this column
  ggplot( aes(x=class, y=ISE1, fill=class, color=class, alpha=class)) + 
  labs(title = "Lambda 1") +
  theme_bw() +
  theme(text = element_text(size = 9.5), plot.title = element_text(hjust = 0.5)) +
  geom_boxplot() +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  scale_alpha_manual(values=rep(0.7,18)) +
  scale_fill_manual(values=c("chocolate", "chocolate", "chocolate", "chocolate", "chocolate",
                             "darkolivegreen4", "darkolivegreen1", 
                             "darkgoldenrod1", "darkgoldenrod1", 
                             "darkslategray4", "darkslategray1", "darkslategray4", "darkslategray1",
                             "darkorchid4", "darkorchid1", 
                             "darkorchid4", "darkorchid1")) +
  scale_color_manual(values = c(rep("black",10),"red",rep("black",6)), guide = "none") +
  theme(legend.position = "none") +
  xlab("") +
  ylab("ISE") +
  ylim(0, 0.1)



ISE2 <- c(ISE.2b.diggle,
          ISE.2b.scott,
          ISE.2b.mat.lscv,
          ISE.2b.diag.lscv,
          ISE.2b.lcv,
          ISE.mat.2b.scv,
          ISE.mat.2b.pi ,
          ISE.2b.oracle,
          ISE.2b.empboot,
          ISE.mat.2b.scv.FS,
          ISE.diag.2b.scv.FS,
          ISE.mat.2b.pi.FS,
          ISE.diag.2b.pi.FS,
          ISE.mat.2b.SB.scv.1 ,
          ISE.mat.2b.SB.pi.1 ,
          ISE.mat.2b.SB.scv.5 ,
          ISE.mat.2b.SB.pi.5)


MISE2 <- c(mean(ISE.2b.diggle),
           mean(ISE.2b.scott),
           mean(ISE.2b.mat.lscv),
           mean(ISE.2b.diag.lscv),
           mean(ISE.2b.lcv),
           mean(ISE.mat.2b.scv),
           mean(ISE.mat.2b.pi),
           mean(ISE.2b.oracle),
           mean(ISE.2b.empboot),
           mean(ISE.mat.2b.scv.FS),
           mean(ISE.diag.2b.scv.FS),
           mean(ISE.mat.2b.pi.FS),
           mean(ISE.diag.2b.pi.FS),
           mean(ISE.mat.2b.SB.scv.1),
           mean(ISE.mat.2b.SB.pi.1),
           mean(ISE.mat.2b.SB.scv.5),
           mean(ISE.mat.2b.SB.pi.5))




ISE2df = data.frame(ISE2, class)

ISE2df$class <- factor(ISE2df$class, levels = c("Diggle", "Scott", "Flscv", "lscv",  "lcv",
                                                "F, SCV, Ori", "F, PI, Ori", 
                                                "Oracle", "Empirical",
                                                "F, SCV, F-S", "F, PI, F-S",
                                                "D, SCV, F-S", "D, PI, F-S",
                                                "F, SCV, 1", "F, PI, 1",
                                                "F, SCV, 5", "F, PI, 5"))

which(MISE2[-8] == min(MISE2[-8]))


f2 = ISE2df %>% 
  
  # Build the boxplot. In the 'fill' argument, give this column
  ggplot( aes(x=class, y=ISE2, fill=class, color=class, alpha=class)) + 
  labs(title = "Lambda 2") +
  theme_bw() +
  theme(text = element_text(size = 9.5), plot.title = element_text(hjust = 0.5)) +
  geom_boxplot() +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  scale_alpha_manual(values=rep(0.7,18)) +
  scale_fill_manual(values=c("chocolate", "chocolate", "chocolate", "chocolate", "chocolate",
                             "darkolivegreen4", "darkolivegreen1", 
                             "darkgoldenrod1", "darkgoldenrod1", 
                             "darkslategray4", "darkslategray1", "darkslategray4", "darkslategray1",
                             "darkorchid4", "darkorchid1", 
                             "darkorchid4", "darkorchid1")) +
  scale_color_manual(values = c(rep("black",13),"red",rep("black",3)), guide = "none") +
  theme(legend.position = "none") +
  xlab("") +
  ylab("ISE") +
  ylim(0, 0.2)


grid.arrange(f1, f2, nrow=2, ncol=1)



####### Lambda 1, 2(Thomas) #######

ISE1th <- c(ISE.1b.thomas.diggle,
            ISE.1b.thomas.scott,
            ISE.1b.thomas.mat.lscv,
            ISE.1b.thomas.diag.lscv,
            ISE.1b.thomas.lcv,
            ISE.mat.1b.thomas.scv,
            ISE.mat.1b.thomas.pi ,
            ISE.1b.thomas.oracle,
            ISE.1b.thomas.empboot,
            ISE.mat.1b.thomas.scv.FS,
            ISE.diag.1b.thomas.scv.FS,
            ISE.mat.1b.thomas.pi.FS,
            ISE.diag.1b.thomas.pi.FS,
            ISE.mat.1b.thomas.SB.scv.1 ,
            ISE.mat.1b.thomas.SB.pi.1 ,
            ISE.mat.1b.thomas.SB.scv.5 ,
            ISE.mat.1b.thomas.SB.pi.5)


MISE1th <- c(mean(ISE.1b.thomas.diggle),
             mean(ISE.1b.thomas.scott),
             mean(ISE.1b.thomas.mat.lscv),
             mean(ISE.1b.thomas.diag.lscv),
             mean(ISE.1b.thomas.lcv),
             mean(ISE.mat.1b.thomas.scv),
             mean(ISE.mat.1b.thomas.pi),
             mean(ISE.1b.thomas.oracle),
             mean(ISE.1b.thomas.empboot),
             mean(ISE.mat.1b.thomas.scv.FS),
             mean(ISE.diag.1b.thomas.scv.FS),
             mean(ISE.mat.1b.thomas.pi.FS),
             mean(ISE.diag.1b.thomas.pi.FS),
             mean(ISE.mat.1b.thomas.SB.scv.1),
             mean(ISE.mat.1b.thomas.SB.pi.1),
             mean(ISE.mat.1b.thomas.SB.scv.5),
             mean(ISE.mat.1b.thomas.SB.pi.5))



ISE1thdf = data.frame(ISE1th, class)

ISE1thdf$class <- factor(ISE1thdf$class, levels = c("Diggle", "Scott", "Flscv", "lscv",  "lcv",
                                                    "F, SCV, Ori", "F, PI, Ori", 
                                                    "Oracle", "Empirical",
                                                    "F, SCV, F-S", "F, PI, F-S",
                                                    "D, SCV, F-S", "D, PI, F-S",
                                                    "F, SCV, 1", "F, PI, 1",
                                                    "F, SCV, 5", "F, PI, 5"))

which(MISE1th[-8] == min(MISE1th[-8]))


f1th = ISE1thdf %>% 
  
  # Build the boxplot. In the 'fill' argument, give this column
  ggplot( aes(x=class, y=ISE1th, fill=class, color=class, alpha=class)) + 
  labs(title = "Lambda 1(Thomas)") +
  theme_bw() +
  theme(text = element_text(size = 9.5), plot.title = element_text(hjust = 0.5)) +
  geom_boxplot() +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  scale_alpha_manual(values=rep(0.7,18)) +
  scale_fill_manual(values=c("chocolate", "chocolate", "chocolate", "chocolate", "chocolate",
                             "darkolivegreen4", "darkolivegreen1", 
                             "darkgoldenrod1", "darkgoldenrod1", 
                             "darkslategray4", "darkslategray1", "darkslategray4", "darkslategray1",
                             "darkorchid4", "darkorchid1", 
                             "darkorchid4", "darkorchid1")) +
  scale_color_manual(values = c(rep("black",11),"red",rep("black",5)), guide = "none") +
  theme(legend.position = "none") +
  xlab("") +
  ylab("ISE") +
  ylim(0, 0.5)



ISE2th <- c(ISE.2b.thomas.diggle,
            ISE.2b.thomas.scott,
            ISE.2b.thomas.mat.lscv,
            ISE.2b.thomas.diag.lscv,
            ISE.2b.thomas.lcv,
            ISE.mat.2b.thomas.scv,
            ISE.mat.2b.thomas.pi ,
            ISE.2b.thomas.oracle,
            ISE.2b.thomas.empboot,
            ISE.mat.2b.thomas.scv.FS,
            ISE.diag.2b.thomas.scv.FS,
            ISE.mat.2b.thomas.pi.FS,
            ISE.diag.2b.thomas.pi.FS,
            ISE.mat.2b.thomas.SB.scv.1 ,
            ISE.mat.2b.thomas.SB.pi.1 ,
            ISE.mat.2b.thomas.SB.scv.5 ,
            ISE.mat.2b.thomas.SB.pi.5)


MISE2th <- c(mean(ISE.2b.thomas.diggle),
             mean(ISE.2b.thomas.scott),
             mean(ISE.2b.thomas.mat.lscv),
             mean(ISE.2b.thomas.diag.lscv),
             mean(ISE.2b.thomas.lcv),
             mean(ISE.mat.2b.thomas.scv),
             mean(ISE.mat.2b.thomas.pi),
             mean(ISE.2b.thomas.oracle),
             mean(ISE.2b.thomas.empboot),
             mean(ISE.mat.2b.thomas.scv.FS),
             mean(ISE.diag.2b.thomas.scv.FS),
             mean(ISE.mat.2b.thomas.pi.FS),
             mean(ISE.diag.2b.thomas.pi.FS),
             mean(ISE.mat.2b.thomas.SB.scv.1),
             mean(ISE.mat.2b.thomas.SB.pi.1),
             mean(ISE.mat.2b.thomas.SB.scv.5),
             mean(ISE.mat.2b.thomas.SB.pi.5))




ISE2thdf = data.frame(ISE2th, class)

ISE2thdf$class <- factor(ISE2thdf$class, levels = c("Diggle", "Scott", "Flscv", "lscv",  "lcv",
                                                    "F, SCV, Ori", "F, PI, Ori", 
                                                    "Oracle", "Empirical",
                                                    "F, SCV, F-S", "F, PI, F-S",
                                                    "D, SCV, F-S", "D, PI, F-S",
                                                    "F, SCV, 1", "F, PI, 1",
                                                    "F, SCV, 5", "F, PI, 5"))

which(MISE2th[-8] == min(MISE2th[-8]))


f2th = ISE2thdf %>% 
  
  # Build the boxplot. In the 'fill' argument, give this column
  ggplot( aes(x=class, y=ISE2th, fill=class, color=class, alpha=class)) + 
  labs(title = "Lambda 2(Thomas)") +
  theme_bw() +
  theme(text = element_text(size = 9.5), plot.title = element_text(hjust = 0.5)) +
  geom_boxplot() +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  scale_alpha_manual(values=rep(0.7,18)) +
  scale_fill_manual(values=c("chocolate", "chocolate", "chocolate", "chocolate", "chocolate",
                             "darkolivegreen4", "darkolivegreen1", 
                             "darkgoldenrod1", "darkgoldenrod1", 
                             "darkslategray4", "darkslategray1", "darkslategray4", "darkslategray1",
                             "darkorchid4", "darkorchid1", 
                             "darkorchid4", "darkorchid1")) +
  scale_color_manual(values = c(rep("black",13),"red",rep("black",3)), guide = "none") +
  theme(legend.position = "none") +
  xlab("") +
  ylab("ISE") +
  ylim(0, 0.4)


grid.arrange(f1th, f2th, nrow=2, ncol=1)



####### Lambda 3, 4 #######



ISE3 <- c(ISE.3b.diggle,
          ISE.3b.scott,
          ISE.3b.mat.lscv,
          ISE.3b.diag.lscv,
          ISE.3b.lcv,
          ISE.mat.3b.scv,
          ISE.mat.3b.pi ,
          ISE.3b.oracle,
          ISE.3b.empboot,
          ISE.mat.3b.scv.FS,
          ISE.diag.3b.scv.FS,
          ISE.mat.3b.pi.FS,
          ISE.diag.3b.pi.FS,
          ISE.mat.3b.SB.scv.1 ,
          ISE.mat.3b.SB.pi.1 ,
          ISE.mat.3b.SB.scv.5 ,
          ISE.mat.3b.SB.pi.5)


MISE3 <- c(mean(ISE.3b.diggle),
           mean(ISE.3b.scott),
           mean(ISE.3b.mat.lscv),
           mean(ISE.3b.diag.lscv),
           mean(ISE.3b.lcv),
           mean(ISE.mat.3b.scv),
           mean(ISE.mat.3b.pi),
           mean(ISE.3b.oracle),
           mean(ISE.3b.empboot),
           mean(ISE.mat.3b.scv.FS),
           mean(ISE.diag.3b.scv.FS),
           mean(ISE.mat.3b.pi.FS),
           mean(ISE.diag.3b.pi.FS),
           mean(ISE.mat.3b.SB.scv.1),
           mean(ISE.mat.3b.SB.pi.1),
           mean(ISE.mat.3b.SB.scv.5),
           mean(ISE.mat.3b.SB.pi.5))




ISE3df = data.frame(ISE3, class)

ISE3df$class <- factor(ISE3df$class, levels = c("Diggle", "Scott", "Flscv", "lscv",  "lcv",
                                                "F, SCV, Ori", "F, PI, Ori", 
                                                "Oracle", "Empirical",
                                                "F, SCV, F-S", "F, PI, F-S",
                                                "D, SCV, F-S", "D, PI, F-S",
                                                "F, SCV, 1", "F, PI, 1",
                                                "F, SCV, 5", "F, PI, 5"))

which(MISE3[-8] == min(MISE3[-8]))


f3 = ISE3df %>% 
  
  # Build the boxplot. In the 'fill' argument, give this column
  ggplot( aes(x=class, y=ISE3, fill=class, color=class, alpha=class)) + 
  labs(title = "Lambda 3") +
  theme_bw() +
  theme(text = element_text(size = 9.5), plot.title = element_text(hjust = 0.5)) +
  geom_boxplot() +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  scale_alpha_manual(values=rep(0.7,18)) +
  scale_fill_manual(values=c("chocolate", "chocolate", "chocolate", "chocolate", "chocolate",
                             "darkolivegreen4", "darkolivegreen1", 
                             "darkgoldenrod1", "darkgoldenrod1", 
                             "darkslategray4", "darkslategray1", "darkslategray4", "darkslategray1",
                             "darkorchid4", "darkorchid1", 
                             "darkorchid4", "darkorchid1")) +
  scale_color_manual(values = c(rep("black",5),"red",rep("black",11)), guide = "none") +
  theme(legend.position = "none") +
  xlab("") +
  ylab("ISE") +
  ylim(0, 4)



ISE4 <- c(ISE.4b.diggle,
          ISE.4b.scott,
          ISE.4b.mat.lscv,
          ISE.4b.diag.lscv,
          ISE.4b.lcv,
          ISE.mat.4b.scv,
          ISE.mat.4b.pi ,
          ISE.4b.oracle,
          ISE.4b.empboot,
          ISE.mat.4b.scv.FS,
          ISE.diag.4b.scv.FS,
          ISE.mat.4b.pi.FS,
          ISE.diag.4b.pi.FS,
          ISE.mat.4b.SB.scv.1 ,
          ISE.mat.4b.SB.pi.1 ,
          ISE.mat.4b.SB.scv.5 ,
          ISE.mat.4b.SB.pi.5)


MISE4 <- c(mean(ISE.4b.diggle),
           mean(ISE.4b.scott),
           mean(ISE.4b.mat.lscv),
           mean(ISE.4b.diag.lscv),
           mean(ISE.4b.lcv),
           mean(ISE.mat.4b.scv),
           mean(ISE.mat.4b.pi),
           mean(ISE.4b.oracle),
           mean(ISE.4b.empboot),
           mean(ISE.mat.4b.scv.FS),
           mean(ISE.diag.4b.scv.FS),
           mean(ISE.mat.4b.pi.FS),
           mean(ISE.diag.4b.pi.FS),
           mean(ISE.mat.4b.SB.scv.1),
           mean(ISE.mat.4b.SB.pi.1),
           mean(ISE.mat.4b.SB.scv.5),
           mean(ISE.mat.4b.SB.pi.5))




ISE4df = data.frame(ISE4, class)

ISE4df$class <- factor(ISE4df$class, levels = c("Diggle", "Scott", "Flscv", "lscv",  "lcv",
                                                "F, SCV, Ori", "F, PI, Ori", 
                                                "Oracle", "Empirical",
                                                "F, SCV, F-S", "F, PI, F-S",
                                                "D, SCV, F-S", "D, PI, F-S",
                                                "F, SCV, 1", "F, PI, 1",
                                                "F, SCV, 5", "F, PI, 5"))

which(MISE4[-8] == min(MISE4[-8]))


f4 = ISE4df %>% 
  
  # Build the boxplot. In the 'fill' argument, give this column
  ggplot( aes(x=class, y=ISE4, fill=class, color=class, alpha=class)) + 
  labs(title = "Lambda 4") +
  theme_bw() +
  theme(text = element_text(size = 9.5), plot.title = element_text(hjust = 0.5)) +
  geom_boxplot() +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  scale_alpha_manual(values=rep(0.7,18)) +
  scale_fill_manual(values=c("chocolate", "chocolate", "chocolate", "chocolate", "chocolate",
                             "darkolivegreen4", "darkolivegreen1", 
                             "darkgoldenrod1", "darkgoldenrod1", 
                             "darkslategray4", "darkslategray1", "darkslategray4", "darkslategray1",
                             "darkorchid4", "darkorchid1", 
                             "darkorchid4", "darkorchid1")) +
  scale_color_manual(values = c(rep("black",6),"red",rep("black",10)), guide = "none") +
  theme(legend.position = "none") +
  xlab("") +
  ylab("ISE") +
  ylim(0.075, 0.2)


grid.arrange(f3, f4, nrow=2, ncol=1)





####### Lambda 5, 6 #######


ISE5 <- c(ISE.5b.diggle,
          ISE.5b.scott,
          ISE.5b.mat.lscv,
          ISE.5b.diag.lscv,
          ISE.5b.lcv,
          ISE.mat.5b.scv,
          ISE.mat.5b.pi ,
          ISE.5b.oracle,
          ISE.5b.empboot,
          ISE.mat.5b.scv.FS,
          ISE.diag.5b.scv.FS,
          ISE.mat.5b.pi.FS,
          ISE.diag.5b.pi.FS,
          ISE.mat.5b.SB.scv.1 ,
          ISE.mat.5b.SB.pi.1 ,
          ISE.mat.5b.SB.scv.5 ,
          ISE.mat.5b.SB.pi.5)


MISE5 <- c(mean(ISE.5b.diggle),
           mean(ISE.5b.scott),
           mean(ISE.5b.mat.lscv),
           mean(ISE.5b.diag.lscv),
           mean(ISE.5b.lcv),
           mean(ISE.mat.5b.scv),
           mean(ISE.mat.5b.pi),
           mean(ISE.5b.oracle),
           mean(ISE.5b.empboot),
           mean(ISE.mat.5b.scv.FS),
           mean(ISE.diag.5b.scv.FS),
           mean(ISE.mat.5b.pi.FS),
           mean(ISE.diag.5b.pi.FS),
           mean(ISE.mat.5b.SB.scv.1),
           mean(ISE.mat.5b.SB.pi.1),
           mean(ISE.mat.5b.SB.scv.5),
           mean(ISE.mat.5b.SB.pi.5))




ISE5df = data.frame(ISE5, class)

ISE5df$class <- factor(ISE5df$class, levels = c("Diggle", "Scott", "Flscv", "lscv",  "lcv",
                                                "F, SCV, Ori", "F, PI, Ori", 
                                                "Oracle", "Empirical",
                                                "F, SCV, F-S", "F, PI, F-S",
                                                "D, SCV, F-S", "D, PI, F-S",
                                                "F, SCV, 1", "F, PI, 1",
                                                "F, SCV, 5", "F, PI, 5"))

which(MISE5[-8] == min(MISE5[-8]))


f5 = ISE5df %>% 
  
  # Build the boxplot. In the 'fill' argument, give this column
  ggplot( aes(x=class, y=ISE5, fill=class, color=class, alpha=class)) + 
  labs(title = "Lambda 5") +
  theme_bw() +
  theme(text = element_text(size = 9.5), plot.title = element_text(hjust = 0.5)) +
  geom_boxplot() +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  scale_alpha_manual(values=rep(0.7,18)) +
  scale_fill_manual(values=c("chocolate", "chocolate", "chocolate", "chocolate", "chocolate",
                             "darkolivegreen4", "darkolivegreen1", 
                             "darkgoldenrod1", "darkgoldenrod1", 
                             "darkslategray4", "darkslategray1", "darkslategray4", "darkslategray1",
                             "darkorchid4", "darkorchid1", 
                             "darkorchid4", "darkorchid1")) +
  scale_color_manual(values = c("red",rep("black",16)), guide = "none") +
  theme(legend.position = "none") +
  xlab("") +
  ylab("ISE") +
  ylim(0, 0.3)



ISE6 <- c(ISE.6b.diggle,
          ISE.6b.scott,
          ISE.6b.mat.lscv,
          ISE.6b.diag.lscv,
          ISE.6b.lcv,
          ISE.mat.6b.scv,
          ISE.mat.6b.pi ,
          ISE.6b.oracle,
          ISE.6b.empboot,
          ISE.mat.6b.scv.FS,
          ISE.diag.6b.scv.FS,
          ISE.mat.6b.pi.FS,
          ISE.diag.6b.pi.FS,
          ISE.mat.6b.SB.scv.1 ,
          ISE.mat.6b.SB.pi.1 ,
          ISE.mat.6b.SB.scv.5 ,
          ISE.mat.6b.SB.pi.5)


MISE6 <- c(mean(ISE.6b.diggle),
           mean(ISE.6b.scott),
           mean(ISE.6b.mat.lscv),
           mean(ISE.6b.diag.lscv),
           mean(ISE.6b.lcv),
           mean(ISE.mat.6b.scv),
           mean(ISE.mat.6b.pi),
           mean(ISE.6b.oracle),
           mean(ISE.6b.empboot),
           mean(ISE.mat.6b.scv.FS),
           mean(ISE.diag.6b.scv.FS),
           mean(ISE.mat.6b.pi.FS),
           mean(ISE.diag.6b.pi.FS),
           mean(ISE.mat.6b.SB.scv.1),
           mean(ISE.mat.6b.SB.pi.1),
           mean(ISE.mat.6b.SB.scv.5),
           mean(ISE.mat.6b.SB.pi.5))




ISE6df = data.frame(ISE6, class)

ISE6df$class <- factor(ISE6df$class, levels = c("Diggle", "Scott", "Flscv", "lscv",  "lcv",
                                                "F, SCV, Ori", "F, PI, Ori", 
                                                "Oracle", "Empirical",
                                                "F, SCV, F-S", "F, PI, F-S",
                                                "D, SCV, F-S", "D, PI, F-S",
                                                "F, SCV, 1", "F, PI, 1",
                                                "F, SCV, 5", "F, PI, 5"))

which(MISE6[-8] == min(MISE6[-8]))


f6 = ISE6df %>% 
  
  # Build the boxplot. In the 'fill' argument, give this column
  ggplot( aes(x=class, y=ISE6, fill=class, color=class, alpha=class)) + 
  labs(title = "Lambda 6") +
  theme_bw() +
  theme(text = element_text(size = 9.5), plot.title = element_text(hjust = 0.5)) +
  geom_boxplot() +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  scale_alpha_manual(values=rep(0.7,18)) +
  scale_fill_manual(values=c("chocolate", "chocolate", "chocolate", "chocolate", "chocolate",
                             "darkolivegreen4", "darkolivegreen1", 
                             "darkgoldenrod1", "darkgoldenrod1", 
                             "darkslategray4", "darkslategray1", "darkslategray4", "darkslategray1",
                             "darkorchid4", "darkorchid1", 
                             "darkorchid4", "darkorchid1")) +
  scale_color_manual(values = c(rep("black",5),"red",rep("black",11)), guide = "none") +
  theme(legend.position = "none") +
  xlab("") +
  ylab("ISE") +
  ylim(0, 0.2)


grid.arrange(f5, f6, nrow=2, ncol=1)

