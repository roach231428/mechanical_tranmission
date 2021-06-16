remove(list = ls())
library(deSolve)
library(magrittr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(dplyr)


##### Functions #####


## main simulation ODE
dengue_simulation = function(t, x, paras){
  with(as.list(c(paras, x)),{
    S_H = x["S_H"]; E_H = x["E_H"]; I_H = x["I_H"]; R_H = x["R_H"]; C_H = x["C_H"];
    S_V = x["S_V"]; M_V = x["M_V"]; I_V = x["I_V"]; E_V = x["E_V"];
    vec = x["vec"]; mec = x["mec"]
    
    lambda_V = b * p * I_V / N_H
    
    if(mechanical) lambda_M = b_M * p_M * M_V / N_H # flag of whether considering mechanical transmission
    else lambda_M = 0 
    
    dS_H = (N_H - S_H) * mu_H - (lambda_V + lambda_M) * S_H
    dE_H = (lambda_V + lambda_M) * S_H - (sigma_H + mu_H) * E_H
    dI_H = sigma_H * E_H - (gamma + mu_H) * I_H
    dR_H = gamma * I_H - mu_H * R_H
    dC_H = sigma_H * E_H # cumulative I_H
    
    lambda_H = b * q * I_H / N_H
    
    dS_V = (k * N_H * (1 - a * cos(2*pi*t / T)) - S_V) * mu_V - lambda_H * S_V
    dM_V = lambda_H * S_V - (sigma_M + mu_V) * M_V
    
    ## if t > the average latent period in mosquito host (tau_V), use M_V(t - tau_V) by lagvalue() function
    ## otherwise, use 0
    if(t >= tau_V+1){
      dE_V = sigma_M * M_V - sigma_M * lagvalue(t - (tau_V-0.5), 7) * exp(-mu_V * tau_V) - mu_V * E_V
      dI_V = sigma_M * lagvalue(t - (tau_V-0.5), 7) * exp(-mu_V * tau_V) - mu_V * I_V
    }
    else{
      dE_V = sigma_M * M_V - sigma_M * 0 * exp(-mu_V * tau_V) - mu_V * E_V
      dI_V = sigma_M * 0 * exp(-mu_V * tau_V) - mu_V * I_V
    }
    
    dvec = lambda_V * S_H  # Those who are infected via vector transmission
    dmec = lambda_M * S_H  # Those who are infected via mechanical transmission
    
    return(list(c(dS_H, dE_H, dI_H, dR_H, dC_H, dS_V, dM_V, dI_V, dE_V, dvec, dmec)))
  })
}

rootfun <- function (t,x,params) {
  dstate <- unlist(dengue_simulation(t,x,params))["I_H"] 
  return(dstate - 1e-10)
}

## Find peak
findPeak = function(x, getValue = FALSE){
  peakIdx = -1
  max = -Inf
  for(i in 1:(length(x)-2)){
    diffVal = x[i+1] - x[i]
    
    if((x[i+1] - x[i]) > 0){
      if((x[i+2] - x[i+1]) <= 0 & x[i+1] > max){
        peakIdx = i+1
        max = x[i+1]
      }
    }
  }
  return(peakIdx)
}


##### Daily-based simulation #####

k_list = as.character(c(1.83, 2.21, 2.75, 10))
pM_list = as.character(c(0.02075, 0.03, 0.083, 0.12))
vars = expand.grid(b = c(0.3, 0.5), sigma_M = c(24/2, 24/1), n_M = c(1, 3))
vars = vars[c(-1, -2, -7, -8), ]
infect_prop_noMec = array(NA, dim = c(length(k_list), nrow(vars)), 
                          dimnames = list(k_list, 
                                          apply(vars, 1, function(x) paste0("b = ", x[1], ", 1/sigma_M = 1/", x[2], ", n_M = ", x[3]))))

gt = list()

for (i in 1:nrow(vars)) {
  
  mech_prop = peak_diff = infect_diff = array(NA, dim = c(length(k_list), length(pM_list)), 
                                              dimnames = list(k_list, pM_list))
  
  for (k in k_list) {
    for(pM in pM_list){
      
  
      ## Parameters
      paras = list(mu_H = 1/(70*365), N_H = 10000, mu_V = 1/14,
                   sigma_H = 1/5, sigma_V = 1/10, sigma_M = 24/1 , # (24 ~ 1)
                   gamma = 1/6, k = 1.83, b = 0.3, p = 0.38, q = 0.38, a = 0.5,
                   T = 1*365, n_M = 1 , # (1 ~ 3)
                   p_M = 0.12,   # (0.083 ~ 0.12)
                   mechanical = TRUE
                   )
      
      ## For different variable combinations
      paras$k = as.numeric(k)
      paras$p_M = as.numeric(pM)
      paras$b = vars$b[i]
      paras$sigma_M = vars$sigma_M[i]
      paras$n_M = vars$n_M[i]
      
      ## other parameters
      paras$b_M = paras$n_M * paras$sigma_M
      paras$N_V = paras$k * paras$N_H
      paras$tau_V = 1 / paras$sigma_V
      
      ## Initial condition
      xstart= c(S_H = 0, E_H = 0, I_H = 10, R_H = 0, C_H = 0, 
                S_V = 0, M_V = 0, I_V = 100, E_V = 0, vec = 0, mec = 0)
      xstart["S_H"] = paras$N_H - xstart["E_H"] - xstart["I_H"]
      xstart["S_V"] = paras$N_H * paras$k - xstart["M_V"] - xstart["I_V"]
      paras$start = xstart
      
      ## Start solving ODE
      result = dede(xstart, 1:10^5, func = dengue_simulation, paras, method = "lsodar", rootfun = rootfun)
      result = data.frame(result)
      result$mec_ratio = result$mec / (result$mec + result$vec)
      
      ## The condition of no mechanical transmission
      paras$mechanical = FALSE
      result_noMec = dede(xstart, 1:10^5, func = dengue_simulation, paras, method = "lsodar", rootfun = rootfun)
      result_noMec = data.frame(result_noMec)
      result_noMec$mec_ratio = result_noMec$mec / (result_noMec$mec + result_noMec$vec)
      paras$mechanical = TRUE
      
      ## Show some information
      # message(paste0("Final mechanical transmission proportion: ", tail(result$mec_ratio, 1)))
      # message(paste0(round((1-tail(result_noMec$S_H, 1) / paras$start["S_H"]) * 100, 2), 
      #                "% susceptible were infected when simulation finished if no mechanical trasmission."))
      # message(paste0("The difference between peaks of mechanical transmission pandemic and non-mechanical is ",
      #                (findPeak(result$I_H) - findPeak(result_noMec$I_H)), " days."))
      # message(paste0("The difference of culmative infected number between mechanical and non-mechanical is ", 
      #                (tail(result$C_H, 1) - tail(result_noMec$C_H, 1)), "."))
      
      mech_prop[k, pM] = tail(result$mec, 1) / tail(result$C_H, 1)
      peak_diff[k, pM] = findPeak(result$I_H) - findPeak(result_noMec$I_H)
      if(abs(peak_diff[k, pM]) > 180 | peak_diff[k, pM] == 0) peak_diff[k, pM] = NA
      infect_diff[k, pM] = (tail(result$C_H, 1) - tail(result_noMec$C_H, 1)) / paras$N_H
      
      
    }
    
    infect_prop_noMec[k, i] = (1-tail(result_noMec$S_H, 1) / paras$start["S_H"]) * 100
    remove(result, result_noMec)
  }
  
  gt[[i*3-2]] = pheatmap(mech_prop, main = paste0("b = ", paras$b, ", 1/simga_M = ", 24/paras$sigma_M, "/24", ", n_M = ", paras$n_M), angle_col = 0, 
                         color = colorRampPalette(brewer.pal(n = 5, name = "OrRd"))(100),
                         display_numbers = round(mech_prop, 3), fontsize_number = 14, fontsize = 14, 
                         #breaks = breakList, 
                         cluster_rows = FALSE, cluster_cols = FALSE, border_color = "white")[[4]]
  gt[[i*3-1]] = pheatmap(peak_diff, main = paste0("b = ", paras$b, ", 1/simga_M = ", 24/paras$sigma_M, "/24", ", n_M = ", paras$n_M), angle_col = 0, 
                         color = rev(colorRampPalette(brewer.pal(n = 5, name = "OrRd"))(100)),
                         display_numbers = peak_diff, fontsize_number = 14, fontsize = 14, 
                         #breaks = breakList, 
                         cluster_rows = FALSE, cluster_cols = FALSE, border_color = "white")[[4]]
  gt[[i*3]] = pheatmap(infect_diff, main = paste0("b = ", paras$b, ", 1/simga_M = ", 24/paras$sigma_M, "/24", ", n_M = ", paras$n_M), angle_col = 0, 
                       color = colorRampPalette(brewer.pal(n = 5, name = "OrRd"))(100),
                       display_numbers = round(infect_diff, 3), fontsize_number = 14, fontsize = 14, 
                       #breaks = breakList, 
                       cluster_rows = FALSE, cluster_cols = FALSE, border_color = "white")[[4]]
  
  remove(peak_diff, mech_prop, infect_diff, paras)
}
remove(k, pM, i)
png("combinations.png", width = 3000, height = 4000, res = 200)
grid.arrange(arrangeGrob(grobs= gt, ncol=3), 
             top = paste("Mechanical transmission proportion", "Peak day delay compare to non-mech", "Total infected proportion difference", 
                         sep = "                                                              "),
             bottom = paste("p_M", "p_M", "p_M", 
                            sep = "                                                                                                        "))
dev.off()


##### Plots #####

## plot every variables
layout(matrix(c(1, 6, 10, 2, 7, 11, 3, 8, 12, 4, 9, 0, 5, 0, 0), nrow = 3, ncol = 5))
for (type in colnames(result)[2:ncol(result)]) {
  plot(result[,type], type = "l", main = type)
}
remove(type)

layout(matrix(c(1), nrow = 1, ncol = 1))

## suspect left - k plot
# df = data.frame(k = seq(5, 1, -0.01))
# df$susLeft = rep(NA, nrow(df))
# for(i in 1:nrow(df)){
#   paras$k = df$k[i]
#   result_noMec = dede(xstart, 1:10^5, func = dengue_simulation, paras, method = "lsodar", rootfun = rootfun)
#   result_noMec = data.frame(result_noMec)
#   df$susLeft[i] = result_noMec$S_H[nrow(result_noMec)]
# }
# plot(df$k, df$susLeft, type = "l", ylab = "Susceptible left", xlab = "k", main = "Biting rate = 0.3")


## daily incidence plot (Figure 4A)

df = data.frame(Time = seq_len(nrow(result)-1),
                `Standard transmission` = result$vec[2:nrow(result)] - result$vec[1:(nrow(result)-1)],
                `Mecanical transmission` = result$mec[2:nrow(result)] - result$mec[1:(nrow(result)-1)])
df$`Total incidence` = df$Standard.transmission + df$Mecanical.transmission
df$`Total incidence w/o mechanical transmission` = result_noMec$vec[2:nrow(result)] - result_noMec$vec[1:(nrow(result)-1)] + 
                                                      result_noMec$mec[2:nrow(result)] - result_noMec$mec[1:(nrow(result)-1)]
                                                     
melt(df, id.vars = "Time", value.name = "Incidence", variable.name = "Type") %>%
  mutate(lineType = sub("Standard.transmission|Mecanical.transmission|Total incidence", "a", Type)) %>% 
  mutate(lineType = sub("a w/o mechanical transmission", "b", lineType)) %>% 
  # mutate(lineType = as.numeric(lineType)) %>% 
  ggplot(aes(x = Time, y = Incidence, group = Type)) +
  geom_line(aes(col = Type, linetype = lineType), size = 2) +
  guides(linetype = FALSE) + 
  theme_light(base_size = 26, base_line_size = 0) +
  theme(legend.position = c(0.8, 0.8),
        legend.key.width = unit(4,"line"),
        legend.key.height = unit(4, "line")) +
  scale_color_manual(labels = c("Standard\ntransmission", "Mechanical\ntransmission", "Total\nincidence", 
                                "Total incidence \nw/o mechanical\n transmission"),
                     values = c("blue", "red", "black", "grey")) +
  labs(x = "Time (days)", y = "Daily incidence", color = "")


## figure 4B
df = data.frame(pM = seq(0, 15, by = 0.05),
                nM1 = rep(NA, length(seq(0, 15, by = 0.05))),
                nM2 = rep(NA, length(seq(0, 15, by = 0.05))),
                nM3 = rep(NA, length(seq(0, 15, by = 0.05))))
for (i in 1:nrow(df)) {
  paras$p_M = df$pM[i] / 100
  for (nM in 1:3) {
    paras$b_M = nM * paras$sigma_M
    df[i, (nM+1)] = dede(xstart, 1:10^5, func = dengue_simulation, paras, method = "lsodar", rootfun = rootfun) %>%
                    as.data.frame() %>%
                    {sum(.[, "mec"]) / (sum(.[, "mec"]) + sum(.[, "vec"]))}
  }
}
remove(nM, i)

df[, c(2:4)] = df[, c(2:4)] * 100
melt(df, id.vars = "pM", value.name = "ratio", variable.name = "nM") %>%
  ggplot(aes(x = pM, y = ratio, group = nM)) +
  geom_line(aes(col = nM), size = 2) +
  theme_light(base_size = 26, base_line_size = 0) +
  theme(legend.position = c(0.2, 0.8),
        legend.key.width = unit(4,"line"),
        legend.key.height = unit(4, "line")) +
  scale_color_manual(labels = c(bquote(~italic(n)[italic(M)]~"=1"), bquote(~italic(n)[italic(M)]~"=2"), bquote(~italic(n)[italic(M)]~"=3")),
                     values = c("red", "blue", "green")) +
  labs(x = "Percentage of cases transmitted mechanically (%)",
       y = "Case transmitted mechanically (%)", color = "")