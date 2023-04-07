rm(list=ls())

load("models_all_chosen_eps_10_reps.Rdata")
library(tidyverse)
library(spls)
data("yeast")
source("cross_validation.R")

pred <- yeast$x
resp <- yeast$y

Y    <- t(sweep(resp, 2, apply(resp, 2, mean, na.rm = TRUE), FUN =  "-"))
X    <- t(scale(pred, center=TRUE, scale=TRUE))
n    <- ncol(Y); p <- nrow(X); q <- nrow(Y); 


# Result for Table 5

model_char_hist <- NULL
for (rep.num in seq_along(CV.models)){
  cvmodel.cur.rep <- CV.models[[rep.num]]
  for (j in 1:length(cvmodel.cur.rep)){
    mf   <- cvmodel.cur.rep[[j]]$ModelFreq
    eps  <- cvmodel.cur.rep[[j]]$eps
    df   <- bind_rows(lapply(seq_along(mf), function(i) { nm   <- names(mf); 
    vars <- as.numeric(unlist(str_split(nm[i], " ")));
    var_names <- str_remove(colnames(pred)[vars], "_YPD")
    TF_chars <- paste0(var_names, collapse = " ") 
    freq <- unname(mf[i]);
    dat  <- tibble("model_len" = length(vars), "model_char" = nm[i], 
                   "models" = list(vars),
                   "TF_chars" = TF_chars,
                   "freq" = freq, 
                   "epsilon" = eps,  "rep" = rep.num);
    dat 
    }))
    model_char_hist <- bind_rows(model_char_hist, df)
  }
}

model.char.prob <- model_char_hist %>% group_by(epsilon, 
                                                 model_char, 
                                                 TF_chars, models, model_len) %>% 
  dplyr::summarise(Total_freq = sum(freq), .groups = "drop") %>%
  mutate(probs = Total_freq/(10*5000)) %>% arrange(epsilon, desc(probs))


model.char.prob_gt10 <- model.char.prob %>% filter(epsilon==0.09) %>%
  filter(probs >= 0.1) %>%
  rowwise() %>%
  mutate(CV.score = list(rowMeans(replicate(1000, cross_validation(Y = t(resp),
                                                                   X = t(pred), 
                                                                   model = unlist(models),
                                                                   nfold = 10)))), 
         MSE = unlist(CV.score)[1], MAD = unlist(CV.score)[2]) %>%
  select(TF_chars, model_len, probs, MSE, MAD)

print(model.char.prob_gt10)

# Result for Table 6

model_hist <- NULL
for (rep.num in seq_along(CV.models)){
  cvmodel.cur.rep <- CV.models[[rep.num]]
  for (j in 1:length(cvmodel.cur.rep)){
    mf   <- cvmodel.cur.rep[[j]]$ModelFreq
    eps  <- cvmodel.cur.rep[[j]]$eps
    print(eps)
    df   <- bind_rows(lapply(seq_along(mf), function(i) { nm   <- names(mf); 
    vars <- as.numeric(unlist(str_split(nm[i], " ")));
    freq <- unname(mf[i]);
    dat  <- data.frame("vars" = vars, "freq" = rep(freq, length(vars)), 
                       "epsilon" = rep(eps, length(vars)),  "rep" = rep(rep.num, length(vars)) );
    dat 
    }))
    model_hist <- bind_rows(model_hist, df)
  }
}

data.genes <- data.frame("vars" = 1:ncol(pred), "TF" = colnames(pred))
model.prob <- model_hist %>% group_by(vars) %>% dplyr::summarise(Total_freq = sum(freq), .groups = "drop") %>%
  mutate(probs = Total_freq/(10*5000*4)) %>% arrange(desc(probs)) %>%
  inner_join(data.genes)

print(model.prob[1:13,c(4,3)])
