library(deeptrafo)
set.seed(-1)
options(width = 60)

# Params ------------------------------------------------------------------

bpath <- "./inst/illustrations"
nr_words <- 10000
embedding_size <- 100
maxlen <- 100
order_bsp <- 25
repetitions <- 2

# Loading the data --------------------------------------------------------

data_list <- readRDS(file.path(bpath, "data_splitted.RDS"))[[1]]
train <- data_list[[1]]
test <- data_list[[2]]

tokenizer <- text_tokenizer(num_words = nr_words)
tokenizer %>% fit_text_tokenizer(readRDS(file.path(bpath, "mov_ov.RDS")))

# Formula interface -------------------------------------------------------

fm <- vote_count | genreAction ~ s(budget) + popularity

# Setting up DCTMs --------------------------------------------------------

opt <- optimizer_adam(learning_rate = 0.1, decay = 4e-4)
m_fm <- cotramNN(formula = fm, data = train, optimizer = opt)
m_fm

# Fitting DCTMs -----------------------------------------------------------

hist <- fit(m_fm, epochs = 200, validation_split = 0.1, batch_size = 64)
unlist(coef(m_fm, which = "shifting"))

# Working with neural networks --------------------------------------------

embd_mod <- function(x) x %>%
  layer_embedding(input_dim = tokenizer$num_words,
                  output_dim = embedding_size) %>%
  layer_lstm(units = 10, return_sequences = TRUE) %>%
  layer_dropout(rate = 0.3) %>% layer_flatten() %>%
  layer_dense(25) %>% layer_dropout(rate = 0.1) %>%
  layer_dense(5) %>% layer_dropout(rate = 0.1) %>%
  layer_dense(1)

fm_deep <- update(fm, . ~ . + deep(texts))
m_deep <- deeptrafo(fm_deep, data = train,
                    list_of_deep_models = list(deep = embd_mod))

dhist <- fit(m_deep, epochs = 5, validation_split = 0.1, batch_size = 32)
# TODO: Show embedding

# Ensembling DCTMs --------------------------------------------------------

ens_deep <- ensemble(m_deep, n_ensemble = 5, epochs = 10, batch_size = 64)

# Cross-validating DCTMs --------------------------------------------------

cv_deep <- cv(m_deep, epochs = 10, cv_folds = 5, batch_size = 64)

# Methods overview --------------------------------------------------------

### "deeptrafo"

coef(mod, which = "interacting")
coef(mod, which = "shifting")

# In-sample
predict(mod, type = "pdf")
fitted(mod)
logLik(mod)

# <FIXME> Out-of-sample </FIXME>
logLik(mod, newdata = test)

### "dtEnsemble"
coef(dctm_ens, which = "interacting")
coef(dctm_ens, which = "shifting")
