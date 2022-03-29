library(devtools)
library(parallel)

##########

bpath <- "./inst/illustrations"
nr_words <- 10000
embedding_size <- 100
maxlen <- 100
order_bsp <- 25

if (!file.exists(file.path(bpath, "data_splitted.RDS"))) {

  library(keras)
  library(tidyverse)
  library(jsonlite)
  library(lubridate)
  library(ggplot2)
  library(ggsci)
  library(tidytext)
  library(tm)
  library(data.table)
  theme_set(theme_bw())

  # read in data
  movies <- read_csv(file.path(bpath, "movies.csv"))

  # get genres
  movies <- movies %>%
    filter(original_language == "en" & status == "Released") %>%
    # filter(nchar(genres)>2) %>%
    mutate(genres = lapply(genres, function(x)
      fromJSON(x)$name)) %>%
    select(
      genres,
      budget,
      overview,
      popularity,
      production_countries,
      release_date,
      revenue,
      runtime,
      vote_average,
      vote_count
    ) %>%
    filter(vote_count > 0) %>%
    mutate(release_date = as.numeric(as.Date("2020-01-01") - as.Date(release_date)))

  genres <-
    movies %>% unnest(c(genres), names_repair = "unique") %>% select(genres)
  table(genres)

  movies <- movies[unlist(sapply(movies$genres, function(x)
    length(x) > 1 | (!("TV Movie" %in% x) & !("Foreign" %in% x)))), ]

  # vote_average
  ggplot(
    movies %>%  unnest(c(genres), names_repair = "unique") %>%
      filter(!genres %in% c(NA, "TV Movie", "Foreign")),
    aes(x = log(revenue), fill = genres)
  ) +
    geom_density(alpha = 0.3) + xlab("log-revenue") +
    theme_bw() +
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      legend.title = element_blank(),
      text = element_text(size = 14),
      legend.position = "bottom"
    ) +
    scale_color_jco() +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE))

  all_genres <- sort(unique(unlist(movies$genres)))
  genres_wide <- t(sapply(movies$genres,
                          function(x)
                            colSums(
                              model.matrix( ~ -1 + genre,
                                            data = data.frame(genre =
                                                                factor(x, levels = all_genres)))
                            )))

  colnames(genres_wide) <- gsub(" ", "_", colnames(genres_wide))

  movies <- cbind(movies %>% select(-genres), genres_wide)


  # init tokenizer
  tokenizer <- text_tokenizer(num_words = nr_words)

  # remove stopwords
  data("stop_words")
  stopwords_regex <- paste(c(stopwords('en'), stop_words$word),
                           collapse = '\\b|\\b')
  stopwords_regex <- paste0('\\b', stopwords_regex, '\\b')
  movies$overview <- tolower(movies$overview)
  movies$overview <-
    stringr::str_replace_all(movies$overview, stopwords_regex, '')
  movies$overview <- gsub('[[:punct:] ]+', ' ', movies$overview)

  saveRDS(movies$overview, file = file.path(bpath, "mov_ov.RDS"))

  tokenizer %>% fit_text_tokenizer(movies$overview)

  # text to sequence
  text_seqs <- texts_to_sequences(tokenizer, movies$overview)

  # pad text sequences
  text_padded <- text_seqs %>%
    pad_sequences(maxlen = maxlen)

  # save words for later
  words <- tibble(word = names(tokenizer$word_index),
                  id = as.integer(unlist(tokenizer$word_index)))

  words <- words %>%
    filter(id <= tokenizer$num_words) %>%
    arrange(id)

  saveRDS(words, file = file.path(bpath, "words.RDS"))
  rm(words)
  gc()

  # text sequences as list of one array
  text_embd <-
    list(texts = array(text_padded, dim = c(NROW(movies), maxlen)))

  # create input list
  data <- append(movies, text_embd)

  rm(movies, text_embd)
  gc()

  repetitions <- 10

  data_list <- list()

  # repeat analysis 10 times
  for (repl in 1:repetitions) {
    set.seed(41 + repl)

    train_ind <-
      sample(1:NROW(data$runtime), round(0.8 * NROW(data$runtime)))
    test_ind <- setdiff(1:NROW(data$runtime), train_ind)

    train <-
      lapply(data, function(x)
        if (length(dim(x)) == 2)
          x[train_ind, ]
        else
          x[train_ind])
    test <-
      lapply(data, function(x)
        if (length(dim(x)) == 2)
          x[test_ind, ]
        else
          x[test_ind])

    data_list[[repl]] <- list(train = train, test = test)
  }

  saveRDS(data_list, file.path(bpath, "data_splitted.RDS"))
} else {
  data_list <- readRDS(file.path(bpath, "data_splitted.RDS"))
}

res <- list()
res_list <- list()
epochs <- 1e3
nrCores <- 2 # length(data_list)

# repeat analysis 10 times
res_list <- mclapply(1:length(data_list), function(repl) {
  library(lubridate)

  # load software
  load_all("../deepregression")
  load_all(".")

  # reinitialize tokenizer
  tokenizer <- text_tokenizer(num_words = nr_words)
  tokenizer %>% fit_text_tokenizer(readRDS(file.path(bpath, "mov_ov.RDS")))

  train <- data_list[[repl]]$train
  test <- data_list[[repl]]$test

  # optimizer
  optimizer <- optimizer_adadelta(learning_rate = 0.1)

  # Models ------------------------------------------------------------------

  # fml <- vote_count ~ s(budget) + s(popularity) + s(runtime) + s(revenue) + genreAction
  train$vote_count <- as.integer(train$vote_count)
  fml_deep <- vote_count ~ s(budget) + s(popularity) + s(runtime) + s(revenue) +
    deep(texts) + genreAction

  embd_mod <- function(x) x %>%
    layer_embedding(input_dim = tokenizer$num_words,
                    output_dim = embedding_size) %>%
    layer_lstm(units = 10, return_sequences = TRUE) %>%
    layer_dropout(rate = 0.3) %>%
    layer_flatten() %>%
    layer_dense(25) %>%
    layer_dropout(rate = 0.1) %>%
    layer_dense(5) %>%
    layer_dropout(rate = 0.1) %>%
    layer_dense(1)

  mod <- deeptrafo(
    fml_deep,
    list_of_deep_models = list(deep = embd_mod),
    data = train,
    optimizer = optimizer,
    order = order_bsp
  )

  mod %>% fit(
    epochs = epochs,
    validation_split = 0.2,
    batch_size = 2^8,
    view_metrics = FALSE,
    verbose = TRUE,
    callbacks =
      list(
        callback_early_stopping(patience = 60),
        callback_reduce_lr_on_plateau(patience = 20)
      )
  )

  coef(mod, which_param = "interacting")
  coef(mod, which_param = "shifting")
  plot(mod, which = "s(budget)")
  plot(mod, which = "s(popularity)")
  plot(mod, which = "s(runtime)")
  - logLik(mod) / length(train$genreAction)

})
