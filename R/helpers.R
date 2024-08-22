
#' Balance a data.frame
#'
#' Balance with `n` individuals in each group, taken from the first column
#'
#' @param x data.frame with first column as grouping column
#' @param n integer number of individuals in each group. Smallest sample size by default.
#'
#' @details Note that [sample]() is done with `replace=TRUE`
#'
#' @return balanced data.frame
#'
#' @examples
#' # untouched dataset is not balanced
#' table(pig$sp)
#' # balance helps
#' table(balance(pig)$sp)
#' @export
balance <- function(x, n=min(table(x[[1]]))) {
  ids <- unlist(lapply(split(seq_along(x[[1]]), x[[1]]), sample, n, replace=TRUE))
  x[ids, ]
}

#' Randomize a data.frame
#'
#' Randomize a data.frame by shuffling the first column
#'
#' @param x data.frame with first column as grouping column
#'
#' @return shuffled data.frame
#'
#' @examples
#' # untouched dataset
#' head(pig$sp, 10)
#' # balance helps
#' head(randomize(pig)$sp, 10)
#' @export
randomize <- function(x){
  x[[1]] <- sample(x[[1]])
  x
}

#' Retain some variables from a data.frame
#'
#' Retain only the `n` columns, the first let apart.
#'
#' @param x data.frame with first column as grouping column
#' @param n integer number of columns to retain. All by default.
#'
#' @return data.frame
#'
#' @examples
#' # untouched dataset with all its columns
#' ncol(pig)
#' # after retaining only 12 (plus the 1st as the grouping)
#' ncol(retain(pig, 12))
#' @export
retain <- function(x, n=ncol(x)-1){
  x[, 1:(n+1)]
}

## core functions -----
#' Linear discricriminant analysis
#'
#' A thin wrapper around `MASS::lda` that returns a (named) confusion matrix
#'
#' @param x data.frame
#'
#' @return confusion matrix (of class `table`)
#'
#' @examples
#' pig %>% lda1()
#' @export
lda1 <- function(x){
  actual <- x[[1]]
  pred <- MASS::lda(x[, -1], actual, CV=TRUE)$class
  table(actual, pred, dnn=c("actual", "predicted"))
}

# # confusion matrix helpers -----
#' Confusion matrix metrics
#'
#' Only accuracy (global and by class) so far but more may come.
#'
#' @param x table typically a confusion matrix
#'
#' @return a named vector
#'
#' @rdname metrics
#' @examples
#' x <- iris %>% lda1
#' x %>% acc()
#' x %>% acc_classes()
#' x %>% acc_all()
#
#' @export
acc <- function(x){
  c("acc"=sum(diag(x))/sum(x))
}

#' @rdname metrics
#' @export
acc_classes <- function(x){
  diag(x)/rowSums(x)
}

#' @rdname metrics
#' @export
acc_all <- function(x){
  c(acc(x), acc_classes(x))
}

# # cross-validation ----
# #    - original: the dataset passed, untouched
# #    - random:   original but shuffled on first (ie grouping) column
# #    - balanced: original but balanced on first column
# #    - balanced_random: balanced but randomized on first column
#
# partition_names <- c("original", "random", "balanced", "balanced_random")
# partition <- function(x){
#   list(original=x,
#        random=x %>% randomize(),
#        balanced=x %>% balance(),
#        balanced_random=x %>% balance() %>% randomize()
#   )
# }
#
# # hooks ----
# hook_bypass <- function(x) x
# hook_pca <- function(x, ...) {
#   scores <- stats::prcomp(x[, -1], ...)$x %>% tibble::as_tibble()
#   bind_cols(dplyr::select(x, 1), scores)
# }
#
#
# iris %>% partition()
# iris %>% partition() %>% map(hook_bypass)
# iris %>% partition() %>% map(hook_pca)
# iris %>% partition() %>%
#   map(lda1) %>%
#   map_dfr(~.x %>% acc_all) %>%
#   dplyr::mutate(partition=partition_names, .before=1)
#
# cv1 <- function(x, retain){
#   x %>%
#     retain(retain) %>%
#     partition() %>%
#     purrr::map_dfr(~.x %>% lda1 %>% acc_all) %>%
#     dplyr::mutate(retain=retain, partition=partition_names, .before=1)
# }
# cv12 <- function(x, retain){
#   x %>%
#     retain(retain) %>%
#     partition() %>%
#     map_dfr(~.x %>% lda1 %>% acc_all) -> x
#   x$partition <- x$partition_names
#   x
# }
#
# microbenchmark::microbenchmark(cv1(pig2, 3), cv12(pig2, 3))
#
# cv1(pig2, 3)
#
#
#
#
# cv2 <- function(x, max_retain=30){
#   map_dfr(seq_len(max_retain), ~x %>% cv1(retain=.x))
# }
#
# cv <- function(x, max_retain=30, iter=5){
#   map_dfr(seq_len(iter), ~x %>% cv2(max_retain) %>% mutate(iter=.x, .before=1))
# }
#
#
# pig2 <- pig$mat %>% as_tibble() %>% mutate(sp=pig$gp, .before=1) %>% retain(30)
#
# pig2 %>% cv1(2)
#
# z <- map_dfr(seq_len(20), ~pig2 %>% cv2(30) %>% mutate(iter=.x, .before=1))
#
# partition_cols <- c("original"="#0392cf", balanced="#7bc043",
#                     "balanced_random"="#f37736", "random"="#ee4035")
# z %>%
#   pivot_longer(acc:last_col(), names_to = "metric") %>%
#   mutate(metric=ifelse(metric=="acc", "global", metric)) %>%
#   mutate(partition=factor(partition, levels=names(partition_cols), ordered=TRUE)) %>%
#   mutate(metric=factor(metric, c("global", "DP", "WB"), ordered = TRUE)) %>%
#   ggplot() +
#   aes(x=retain, y=value, col=partition) +
#   geom_point(size=0.01, alpha=0.1) +
#   # geom_density2d() +
#   geom_smooth(linewidth=0.5) +
#   scale_color_manual(values=partition_cols) +
#   facet_grid(~metric) +
#   ggtitle("mevolCVP rewrite") +
#   xlab("nb of components") +
#   ylab("accuracies") +
#   theme_minimal() +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_x_continuous(expand = c(0, 0)) +
#   theme(panel.grid.minor = element_blank(), legend.position = "bottom",
#         panel.background = element_rect(colour="grey50"))
# ggsave("mevolCVP2.pdf", width=9, height=5)
#
# z %>%
#   pivot_longer(acc:last_col(), names_to = "metric") %>%
#   mutate(metric=ifelse(metric=="acc", "global", metric)) %>%
#   mutate(partition=factor(partition, levels=names(partition_cols), ordered=TRUE)) %>%
#   mutate(metric=factor(metric, c("global", "DP", "WB"), ordered = TRUE)) %>%
#   ggplot() +
#   aes(x=retain, y=value, col=partition) +
#   geom_point(size=0.01, alpha=0.1) +
#   # geom_density2d() +
#   geom_smooth(formula="y~x", method="loess", linewidth=0.5, level=0.95) +
#   scale_color_manual(values=partition_cols) +
#   facet_grid(~metric) +
#   ggtitle("mevolCVP rewrite") +
#   xlab("nb of components") +
#   ylab("accuracies") +
#   theme_minimal() +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_x_continuous(expand = c(0, 0)) +
#   theme(panel.grid.minor = element_blank(), legend.position = "bottom",
#         panel.background = element_rect(colour="grey50"))
#
#
# table(g)
#
#
#
# expand_grid(iter=1:5, retain=1:10, partition=partition_names)
#
#
#
#
#
#
#
#
#
#
#
