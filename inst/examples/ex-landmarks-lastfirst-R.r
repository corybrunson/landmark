set.seed(4)
# small circle sample
X <- tdaunif::sample_circle(n = 6)
# random seed index
l <- landmarks_lastfirst_R(X, seed_index = "random")
# plot landmark order at point positions
plot(X, asp = 1, pch = NA)
text(X, labels = order(l))
# add duplicate points
X <- rbind(X, X[sample(nrow(X), 6, replace = TRUE), , drop = FALSE])
# firstlast seed index
l <- landmarks_lastfirst_R(X, seed_index = "firstlast")
# plot landmark order at point positions
plot(X, asp = 1, pch = NA)
text(X + .1 * cbind(floor((1:12 - 1) / 6) - .5, 0), labels = order(l))
