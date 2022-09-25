set.seed(3)
# small circle sample
X <- tdaunif::sample_circle(n = 12L)
# random seed index
(l <- landmarks_maxmin(X, seed_index = "random"))
# plot landmark order at point positions
plot(X, asp = 1, pch = NA)
text(X, labels = order(l$landmark))
# minmax seed index
(l <- landmarks_maxmin(X, seed_index = "minmax"))
# plot landmark order at point positions
plot(X, asp = 1, pch = NA)
text(X, labels = order(l$landmark))
