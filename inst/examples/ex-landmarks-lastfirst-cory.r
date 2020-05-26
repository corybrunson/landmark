# small circle sample
X <- tdaunif::sample_circle(n = 12)
# random seed index
s <- sample(nrow(X), size = 1)
l <- landmarks_lastfirst_cory(X, seed_index = s)
# plot landmark order at point positions
plot(X, asp = 1, pch = NA)
text(X, labels = order(l))

# add duplicate points
X <- rbind(X, X[sample(nrow(X), 6, replace = TRUE), , drop = FALSE])
# firstlast seed index
f <- chebycenter_firstlast_cory(X)
l <- landmarks_lastfirst_cory(X, seed_index = f)
# plot landmark order at point positions
plot(X, asp = 1, pch = NA)
text(X, labels = order(l))
