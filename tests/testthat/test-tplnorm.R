
# The goal of this file is to test the computation of the sum:
#
#   x^(-a)exp(-bx) for x from 0 to infinity
#
# which is done by tplinfsum in spatialwarnings


# very small to 10, log-separated
res <- 8
grid <- expand.grid(a = seq(-.5, 10, l = res),
                    b = seq(.Machine$double.eps, 10, l = res))

reltol <- 1e-8
maxit <- 1e7
chunksize <- 1e6

results <- plyr::ldply(seq.int(nrow(grid)), function(gridi) {

  a <- grid[gridi, "a"]
  b <- grid[gridi, "b"]
  exps <- pows <- prod <- 0

  i <- inext <- 1
  while ( i < maxit ) {
    inext <- min(i + chunksize + 1, maxit)
    xs <- seq(i, inext, by = 1)

    exps <- exps + sum( exp(-b * xs ) )
    pows <- pows + sum( xs^(-a)       )
    prod <- prod + sum( exp(-b * xs ) * xs^(-a) )
    i <- inext
  }

  ref <- tplinfsum(a, b, xmin = 1, maxit = maxit, reltol = reltol)

  data.frame(grid[gridi, ],
             imax = inext,
             maxit = maxit,
             expsum = exps,
             powsum = pows,
             totsum = prod,
             ref = ref)
}, .progress = "time")

with(results, plot(log10(totsum), log10( abs(totsum - ref) / totsum )))


ggplot(results, aes(x = a, y = b)) +
  geom_raster(aes(fill = log10(totsum - ref))) +
  scale_fill_distiller(palette = "Spectral")



