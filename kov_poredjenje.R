# Prethodno definisana funkcija za generisanje geometrijskog Braunovog kretanja.
GB_v <- function(S_0, r, sigma, T, n, N) {
  delta_t <- T / n
  c_1 <- (r - (sigma ^ 2 / 2)) * delta_t
  c_2 <- sqrt(delta_t) * sigma
  norm <- matrix(rnorm(n * N, c_1, c_2), nrow = N, ncol = n)
  log_prirastaj <- cbind(log(S_0), norm)
  log_putanja <- t(apply(log_prirastaj, 1, cumsum))
  return(exp(log_putanja))
}

AzijskaOpcijaKovPoredjenje <- function(S_0, K, T, r, sigma, n, N) {
  # vektor cijena opcije dobijenih naivnom Monte-Karlo metodom
  prosjecne_cijene <- rep(0:N)
  
  # vektor cijena dobijen na osnovu prve kontrolne promjenljive
  aritmeticke_cijene <- rep(0:N)
  # vektor cijena evropske kol opcije
  evropske_cijene <- rep(0:N)
  
  # vektor cijena opcija u slucaju da se prosek cijena akcija racuna preko 
  # geometrijske sredine
  geometrijske_cijene <- rep(0:N)

  faktor_diskontovanja <- exp(-r * T)
  # generisemo N simulacija i u svakoj od njih racunamo odgovarajucu cenu opcije
  for (i in 1:N) {
    putanja <- GB_v(S_0, r, sigma, T, n, 1)
    prosjecne_cijene[i] <- mean(putanja[-1])
    evropske_cijene[i] <-
      faktor_diskontovanja * max(putanja[length(putanja)] - K, 0)
    geometrijske_cijene[i] <-
      faktor_diskontovanja * max((prod(putanja[-1])) ^ (1 / n) - K, 0)
    aritmeticke_cijene[i] <- faktor_diskontovanja * max(mean(putanja[-1]) - K, 0)
  }
  #racunamo 3 korelacije (trazimo onu kombinaciju koja ima najvecu korelaciju)
  x <- cor(prosjecne_cijene, aritmeticke_cijene)
  y <- cor(evropske_cijene, aritmeticke_cijene)
  z <- cor(geometrijske_cijene, aritmeticke_cijene)
  cat("Korelacija izmedju cena opcija dobijenih iz Naive Monte Carlo metode
      (aritmeticka sredina) i prosecnih cena akcija je: ",x,"\n\n")
  cat(
    "Korelacija izmedju cena opcija dobijenih iz Naive Monte Carlo metode
     (aritmeticka sredina) i cena evropske kol opcije je: ", y,"\n\n")
  cat(
    "Korelacija izmedju cena opcija dobijenih iz Naive Monte Carlo metode
    (aritmeticka sredina i cena azijske opcije (geometrijska sredina): ",z)
}
set.seed(3333)
AzijskaOpcijaKovPoredjenje(100, 80, 1, 0.05, 0.4, 12, 500)

