#' @title Funkcija za generisanje geometrijskog Braunovog kretanja koja je 
#' prilagođena metodi antitetičkog uzorkovanja za redukciju disperzije 
#' (Antithetic Variates).
#'
#' @description Izbjegava se korišćenje ugnježđenih for petlji,
#' funkcija je vektorizovana. U ovoj verziji funkcije se odvojeno generišu 
#' putanje za parove antitetičkihslučajnih veličina. 
#' Prilagođena je funkciji AzijskaOpcija_AV.
#' 
#' Generisanje uzorka iz normalne raspodjele se dijeli na dva dijela.
#' Ideja je da se generiše uzorak tako da je prva polovina uzorka negativno 
#' korelisana sa drugom polovinom i da su obje polovine uzorka iz normalne 
#' N(c_1, c_2) raspodjele.
#' 
#' @param S_0 numeric. Početna cijena akcije
#' @param r numeric. Bezrizična kamatna stopa
#' @param sigma numeric. Volatilnost akcije
#' @param n numeric. Broj vremenskih intervala kojima diskretizujemo proces,
#' odnosno broj koraka
#' @param N numeric. Broj trajektorija geometrijskog Braunovog kretanja 
#' @param idx numeric. Uzima vrijednosti 0 ili 1. Ako je vrijednost 1, to je 
#' indikator da funkcija treba da vrati prve elemente parova antitetičkog 
#' uzorka, odnosno prvi element svakog para. 
#' Ako je vrijednost 2, to je indikator da funkcija treba da vrati druge 
#' elemente parova, odnosno drugi element svakog para.
#' @return Matrica trajektorija
#'
#' @examples
#' putanja <- AV_GBM_v(100, 0.05, 0.3, 1, 100, 1)
AV_GBM_v <- function(S_0, r, sigma, T, n, N, idx) {
  # podjela vremenskog intervala
  delta_t <- T / n
  
  # izračunavamo unaprijed poznate konstante
  c_1 <- (r - (sigma ^ 2 / 2)) * delta_t
  c_2 <- sqrt(delta_t) * sigma
  
  # u matricu Z1 sa N/2 vrsta smještamo slučajne veličine iz normalne raspodjele
  Z1 <- matrix(rnorm(n * ceiling(N / 2), c_1, c_2),
               nrow = ceiling(N / 2),
               ncol = n)
  # odgovarajući antitetički par je definisan sa Z2=2*c_1-Z1
  # Z2 takođe ima normalnu N(c_1, c_2) raspodjelu
  Z2 <- 2 * c_1 - Z1
  # ovako definisane Z1 i Z2 su negativno korelisane
  
  # u matrice av_log_prirastaj_i, i=1,2 dodata je kolona početnih vrijednosti 
  # cijena akcije
  av_log_prirastaj_1 <- cbind(log(S_0), Z1)
  av_log_prirastaj_2 <- cbind(log(S_0), Z2)
  
  # log_putanja_i, i=1,2, sadrži trajektorije geometrijskog Braunovog kretanja
  # dobijene preko prvih (drugih) članova niza antitetičkog para (Z1, Z2)
  log_putanja_1 <- t(apply(av_log_prirastaj_1, 1, cumsum))
  log_putanja_2 <- t(apply(av_log_prirastaj_2, 1, cumsum))
  
  # biramo koje putanje treba da se vrate u zavisnosti od izbora vrijednosti
  # parametra idx - upotreba će biti jasnija u sledećoj funkciji
  if (idx == 1)
    return(exp(log_putanja_1))
  else if (idx == 2)
    return(exp(log_putanja_2))
}   

#' @title Monte-Karlo metod za računanje cijene aritmetičke azijske opcije
#' koristeći metodu antitetičkog uzorkovanja za redukciju disperzije.
#' 
#' @description Za vrednovanje opcije koristi se princip vrednovanja 
#' neutralan od rizika.
#' Cijena opcije je diskontovana očekivana vrijednost dobiti opcije pri 
#' bezrizičnoj kamatnoj stopi.
#' Monte-Karlo ocjena cijene je prosjek cijena opcije dobijenih u simulacijama.
#' 
#' @param K numeric. Ugovorena cijena akcije na koju se opcija odnosi
#' @param S_0 numeric. Početna cijena akcije
#' @param r numeric. Bezrizična kamatna stopa
#' @param sigma numeric. Volatilnost akcije
#' @param n numeric. Broj vremenskih intervala kojima diskretizujemo proces,
#' odnosno broj koraka
#' @param N numeric. Broj Monte-Karlo simulacija (ujedno i broj 
#' generisanih trajektorija)
#' @param beta numeric. Nivo povjerenja za konstrukciju intervala povjerenja
#'
#' @examples
#' AzijskaOpcija_AV(70, 100, 0.05, 0.3, 1, 100, 5000, 0.95)

AzijskaOpcija_AV <- function(K, S_0, r, sigma, T, n, N, beta) {
  cijena_opcije <- 0
  
  # cijene_i, i=1,2, su vektori u koje smještamo cijene opcije dobijene 
  # martingalskim pristupom
  # u metodi antitetičkog uzorkovanja se računa cijena opcije dobijena na osnovu
  # prosjeka prvih elemenata antitetičkih parova, a potom se računa cijena na
  # osnovu drugih elemenata antitetičkih parova
  # ovo se ponavlja za svaku simulaciju, odnosno trajektoriju
  cijene1 <- rep(0, ceiling(N / 2))
  cijene2 <- rep(0, ceiling(N / 2))
  
  # konačni vektor cijena se dobija uprosječavanjem cijena dobijenih po parovima
  # ključno je da iako je unutar svakog para dozvoljena korelacija, parovi su
  # među sobom nezavisni, pa ako posmatramo niz dobijen na osnovu prosjeka 
  # unutar svakog para, onda dobijamo niz nezavisnih slučajnih veličina i možemo 
  # konsutrisati intervale povjerenja na uobičajen način
  cijene <- rep(0, 2 * ceiling(N / 2))
  
  # faktor diskontovanja je konstantna vrijednost
  faktor_diskontovanja <- exp(-r * T)
  
  for (i in 1:N) {
    # jedna trajektorija geometrijskog Braunovog kretanja dobijena AV metodom
    # putanja_i, i=1,2, je trajektorija dobijena na osnovu i_tih, i=1,2
    # elemenata antitetičih parova
    putanja1 <- AV_GBM_v(S_0, r, sigma, T, n, 1, 1)
    putanja2 <- AV_GBM_v(S_0, r, sigma, T, n, 1, 2)
    
    # cijena akcije pri valuaciji neutralnoj od rizika
    cijene1[i] <-
      faktor_diskontovanja * max(0, mean(putanja1[-1]) - K)
    cijene2[i] <-
      faktor_diskontovanja * max(0, mean(putanja2[-1]) - K)
    # konačni vektor niz cijena na osnovu kog se računa MK ocjena cijene opcije
    cijene[i] <- (cijene1[i] + cijene2[i]) / 2
  }
  # standardna devijacija niza dobijenih cijena opcije
  std <- sd(cijene)
  
  # Monte-Karlo ocjena cijene opcije
  cijena_opcije <- mean(cijene)
  
  # primjenjujemo studentov test na konačni niz cijena da bismo dobili 
  studentov_test <- t.test(cijene, conf.level = beta)
  # interval povjerenja
  interval_povjerenja <- studentov_test$conf.int
  # računamo grešku Monte-Karlo ocjene
  greska <-
    100 * (1 / (as.numeric(studentov_test$estimate))) *
    ((studentov_test$conf.int[2] - studentov_test$conf.int[1]) / 2)
  # štampamo vrijednost ocjene, standardnu devijaciju, interval povjerenja i 
  # samu grešku
  cat("Tačkasta ocena za cenu azijske opcije je: ", cijena_opcije, "\n")
  cat("Standardna devijacija je: ", std, "\n")
  cat(beta * 100, "% interval povjerenja je: ", interval_povjerenja, "\n")
  cat("Greška ocjene izražena u procentima je: ", greska)
}

# najpreciznije bi bilo uzimati za N paran broj, jer se uzorak u kodu 
# razdvaja na dva dijela koji se odvojeno generisu
# razlika zanemarljiva u slučaju da je N neparno
set.seed(777)
AzijskaOpcija_AV(50, 50, 0.05, 0.6, 1, 12, 5000, 0.95)
