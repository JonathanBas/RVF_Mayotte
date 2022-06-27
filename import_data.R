rm(list=ls())

library("dplyr")
library("data.table")
library("ggplot2")

dat0 <- as.data.frame(fread(file="./data/res_mayotte_rvf_2_CSV.csv",
                            select = c("identifiant", "IPSOS_COMMUNE", "P15_AGE", "p15_INF14", "RVF IgG"),
                            header=T, sep=";", stringsAsFactors=F,na.strings = "", dec=",", encoding = "Latin-1"))
dat0 <- as.data.frame(dat0 %>% rename("id"="identifiant", "comm"="IPSOS_COMMUNE", "age"="P15_AGE", "date_prev"="p15_INF14", "rvf_igg"="RVF IgG")
                      %>% filter(!is.na(rvf_igg)))

head(dat0)

dat0$date_prev <- as.Date(dat0$date_prev, origin="1960-01-01")
summary(dat0$date_prev)

dat0$inf <- "anc"
dat0$inf[dat0$rvf_igg <= 3] <- "neg"

dat0$week <- format(dat0$date_prev, "%Y-%V")
dat0$week <- gsub("2018-01", "2019-01", dat0$week)
dat0$comm <- tolower(dat0$comm)
dat0$type <- "sero"

cas_dec <- as.data.frame(fread(file="./data/cas_fvr_mayotte.csv",
                               select = c("Num", "date de la declaration", "semaine", "date de naissance", "COMMUNE", "Date debut des signes"),
                               header=T, sep=";", stringsAsFactors=F,na.strings = "", dec=",", encoding = "Latin-1"))

cas_dec <- as.data.frame(cas_dec %>% rename("id"="Num", "date_dec"="date de la declaration", "raw_wk"="semaine", "DOB"="date de naissance", "comm"="COMMUNE", "DSO"="Date debut des signes"))

# Dates: if available, we take the date of symptoms onset (DSO). If not, we take the week in database
cas_dec$date_dec <- gsub("02/2018", "02/2019", cas_dec$date_dec)
cas_dec$date_dec <- gsub("01/2018", "01/2019", cas_dec$date_dec)
cas_dec$date_dec <- as.Date(cas_dec$date_dec, format="%d/%m/%Y")
cas_dec$DSO <- as.Date(cas_dec$DSO, format="%d/%m/%Y", na.)
cas_dec$week <- format(cas_dec$DSO, "%Y-%V")
cas_dec$week <- ifelse(is.na(cas_dec$week),
                          ifelse(cas_dec$raw_wk >= 40,
                                 paste0("2018-", cas_dec$raw_wk),
                                 ifelse(cas_dec$raw_wk >= 10,
                                        paste0("2019-", cas_dec$raw_wk),
                                        paste0("2019-0", cas_dec$raw_wk))),
                          cas_dec$week)
cas_dec$comm <- tolower(cas_dec$comm)
cas_dec$type <- "cas"
# Age of cases:
cas_dec$DOB <- as.Date(cas_dec$DOB, format="%d/%m/%Y")
cas_dec$age <- 2019 - as.numeric(format(cas_dec$DOB, "%Y"))
print(paste(sum(!is.na(cas_dec$age) & (cas_dec$age>=15)), "cases included"))

dat1 <- merge(x = dat0,
              y = filter(cas_dec, !is.na(age), age>=15),
              by = c("id", "comm", "week", "type", "age"), all=T)

dat1$aggr_geo1 <- "Mayotte"

dat1$aggr_geo3 <- "Outer communes"
dat1$aggr_geo3[is.na(dat1$comm)] <- NA
dat1$aggr_geo3[dat1$comm %in% c("mamoudzou", "dembeni","ouangani","tsingoni","sada")] <- "Central communes"
