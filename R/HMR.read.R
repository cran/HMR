.HMR.read<-function(filename,dec,sep)
{
  ## Input
  ## -----
  ## filename: En tekststreng indeholdende filnavnet. Det forudsaettes, at datamappen i forvejen er sat med 'setwd'.
  ##           Datafilen skal vaere organiseret i fem kolonner og opfylde flg. betingelser:
  ##
  ##           1. Indholdet af 1. raekke ignoreres - kan evt. bruges til variabelnavne eller anden hjaelpetekst.
  ##           2. Kolonnerne indeholder (1.-5.):
  ##                - Serienavn (tekststreng)
  ##                - Kammerets volumen
  ##                - Kammerets tvaersnitsareal
  ##                - Proevetagningstidspunkterne
  ##                - Maalte koncentrationer
  ##           3. Serienavnet identificerer den enkelte dataserie, saa forskellige dataserier SKAL have forskellige
  ##              serienavne. Volumen og areal skal indeholde positive tal, der er konstante indenfor den enkelte
  ##              dataserie. Indholdet af kolonne 1-3 er altsaa konstant indenfor serienavn. Proevetagningstidspunkterne
  ##              skal vaere (strengt) stigende tal, og det foerste tidspunkt skal vaere nul. De maalte koncentrationer skal
  ##              vaere positive tal, og manglende vaerdier er ikke tilladte.
  ##           4. Der stilles ingen krav til fysiske enheder.
  ##
  ## dec     : Decimaltegn.
  ## sep     : Soejleseparator.

  ## Den basale indlaesning vha. 'read.table'
  ## ---------------------------------------
  Rdata<-try(read.table(file=filename,header=FALSE,skip=1,blank.lines.skip=TRUE,
  dec=dec,sep=sep,col.names=c('serie','V','A','tid','konc'),
  colClasses=c('character','numeric','numeric','numeric','numeric')),silent=TRUE)

  ## Yderligere kontrol af data  
  ## --------------------------
  ##
  ## ANTAGELSE   : 'serie' definerer dataserierne
  ##
  ## Kontrollerer:
  ##
  ##   1. 'V' og 'A' er konstante indenfor 'serie'
  ##   2. Mindst 3 obs. af 'V', 'A', 'tid' og 'konc' indenfor 'serie'
  ##   3. Ingen 'NA' i 'V', 'A', 'tid' eller 'konc'
  ##   4. Ingen 'Inf' eller '-Inf' i 'V', 'A', 'tid' eller 'konc'
  ##   5. Ingen negative vaerdier i 'V', 'A' eller 'konc'
  ##   6. Vaerdierne af 'tid' skal vaere positive og stigende
  FATAL<-FALSE
  if (inherits(Rdata,'try-error'))
  {
    # Elementaere datafejl
    FATAL<-TRUE; HMRdata<-NA; nserie<-NA
  } else
  {
    # Ingen elementaere datafejl
    serie<-Rdata$serie; userie<-unique(serie); nserie<-length(userie)
    V<-Rdata$V
    A<-Rdata$A
    tid<-Rdata$tid
    konc<-Rdata$konc
    nfejlserie<-0
    HMRdata<-vector(mode='list',length=nserie)
    for (i in 1:nserie)
    {
      # Lokale variable
      itid<-tid[serie==userie[i]]
      ikonc<-konc[serie==userie[i]]
      iV<-V[serie==userie[i]]
      iA<-A[serie==userie[i]]
      serieOK<-TRUE
      # Mere end eet 'V' eller 'A'
      if ((length(unique(iV))>1)|(length(unique(iA))>1)) {serieOK<-FALSE} else
      {
        # Mindre end 3 obs.
        if ((length(itid)<3)|(length(ikonc)<3)|(length(iV)<3)|(length(iA)<3)) {serieOK<-FALSE} else
        {
          # 'NA' i numeriske variable
          if ((sum(is.na(itid))>0)|(sum(is.na(ikonc))>0)|(sum(is.na(iV))>0)|(sum(is.na(iA))>0)) {serieOK<-FALSE} else
          {
            # 'Inf' eller '-Inf' i numeriske variable
            if ((max(abs(itid))==Inf)|(max(abs(ikonc))==Inf)|(max(abs(iV))==Inf)|(max(abs(iA))==Inf)) {serieOK<-FALSE} else
            {
              # Ikke-positive vaerdier af 'V', 'A', 'konc' eller 'tid'
              if ((min(iV)<=0)|(min(iA)<=0)|(min(ikonc)<=0)|(min(itid)<0)) {serieOK<-FALSE} else
              {
                # 'tid' er ikke-stigende
                if (min(itid[-1]-itid[-length(itid)])<=0) {serieOK<-FALSE} else
                {
                  # Data 'userie[i]' er OK!
                  HMRdata[[i]]<-list(serie=userie[i],V=unique(iV),A=unique(iA),tid=itid,konc=ikonc,status=1) # status=0: fejl, status=1: OK
                }
              }
            }
          }
        }
      }
      # Hvis 'userie[i]' indeholder fejl
      if (!serieOK) {nfejlserie<-nfejlserie+1; HMRdata[[i]]<-list(serie=userie[i],status=0)}
    }
    # Hvis fejl i alle serier
    if (nfejlserie==nserie) {FATAL<-TRUE; HMRdata<-NA}
  }

  ## Output
  list(nserier=nserie,HMRdata=HMRdata,FATAL=FATAL)
}
