library(ggplot2)
library(data.table)
library(readstata13)
library(heatmapEco)

destring <- function(x, ignore) as.numeric(gsub(ignore, "", x))

genHOColumns <- function(data) {
  for (i in 2010:2014) {
  istr <- as.character(i)
  ids <- c(paste0("ho", istr), paste0("owners", istr), paste0("pop", istr))
  data[,(ids[1]) := get(ids[2])/get(ids[3])]
  }
  for (i in 2011:2014) {
  istr <- as.character(i)
  data[,(paste0("hodiff", istr)) := log(get(paste0("ho", istr)))
         - log(ho2010)]
  }
}
hoDemean <- function(data1, colname, val) data1[,c(colname), with=FALSE] - val

ageShiftView <- function(data, binvar, pol_yr=2009) {
  tryCatch(data[, head(tax_yr), by=tax_yr],
           error = function(x) stop("Cannot find year variable"))
  # Plot of all distributions (like in first paper)
  print(qplot(x=get(binvar), y=val_share, group=factor(tax_yr),
        data=data[tax_yr != pol_yr], geom="line") +
        geom_line(data=data[tax_yr == pol_yr], colour="red"))
  # Counterfactual distributions (mean of all non-policy yrs)
  val_cf <- lm(val_share ~ factor(get(binvar)), data=data[tax_yr != pol_yr])
  dataSums <- summary(data[[binvar]])
  limits <- unique(data[[binvar]])
  newD <- data.frame(limits)
  colnames(newD) <- binvar
  predicted <- data.frame(xvar=limits,
                          val_hat=predict(val_cf, newdata=newD))
  # Plot of counterfactual vs. policy distributions
  print(qplot(x=xvar, y=val_hat, data=predicted, geom="line") + 
        geom_line(data=data[tax_yr == pol_yr], aes(x=get(binvar), y=val_share),
                  colour="red"))
  predicted
}

# HO RATE BY EXPOSURE BIN STUFF
setwd("~/dropbox/fthb_model/data")
crosswalk <- data.table::data.table(read.dta13("fthb_zipbins_master.dta"))
ho_zip <- data.table::data.table(read.dta13("horate_by_zip.dta"))
setkey(ho_zip, zip)
setkey(crosswalk, zip)
ho_merge <- crosswalk[ho_zip, nomatch=0]
ho_collapseCols <- names(ho_merge) %like% "pop|owners"

ho_agg <- ho_merge[,lapply(.SD, function(x) sum(destring(x, ","))),.SDcols=ho_collapseCols]
genHOColumns(ho_agg)

ho_collapsed <- ho_merge[!is.na(treat_decile),
                         lapply(.SD, function(x) sum(destring(x, ","))),
                         .SDcols=ho_collapseCols, by=.(treat_decile)]
genHOColumns(ho_collapsed)
for (i in 2010:2014) {
  colname <- paste0("ho", as.character(i))
  ho_collapsed[,(colname) := hoDemean(
               ho_collapsed, colname,
               as.numeric(ho_agg[,c(colname), with=FALSE]))]
}
qplot(data = melt(ho_collapsed[, .(treat_decile, ho2010, ho2011,
                                   ho2012, ho2013, ho2014)],
                  c("treat_decile")), x=variable, y=value,
      group=treat_decile, colour=factor(treat_decile), geom="line")

ho_collapsed <- ho_merge[!is.na(incdrop_decile)][!is.na(treat_decile),
                                                lapply(.SD, function(x) sum(destring(x, ","))),
                                                .SDcols=ho_collapseCols,
                                                by=.(treat_decile, incdrop_decile)]
genHOColumns(ho_collapsed)
heatmapEco(hodiff2014 ~C(treat_decile,treat_decile):incdrop_decile,
           data=ho_collapsed, factor.ax=TRUE)

ho_collapsed <- ho_merge[!is.na(polinc_decile)][!is.na(treat_decile),
                                                lapply(.SD, function(x) sum(destring(x, ","))),
                                                .SDcols=ho_collapseCols,
                                                by=.(treat_decile, polinc_decile)]
genHOColumns(ho_collapsed)
heatmapEco(hodiff2014 ~C(treat_decile,treat_decile):polinc_decile,
           data=ho_collapsed, factor.ax=TRUE)


# STATE-LEVEL FHA VARIATION STUFF
hmda <- fread('hmda_type_lar.csv')

hmda$count <- as.numeric(hmda$count)
setkey(hmda, as_of_year, state_name)


hmda[,denom:=sum(count), by=.(as_of_year, state_name)]
hmda[, prop := count/denom*100]
hmda[,logprop:=log(prop)]
hmdapre <- hmda[as_of_year==2008][loan_type_name=="FHA-insured", .(state_name, loan_type_name, logprop)]
hmdapost <- hmda[as_of_year==2009][loan_type_name=="FHA-insured", .(state_name, loan_type_name, logprop)]
hmdapre[,logdiff := hmdapost$logprop - logprop]
setkey(hmdapre, logprop)

ggplot(hmda[loan_type_name=="FHA-insured"][
            grep("^[A-F]", state_name)],
       aes(y=prop, x=as_of_year, group=state_name, color=state_name)) +
       geom_line()
ggplot(hmda[loan_type_name=="FHA-insured"][
            grep("^[G-L]", state_name)],
       aes(y=prop, x=as_of_year, group=state_name, color=state_name))  +
       geom_line()


# FTHB DISTRIBUTION STUFF
cross <- fread("IRS/states.csv")
colnames(cross) <- c("state_name", "state")
setkey(cross, state)

ages <- data.table::data.table(read.dta13("IRS/AGE_YR.dta"))
ages[,val_tot := sum(count), by=tax_yr]
ages[, val_share := count/val_tot]
ageShiftView(ages, "age_primary")


ages_state <- data.table::data.table(read.dta13("IRS/AGE_YR_state.dta"))
ages_state <- ages_state[!which(ages_state$state %in% c("GU", "PR")),]
ages_state[!is.na(count),val_tot := sum(count), by=.(tax_yr, state)]
ages_state[, val_share := count/val_tot]
stateExcess <- data.table(state=c(""), excess=c(NA))
for (STATE in unique(ages_state$state[1:51])) {
  final <- ageShiftView(ages_state[state==STATE][!is.na(count)],
                        "age_primary")
  polDist <- (ages_state[state==STATE][age_primary <= 30]
              [tax_yr == 2009, val_share])
  stateExcess <- rbind(stateExcess,
                       list(state=STATE, excess=(
                       sum(polDist) - sum(final[final$xvar <= 30,c("val_hat")])
                       )))
}
setkey(stateExcess, state)
stateNew <- stateExcess[cross]
setkey(hmdapre, state_name)
setkey(stateNew, state_name)
state_final <- hmdapre[stateNew]
qplot(x=logprop, y=excess, data=state_final) +stat_smooth(method="lm")

ages_pct <- data.table::data.table(read.dta13("IRS/AGE_YR_TPCT.dta"))
ages_pct[!is.na(count),val_tot := sum(count), by=.(tax_yr, treat_pct)]
ages_pct[, val_share := count/val_tot]
pctExcess <- data.frame(trear_pct=c(), excess=c(), stringsAsFactors = FALSE)
for (i in 1:100) {
  final <- ageShiftView(ages_pct[treat_pct==i][!is.na(count)],
                        "age_primary")
  polDist <- (ages_pct[treat_pct==i][age_primary <= 30]
              [tax_yr == 2009, val_share])
  pctExcess <- rbind(pctExcess,
                     list(treat_pct=i, excess=(
                     sum(polDist) - sum(final[final$xvar <= 30,c("val_hat")])
                       )), stringsAsFactors = FALSE)
}

incs_pct <- data.table::data.table(read.dta13("IRS/INC_YR.dta"))
incs_pct <- incs_pct[inc_bin < 2.5e5][inc_bin > 0]
incs_pct[!is.na(count),val_tot := sum(count), by=tax_yr]
incs_pct[,val_share := count/val_tot]
incs_pct[,inc_bin := inc_bin/1000]
incs_pct2 <- incs_pct[!which(incs_pct$tax_yr %in% c(2004, 2005, 2006))]
incs_pct3 <- incs_pct[!which(incs_pct$tax_yr %in% c(2004, 2006, 2009))]
ageShiftView(incs_pct, "inc_bin")
ageShiftView(incs_pct2, "inc_bin")
ageShiftView(incs_pct3, "inc_bin", pol_yr=2005)

incs_state <- data.table::data.table(read.dta13("IRS/INC_YR_STATE.dta"))
incs_state <- incs_state[inc_bin5 < 2.5e5][inc_bin5 > 0]
incs_state <- incs_state[!which(
    incs_state$state %in% c("AE", "AP", "GU", "PR", "VI"))][!is.na(count)]
incs_state[,val_tot := sum(count), by=tax_yr]
incs_state[,val_share := count/val_tot]
incs_state[,inc_bin5 := inc_bin5/1000]
stateExcess <- data.frame(state=c(), excess=c(), stringsAsFactors = FALSE)
for (STATE in unique(incs_state$state)) {
  tryCatch({
  final <- ageShiftView(incs_state[state==STATE][!is.na(count)],
                        "inc_bin5")
  polDist <- (incs_state[state==STATE][inc_bin5 <= 6e4]
              [tax_yr == 2009, val_share])
  print(sum(polDist) - sum(final[final$xvar <= 6e4,c("val_hat")]))
  stateExcess <- rbind(stateExcess,
                       list(state=STATE, excess=(
                         sum(polDist) - sum(final[final$xvar <= 6e4,c("val_hat")])
                       )), stringsAsFactors = FALSE)
  }, error = function(x) {
    stateExcess <- rbind(stateExcess,
                         list(state=STATE, excess=NA),
                         stringsAsFactors = FALSE)
  })
}
