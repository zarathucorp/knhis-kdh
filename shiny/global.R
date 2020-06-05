##### global.R: 전처리

## Load packages
library(data.table)
library(jstable)
library(lubridate)

## Read data
#setwd("/home/jinseob2kim/Desktop")
a <- fread("data_draft.csv")[AGE_GROUP >= 5]

a$AGE <- a$AGE_GROUP                               ## Use AGE_GROUP as continuous variable.
a$PERSON_ID <- as.integer(as.vector(a$PERSON_ID))  ## If class(a$PERSON_ID) != "integer"


## Assign period < 30 as non-User
a$group_ppi2[a$period.ppi <= 30] <- 0
a$group_ppi3[a$period.ppi <= 30] <- 0
a$group_h2ra2[a$period.h2ra <= 30] <- 0
a$group_h2ra3[a$period.h2ra <= 30] <- 0

## Gastric cancer type: Cardia vs nonCardia
a[, type_gc := ifelse(type_gc %in% c("C160", "C1600", "C1601", "C1609"), "cardia", ifelse(type_gc == "", "nongc", "noncardia"))]

## New group : "Exact 0 day" use VS "> 0 day"
a[, group_ppi2_new := factor(as.integer(period.ppi > 0))]
a[, group_h2ra2_new := factor(as.integer(period.h2ra > 0))]

## Day from PPI(H2RA) to Gastric cancer
a[, ppi_to_gc := day_gc - as.integer(as_date(start.ppi) - as_date(date))]
a[, h2ra_to_gc := day_gc - as.integer(as_date(start.h2ra) - as_date(date))]


## Make variable lists to use
varlist <- list(
  Base = c("PERSON_ID", "SEX", "AGE", "AGE_GROUP", grep("Pre_", names(a), value = T)),
  Group = c("group_ppi2", "group_ppi2_new", "group_ppi3", "period.ppi4", "group_h2ra2", "group_h2ra2_new", "group_h2ra3", "period.h2ra4"),
  Event = c("event_gc", "type_gc", "day_gc", "period.ppi", "period.h2ra", "ppi_to_gc")
)

## Factor variables
factor_vars <- setdiff(unlist(varlist), c("AGE", "day_gc","period.ppi", "period.h2ra", "ppi_to_gc"))
a[, (factor_vars) := lapply(.SD, factor), .SDcols = factor_vars]


## Data to analysis : Exclude GC event within 1 year
out <- a[(event_gc == 0) | (day_gc > 365), .SD, .SDcols = unlist(varlist)]
#out.ppi <- a[(event_gc == 0) |  ((day_gc > 365) & ((start.ppi == "") | (ppi_to_gc > 180))), .SD, .SDcols = c(unlist(varlist), "start.ppi")]
#out.h2ra <- a[(event_gc == 0)| ((day_gc > 365) & ((start.h2ra == "") |(h2ra_to_gc > 180))), .SD, .SDcols = unlist(varlist)]
#out.ppih2ra <- out.ppi[group_ppi2 != group_h2ra2][]


## Make label information: use jstable::mk.lev
out.label <- mk.lev(out)
out.label[variable == "SEX", val_label := c("Male", "Female")]
out.label[variable == "AGE_GROUP", ":="(var_label = "AGE group", val_label = c("20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                                                                                "60-64", "65-69", "70-74", "75-79", "80-84", "> 85"))]

for (v in grep("Pre_", names(a), value = T)){
  out.label[variable == v, val_label := c("No", "Yes")]
}

out.label[variable == "group_ppi2", ":="(var_label = "PPI User", val_label = c("Non-PPI(< weekly use or < 30 days)" ,"PPI"))]
out.label[variable == "group_ppi3", ":="(var_label = "PPI User", val_label = c("Non-PPI(< weekly use or < 30 days)" ,"Weekly to <daily", "Daily"))]
out.label[variable == "group_ppi2_new", ":="(var_label = "PPI User", val_label = c("Non-PPI(exact 0 day)" ,"PPI(> 0 days)"))]
out.label[variable == "period.ppi4", ":="(var_label = "PPI period", val_label = c("1 month - 1 yrs" ,"1-2 yrs", "2-3 yrs", "> 3 yrs"))]
out.label[variable == "group_h2ra2", ":="(var_label = "H2RA User", val_label = c("Non-H2RA(< weekly use or < 30 days)" ,"H2RA"))]
out.label[variable == "group_h2ra3", ":="(var_label = "H2RA User", val_label = c("Non-H2RA(< weekly use or < 30 days)" ,"Weekly to <daily", "Daily"))]
out.label[variable == "group_h2ra2_new", ":="(var_label = "H2RA User", val_label = c("Non-H2RA(exact 0 days)" ,"H2RA(> 0 days)"))]
out.label[variable == "period.h2ra4", ":="(var_label = "H2RA period", val_label = c("1 month - 1 yrs" ,"1-2 yrs", "2-3 yrs", "> 3 yrs"))]
out.label[variable == "event_gc", ":="(var_label = "Gastric cancer", val_label = c("No", "Yes"))]
out.label[variable == "day_gc", ":="(var_label = "FU days")]
