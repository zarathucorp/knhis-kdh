## Load global.R
source("global.R")

## Load packages
library(ggplot2) 
library(DT)                   ## datatable
library(survival)             ## Survival analysis
library(jskm)                 ## Survival curve: my package
library(shinycustomloader)    ## Loading image in shiny
library(forestplot)           ## Forestplot in subgroup analysis


## Matching, IPTW weighting
library(MatchIt);library(survey)


## Own function: kaplan-meier survival curve
Makekaplan <- function(var.day = "day_gc", var.event = "event_gc", var.group, data = out, data.label = out.label, yr.event = 10, cut.yr = F){
  
  form <- as.formula(paste("Surv(", var.day, ",", var.event, ") ~ ", var.group, sep = ""))
  if ("survey.design" %in% class(data)){
    if (cut.yr == T){
      data$variables[[var.event]] <- ifelse(data[[var.day]] > 365 * yr.event & data[[var.event]] == "1", 0,  as.numeric(as.vector(data[[var.event]])))
      data$variables[[var.day]]  <- ifelse(data[[var.day]] > 365 * yr.event, 365 * yr.event, data[[var.day]])
    }
    data$variables[[var.event]] <- as.numeric(as.vector(data$variables[[var.event]]))
    data$variables[[var.day]] <- data$variables[[var.day]]/365
    res.kap <- survey::svykm(form, design = data)  
    p <- svyjskm(res.kap, xlabs = "Years", ylab = "Cumulative Incidence", cumhaz = T, xlims = c(0, yr.event), ylims = c(0, 0.05),
                 pval.coord = c(2.5, 0.01), legendposition = c(0.3, 0.8),
                 ystrataname = data.label[variable == var.event, var_label][1], ystratalabs = data.label[variable == var.group][level %in% levels(data$variables[[var.group]]), val_label],
                 surv.scale = "percent", mark = F, pval = T, table = T, design = data, pval.testname = F)
  } else{
    if (cut.yr == T){
      data[[var.event]] <- ifelse(data[[var.day]] > 365 * yr.event & data[[var.event]] == "1", 0,  as.numeric(as.vector(data[[var.event]])))
      data[[var.day]] <- ifelse(data[[var.day]] > 365 * yr.event, 365 * yr.event, data[[var.day]])
    }
    
    data[[var.event]] <- as.numeric(as.vector(data[[var.event]]))
    data[[var.day]] <- data[[var.day]]/365
    
    res.kap <- survfit(form, data = data)
    res.kap$call$formula <- form
    p <- jskm(res.kap, xlabs = "Years", ylab = "Cumulative Incidence", cumhaz = T, xlims = c(0, yr.event), ylims = c(0, 0.05), 
              ystrataname = data.label[variable == var.event, var_label][1], ystratalabs = data.label[variable == var.group][level %in% levels(data[[var.group]]), val_label],
              pval.coord = c(2.5, 0.01), legendposition = c(0.3, 0.8), timeby = 1,
              surv.scale = "percent", mark = F, pval = T, table = T, data = data)
  }
  
  
  return(p)
  
  
  
}




## UI
ui <- navbarPage("PPI NHIS",
                 tabPanel("Table 1",
                          sidebarLayout(
                            sidebarPanel(
                              radioButtons("data_tb1", "Study ", choices = c("PPI", "H2RA", "PPI vs H2RA"), selected = "PPI", inline = T),
                              radioButtons("group_tb1", "Group by", c("non-Use vs Use", "Dose", "Exact 0 days vs at least 1 days"), "non-Use vs Use", inline = T),
                              selectInput("varmat_ct", "Variables to match", choices = varlist$Base[-c(1, 3)], selected = varlist$Base[-c(1, 3)], multiple = T),
                              radioButtons("ratio_mat", "Matching ratio", choices = c("1:1", "1:2"), selected = "1:1", inline = T),
                              sliderInput("caliper", "Caliper(0: no caliper)", min = 0, max = 2, val =0, step = 0.1),
                              selectInput("nonnormal", "Variables summarized with Median[IQR]", choices = c("day_gc", "period.ppi", "period.h2ra", "ppi_to_gc"), selected = c("day_gc", "period.ppi", "period.h2ra", "ppi_to_gc"), multiple = T),
                              h4("Change exclusion criteria"),
                              sliderInput("cr_hp","제균 치료 후 XX 일 이전 위암 발생 제외", min = 30, max = 500, value = 365, step =5),
                              sliderInput("cr_ppi","PPI(H2RA) 복용 XX 일 이전  위암 발생 제외", min = 30, max = 500, value = 180, step =5),
                              h4("19세 이하 제외"),
                              h4("Smoking:COPD, 담배코드, Alcohol: 알콜중독, 알콜간질환 코드"),
                              h4(textOutput("excl"))
                            ),
                            mainPanel(
                              tabsetPanel(type = "pills",
                                          tabPanel("Original", withLoader(DTOutput("table1"), type="html", loader="loader6")),
                                          tabPanel("PS matching", withLoader(DTOutput("table1_ps"), type="html", loader="loader6")),
                                          tabPanel("IPTW", withLoader(DTOutput("table1_IPTW"), type="html", loader="loader6"))
                              )
                              
                            )
                          )
                          
                 ),
                 tabPanel("Cox model",
                          sidebarLayout(
                            sidebarPanel(
                              radioButtons("data_cox", "Study ", choices = c("PPI", "H2RA", "PPI vs H2RA"), selected = "PPI", inline = T),
                              radioButtons("sub_cox", "Subgroup", choices = c("All", "1 month - 1 yrs", "1-2 yrs", "2-3 yrs", " > 3yrs", "1 month - 3yrs", "> 1 yrs", "> 2 yrs", "0-2 yrs"), selected = "All", inline = T),
                              selectInput("dep_cox", "Outcome", choices = "Gastric cancer", selected = "Gastric cancer"),
                              uiOutput("gcox")
                            ),
                            mainPanel(
                              tabsetPanel(type = "pills",
                                          tabPanel("Original", withLoader(DTOutput("tablecox"), type="html", loader="loader6")),
                                          tabPanel("PS matching", withLoader(DTOutput("tablecox_ps"), type="html", loader="loader6")),
                                          tabPanel("IPTW", withLoader(DTOutput("tablecox_iptw"), type="html", loader="loader6"))
                              )
                            )
                          )
                          
                 ),
                 
                 tabPanel("Kaplan-meier plot",
                          sidebarLayout(
                            sidebarPanel(
                              radioButtons("data_kap", "Study ", choices = c("PPI", "H2RA", "PPI vs H2RA"), selected = "PPI", inline = T),
                              radioButtons("sub_kap", "Subgroup", choices = c("All", "1 month - 1 yrs", "1-2 yrs", "2-3 yrs", " > 3yrs", "1 month - 3yrs", "> 1 yrs", "> 2 yrs", "0-2 yrs"), selected = "All", inline = T),
                              selectInput("event_kap", "Event", choices = "Gastric cancer", selected = "Gastric cancer"),
                              uiOutput("gkap"),
                              uiOutput("cutconti")
                            ),
                            mainPanel(
                              tabsetPanel(type = "pills",
                                          tabPanel("Original", 
                                                   withLoader(plotOutput("kap"), type="html", loader="loader6"),
                                                   h3("Download options"),
                                                   wellPanel(
                                                     uiOutput("downloadControls_kap"),
                                                     downloadButton("downloadButton_kap", label = "Download the plot")
                                                   )),
                                          tabPanel("PS matching", 
                                                   withLoader(plotOutput("kap_ps"), type="html", loader="loader6"),
                                                   h3("Download options"),
                                                   wellPanel(
                                                     uiOutput("downloadControls_kap_ps"),
                                                     downloadButton("downloadButton_kap_ps", label = "Download the plot")
                                                   )),
                                          tabPanel("IPTW", 
                                                   withLoader(plotOutput("kap_iptw"), type="html", loader="loader6"),
                                                   h3("Download options"),
                                                   wellPanel(
                                                     uiOutput("downloadControls_kap_iptw"),
                                                     downloadButton("downloadButton_kap_iptw", label = "Download the plot")
                                                   ))
                              )
                            )
                          )
                          
                 ),
                 tabPanel("Subgroup analysis",
                          sidebarLayout(
                            sidebarPanel(
                              radioButtons("data_tbsub", "Study ", choices = c("PPI", "H2RA", "PPI vs H2RA"), selected = "PPI", inline = T),
                              radioButtons("sub_tbsub", "Subgroup", choices = c("All", "1 month - 1 yrs", "1-2 yrs", "2-3 yrs", " > 3yrs", "1 month - 3yrs", "> 1 yrs", "> 2 yrs", "0-2 yrs"), selected = "All", inline = T),
                              selectInput("dep_tbsub", "Event", choices = "Gastric cancer", selected = "Gastric cancer"),
                              uiOutput("gsub"),
                              sliderInput("year_tbsub", "Cut month", min = 0 , max = 120, value = c(0, 120)),
                              selectInput("subvar_tbsub", "Subgroup to include", choices = grep("Pre_", names(out), value = T), selected = grep("Pre_", names(out), value = T), multiple = T),
                              selectInput("cov_tbsub", "Additional covariates", c("SEX", "AGE", grep("Pre_", names(out), value = T)), selected = NULL, multiple = T),
                              downloadButton("forest", "Download forest plot")
                            ),
                            mainPanel(
                              radioButtons("type_tbsub", "Analysis", c("Original", "PS matching", "IPTW"), inline = T),
                              withLoader(DTOutput("tablesub"), type="html", loader="loader6")
                            )
                          )
                          
                 )
)


## Server
server <- function(input, output, session) {
  
  ## Dataset 
  data.tb1 <- reactive({
    if (input$data_tb1 == "PPI") {
      a[(event_gc == 0) |  ((day_gc > input$cr_hp) & ((start.ppi == "") | (ppi_to_gc > input$cr_ppi))), .SD, .SDcols = c(unlist(varlist), "start.ppi")]
      
      #out.ppi
    } else if (input$data_tb1 == "H2RA"){
      a[(event_gc == 0)| ((day_gc > input$cr_hp) & ((start.h2ra == "") |(h2ra_to_gc > input$cr_ppi))), .SD, .SDcols = unlist(varlist)]
      #out.h2ra
    } else{
      a[(event_gc == 0) |  ((day_gc > input$cr_hp) & ((start.ppi == "") | (ppi_to_gc > input$cr_ppi))), .SD, .SDcols = c(unlist(varlist), "start.ppi")][group_ppi2 != group_h2ra2][]
      #out.ppih2ra
    }
  })
  
  ## Exclusion criteria
  output$excl <- renderText({
    if (input$data_tb1 == "PPI") {
      ex.hp <- nrow(a) - nrow(a[(event_gc == 0) |  ((day_gc > input$cr_hp))])
      ex.ppi <- nrow(a[(event_gc == 0) |  ((day_gc > input$cr_hp))]) - nrow(a[(event_gc == 0) |  ((day_gc > input$cr_hp) & ((start.ppi == "") | (ppi_to_gc > input$cr_ppi)))])
      
      paste0("제균 치료 후 ", input$cr_hp, "일 이전 위암 발생: ", ex.hp, " 명 제외, ", "PPI 복용 후 ", input$cr_ppi, "일 이전 위암 발생: ", ex.ppi, " 명 제외")
      #out.ppi
    } else if (input$data_tb1 == "H2RA"){
      ex.hp <- nrow(a) - nrow(a[(event_gc == 0) |  ((day_gc > input$cr_hp))])
      ex.h2ra <- nrow(a[(event_gc == 0) |  ((day_gc > input$cr_hp))]) - nrow(a[(event_gc == 0) |  ((day_gc > input$cr_hp) & ((start.h2ra == "") | (h2ra_to_gc > input$cr_ppi)))])
      paste0("제균 치료 후 ", input$cr_hp, "일 이전 위암 발생: ", ex.hp, " 명 제외, ", "H2RA 복용 후 ", input$cr_ppi, "일 이전 위암 발생: " ,ex.h2ra, " 명 제외")
      #out.h2ra
    } else{
      ex.hp <- nrow(a) - nrow(a[(event_gc == 0) |  ((day_gc > input$cr_hp))])
      ex.ppi <- nrow(a[(event_gc == 0) |  ((day_gc > input$cr_hp))]) - nrow(a[(event_gc == 0) |  ((day_gc > input$cr_hp) & ((start.ppi == "") | (ppi_to_gc > input$cr_ppi)))])
      ex.dupli <- nrow(a[(event_gc == 0) |  ((day_gc > input$cr_hp) & ((start.ppi == "") | (ppi_to_gc > input$cr_ppi))), .SD, .SDcols = c(unlist(varlist), "start.ppi")][group_ppi2 == group_h2ra2]) 
      
      paste0("제균 치료 후 ", input$cr_hp, "일 이전 위암 발생: ", ex.hp, " 명 제외", "PPI 복용 후 ", input$cr_ppi, "일 이전 위암 발생: ", ex.ppi, " 명 제외, PPI와 H2RA 모두 복용력 있는 사람: ", ex.dupli, " 명 제외")
      
      #a[(event_gc == 0) |  ((day_gc > input$cr_hp) & ((start.ppi == "") | (ppi_to_gc > input$cr_ppi))), .SD, .SDcols = c(unlist(varlist), "start.ppi")][group_ppi2 != group_h2ra2][]
      #out.ppih2ra
    }
  })
  
  
  ## Table 1 group variable
  group.tb1 <- reactive({
    if (input$data_tb1 %in% c("PPI", "PPI vs H2RA") & input$group_tb1 == "non-Use vs Use"){
      "group_ppi2"
    } else if (input$data_tb1 %in% c("PPI", "PPI vs H2RA") & input$group_tb1 == "Dose"){
      "group_ppi3"
    } else if (input$data_tb1 %in% c("PPI", "PPI vs H2RA") & input$group_tb1 == "Exact 0 days vs at least 1 days"){
      "group_ppi2_new"
    } else if (input$data_tb1 == "H2RA" & input$group_tb1 == "non-Use vs Use"){
      "group_h2ra2"
    } else if (input$data_tb1 == "H2RA" & input$group_tb1 == "Dose"){
      "group_h2ra3"
    } else if (input$data_tb1 == "H2RA" & input$group_tb1 == "Exact 0 days vs at least 1 days"){
      "group_h2ra2_new"
    }
  })
  
  
  ## Propensity score matching, IPTW
  matct <- reactive({
    req(input$group_tb1 %in%  c("non-Use vs Use", "Exact 0 days vs at least 1 days"))
    req(group.tb1())
    zz <- data.tb1()
    set.seed(1234)
    m.form <- paste0(group.tb1(), "~ ", paste(input$varmat_ct, collapse = " + "))
    #col.mat <- c(group.tb1(), input$varmat_ct, "PERSON_ID")
    if (diff(table(zz[[group.tb1()]])) > 0){
      zz$depvar <- as.integer(zz[[group.tb1()]] == levels(zz[[group.tb1()]])[1])
      m.form <- paste0("depvar ~ ", paste(input$varmat_ct, collapse = " + "))
      #col.mat <- c(group.tb1(), input$varmat_ct, "PERSON_ID", "depvar")
    }
    
    ratio <- ifelse(input$ratio_mat == "1:1", 1, 2)
  
    m.out <- matchit(as.formula(m.form), data = zz[, -c("ppi_to_gc", "start.ppi")], ratio = ratio, caliper= input$caliper)
    zz$wt.iptw.gacheon <- ifelse(as.character(m.out$treat) == "1", 1/m.out$distance, 1/(1-m.out$distance))
    zz[, wt.iptw.gacheon := ifelse(wt.iptw.gacheon > 10, 10, wt.iptw.gacheon)]
    out.iptw <- svydesign(id = ~1, weights = ~wt.iptw.gacheon, data = zz )
    #out.ps <- match.data(m.out)
    out.ps <- zz[PERSON_ID %in% match.data(m.out)$PERSON_ID]
    
    
    return(list(ps = out.ps, iptw = out.iptw))
  }) 
  
  
  
  
  ## Table 1
  output$table1 <- renderDT({
    data <- data.tb1()
    label <- out.label[, .SD]
    if (input$data_tb1 == "PPI vs H2RA"){
      label[variable == "group_ppi2", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI"))]
      label[variable == "group_ppi3", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI(Weekly to <daily)", "PPI(Daily)"))]
    }
    
    tb1 <- CreateTableOneJS(c(varlist$Base[-c(1, 3)], varlist$Event), strata = group.tb1(), data = data, Labels = T, labeldata = label, showAllLevels = F, smd=T, nonnormal = input$nonnormal)$table
    datatable(tb1, caption = paste0(input$data_tb1, "_", group.tb1()), rownames = T, extensions= "Buttons",
              options = c(opt.tb1("tb1"),
                          list(columnDefs = list(list(visible=FALSE, targets= which(colnames(tb1) %in% c("test","sig"))))
                          ),
                          list(scrollX = TRUE)
              )) %>% formatStyle("sig", target = 'row' ,backgroundColor = styleEqual("**", 'yellow'))
  })
  
  output$table1_ps <- renderDT({
    data <- matct()$ps
    label <- out.label[, .SD]
    if (input$data_tb1 == "PPI vs H2RA"){
      label[variable == "group_ppi2", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI"))]
      label[variable == "group_ppi3", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI(Weekly to <daily)", "PPI(Daily)"))]
    }
    
    tb1 <- CreateTableOneJS(c(varlist$Base[-c(1,3)], varlist$Event), strata = group.tb1(), data = data, Labels = T, labeldata = label,  showAllLevels = F, smd=T, nonnormal = input$nonnormal)$table
    datatable(tb1, caption = paste0(input$data_tb1, "_", group.tb1()), rownames = T, extensions= "Buttons",
              options = c(opt.tb1("tb1"),
                          list(columnDefs = list(list(visible=FALSE, targets= which(colnames(tb1) %in% c("test","sig"))))
                          ),
                          list(scrollX = TRUE)
              )) %>% formatStyle("sig", target = 'row' ,backgroundColor = styleEqual("**", 'yellow'))
  })
  
  
  
  output$table1_IPTW <- renderDT({
    data <- matct()$iptw
    label <- out.label[, .SD]
    if (input$data_tb1 == "PPI vs H2RA"){
      label[variable == "group_ppi2", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI"))]
      label[variable == "group_ppi3", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI(Weekly to <daily)", "PPI(Daily)"))]
    }
    
    tb1 <- svyCreateTableOneJS(c(varlist$Base[-c(1,3)], varlist$Event), strata = group.tb1(), data = data, Labels = T, labeldata = label,  showAllLevels = F, smd=T, nonnormal = input$nonnormal)$table
    datatable(tb1, caption = paste0(input$data_tb1, "_", group.tb1()), rownames = T, extensions= "Buttons",
              options = c(opt.tb1("tb1"),
                          list(columnDefs = list(list(visible=FALSE, targets= which(colnames(tb1) %in% c("test","sig"))))
                          ),
                          list(scrollX = TRUE)
              )) %>% formatStyle("sig", target = 'row' ,backgroundColor = styleEqual("**", 'yellow'))
  })
  
 
  
  
  ## Kaplan-meier plot
  data.kap <- reactive({
    if (input$data_kap == "PPI") {
      a[(event_gc == 0) |  ((day_gc > input$cr_hp) & ((start.ppi == "") | (ppi_to_gc > input$cr_ppi))), .SD, .SDcols = c(unlist(varlist), "start.ppi")]
      
      #out.ppi
    } else if (input$data_kap == "H2RA"){
      a[(event_gc == 0)| ((day_gc > input$cr_hp) & ((start.h2ra == "") |(h2ra_to_gc > input$cr_ppi))), .SD, .SDcols = unlist(varlist)]
      #out.h2ra
    } else{
      a[(event_gc == 0) |  ((day_gc > input$cr_hp) & ((start.ppi == "") | (ppi_to_gc > input$cr_ppi))), .SD, .SDcols = c(unlist(varlist), "start.ppi")][group_ppi2 != group_h2ra2][]
      #out.ppih2ra
    }
  })
  
  subvar.kap <- reactive({
    if (input$data_kap %in% c("PPI", "PPI vs H2RA")) {
      "period.ppi4"
    } else{
      "period.h2ra4"
    }
  })
  
  output$gkap <- renderUI({
    zz <- c("group_ppi2", "group_ppi2_new", "group_ppi3")
    if (input$data_kap == "H2RA"){
      zz <- c("group_h2ra2", "group_h2ra2_new", "group_h2ra3")
    }
    selectInput("group_kap", "Group by", choices = c(zz, varlist$Base[-1]), selected = zz[1])
  })
  
  
  
  output$cutconti <- renderUI({
    req(input$group_kap)
    if (class(data.kap()[[input$group_kap]]) %in% c("integer", "numeric")){
      data.km <- switch(input$sub_kap,
                        "All" = data.kap(), 
                        "1 month - 1 yrs" = data.kap()[get(subvar.kap()) == 0], 
                        "1-2 yrs" = data.kap()[get(subvar.kap()) == 1], 
                        "2-3 yrs" = data.kap()[get(subvar.kap()) == 2], 
                        " > 3yrs" = data.kap()[get(subvar.kap()) == 3],
                        "1 month - 3yrs" = data.kap()[get(subvar.kap()) %in% 0:2],
                        "> 1 yrs" = data.kap()[get(subvar.kap()) %in% 1:3], 
                        "> 2 yrs" = data.kap()[get(subvar.kap()) %in% 2:3], 
                        "0-2 yrs" = data.kap()[get(subvar.kap()) %in% 0:1]
                        )
      var.event <- 'event_gc'
      var.day <- "day_gc"
      data.km[[var.event]] <- as.numeric(as.vector(data.km[[var.event]]))
      mstat <- maxstat::maxstat.test(as.formula(paste("survival::Surv(",var.day,",", var.event,") ~ ", input$group_kap, sep="")), data= data.km, smethod="LogRank", pmethod="condMC", B=999)
      cut5 <- mstat$cuts[order(-mstat$stats)][1:5]
      
      vec <- data.km[[input$group_kap]][!is.na(data.km[[input$group_kap]])]
      numericInput("cut_conti", "Cut-off", value = cut5[1], min = quantile(vec, 0.05), max = quantile(vec, 0.95))
    }
  })
  
  obj.km <- reactive({
    req(input$data_kap)
    req(input$group_kap)
    req(input$event_kap)
    
    data <- switch(input$sub_kap,
                      "All" = data.kap(), 
                      "1 month - 1 yrs" = data.kap()[get(subvar.kap()) == 0], 
                      "1-2 yrs" = data.kap()[get(subvar.kap()) == 1], 
                      "2-3 yrs" = data.kap()[get(subvar.kap()) == 2], 
                      " > 3yrs" = data.kap()[get(subvar.kap()) == 3],
                      "1 month - 3yrs" = data.kap()[get(subvar.kap()) %in% 0:2],
                      "> 1 yrs" = data.kap()[get(subvar.kap()) %in% 1:3], 
                      "> 2 yrs" = data.kap()[get(subvar.kap()) %in% 2:3], 
                      "0-2 yrs" = data.kap()[get(subvar.kap()) %in% 0:1]
    )
    label <- out.label[, .SD]
    if (input$data_kap == "PPI vs H2RA"){
      label[variable == "group_ppi2", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI"))]
      label[variable == "group_ppi3", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI(Weekly to <daily)", "PPI(Daily)"))]
    }
    
    
    if (class(data[[input$group_kap]]) %in% c("numeric", "integer")){
      req(input$cut_conti)
      data$xcat <- factor(as.integer(data[[input$group_kap]] > input$cut_conti))
      gvar <- "xcat"
      addlabel <- mk.lev(data[, .SD, .SDcols = "xcat"])
      addlabel[, var_label := paste(label[variable == input$group_kap, var_label][1], "group")]
      addlabel[, val_label := paste(label[variable == input$group_kap, var_label][1], paste(c("\u2264", ">"), input$cut_conti, sep=""))]
      label <- rbind(label, addlabel)
      setkey(label, variable)
    } else{
      gvar <- input$group_kap
    }
    
    Makekaplan(var.group = gvar, data = data, data.label = label) 
    
  })
  
  output$kap <- renderPlot({
    print(obj.km())
  })
  
  output$downloadControls_kap <- renderUI({
    fluidRow(
      column(4,
             selectizeInput("kap_file_ext", "File extension (dpi = 300)", 
                            choices = c("jpg","pdf", "tiff", "svg", "emf"), multiple = F, 
                            selected = "jpg"
             )
      ),
      column(4,
             sliderInput("fig_width_kap", "Width (in):",
                         min = 5, max = 20, value = 8
             )
      ),
      column(4,
             sliderInput("fig_height_kap", "Height (in):",
                         min = 5, max = 20, value = 6
             )
      )
    )
  })
  
  output$downloadButton_kap <- downloadHandler(
    filename =  function() {
      paste(input$event_kap, "_", input$data_kap, "_", input$group_kap , "_plot.", input$kap_file_ext ,sep="")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      withProgress(message = 'Download in progress',
                   detail = 'This may take a while...', value = 0, {
                     for (i in 1:15) {
                       incProgress(1/15)
                       Sys.sleep(0.01)
                     }
                     
                     if (input$kap_file_ext == "emf"){
                       devEMF::emf(file, width = input$fig_width_kap, height =input$fig_height_kap, coordDPI = 300, emfPlus = F)
                       plot(obj.km())
                       dev.off()
                       
                     } else{
                       ggsave(file, obj.km(), dpi = 300, units = "in", width = input$fig_width_kap, height =input$fig_height_kap)
                     }
                     
                   })
      
    })
  
  
  obj.kmps <- reactive({
    req(input$data_kap)
    req(input$group_kap)
    req(input$event_kap)
    
    data <- switch(input$sub_kap,
                   "All" = matct()$ps, 
                   "1 month - 1 yrs" = matct()$ps[get(subvar.kap) == 0], 
                   "1-2 yrs" = matct()$ps[get(subvar.kap()) == 1], 
                   "2-3 yrs" = matct()$ps[get(subvar.kap()) == 2], 
                   " > 3yrs" = matct()$ps[get(subvar.kap()) == 3],
                   "1 month - 3yrs" = matct()$ps[get(subvar.kap()) %in% 0:2],
                   "> 1 yrs" = matct()$ps[get(subvar.kap()) %in% 1:3], 
                   "> 2 yrs" = matct()$ps[get(subvar.kap()) %in% 2:3], 
                   "0-2 yrs" = matct()$ps[get(subvar.kap()) %in% 0:1]
    )
    
    label <- out.label[, .SD]
    if (input$data_kap == "PPI vs H2RA"){
      label[variable == "group_ppi2", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI"))]
      label[variable == "group_ppi3", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI(Weekly to <daily)", "PPI(Daily)"))]
    }
    
    if (class(data[[input$group_kap]]) %in% c("numeric", "integer")){
      req(input$cut_conti)
      data$xcat <- factor(as.integer(data[[input$group_kap]] > input$cut_conti))
      gvar <- "xcat"
      addlabel <- mk.lev(data[, .SD, .SDcols = "xcat"])
      addlabel[, var_label := paste(label[variable == input$group_kap, var_label][1], "group")]
      addlabel[, val_label := paste(label[variable == input$group_kap, var_label][1], paste(c("\u2264", ">"), input$cut_conti, sep=""))]
      label <- rbind(label, addlabel)
      setkey(label, variable)
    } else{
      gvar <- input$group_kap
    }
    
    
    Makekaplan(var.group = gvar, data = data, data.label = label) 
    
  })
  
  output$kap_ps <- renderPlot({
    print(obj.kmps())
  })
  
  output$downloadControls_kap_ps <- renderUI({
    fluidRow(
      column(4,
             selectizeInput("kap_ps_file_ext", "File extension (dpi = 300)", 
                            choices = c("jpg","pdf", "tiff", "svg", "emf"), multiple = F, 
                            selected = "jpg"
             )
      ),
      column(4,
             sliderInput("fig_width_kap_ps", "Width (in):",
                         min = 5, max = 20, value = 8
             )
      ),
      column(4,
             sliderInput("fig_height_kap_ps", "Height (in):",
                         min = 5, max = 20, value = 6
             )
      )
    )
  })
  
  output$downloadButton_kap_ps <- downloadHandler(
    filename =  function() {
      paste(input$event_kap, "_", input$data_kap, "_", input$group_kap , "_plot_matching.", input$kap_ps_file_ext ,sep="")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      withProgress(message = 'Download in progress',
                   detail = 'This may take a while...', value = 0, {
                     for (i in 1:15) {
                       incProgress(1/15)
                       Sys.sleep(0.01)
                     }
                     
                     if (input$kap_ps_file_ext == "emf"){
                       devEMF::emf(file, width = input$fig_width_kap_ps, height =input$fig_height_kap_ps, coordDPI = 300, emfPlus = F)
                       plot(obj.kmps())
                       dev.off()
                       
                     } else{
                       ggsave(file, obj.kmps(), dpi = 300, units = "in", width = input$fig_width_kap_ps, height =input$fig_height_kap_ps)
                     }
                     
                   })
      
    })
  
  
  obj.kmiptw <- reactive({
    req(input$data_kap)
    req(input$group_kap)
    req(input$event_kap)
    
    data <- switch(input$sub_kap,
                   "All" = matct()$iptw, 
                   "1 month - 1 yrs" = subset(matct()$iptw, get(subvar.kap()) == 0), 
                   "1-2 yrs" = subset(matct()$iptw, get(subvar.kap()) == 1), 
                   "2-3 yrs" = subset(matct()$iptw, get(subvar.kap()) == 2), 
                   " > 3yrs" = subset(matct()$iptw, get(subvar.kap()) == 3),
                   "1 month - 3yrs" = subset(matct()$iptw, get(subvar.kap()) %in% 0:2),
                   "> 1 yrs" = subset(matct()$iptw, get(subvar.kap()) %in% 1:3), 
                   "> 2 yrs" = subset(matct()$iptw, get(subvar.kap()) %in% 2:3), 
                   "0-2 yrs" = subset(matct()$iptw, get(subvar.kap()) %in% 0:1)
    )
    label <- out.label[, .SD]
    if (input$data_kap == "PPI vs H2RA"){
      label[variable == "group_ppi2", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI"))]
      label[variable == "group_ppi3", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI(Weekly to <daily)", "PPI(Daily)"))]
    }
    
    if (class(data$variables[[input$group_kap]]) %in% c("numeric", "integer")){
      req(input$cut_conti)
      data$variables$xcat <- factor(as.integer(data$variables[[input$group_kap]] > input$cut_conti))
      gvar <- "xcat"
      
      addlabel <- mk.lev(data$variables[, .SD, .SDcols = "xcat"])
      addlabel[, var_label := paste(label[variable == input$group_kap, var_label][1], "group")]
      addlabel[, val_label := paste(label[variable == input$group_kap, var_label][1], paste(c("\u2264", ">"), input$cut_conti, sep=""))]
      
      label <- rbind(label, addlabel)
      setkey(label, variable)
    } else{
      gvar <- input$group_kap
    }
    
    
    Makekaplan(var.group = gvar, data = data, data.label = label) 
    
  })
  
  output$kap_iptw <- renderPlot({
    print(obj.kmiptw())
  })
  
  output$downloadControls_kap_iptw <- renderUI({
    fluidRow(
      column(4,
             selectizeInput("kap_iptw_file_ext", "File extension (dpi = 300)", 
                            choices = c("jpg","pdf", "tiff", "svg", "emf"), multiple = F, 
                            selected = "jpg"
             )
      ),
      column(4,
             sliderInput("fig_width_kap_iptw", "Width (in):",
                         min = 5, max = 20, value = 8
             )
      ),
      column(4,
             sliderInput("fig_height_kap_iptw", "Height (in):",
                         min = 5, max = 20, value = 6
             )
      )
    )
  })
  
  output$downloadButton_kap_iptw <- downloadHandler(
    filename =  function() {
      paste(input$event_kap, "_", input$data_kap, "_", input$group_kap , "_plot_IPTW.", input$kap_iptw_file_ext ,sep="")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      withProgress(message = 'Download in progress',
                   detail = 'This may take a while...', value = 0, {
                     for (i in 1:15) {
                       incProgress(1/15)
                       Sys.sleep(0.01)
                     }
                     
                     if (input$kap_iptw_file_ext == "emf"){
                       devEMF::emf(file, width = input$fig_width_kap_iptw, height =input$fig_height_kap_iptw, coordDPI = 300, emfPlus = F)
                       plot(obj.kmiptw())
                       dev.off()
                       
                     } else{
                       ggsave(file, obj.kmiptw(), dpi = 300, units = "in", width = input$fig_width_kap_iptw, height =input$fig_height_kap_iptw)
                     }
                     
                   })
      
    })
  
  
  ## Cox analysis
  
  data.cox <- reactive({
    if (input$data_cox == "PPI") {
      a[(event_gc == 0) |  ((day_gc > input$cr_hp) & ((start.ppi == "") | (ppi_to_gc > input$cr_ppi))), .SD, .SDcols = c(unlist(varlist), "start.ppi")]
      
      #out.ppi
    } else if (input$data_cox == "H2RA"){
      a[(event_gc == 0)| ((day_gc > input$cr_hp) & ((start.h2ra == "") |(h2ra_to_gc > input$cr_ppi))), .SD, .SDcols = unlist(varlist)]
      #out.h2ra
    } else{
      a[(event_gc == 0) |  ((day_gc > input$cr_hp) & ((start.ppi == "") | (ppi_to_gc > input$cr_ppi))), .SD, .SDcols = c(unlist(varlist), "start.ppi")][group_ppi2 != group_h2ra2][]
      #out.ppih2ra
    }
  })
  
  subvar.cox <- reactive({
    if (input$data_cox %in% c("PPI", "PPI vs H2RA")) {
      "period.ppi4"
    } else{
      "period.h2ra4"
    }
  })
  
  output$gcox <- renderUI({
    zz <- c("group_ppi2", "group_ppi2_new", "group_ppi3")
    if (input$data_cox == "H2RA"){
      zz <- c("group_h2ra2", "group_h2ra2_new", "group_h2ra3")
    }
    selectInput("cov_cox", "Covariates", choices = c(zz, varlist$Base[-1]), selected = zz[1], multiple = T)
  })

  

  

  
  output$tablecox <- renderDT({
    validate(
      need(!is.null(input$cov_cox), "Please select at least 1 independent variable.")
    )
    #data <- data.cox()
    data <- switch(input$sub_cox,
                   "All" = data.cox(), 
                   "1 month - 1 yrs" = data.cox()[get(subvar.cox()) == 0], 
                   "1-2 yrs" = data.cox()[get(subvar.cox()) == 1], 
                   "2-3 yrs" = data.cox()[get(subvar.cox()) == 2], 
                   " > 3yrs" = data.cox()[get(subvar.cox()) == 3],
                   "1 month - 3yrs" = data.cox()[get(subvar.cox()) %in% 0:2],
                   "> 1 yrs" = data.cox()[get(subvar.cox()) %in% 1:3], 
                   "> 2 yrs" = data.cox()[get(subvar.cox()) %in% 2:3], 
                   "0-2 yrs" = data.cox()[get(subvar.cox()) %in% 0:1]
                   
                   
    )
    label <- out.label[, .SD]
    if (input$data_cox == "PPI vs H2RA"){
      label[variable == "group_ppi2", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI"))]
      label[variable == "group_ppi3", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI(Weekly to <daily)", "PPI(Daily)"))]
    }
    
    
    var.event <- "event_gc"
    var.day <- "day_gc"
    data[[var.event]] <- as.numeric(as.vector(data[[var.event]]))
    
    forms.cox <- as.formula(paste("Surv(", var.day,",", var.event,") ~ ", paste(input$cov_cox, collapse="+"), sep=""))
    
    cc <- substitute(survival::coxph(.form, data= data, model = T), list(.form= forms.cox))
    res.cox <- eval(cc)
    tb.cox <- jstable::cox2.display(res.cox, dec = 2)
    tb.cox <- jstable::LabeljsCox(tb.cox, ref = label)
    out.cox <- rbind(tb.cox$table, tb.cox$metric)
    sig <- out.cox[, ncol(out.cox)]
    sig <- gsub("< ", "", sig)
    sig <- ifelse(as.numeric(as.vector(sig)) <= 0.05, "**", NA)
    out.cox <- cbind(out.cox, sig)
    cap.cox <- paste("Cox's proportional hazard model on time ('", label[variable == var.day, var_label][1] , "') to event ('", label[variable == var.event, var_label][1], "')", sep="")
    
    hide <- which(colnames(out.cox) == c("sig"))
    datatable(out.cox, rownames=T, extensions= "Buttons", caption = cap.cox,
              options = c(opt.tbreg(cap.cox),
                          list(columnDefs = list(list(visible=FALSE, targets= hide))
                          )
              )
    )  %>% formatStyle("sig", target = 'row',backgroundColor = styleEqual("**", 'yellow'))
    
    
  })
  
  
  output$tablecox_ps <- renderDT({
    validate(
      need(!is.null(input$cov_cox), "Please select at least 1 independent variable.")
    )
    data <- switch(input$sub_cox,
                   "All" = matct()$ps, 
                   "1 month - 1 yrs" = matct()$ps[get(subvar.cox()) == 0], 
                   "1-2 yrs" = matct()$ps[get(subvar.cox()) == 1], 
                   "2-3 yrs" = matct()$ps[get(subvar.cox()) == 2], 
                   " > 3yrs" = matct()$ps[get(subvar.cox()) == 3],
                   "1 month - 3yrs" = matct()$ps[get(subvar.cox()) %in% 0:2],
                   "> 1 yrs" = matct()$ps[get(subvar.cox()) %in% 1:3], 
                   "> 2 yrs" = matct()$ps[get(subvar.cox()) %in% 2:3], 
                   "0-2 yrs" = matct()$ps[get(subvar.cox()) %in% 0:1]
    )
    label <- out.label[, .SD]
    if (input$data_cox == "PPI vs H2RA"){
      label[variable == "group_ppi2", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI"))]
      label[variable == "group_ppi3", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI(Weekly to <daily)", "PPI(Daily)"))]
    }
    
    
    var.event <- "event_gc"
    var.day <- "day_gc"

    data[[var.event]] <- as.numeric(as.vector(data[[var.event]]))
    
    forms.cox <- as.formula(paste("Surv(", var.day,",", var.event,") ~ ", paste(input$cov_cox, collapse="+"), sep=""))
    
    cc <- substitute(coxph(.form, data= data, model = T), list(.form= forms.cox))
    res.cox <- eval(cc)
    tb.cox <- jstable::cox2.display(res.cox, dec = 2)
    tb.cox <- jstable::LabeljsCox(tb.cox, ref = label)
    out.cox <- rbind(tb.cox$table, tb.cox$metric)
    sig <- out.cox[, ncol(out.cox)]
    sig <- gsub("< ", "", sig)
    sig <- ifelse(as.numeric(as.vector(sig)) <= 0.05, "**", NA)
    out.cox <- cbind(out.cox, sig)
    cap.cox <- paste("Cox's proportional hazard model on time ('", label[variable == var.day, var_label][1] , "') to event ('", label[variable == var.event, var_label][1], "') : PS matching", sep="")
    
    hide <- which(colnames(out.cox) == c("sig"))
    datatable(out.cox, rownames=T, extensions= "Buttons", caption = cap.cox,
              options = c(opt.tbreg(cap.cox),
                          list(columnDefs = list(list(visible=FALSE, targets= hide))
                          )
              )
    )  %>% formatStyle("sig", target = 'row',backgroundColor = styleEqual("**", 'yellow'))
    
    
  })
  
  
  
  output$tablecox_iptw <- renderDT({
    validate(
      need(!is.null(input$cov_cox), "Please select at least 1 independent variable.")
    )
    data <- switch(input$sub_cox,
                   "All" = matct()$iptw, 
                   "1 month - 1 yrs" = subset(matct()$iptw, get(subvar.cox()) == 0), 
                   "1-2 yrs" = subset(matct()$iptw, get(subvar.cox()) == 1), 
                   "2-3 yrs" = subset(matct()$iptw, get(subvar.cox()) == 2), 
                   " > 3yrs" = subset(matct()$iptw, get(subvar.cox()) == 3),
                   "1 month - 3yrs" = subset(matct()$iptw, get(subvar.cox()) %in% 0:2),
                   "> 1 yrs" = subset(matct()$iptw, get(subvar.cox()) %in% 1:3), 
                   "> 2 yrs" = subset(matct()$iptw, get(subvar.cox()) %in% 2:3), 
                   "0-2 yrs" = subset(matct()$iptw, get(subvar.cox()) %in% 0:1)
    )
    label <- out.label[, .SD]
    if (input$data_cox == "PPI vs H2RA"){
      label[variable == "group_ppi2", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI"))]
      label[variable == "group_ppi3", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI(Weekly to <daily)", "PPI(Daily)"))]
    }
    
    
    var.event <- "event_gc"
    var.day <- "day_gc"
    data$variables[[var.event]] <- as.numeric(as.vector(data$variables[[var.event]]))
    
    forms.cox <- as.formula(paste("Surv(", var.day,",", var.event,") ~ ", paste(input$cov_cox, collapse="+"), sep=""))
    
    cc <- substitute(survey::svycoxph(.form, design= data, model = T), list(.form= forms.cox))
    res.cox <- eval(cc)
    tb.cox <- jstable::svycox.display(res.cox, dec = 2)
    tb.cox <- jstable::LabeljsCox(tb.cox, ref = label)
    out.cox <- rbind(tb.cox$table, tb.cox$metric)
    sig <- out.cox[, ncol(out.cox)]
    sig <- gsub("< ", "", sig)
    sig <- ifelse(as.numeric(as.vector(sig)) <= 0.05, "**", NA)
    out.cox <- cbind(out.cox, sig)
    cap.cox <- paste("Cox's proportional hazard model on time ('", label[variable == var.day, var_label][1] , "') to event ('", label[variable == var.event, var_label][1], "') : IPTW", sep="")
    
    hide <- which(colnames(out.cox) == c("sig"))
    datatable(out.cox, rownames=T, extensions= "Buttons", caption = cap.cox,
              options = c(opt.tbreg(cap.cox),
                          list(columnDefs = list(list(visible=FALSE, targets= hide))
                          )
              )
    )  %>% formatStyle("sig", target = 'row',backgroundColor = styleEqual("**", 'yellow'))
    
    
  })
  
  
  ## Subgroup analysis
  data.tbsub <- reactive({
    if (input$data_tbsub == "PPI") {
      a[(event_gc == 0) |  ((day_gc > input$cr_hp) & ((start.ppi == "") | (ppi_to_gc > input$cr_ppi))), .SD, .SDcols = c(unlist(varlist), "start.ppi")]
      
      #out.ppi
    } else if (input$data_tbsub == "H2RA"){
      a[(event_gc == 0)| ((day_gc > input$cr_hp) & ((start.h2ra == "") |(h2ra_to_gc > input$cr_ppi))), .SD, .SDcols = unlist(varlist)]
      #out.h2ra
    } else{
      a[(event_gc == 0) |  ((day_gc > input$cr_hp) & ((start.ppi == "") | (ppi_to_gc > input$cr_ppi))), .SD, .SDcols = c(unlist(varlist), "start.ppi")][group_ppi2 != group_h2ra2][]
      #out.ppih2ra
    }
  })
  
  subvar.tbsub <- reactive({
    if (input$data_tbsub %in% c("PPI", "PPI vs H2RA")) {
      "period.ppi4"
    } else{
      "period.h2ra4"
    }
  })
  
  output$gsub <- renderUI({
    zz <- c("group_ppi2", "group_ppi2_new", "group_ppi3")
    if (input$data_tbsub == "H2RA"){
      zz <- c("group_h2ra2", "group_h2ra2_new", "group_h2ra3")
    }
    radioButtons("group_tbsub", "Variable for HR", choices = zz[1:2], selected = zz[1], inline = T)
  })
  
  group.tbsub <- reactive(
    input$group_tbsub
  )

  
  tbsub <- reactive({
    data <- switch(input$type_tbsub,
                   "Original" = data.tbsub(),
                   "PS matching" = matct()$ps,
                   "IPTW" = matct()$iptw)
    
    data <- switch(input$sub_tbsub,
                   "All" = data, 
                   "1 month - 1 yrs" = subset(data, get(subvar.tbsub()) == 0), 
                   "1-2 yrs" = subset(data, get(subvar.tbsub()) == 1), 
                   "2-3 yrs" = subset(data, get(subvar.tbsub()) == 2), 
                   " > 3yrs" = subset(data, get(subvar.tbsub()) == 3),
                   "1 month - 3yrs" = subset(data, get(subvar.tbsub()) %in% 0:2),
                   "> 1 yrs" = subset(data, get(subvar.tbsub()) %in% 1:3), 
                   "> 2 yrs" = subset(data, get(subvar.tbsub()) %in% 2:3), 
                   "0-2 yrs" = subset(data, get(subvar.tbsub()) %in% 0:1)
                   
    )
    #data <- data.tbsub()
    label <- out.label[, .SD]
    if (input$data_tbsub == "PPI vs H2RA"){
      label[variable == "group_ppi2", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI"))]
      label[variable == "group_ppi3", ":="(var_label = "PPI User", val_label = c("H2RA" ,"PPI(Weekly to <daily)", "PPI(Daily)"))]
    }
    
    
    var.event <- "event_gc"
    var.day <- "day_gc"
    
    if (input$type_tbsub != "IPTW"){
      data <- data[!( get(var.day) < 365 * input$year_tbsub[1]/12)]
      data[[var.event]] <- ifelse(data[[var.day]] >= 365 * input$year_tbsub[2]/12 & data[[var.event]] == "1", 0,  as.numeric(as.vector(data[[var.event]])))
      data[[var.day]] <- ifelse(data[[var.day]] >= 365 * input$year_tbsub[2]/12, 365 * input$year_tbsub[2]/12, data[[var.day]])
      
      data[[var.event]] <- as.numeric(as.vector(data[[var.event]]))
    } else{
      data <- subset(data, !(get(var.day) < 365 * input$year_tbsub[1]/12))
      data$variables[[var.event]] <- ifelse(data$variables[[var.day]] >= 365 * input$year_tbsub[2]/12 & data$variables[[var.event]] == "1", 0,  as.numeric(as.vector(data$variables[[var.event]])))
      data$variables[[var.day]] <- ifelse(data$variables[[var.day]] >= 365 * input$year_tbsub[2]/12, 365 * input$year_tbsub[2]/12, data$variables[[var.day]])
      
      
      data$variables[[var.event]] <- as.numeric(as.vector(data$variables[[var.event]]))
    }
    
    
    
    form <- as.formula(paste("Surv(", var.day, ",", var.event, ") ~ ", group.tbsub(), sep = ""))
    vs <- input$subvar_tbsub
    tbsub <-  TableSubgroupMultiCox(form, var_subgroups = vs, var_cov = setdiff(input$cov_tbsub, vs),
                                    data=data, time_eventrate = 365 * as.numeric(input$year_tbsub[2])/12, line = F)
    #data[[var.event]] <- ifelse(data[[var.day]] > 365 * 5 & data[[var.event]] == 1, 0,  as.numeric(as.vector(data[[var.event]])))
    
    if (input$type_tbsub != "IPTW"){
      lapply(vs, 
             function(x){
               cc <- data.table(t(c(x, NA, NA)))
               
               dd <- lapply(levels(data[[group.tbsub()]]),
                            function(y){
                              ev <- data[!is.na(get(x)) & get(group.tbsub()) == y, sum(as.numeric(as.vector(get(var.event)))), keyby = get(x)]
                              nn <- data[!is.na(get(x)) & get(group.tbsub()) == y, .N, keyby = get(x)]
                              vv <- paste0(ev[, V1], "/", nn[, N], " (", round(ev[, V1]/ nn[, N] * 100, 1), ")")
                              data.table(ev[, 1], vv)
                            })
               dd.bind <- cbind(dd[[1]], dd[[2]][, -1])
               names(cc) <- names(dd.bind)
               rbind(cc, dd.bind)
             }) %>% rbindlist -> ll
      
      ev.ov <- data[, sum(as.numeric(as.vector(get(var.event)))), keyby = get(group.tbsub())][, V1]
      nn.ov <- data[, .N, keyby = get(group.tbsub())][, N]
      
      ov <- data.table(t(c("OverAll", paste0(ev.ov, "/", nn.ov, " (", round(ev.ov/nn.ov * 100, 1), ")"))))
      names(ov) <- names(ll)
      cn <- rbind(ov, ll)
    } else{
      lapply(vs, 
             function(x){
               cc <- data.table(t(c(x, NA, NA)))
               
               dd <- lapply(levels(data$variables[[group.tbsub()]]),
                            function(y){
                              ev <- round(svytable(as.formula(paste0("~", var.event, "+", x)), design = subset(data, !is.na(get(x)) & get(group.tbsub()) == y))[2, ], 1)
                              nn <- round(svytable(as.formula(paste0("~", x)), design = subset(data, !is.na(get(x)) & get(group.tbsub()) == y)), 1)
                              vv <- paste0(ev, "/", nn, " (", round(ev/ nn * 100, 1), ")")
                              data.table(ev, vv)
                            })
               dd.bind <- cbind(dd[[1]], dd[[2]][, -1])
               names(cc) <- names(dd.bind)
               rbind(cc, dd.bind)
             }) %>% rbindlist -> ll
      
      #ev.ov <- data$variables[, sum(as.numeric(as.vector(get(var.event)))), keyby = get(group.tbsub())][, V1]
      ev.ov <- round(svytable(as.formula(paste0("~", var.event, "+", group.tbsub())), design = data)[2, ], 1)
      nn.ov <- round(svytable(as.formula(paste0("~", group.tbsub())), design = data), 1)
      
      ov <- data.table(t(c("OverAll", paste0(ev.ov, "/", nn.ov, " (", round(ev.ov/nn.ov * 100, 1), ")"))))
      names(ov) <- names(ll)
      cn <- rbind(ov, ll)
    }
    
    
    names(cn)[-1] <- label[variable == group.tbsub(), val_label]
    tbsub <- cbind(Variable = tbsub[, 1], cn[, -1], tbsub[, c(7, 8, 4, 5, 6, 9, 10)])
    
    tbsub[-1, 1] <- unlist(lapply(vs, function(x){c(label[variable == x, var_label][1], paste0("     ", label[variable == x, val_label]))}))
    colnames(tbsub)[1:6] <- c("Subgroup", paste0("N(%): ", label[variable == group.tbsub(), val_label]), paste0(paste(input$year_tbsub, collapse = "~"), "-month KM rate(%): ", label[variable == group.tbsub(), val_label]), "HR")
    #colnames(tbsub)[c(4)] <- c("HR")
    #tbsub[27:30, Subgroup := c("EES", "ZES", "BES", "Others")]
    
    return(tbsub)
  })
  
  output$tablesub <- renderDT({
    
    datatable(tbsub(), caption = paste0(input$dep_tbsub, " subgroup analysis: ", input$type_tbsub), rownames = F, extensions= "Buttons",
              options = c(opt.tb1(paste0("tbsub_original_", input$type_tbsub)),
                          list(scrollX = TRUE, columnDefs = list(list(className = 'dt-right', targets = 0)))
              )) 
    
  })
  
  
  
  
  
  output$forest <- downloadHandler(
    filename =  function() {
      paste(input$dep_tbsub,"_forestplot.emf", sep="")
      
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      withProgress(message = 'Download in progress',
                   detail = 'This may take a while...', value = 0, {
                     for (i in 1:15) {
                       incProgress(1/15)
                       Sys.sleep(0.01)
                     }
                     
                     label <- out.label
                     devEMF::emf(file, width = 15, height = 20, coordDPI = 300, emfPlus = F)
                     data <- tbsub()
                     
                     
                     tabletext <- cbind(c("Subgroup","\n",data$Subgroup),
                                        c("N(%)\nNo imaging", "\n" , data[[paste0("N(%): ", label[variable == group.tbsub(), val_label][1])]]),
                                        c("N(%)\nImaging", "\n", data[[paste0("N(%): ", label[variable == group.tbsub(), val_label][2])]]),
                                        c("5-yr KM rate(%)\nNo Imaging", "\n", data[[paste0(paste(input$year_tbsub, collapse = "~"), "-month KM rate(%): ", label[variable == group.tbsub(), val_label][1])]]),
                                        c("5-yr KM rate(%)\nImaging", "\n", data[[paste0(paste(input$year_tbsub, collapse = "~"), "-month KM rate(%): ", label[variable == group.tbsub(), val_label][2])]]),
                                        c("P Value","\n",data$`P value`),
                                        c("P for interaction","\n",data$`P for interaction`))
                     
                     tabletext <- tabletext[, c(1, 6, 7)]
                     ## Save as tiff 
                     forestplot::forestplot(labeltext=tabletext, graph.pos=2, xticks = c(0.1, 0.5, 1, 2, 10), xlog= T, align = c("r", rep("c", ncol(tabletext) - 1)),                          ## graph.pos- column number
                                            mean=c(NA,NA,as.numeric(data$HR)), 
                                            lower=c(NA,NA,as.numeric(data$Lower)), upper=c(NA,NA,as.numeric(data$Upper)),
                                            title="Hazard Ratio",
                                            xlab="<---No imaging risky ---    ---Imaging risky --->",    ## You cas modify this.
                                            hrzl_lines=list("3" = gpar(lwd=1, col="#99999922")
                                            ),
                                            
                                            txt_gp=fpTxtGp(label=gpar(cex=1.25),
                                                           ticks=gpar(cex=1.1),
                                                           xlab=gpar(cex = 1.2),
                                                           title=gpar(cex = 1.2)),
                                            col=fpColors(box="black", lines="black", zero = "gray50"),
                                            zero=1, cex=0.9, lineheight = "auto", boxsize=0.2, colgap=unit(6,"mm"),
                                            lwd.ci=2, ci.vertices=TRUE, ci.vertices.height = 0.2) 
                     grDevices::dev.off()
                   })
      
    }
  )
  
  
}


shinyApp(ui, server)