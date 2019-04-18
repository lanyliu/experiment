# (1) First time install
# R
# > install.packages(c("shiny", "shinydashboard", "pwr"))
#   select CRAN mirrors 71
# > q()

# (2) Run app
# Rscript shinyApp.R
# --------------------------------------------------------------------------

# import library
library(shiny)
library(shinydashboard)
library(pwr)


t.test.func <- function(m0, m1, s0, s1, n0, n1, dm=0, equal.variance=FALSE)
{ # m0, m1: sample means
  # s0, s1: sample standard deviations
  # n0, n1: sample sizes
  # dm: expected difference in means to be tested for. Default 0
  # equal.variance: whether or not to assume equal variance. Default FALSE
  if(equal.variance==FALSE) { # welch-satterthwaite df
    se <- sqrt((s0^2/n0) + (s1^2/n1))
    df <- ((s0^2/n0 + s1^2/n1)^2) / ((s0^2/n0)^2/(n0-1) + (s1^2/n1)^2/(n1-1))
  } else { # pooled standard deviation, scaled by the sample sizes
    se <- sqrt((1/n0+1/n1) * ((n0-1)*s0^2 + (n1-1)*s1^2)/(n0+n1-2))
    df <- n0+n1-2}
  
  t <- (m0-m1-dm)/se
  dat <- c(m0-m1, se, t, 2*pt(-abs(t), df))
  
  cat("Delta of means: ", dat[1], "\n")
  cat("Std error: ", dat[2], "\n")
  cat("t: ", dat[3], "\n")
  cat("p-value: ", dat[4])
}

ui <- dashboardPage(
  
  dashboardHeader(title = "Experiment tools"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Idea sample", tabName = "sample", icon = icon("users")),
      menuItem("Test proportions", tabName = "psig", icon = icon("th")),
      menuItem("Test means", tabName = "sig", icon = icon("th"))
    )),
  
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "sample",
              h2("Idea sample"),
              fluidRow(
                box(
                  title = "Idea sample size",
                  verbatimTextOutput("pwr")
                  ),
                
                box(
                  title = "Settings",
                  numericInput("sample.std", "Current standard deviation:", 10),
                  numericInput("sample.m0", "Current mean value:", 10),
                  numericInput("sample.m1", "Expect mean value:", 11),
                  numericInput("sample.power", "Power:", 0.8),
                  numericInput("sample.sig", "Significant level:", 0.05),
                  
                  # Include clarifying text ----
                  helpText("Note: recommend Power and Significant level as default value.")
                  #actionButton("calculate", "Calculate")
                ))),
      
      # Second tab content
      tabItem(tabName = "psig",
              h2("Significant test on proportions"),
              fluidRow(
                box(
                  title = "Significant test",
                  verbatimTextOutput("p.sig")
                ),
                
                box(
                  title = "Settings",
                  numericInput("psig.n0", "Sample size of control:", 8700),
                  numericInput("psig.n1", "Sample size of variant:", 963),
                  numericInput("psig.e0", "Response size of control:", 375),
                  numericInput("psig.e1", "Response size of variant:", 34),
                  numericInput("psig.ci", "Confidence level", 0.95),
                  
                  # Include clarifying text ----
                  helpText("Note: Test for equality of proportions between 2 groups")
                  #actionButton("calculate", "Calculate")
                ))),
              
        # 3rd tab content
        tabItem(tabName = "sig",
                h2("Significant test on means"),
                fluidRow(
                  box(
                    title = "Control group Settings",
                    numericInput("sig.n0", "Sample size of control:", 5000),
                    numericInput("sig.m0", "Mean value of control:", 0.42),
                    numericInput("sig.s0", "Standard deviation of control:", 0.43)
                  ),
                  
                  box(
                    title = "Variant group Settings",
                    numericInput("sig.n1", "Sample size of variant:", 5000),
                    numericInput("sig.m1", "Mean value of variant:", 0.42),
                    numericInput("sig.s1", "Standard deviation of variant:", 0.43),
                    numericInput("sig.dm", "H0 expected difference in means:", 0),
                    numericInput("sig.ci", "Confidence level", 0.95),
                    radioButtons("sig.equal.variance", label = "Is equal variance:",
                                 choices = list("True" = TRUE, "False" = FALSE), 
                                 selected = FALSE),
                    
                    # Include clarifying text ----
                    helpText("Note: Test for equality of means between 2 groups")
                    #actionButton("calculate", "Calculate")
                  ),
                  
                  box(
                    title = "Significant test",
                    verbatimTextOutput("sig")
                  )))
  ))
)



# Define server logic to summarize and view selected dataset ----
server <- function(input, output) {
  
  output$pwr <- renderPrint({
    # https://stats.idre.ucla.edu/r/dae/power-analysis-for-two-group-independent-sample-t-test/
    d = abs(input$sample.m1 - input$sample.m0)/input$sample.std
    pwr.t.test(d = d, power = input$sample.power, sig.level = input$sample.sig, type="two.sample")
  })
  
  output$p.sig <- renderPrint({
    
    ##' method 1: fisher exact test
    print(
      fisher.test(rbind(c(input$psig.e0, input$psig.n0-input$psig.e0), c(input$psig.e1, input$psig.n1-input$psig.e1)), conf.level = input$psig.ci)
      )
    
    ##' method 2: Normal Approximation to Binomial
    #p0.hat <- input$e0/input$n0
    #p1.hat <- input$e1/input$n1
    #p.hat <- (input$n0*p0.hat + input$n1*p1.hat)/(input$n0 + input$n1)
    #z <- (p0.hat - p1.hat)/sqrt(p.hat*(1-p.hat)*(1/input$n0 + 1/input$n1))
    #pnorm(z, lower.tail = F)
    
    ##' method 3: 2-sample test for equality of proportions
    print(
      prop.test(c(input$psig.e0, input$psig.e1), c(input$psig.n0, input$psig.n1), correct=FALSE, conf.level = 0.95)
    )
    
    ##' method 4: Exact binomial test
    # print(binom.test(x=input$e1, n=input$n1, p=input$e0/input$e0, alternative = "greater", conf.level = 0.95))
  })
  
  output$sig <- renderPrint({
    # 2 sample t test
    t.test.func(input$sig.m0, input$sig.m1, input$sig.s0, input$sig.s1, input$sig.n0, input$sig.n1, input$sig.dm, input$sig.equal.variance)
  })
  
}

# Create Shiny app ----
shinyApp(ui, server)
