library(shiny)
library(bslib)

# Model exported from index.qmd (Step 8)
model <- readRDS("model.rds")
spec  <- model$var_spec

fmt_num <- function(x) format(round(x, 2), nsmall = 0, trim = TRUE)

# ---- Input widgets ---------------------------------------------------------

make_input <- function(s) {
  id <- paste0("in_", s$var)
  widget <- switch(
    s$type,
    numeric = numericInput(id, s$label, value = s$median, min = NA, max = NA, step = 1),
    binary  = radioButtons(id, s$label, choices = c("No" = 0, "Yes" = 1),
                           selected = 0, inline = TRUE),
    factor  = {
      ch <- setNames(s$levels, ifelse(
        s$levels == "Other",
        paste0("Other (", paste(s$collapsed, collapse = ", "), ", or any other code)"),
        s$levels))
      selectInput(id, s$label, choices = ch)
    }
  )
  hint <- if (s$type == "numeric") {
    helpText(sprintf("Development data range: %s to %s (typical: %s to %s)",
                     fmt_num(s$min), fmt_num(s$max), fmt_num(s$p025), fmt_num(s$p975)))
  }
  div(widget, hint, uiOutput(paste0("flag_", s$var)), hr(style = "margin: 0.6rem 0;"))
}

badge <- function(text, colour) {
  span(class = paste0("badge bg-", colour), text)
}

# ---- UI --------------------------------------------------------------------

ui <- page_sidebar(
  title = "Rehabilitation length of stay calculator",
  theme = bs_theme(version = 5),
  sidebar = do.call(sidebar, c(
    list(width = 420,
         h5("Patient characteristics"),
         helpText("Only predictors retained by the LASSO model are required.")),
    unname(lapply(spec, make_input))
  )),
  uiOutput("warning_banner"),
  layout_columns(
    col_widths = c(6, 6),
    value_box(title = "Predicted mean length of stay",
              value = textOutput("pred_mean"), theme = "primary"),
    value_box(title = "Typical range for an individual patient (5th to 95th percentile)",
              value = textOutput("pred_range"), theme = "secondary")
  ),
  card(
    card_header("Input check"),
    tableOutput("check_table")
  ),
  card(
    card_header("About this tool"),
    p(sprintf(paste0(
      "Negative binomial regression with LASSO penalisation, developed in %d patients. ",
      "Predicted mean LOS = exp(intercept + sum of coefficient x value). ",
      "The range shown is the 5th to 95th percentile of the negative binomial distribution ",
      "(dispersion theta = %.2f) around the predicted mean; it does not include uncertainty ",
      "in the model coefficients."), model$n_dev, model$theta)),
    p("Inputs outside the range of the development data are extrapolations and predictions ",
      "for such patients should not be relied upon. This tool is intended to support, not ",
      "replace, clinical judgement.")
  )
)

# ---- Server ----------------------------------------------------------------

check_numeric <- function(x, s) {
  if (is.null(x) || is.na(x)) return("missing")
  if (x < s$min || x > s$max) return("out")
  if (x < s$p025 || x > s$p975) return("edge")
  "ok"
}

server <- function(input, output, session) {

  status <- reactive({
    vapply(spec, function(s) {
      v <- input[[paste0("in_", s$var)]]
      switch(s$type,
             numeric = check_numeric(suppressWarnings(as.numeric(v)), s),
             binary  = if (is.null(v)) "missing" else "ok",
             factor  = if (is.null(v) || v == "") "missing" else "ok")
    }, character(1))
  })

  # Per-variable markers
  lapply(spec, function(s) {
    output[[paste0("flag_", s$var)]] <- renderUI({
      st <- status()[[s$var]]
      switch(st,
             ok      = badge("In bounds", "success"),
             edge    = tagList(badge("In bounds, but unusual", "warning"),
                               helpText("Within the observed range but in the extreme 5% of the development data.")),
             out     = tagList(badge("OUT OF BOUNDS", "danger"),
                               helpText(sprintf("Outside the development data range (%s to %s). Prediction is an extrapolation.",
                                                fmt_num(s$min), fmt_num(s$max)))),
             missing = badge("Enter a value", "danger"))
    })
  })

  linear_predictor <- reactive({
    lp <- model$intercept
    for (s in spec) {
      v <- input[[paste0("in_", s$var)]]
      if (s$type == "factor") {
        nm <- paste0(s$var, v)
        if (nm %in% names(model$coefs)) lp <- lp + model$coefs[[nm]]
      } else {
        lp <- lp + model$coefs[[s$var]] * as.numeric(v)
      }
    }
    lp
  })

  complete <- reactive(!any(status() == "missing"))
  mu <- reactive({
    req(complete())
    exp(linear_predictor())
  })

  output$pred_mean <- renderText({
    if (!complete()) return("Enter all inputs")
    sprintf("%.1f days", mu())
  })

  output$pred_range <- renderText({
    if (!complete()) return("-")
    q <- qnbinom(c(0.05, 0.95), mu = mu(), size = model$theta)
    sprintf("%d to %d days", q[1], q[2])
  })

  output$warning_banner <- renderUI({
    st <- status()
    if (any(st == "out")) {
      bad <- vapply(spec[st == "out"], `[[`, character(1), "label")
      div(class = "alert alert-danger",
          strong("Warning: out-of-bounds input. "),
          "The prediction is an extrapolation beyond the development data for: ",
          paste(bad, collapse = "; "), ".")
    } else if (any(st == "edge")) {
      div(class = "alert alert-warning",
          strong("Note: "), "one or more inputs are in the extreme 5% of the development data, ",
          "so the prediction is less certain.")
    }
  })

  output$check_table <- renderTable({
    st <- status()
    data.frame(
      Variable = vapply(spec, `[[`, character(1), "label"),
      Status = unname(c(ok = "In bounds", edge = "Unusual", out = "OUT OF BOUNDS",
                        missing = "Missing")[st]),
      check.names = FALSE
    )
  })
}

shinyApp(ui, server)
