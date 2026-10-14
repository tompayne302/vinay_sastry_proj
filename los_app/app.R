library(shiny)
library(bslib)
library(ggplot2)

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
    col_widths = c(4, 4, 4),
    value_box(title = "Predicted mean length of stay",
              value = textOutput("pred_mean"), theme = "primary"),
    value_box(title = "90% confidence interval for the mean",
              value = textOutput("ci_mean"), theme = "info"),
    value_box(title = "90% prediction interval for an individual patient",
              value = textOutput("pred_range"), theme = "secondary")
  ),
  card(
    card_header("Likely length of stay for this patient"),
    plotOutput("pred_plot", height = "320px")
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
      "The interval is the 5th to 95th percentile of the negative binomial distribution ",
      "(dispersion theta = %.2f) around the predicted mean, i.e. the range within which ",
      "90%% of individual patients with these characteristics are expected to fall. ",
      "The confidence interval for the mean is the 5th to 95th percentile of the predicted ",
      "mean across bootstrap refits of the model; it reflects sampling and imputation ",
      "uncertainty, conditional on the selected predictors, and describes the average ",
      "patient with these characteristics, not an individual."),
      model$n_dev, model$theta)),
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

  # Design vector aligned with c("(Intercept)", names(coefs))
  x_vec <- reactive({
    x <- setNames(numeric(length(model$coefs) + 1), c("(Intercept)", names(model$coefs)))
    x[["(Intercept)"]] <- 1
    for (s in spec) {
      v <- input[[paste0("in_", s$var)]]
      if (s$type == "factor") {
        nm <- paste0(s$var, v)
        if (nm %in% names(x)) x[[nm]] <- 1
      } else {
        x[[s$var]] <- as.numeric(v)
      }
    }
    x
  })

  linear_predictor <- reactive({
    x <- x_vec()
    sum(x * c("(Intercept)" = model$intercept, model$coefs)[names(x)])
  })

  # Predicted mean under each bootstrap model
  boot_mu <- reactive({
    x <- x_vec()
    as.numeric(exp(model$boot_coefs[, names(x), drop = FALSE] %*% x))
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

  output$ci_mean <- renderText({
    if (!complete()) return("-")
    ci <- quantile(boot_mu(), c(0.05, 0.95), names = FALSE)
    sprintf("%.1f to %.1f days", ci[1], ci[2])
  })

  output$pred_range <- renderText({
    if (!complete()) return("-")
    q <- qnbinom(c(0.05, 0.95), mu = mu(), size = model$theta)
    sprintf("%d to %d days", q[1], q[2])
  })

  output$pred_plot <- renderPlot({
    req(complete())
    m <- mu()
    q <- qnbinom(c(0.05, 0.95), mu = m, size = model$theta)
    x <- 0:max(qnbinom(0.995, mu = m, size = model$theta), q[2] + 5)
    d <- data.frame(los = x, p = dnbinom(x, mu = m, size = model$theta))
    d$inside <- d$los >= q[1] & d$los <= q[2]
    ci <- quantile(boot_mu(), c(0.05, 0.95), names = FALSE)
    y_ci <- -0.06 * max(d$p)

    ggplot(d, aes(los, p, fill = inside)) +
      geom_col(width = 1) +
      geom_vline(xintercept = m, linewidth = 1) +
      annotate("label", x = m, y = max(d$p), hjust = -0.1, vjust = 1,
               label = sprintf("Mean: %.1f days", m)) +
      annotate("segment", x = ci[1], xend = ci[2], y = y_ci, yend = y_ci,
               linewidth = 2, colour = "black") +
      annotate("point", x = m, y = y_ci, size = 4, colour = "black") +
      annotate("text", x = ci[2], y = y_ci, hjust = -0.1, size = 4.5,
               label = "90% CI for the mean") +
      scale_fill_manual(values = c("TRUE" = "#2c7fb8", "FALSE" = "grey80"),
                        labels = c("TRUE" = "Within 90% interval", "FALSE" = "Outside interval"),
                        name = NULL) +
      scale_y_continuous(breaks = function(l) pretty(c(0, l[2])),
                         labels = scales::percent_format(accuracy = 0.1)) +
      labs(x = "Length of stay (days)", y = "Probability of each exact LOS") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top", panel.grid.minor = element_blank())
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
