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

# ---- Prediction plot -------------------------------------------------------

plot_prediction <- function(m, ci, theta) {
  q <- qnbinom(c(0.05, 0.95), mu = m, size = theta)
  x_max <- max(qnbinom(0.99, mu = m, size = theta), ci[2]) + 3
  d <- data.frame(los = 0:x_max, p = dnbinom(0:x_max, mu = m, size = theta))
  d$band <- factor(ifelse(d$los >= q[1] & d$los <= q[2],
                          sprintf("Inside 90%% prediction interval (%d\u2013%d days)", q[1], q[2]),
                          "Outside 90% prediction interval"))
  p_max <- max(d$p)
  y_ci <- p_max * 1.15

  ggplot(d, aes(los, p)) +
    geom_col(aes(fill = band), width = 0.9) +
    geom_vline(xintercept = m, linetype = "dashed", linewidth = 0.6) +
    # CI for the mean, drawn as a bar with end ticks above the distribution
    annotate("segment", x = ci[1], xend = ci[2], y = y_ci, yend = y_ci, linewidth = 1) +
    annotate("segment", x = ci, xend = ci, y = y_ci - p_max * 0.03,
             yend = y_ci + p_max * 0.03, linewidth = 1) +
    annotate("point", x = m, y = y_ci, size = 3.5) +
    annotate("label", x = m, y = y_ci + p_max * 0.06, vjust = 0, size = 4.5,
             fill = "white",
             label = sprintf("Mean %.1f days (90%% CI %.1f\u2013%.1f)", m, ci[1], ci[2])) +
    scale_fill_manual(values = c("#2c7fb8", "grey80"), name = NULL,
                      breaks = levels(d$band)[c(grep("^Inside", levels(d$band)),
                                                grep("^Outside", levels(d$band)))]) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 8),
                       expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 0.1),
                       expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(ylim = c(0, p_max * 1.35), clip = "off") +
    labs(x = "Length of stay (days)", y = "Probability") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 12),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.title = element_text(face = "bold"))
}

# ---- UI --------------------------------------------------------------------

app_theme <- bs_theme(version = 5) |>
  bs_add_rules("
    .value-box-value { font-size: 1.9rem !important; line-height: 1.2; white-space: nowrap; }
    .value-box-title { font-size: 1rem !important; font-weight: 600; }
  ")

result_box <- function(title, output_id, caption, theme) {
  value_box(
    title = title,
    value = textOutput(output_id, inline = TRUE),
    p(caption, class = "small mb-0"),
    theme = theme,
    min_height = "160px"
  )
}

ui <- page_sidebar(
  title = "Rehabilitation length of stay calculator",
  theme = app_theme,
  # Not fillable: lets cards keep their natural height and the page scroll
  fillable = FALSE,
  sidebar = do.call(sidebar, c(
    list(width = 420,
         h5("Patient characteristics"),
         helpText("Only predictors retained by the LASSO model are required.")),
    unname(lapply(spec, make_input))
  )),
  uiOutput("warning_banner"),
  layout_columns(
    col_widths = breakpoints(sm = 12, lg = 4),
    result_box("Predicted mean LOS", "pred_mean",
               "Expected average stay for patients like this", "primary"),
    result_box("90% CI for the mean", "ci_mean",
               "Uncertainty in the average (bootstrap)", "info"),
    result_box("90% prediction interval", "pred_range",
               "Likely range for an individual patient", "secondary")
  ),
  card(
    full_screen = TRUE,
    card_header("Likely length of stay for this patient"),
    plotOutput("pred_plot", height = "480px"),
    card_footer(class = "small text-muted",
      "Bars: model-predicted probability of each length of stay. Blue bars fall inside the ",
      "90% prediction interval. The black point and line at the top show the predicted mean ",
      "and its 90% bootstrap confidence interval.")
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
    sprintf("%.1f\u2013%.1f days", ci[1], ci[2])
  })

  output$pred_range <- renderText({
    if (!complete()) return("-")
    q <- qnbinom(c(0.05, 0.95), mu = mu(), size = model$theta)
    sprintf("%d\u2013%d days", q[1], q[2])
  })

  output$pred_plot <- renderPlot({
    req(complete())
    plot_prediction(mu(), quantile(boot_mu(), c(0.05, 0.95), names = FALSE), model$theta)
  }, res = 110)

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
