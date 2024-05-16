#'
#' plot_RMSE_LPD
#'
#' Reads in fitting and prediction metrics from CV and plots RMSE and relative LPD
#'
#' @param df Tibble. Tibble of performance metrics of all different models.
#' @param colour_of_train String. Colour of boxplots showing RMSE on training data. 
#' @param colour_of_test String. Colour of boxplots showing Relative LPD and RMSE on testing data. 
#' @param og_data Tibble. Tibble of all data (training + testing) models were fit to.
#'
#' @return ggplot Boxplot of RMSE and relative LPD of models
#'
plot_RMSE_LPD <- function(df, colour_of_train, colour_of_test, og_data){
  df %>%
    filter(is.finite(Mean)) %>%
    filter(Metric%in%c("LPD", "RMSE")) %>%
    group_by(Metric, Training) %>%
    mutate(Diff = max(Mean)-Mean + 1) %>%
    ungroup() %>%
    mutate(reference = factor(reference, levels=c("Our Model", "Reference Logistic Model \nwithout Calibration", "Reference Models \nwithout Calibration", "Reference Models \nwith Calibration", "Extensions"), 
                              labels = c(expression(paste("Our Model", sep="")), expression(paste(atop(textstyle('Reference Logistic Model'), textstyle('without Calibration')), sep="")), expression(paste(atop(textstyle('Reference Models'), textstyle('without Calibration')), sep="")), expression(paste(atop(textstyle('Reference Models'), textstyle('with Calibration')), sep="")), "Extensions"))) %>%
    mutate(values = case_when(Metric=="LPD" ~ Diff, Metric=="RMSE" ~ Mean)) %>%
    filter(case_when(Metric=="LPD" ~ !Training, Metric=="RMSE" ~ Training | !Training)) %>%
    mutate(Metric_for_plot = factor(Metric, levels = c("RMSE", "LPD"), labels=c(expression(log[10]~"RMSE"), expression(log[10]~"Relative LPD")))) -> df_for_plot
  
  df_for_plot %>%
    ggplot(aes(y = Model, x = values, color=factor(Training, levels = c(TRUE, FALSE)))) +
    geom_boxplot(linewidth=0.3, outlier.size=0.3) +
    geom_vline(data = filter(df_for_plot, Metric=="RMSE"), aes(xintercept = accuracy(og_data$OD, rep(mean(og_data$OD), length(og_data$OD)))[2]), color = "darkgray", linetype="dashed", alpha=0.9, size=1) +
    facet_grid(reference~Metric_for_plot, scales = "free", space = "free_y", switch="x", labeller="label_parsed") +
    theme_bw(base_size=11) +
    theme(strip.placement = "outside", strip.background = element_blank()) +
    xlab("") +
    ylab("Model") +
    scale_y_discrete(breaks=c("RW", "mixed_logistic_OD_calibration", "gompertz", "exponential", "gompertz_OD_calibration", "GP_OD_calibration", "GP", "logistic","exponential_OD_calibration", "logistic_OD_calibration", "logistic_OD_calibration_no_delay"),
                     labels=c("Random Walk", "Mixed logistic-\nOD-calibration", "Gompertz-OD", "Exponential-OD", "Gompertz-OD\n-calibration", "GP-OD\n-calibration", "GP-OD", "Logistic-OD", "Exponential-OD\n-calibration", "Logistic-OD\n-calibration", "No-delay logistic-\nOD-calibration")) +
    theme(axis.text.x = element_text(angle=0), strip.text.x=element_text(size=11)) +
    scale_color_manual(name = "Data Split", breaks = c(TRUE, FALSE), labels = c("Train", "Test"), values = c(colour_of_train, colour_of_test)) +
    theme(strip.text.y=element_text(size=9, angle=0), text = element_text(size=11), legend.position = "bottom") -> p
  return(p)
}

