library(lazyeval)

theme_set(theme_classic())
colors = scale_color_manual("", values=c("black", "red"))
fills = scale_fill_manual("", values=c("black", "red"))

plot_syl_transition_matrices_by_date = function(info, labels=NULL) {
  d = info %>% group_by(date) %>% do(process_syllable_matrix(., labels))
  gg = ggplot(d, aes(To, From, fill=value)) + geom_tile() + facet_wrap(~date)
  print(gg)
  return(d)
}

plot_feature_vs_time = function(data, value=NULL, reference_date=NULL, ylab=NULL, xlab=NULL, ylim=NULL, add_baseline=FALSE, grouping_factor=NULL) {
  if (is.null(value)) 
    stop("Specify value to plot")
  
  value_obj = as.symbol(value)
  ref_date = as.symbol(reference_date)
  m = data %>% group_by(date) %>% summarize(mean=mean(value_obj), sd=sd(value_obj))
  m$date = paste(m$date, "12:00:00", sep=" ")
  m$date = as.POSIXct(m$date)
  m$rel_date = with(m, as.numeric(difftime(date, as.POSIXct(reference_date), units="days")))
  
  
  theme_set(theme_classic())
  gg =ggplot(data, aes_string("rel_mtime", value)) + geom_point(alpha=I(1/50)) 
  gg = gg + geom_pointrange(data=m, aes(x=as.numeric(rel_date), y=mean, ymin=(mean-sd), ymax=(mean+sd))) 
  gg = gg + geom_vline(xintercept=0, linetype=2)
  

  if (!is.null(grouping_factor)) {
    grouping_obj = as.symbol(grouping_factor)
    gg = gg + facet_wrap(~grouping_obj)
  }
  
  # Define labs
  if (!is.null(ylab) & is.null(xlab)) {
    gg = gg +  labs(y=ylab)
  } else if (!is.null(xlab) & is.null(ylab)) {
    gg = gg + labs(x=xlab)
  } else if (!is.null(xlab) & !is.null(ylab)) {
    gg = gg + labs(x=xlab, y=ylab)
  }
  
  # Define axis limits
  if (!is.null(ylim)) {
    gg = gg + ylim(ylim)
  }
  
  # Plot baseline statistic
  if (add_baseline) {
    if (is.null(grouping_factor))
      baseline_data = data %>% filter(rel_mtime<ref_date) %>% summarize(baseline_stat = mean(value_obj))
    gg = gg + geom_hline(data=baseline_data, aes(yintercept=baseline_stat), linetype=2)
  }
  gg
}

plot_feature_vs_time_mult = function(data, 
                                     value=NULL, 
                                     reference_date=NULL, 
                                     ylab=NULL, 
                                     xlab=NULL, 
                                     ylim=NULL, 
                                     add_baseline=FALSE, 
                                     grouping_factor=NULL) {
  
  #value = "tempo"
  #ylab = "Syllables / second"
  #xlab = xaxis
  value_obj = as.symbol(value)
  ref_obj = as.symbol(reference_date)
  group_obj = as.symbol(grouping_factor)
  group_dots = c("date", reference_date, grouping_factor)
  m = data %>% group_by_(.dots=group_dots) %>% summarize(mean=mean(value_obj), sd=sd(value_obj))
  m$date = paste(m$date, "12:00:00", sep=" ")
  m$date = as.POSIXct(m$date)
  m = m %>% mutate_(.dots=setNames(list(interp(quote(as.numeric(difftime(date, 
                                                                as.POSIXct(ref_obj), 
                                                                units="days"))),
                           ref_obj=as.name(reference_date))),
                           c("rel_date")))
  
  theme_set(theme_classic())
  gg =ggplot(data, aes_string("rel_mtime", value)) + geom_point(alpha=I(1/100)) 
  gg = gg + geom_pointrange(data=m, aes(x=as.numeric(rel_date), y=mean, ymin=(mean-sd), ymax=(mean+sd))) 
  gg = gg + geom_vline(xintercept=0, linetype=2)
  
  
  #if (!is.null(grouping_factor)) {
  #  grouping_obj = as.symbol(grouping_factor)
  gg = gg + facet_wrap(as.formula(paste("~", grouping_factor, sep="")))
  #}
  
  # Define labs
  if (!is.null(ylab) & is.null(xlab)) {
    gg = gg +  labs(y=ylab)
  } else if (!is.null(xlab) & is.null(ylab)) {
    gg = gg + labs(x=xlab)
  } else if (!is.null(xlab) & !is.null(ylab)) {
    gg = gg + labs(x=xlab, y=ylab)
  }
  
  # Define axis limits
  if (!is.null(ylim)) {
    gg = gg + ylim(ylim)
  }
  
  # Plot baseline statistic
  #if (add_baseline) {
  #if (is.null(grouping_factor))
  baseline_data = data %>% 
    filter_(.dots=interp(quote(rel_mtime<ref), ref=reference_date)) %>% 
    group_by_(.dots=grouping_factor) %>% 
    summarize(baseline_stat = mean(value_obj))
  gg = gg + geom_hline(data=baseline_data, aes(yintercept=baseline_stat), linetype=2)
  #}
  gg
}
