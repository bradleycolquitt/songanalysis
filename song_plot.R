library(lazyeval)
library(scales)

#theme_set(theme_classic())
colors = scale_color_manual("", values=c("black", "red"))
fills = scale_fill_manual("", values=c("black", "red"))
# 
# plot_syl_transition_matrices_by_date = function(info, labels=NULL) {
#   d = info %>% group_by(date) %>% do(process_syllable_matrix(., labels))
#   gg = ggplot(d, aes(To, From, fill=value)) + geom_tile() + facet_wrap(~date)
#   print(gg)
#   return(d)
# }

plot_feature_vs_time = function(data, value=NULL, reference_datetime=NULL, ylab=NULL, xlab=NULL, ylim=NULL, add_baseline=FALSE, grouping_factor=NULL) {
  if (is.null(value)) 
    stop("Specify value to plot")
  
  value_obj = as.symbol(value)
  ref_date = as.symbol(reference_date)
  group_dots = c(grouping_factor, reference_datetime)
  #dots = list(lazyeval::interp(~median(v), v = as.name(value)), 
  #            lazyeval::interp(~sd(v), v = as.name(value)))
  sum_dots = list(lazyeval::interp(~mean(v, na.rm=T), v = as.name(value)),
                  lazyeval::interp(~sd(v, na.rm=T), v = as.name(value)))
  m = data %>% group_by_(.dots=group_dots)%>% summarize_(.dots=setNames(sum_dots, c("mean", "sd")))
  m$date = paste(m$date, "12:00:00", sep=" ")
  m$date = as.POSIXct(m$date)
  m = m %>% mutate_(.dots=setNames(list(lazyeval::interp(quote(as.numeric(difftime(date, 
                                                                                   as.POSIXct(ref_obj), 
                                                                                   units="days"))),
                                                         ref_obj=as.name(reference_date))),
                                   c("rel_date")))
  m = m %>% mutate(upper=mean+sd,
                   lower=mean-sd)
  #m = data %>% group_by(date) %>% summarize_(.dots=setNames(dots, c("mean", "sd")))
  #m$date = paste(m$date, "12:00:00", sep=" ")
  #m$date = as.POSIXct(m$date)
  #m$rel_date = with(m, as.numeric(difftime(date, as.POSIXct(reference_date), units="days")))
  
  
  theme_set(theme_classic())
  gg =ggplot(data, aes_string("rel_mtime", value)) + geom_point(alpha=I(1/50)) 
  gg = gg + geom_pointrange(data=m, aes(x=as.numeric(rel_date), y=mean, ymin=lower, ymax=upper))
 #gg = gg + geom_pointrange(data=m, aes(x=as.numeric(rel_date), y=mean, ymin=(mean-sd), ymax=(mean+sd))) 
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
                                     add_smooth=FALSE,
                                     color_factor=NULL,
                                     grouping_factor=NULL,
                                     grouping_factor2=NULL,
                                     subsample=NULL) {
  
  #value = "tempo"
  #ylab = "Syllables / second"
  #xlab = xaxis
  value_obj = as.symbol(value)
  ref_obj = as.symbol(reference_date)
  group_obj = as.symbol(grouping_factor)
  if (is.null(color_factor)) {
    group_dots = c("date", reference_date, grouping_factor)
  } else {
    group_dots = c("date", reference_date, grouping_factor, color_factor)
  }
  sum_dots = list(lazyeval::interp(~mean(v, na.rm=T), v = as.name(value)),
                  lazyeval::interp(~sd(v, na.rm=T), v = as.name(value)))
  m = data %>% group_by_(.dots=group_dots) %>% summarize_(.dots=setNames(sum_dots, c("mean", "sd")))
  m$date = paste(m$date, "12:00:00", sep=" ")
  m$date = as.POSIXct(m$date)
  m = m %>% mutate_(.dots=setNames(list(lazyeval::interp(quote(as.numeric(difftime(date, 
                                                                as.POSIXct(ref_obj), 
                                                                units="days"))),
                           ref_obj=as.name(reference_date))),
                           c("rel_date")))
  m = m %>% mutate(upper=mean+sd,
                   lower=mean-sd)
  
  theme_set(theme_classic())
  
  if (!is.null(subsample)) {
    data = data %>% group_by_(.dots=group_dots) %>% sample_frac(size = subsample) %>% ungroup(.)
  }
  
  gg =ggplot(data, aes_string("rel_mtime", value)) 
  gg = gg + geom_vline(xintercept=0, linetype=2)
  gg = gg + geom_point(alpha=I(1/50), color=1)
  
  if (add_smooth) {
    k = 201
    sm = data %>% group_by_(.dots=group_dots) %>% do({ 
      d = . 
      d = d %>% group_by_(.dots=c("fname", "rel_mtime", group_dots)) %>% dplyr::summarize_(.dots=setNames(list(lazyeval::interp(quote(mean(x)), 
                                                                                              x=as.name(value))),
                                                                                      "song_mean"))
      tmp = as.data.frame(d)
      #tmp = tmp[order(tmp[,"rel_mtime"]),]
      print(length(tmp[,"song_mean"]))
      if (length(tmp[,"song_mean"]) < k)
        return(data.frame(d, smoothed=NA))
      res = rollmean(tmp[,"song_mean"], k=k, na.pad=T)
      return(data.frame(d, smoothed=res))
    })
    print(summary(sm$smoothed))
    gg = gg + geom_point(data=sm, aes(x=rel_mtime, y=smoothed), color="blue", size=1)
  }
  
  #if (!is.null(grouping_factor)) {
  #  grouping_obj = as.symbol(grouping_factor)
  if (is.null(color_factor)) {
  
    gg = gg + geom_pointrange(data=m, aes(x=as.numeric(rel_date), y=mean, ymin=(mean-sd), ymax=(mean+sd))) 
    gg = gg + facet_wrap(as.formula(paste("~", grouping_factor, sep="")))
  } else {
    gg = gg + geom_pointrange(data=m, aes_string(x="rel_date", y="mean", ymin="lower", ymax="upper", color=color_factor)) 
    gg = gg + facet_wrap(as.formula(paste("~", grouping_factor, sep="")))
  }
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
    filter_(.dots=lazyeval::interp(quote(rel_mtime<ref), ref=reference_date)) %>% 
    group_by_(.dots=grouping_factor) %>% 
    summarize(baseline_stat = mean(value_obj, na.rm=T))
  gg = gg + geom_hline(data=baseline_data, aes(yintercept=baseline_stat), linetype=2)
  #}
  
  gg = gg + theme(panel.grid.major.x=element_line(linetype=2, color="grey"))
  gg
}

generate_by_hour_plot = function(data, value, hour_column="mtime_hour", color_factor=NULL) {
  sum_dots = list(lazyeval::interp(~mean(v, na.rm=T), v = as.name(value)),
                  lazyeval::interp(~sd(v, na.rm=T), v = as.name(value)))
  
  dots = NULL
  if (!is.null(color_factor)) {
    dots = c(hour_column, "labels", color_factor)
  } else {
    dots = c(hour_column, "labels")
  }
  m = data %>% group_by_(.dots=dots) %>% dplyr::summarize_(.dots=setNames(sum_dots, c("mean", "sd")))
  m = m %>% mutate_(.dots=setNames(list(lazyeval::interp(quote(h + lubridate::minutes(30)), 
                                                         h=as.name(hour_column))),
                                   paste(hour_column, "p30", sep="_")))
  m = m %>% mutate(upper=mean+sd,
                   lower=mean-sd)
  return(m)
}
plot_within_day = function(data, value, x_value="mtime", reference_datetime=NULL, group_factor="hour", color_factor=NULL, to_exclude=c(" ", "-" , "n")) {
  data = data %>% dplyr::filter(!(labels %in% to_exclude))
  gg = ggplot(data, aes_string(x_value, value))
  gg = gg + geom_point(alpha=I(1/10))
  gg = gg + facet_wrap(~labels, scales="free_y")
  if (!is.null(reference_datetime)) {
    if (length(reference_datetime)>1) {
      for (i in 1:length(reference_datetime)) {
        gg = gg +  geom_vline(xintercept=as.numeric(reference_datetime[i]), linetype=2)
      }
    } else {
      gg = gg + geom_vline(xintercept=as.numeric(reference_datetime), linetype=2)
    }
  }

  
  m = NULL
  if (group_factor=="hour") {
    data = data %>% mutate(mtime_hour = as.POSIXct(floor_date(mtime, unit="hour")))
    m = generate_by_hour_plot(data, value, color_factor=color_factor)
  }
  if (!is.null(color_factor)) {
    gg = gg + geom_pointrange(data=m, aes_string(x="mtime_hour_p30", 
                                                 y="mean", 
                                                 ymin="lower", 
                                                 ymax="upper", 
                                                 color=color_factor))
    gg = gg + scale_color_manual(values=c("darkblue", "darkred"))
  } else {
    gg = gg + geom_pointrange(data=m, aes(x=mtime_hour_p30, y=mean, ymin=lower, ymax=upper))
  }
  
 # min_time = min(data$mtime)

  #max_time = max(data$mtime)
  #interval = 3600 * 2
  #print(min_time)
  gg = gg + scale_x_datetime(breaks=date_breaks("6 hours"))
  gg = gg + theme(axis.text.x = element_text(angle=45, hjust=1), panel.grid.major.x=element_line(linetype=3, color="grey"))
  gg
}
plot_within_day_split_labels = function(data, value, reference_datetime, group_factor, color_factor, to_exclude=c(" ", "-", "n"), plot_dir=NULL) {
  labels = unique(data$labels)
  labels = labels[!(labels %in% to_exclude)]
  for (l in labels) {
    if (!is.null(plot_dir))
      cairo_pdf(paste(plot_dir, sprintf("%s_%s.pdf", value, l), sep="/"), width=12, height=8)
    
    gg = plot_within_day(data %>% filter(labels==l), value,
                         reference_datetime,
                         group_factor=group_factor, 
                         color_factor=color_factor,
                         to_exclude=to_exclude)
    print(gg)
    if (!is.null(plot_dir))
      dev.off()
  }
}

plot_within_day_grid = function(data, value, reference_datetime=NULL, group_factor="hour", color_factor=NULL, to_exclude=c(" ", "-" , "n")) {
  data = data %>% dplyr::filter(!(labels %in% to_exclude))
  gg = ggplot(data, aes_string("time", value))
  gg = gg + geom_point(alpha=I(1/10))
  gg = gg + facet_grid(labels~date, scales="free")
  if (!is.null(reference_datetime)) {
    if (length(reference_datetime)>1) {
      for (i in 1:length(reference_datetime)) {
        gg = gg +  geom_vline(xintercept=as.numeric(reference_datetime[i]), linetype=2)
      }
    } else {
      gg = gg + geom_vline(xintercept=as.numeric(reference_datetime), linetype=2)
    }
  }
  
  m = NULL
  if (group_factor=="hour") {
    data = data %>% mutate(mtime_hour = as.POSIXct(floor_date(time, unit="hour")))
    m = generate_by_hour_plot(data, value, color_factor=color_factor)
  }
  if (!is.null(color_factor)) {
    gg = gg + geom_pointrange(data=m, aes_string(x="mtime_hour_p30", 
                                                 y="mean", 
                                                 ymin="lower", 
                                                 ymax="upper", 
                                                 color=color_factor))
    gg = gg + scale_color_manual(values=c("darkblue", "darkred"))
  } else {
    gg = gg + geom_pointrange(data=m, aes(x=mtime_hour_p30, y=mean, ymin=lower, ymax=upper))
  }
  
  gg = gg + scale_x_datetime(breaks=date_breaks("6 hours"))
  gg = gg + theme(axis.text.x = element_text(angle=45, hjust=1), panel.grid.major.x=element_line(linetype=3, color="grey"))
  gg
}
