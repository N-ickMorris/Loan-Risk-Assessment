# -----------------------------------------------------------------------------------
# ---- Set Up -----------------------------------------------------------------------
# -----------------------------------------------------------------------------------

# set the path of where the input files are
mywd = "C:/Users/Nick Morris/Downloads/Data-Science/Eilif"

# open up a graphics window
windows()

# -----------------------------------------------------------------------------------
# ---- Packages ---------------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# data handling
require(data.table)
require(tm)
require(gtools)
require(stringdist)
require(zoo)
require(missRanger)
require(TTR)

# plotting
require(ggplot2)
require(gridExtra)
require(GGally)
require(scales)
require(scatterplot3d)
require(gridGraphics)
require(corrplot)
require(VIM)

# modeling
require(fitdistrplus)
require(bestNormalize)
require(fpc)
require(caret)
require(ranger)
require(cluster)
require(car)
require(nortest)
require(neuralnet)
require(h2o)
require(MLmetrics)

# parallel computing
require(foreach)
require(parallel)
require(doSNOW)

# choose how many workers to use when parallel computing
workers = max(1, floor((2/3) * detectCores()))
workers = max(1, detectCores() - 2)

}

# -----------------------------------------------------------------------------------
# ---- Functions --------------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# ---- prints the data types of each column in a data frame -------------------------

types = function(dat)
{
  # make dat into a data.frame
  dat = data.frame(dat)
  
  # get the column names
  column = colnames(dat)
  
  # get the class of the columns
  data.type = sapply(1:ncol(dat), function(i) class(dat[,i]))
  
  # compute the number of levels for each column
  levels = sapply(1:ncol(dat), function(i) ifelse(data.type[i] == "factor", length(levels(droplevels(dat[,i]))), 0))
  
  return(data.frame(column, data.type, levels))
}

# ---- converts all columns to a character data type --------------------------------

tochar = function(dat)
{
  # make dat into a data.frame
  dat = data.frame(dat)
  
  # get the column names
  column = colnames(dat)
  
  # get the values in the columns and convert them to character data types
  values = lapply(1:ncol(dat), function(i) as.character(dat[,i]))
  
  # combine the values back into a data.frame
  dat = data.frame(do.call("cbind", values), stringsAsFactors = FALSE)
  
  # give dat its column names
  colnames(dat) = column
  
  return(dat)
}

# ---- a qualitative color scheme ---------------------------------------------------

mycolors = function(n)
{
  require(grDevices)
  return(colorRampPalette(c("#e41a1c", "#0099ff", "#4daf4a", "#984ea3", "#ff7f00", "#ff96ca", "#a65628"))(n))
}

# ---- emulates the default ggplot2 color scheme ------------------------------------

ggcolor = function(n, alpha = 1)
{
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# ---- plots various cluster plots --------------------------------------------------

plot.clusters = function(dat, cluster.column.name = NULL, distance.matrix = NULL, DC.title = "Discriminant Coordinate Cluster Plot", pairs.title = "Cluster Pairs Plot", silhouette.title = "Silhouette Width", font.size = 20, pairs.plot.font.size = 12, rotate = 0, cor.size = 3)
{
  # load packages we need
  require(data.table)
  require(ggplot2)
  require(GGally)
  require(scatterplot3d)
  require(gridGraphics)
  require(grid)
  require(scales)
  require(cluster)
  require(fpc)
  
  # this function emulates the default ggplot2 color scheme
  ggcolor = function(n)
  {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  # error check
  if(is.null(cluster.column.name))
  {
    print("you must specify a value for the parameter: cluster.column.name")
    
  } else
  {
    # ---- computing discriminant coordiantes -------------------------------
    
    # make dat into a data table
    dat = data.table(dat)
    
    # extract the cluster column from dat
    clusters = as.numeric(unname(unlist(dat[, cluster.column.name, with = FALSE])))
    
    # remove the cluster column from dat
    dat = dat[, !cluster.column.name, with = FALSE]
    
    # compute the number of columns in dat
    numcol = ncol(dat)
    
    # compute the discriminant coordinates, and extract the first 3
    dat.dc = data.table(discrproj(x = dat, 
                                  clvecd = clusters,
                                  method = "dc")$proj[,1:min(c(3, numcol))])
    
    # rename the columns appropriately
    setnames(dat.dc, paste0("DC", 1:min(c(3, numcol))))
    
    # give dat.dc the cluster column, and make it into a factor for plotting purposes
    dat.dc[, Cluster := factor(clusters, levels = sort(unique(clusters)))]
    
    # create an output list to store results
    output = list()
    
    # ---- plotting 2D discriminant coordiantes -----------------------------
    
    if(numcol >= 2)
    {
      # create a 2D cluster plot across the first 2 discriminant coordinates
      plot.2D = ggplot(dat.dc, aes(x = DC1, y = DC2, fill = Cluster)) +
        stat_density_2d(geom = "polygon", color = NA, alpha = 1/3) +
        theme_bw(font.size) +
        ggtitle(DC.title) +
        theme(legend.position = "top", 
              legend.key.size = unit(.25, "in"), 
              plot.title = element_text(hjust = 0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        guides(fill = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1))
      
      # store this plot in output
      output$plot.2D = plot.2D
    }
    
    # ---- plotting 3D discriminant coordiantes -----------------------------
    
    if(numcol >= 3)
    {
      # create a table indicating which cluster should get which ggplot color
      color.dat.dc = data.table(Cluster = levels(dat.dc$Cluster),
                                Color = ggcolor(length(levels(dat.dc$Cluster))))
      
      # set Cluster as the key column in color.dat.dc and dat.dc
      setkey(dat.dc, Cluster)
      setkey(color.dat.dc, Cluster)
      
      # join color.dat.dc onto dat.dc
      dat.dc = color.dat.dc[dat.dc]
      
      # convert Cluster back into a factor data type
      dat.dc[, Cluster := factor(Cluster, levels = sort(unique(Cluster)))]
      
      # here is my default font size for base R plotting
      font.default = 20
      
      # compute the desired adjustment according to the specified value of font.size
      font.adjust = font.size / font.default
      
      # adjust the font of the title, axis titles, axis labels, and legend
      font.title = 2 * font.adjust
      font.axis.title = 1.125 * font.adjust
      font.axis.label = 1.125 * font.adjust
      font.legend = 1.75 * font.adjust
      
      # here are my 4 default angles for viewing a 3D scatterplot
      angles = c(45, 135, 225, 315)
      
      # apply the specified rotation
      angles = angles + rotate
      
      # set up 4 plot ID numbers so each plot angle has a position in the plot window
      plot.id = 2:5
      
      # set up a legend ID number so the legend has a postion across the top of the plot window
      legend.id = c(1, 1)
      
      # set up a matrix that defines the plot layout
      plot.layout = matrix(c(legend.id, plot.id), nrow = 3, ncol = 2, byrow = TRUE)
      
      # create a new plot window
      windows()
      plot.new()
      
      # define plot margins
      par(mar = c(0, 0, 3.5, 0))
      
      # apply the layout to the plot window
      layout(mat = plot.layout, heights = c(1, 1.5, 1.5))
      
      # produce a dummy plot as a place holder for the legend and title
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
      
      # produce the title
      title(main = DC.title, cex.main = font.title)
      
      # produce the legend
      legend("top", inset = 0, bty = "n", cex = font.legend, horiz = TRUE,
             title = "Clusters", legend = levels(dat.dc$Cluster), fill = ggcolor(length(levels(dat.dc$Cluster))))
      
      # create 3D cluster plots across the first 3 discriminant coordinates
      for(angle in angles)
      {
        # build the scatterplot
        scatterplot3d(x = dat.dc$DC1, y = dat.dc$DC2, z = dat.dc$DC3, color = alpha(dat.dc$Color, 1/3), 
                      xlab = "DC1", ylab = "DC2", zlab = "DC3", mar = c(3, 3, 0, 3), cex = 1.5,
                      pch = 16, cex.lab = font.axis.title, cex.axis = font.axis.label, angle = angle)
      }
      
      # save this plot as a grid object
      grid.echo()
      plot.3D = grid.grab()
      
      # add this plot to our output list
      output$plot.3D = plot.3D
      
      # close the graphics window
      graphics.off()
    }
    
    # ---- plotting clusters across variable pairs ---------------------------
    
    # give dat the cluster column, and make it into a factor for plotting purposes
    dat[, Cluster := factor(clusters, levels = sort(unique(clusters)))]
    
    # plot the clusters across all variable pairs
    plot.pairs = ggpairs(dat,
                         mapping = aes(color = Cluster, fill = Cluster),
                         columns = which(names(dat) != "Cluster"),
                         lower = list(continuous = wrap(ggally_points, size = 1.5, alpha = 1/3)), 
                         upper = list(continuous = wrap(ggally_cor, size = cor.size)),
                         diag = list(continuous = wrap(ggally_densityDiag, alpha = 1/3)),
                         title = pairs.title,
                         legend = grab_legend(plot.2D + theme_classic(base_size = pairs.plot.font.size) + theme(legend.position = "top"))) + 
      theme_classic(base_size = pairs.plot.font.size) +
      theme(legend.position = "top", plot.title = element_text(hjust = 0.5))
    
    # remove the cluster column from dat
    dat[, Cluster := NULL]
    
    # add this plot to our output list
    output$plot.pairs = plot.pairs
    
    # ---- plotting silhouette widths ----------------------------------------
    
    # if the user gave a distance matrix then lets compute silhouette widths
    if(!is.null(distance.matrix))
    {
      # compute the silhouette widths
      dat.sil = silhouette(x = clusters,
                           dist = distance.matrix)
      
      # compute the summary of dat.sil
      dat.sil.sum = summary(dat.sil)
      
      # extract the avg widths from dat.sil.sum
      dat.sil.avg = data.table(Cluster = as.numeric(names(dat.sil.sum$clus.avg.widths)),
                               Average_Width = round(unname(dat.sil.sum$clus.avg.widths), 2))
      
      # order dat.sil.avg by Cluster
      dat.sil.avg = dat.sil.avg[order(Cluster)]
      
      # extract the cluster sizes from dat.sil.sum
      dat.sil.size = data.table(Cluster = as.numeric(names(dat.sil.sum$clus.sizes)),
                                Size = as.numeric(unname(dat.sil.sum$clus.sizes)))
      
      # order dat.sil.size by Cluster
      dat.sil.size = dat.sil.size[order(Cluster)]
      
      # combine dat.sil.avg and dat.sil.size into a table that will go on our plot
      dat.sil.tab = cbind(dat.sil.avg, dat.sil.size[,!"Cluster"])
      
      # convert dat.sil into a data table
      dat.sil = data.table(dat.sil[1:nrow(dat.sil), 1:ncol(dat.sil)])
      
      # sort dat.sil by cluster and sil_width
      dat.sil = dat.sil[order(cluster, -sil_width)]
      
      # give dat.sil an ID column for plotting purposes
      dat.sil[, ID := 1:nrow(dat.sil)]
      
      # convert cluster to a factor for plotting purposes
      dat.sil[, cluster := factor(cluster, levels = sort(unique(cluster)))]
      
      # aggregate sil_width by cluster in dat.sil to determine where to place dat.sil.tab in the plot 
      dat.agg = dat.sil[, .(sil_width.min = min(sil_width),
                            sil_width.max = max(sil_width)),
                        by = cluster]
      
      # build the four corners of the dat.sil.tab to place it in the plot
      # find the cluster with the smallest peak and set the peak's sil_width as the ymin
      ymin = as.numeric(min(dat.agg$sil_width.max))
      
      # find the sil_width of the max peak and set it as the ymax
      ymax = as.numeric(max(dat.sil$sil_width))
      
      # extract the cluster with the smallest peak from dat.agg
      small.peak = dat.agg[which.min(sil_width.max), cluster]
      
      # find the first ID number in dat.sil for the cluster with the smallest peak, and set that as xmin
      xmin = min(dat.sil[cluster == small.peak, ID])
      
      # find the last ID number in dat.sil for the cluster with the smallest peak, and set that as xmax
      xmax = max(dat.sil[cluster == small.peak, ID])
      
      # plot the silhouette width and add the dat.sil.tab to it
      plot.sil.width = ggplot(dat.sil, aes(x = ID, y = sil_width, fill = cluster, color = cluster)) + 
        geom_bar(stat = "identity", position = "dodge") +
        # annotation_custom(tableGrob(as.matrix(dat.sil.tab), rows = NULL, 
        #                             theme = ttheme_default(base_size = font.size,
        #                                                    colhead = list(fg_params = list(col = "black"), bg_params = list(fill = "lightgray", col = "black")),
        #                                                    core = list(fg_params = list(hjust = 0.5), bg_params = list(fill = c("white"), col = "black")))),
        #                   xmin = xmin, 
        #                   xmax = xmax, 
        #                   ymin = ymin, 
        #                   ymax = ymax) +
        ggtitle(silhouette.title) +
        labs(x = "Observation", y = "Silhouette Width", fill = "Cluster", color = "Cluster") +
        theme_bw(base_size = font.size) +
        theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) + 
        guides(fill = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1),
               color = guide_legend(nrow = 1))
      
      # add this plot to our output list
      output$plot.sil.width = plot.sil.width
    }
    
    # add the DC data to the output list
    output$dat.dc = dat.dc
    
    # helpful message
    print("use grid.draw() to see the 3D cluster plot, for example: grid.draw(my.plot.clusters$plot.3D)")
    
    return(output)
  }
}

# ---- prints out a dat file object in ampl syntax ----------------------------------

ampl = function(dat, object = "param", name = "c")
{
  # converts all columns to a character data type 
  tochar = function(dat)
  {
    # make dat into a data.frame
    dat = data.frame(dat)
    
    # get the column names
    column = colnames(dat)
    
    # get the values in the columns and convert them to character data types
    values = lapply(1:ncol(dat), function(i) as.character(dat[,i]))
    
    # combine the values back into a data.frame
    dat = data.frame(do.call("cbind", values), stringsAsFactors = FALSE)
    
    # give dat its column names
    colnames(dat) = column
    
    return(dat)
  }
  
  # make sure the data is a data frame object
  dat = tochar(dat)
  
  # every parameter/set object in an ampl dat file must end with a semicolon
  # so set up 1 semicolon to give to dat
  semicolon = c(";", rep(" ", ncol(dat) - 1))
  
  # add this semicolon as the last row of the data frame
  result = data.frame(rbind(dat, semicolon))
  
  # every parameter/set object in an ample dat file must begin with the name of the object and what it equals
  # for example: param c := 
  # so set up a header to give to dat
  header = c(paste(object, name, ":="), rep(" ", ncol(dat) - 1))
  
  # update the column names of dat to be the header we created
  colnames(result) = header
  
  # print out the result without any row names
  # print out the result left adjusted
  # print(result, right = FALSE, row.names = FALSE)
  
  return(result)	
}

# ---- compares the quantiles of emprical data against the quantiles of any statistical distribution 

ggqq = function(x, distribution = "norm", ..., conf = 0.95, probs = c(0.25, 0.75), alpha = 0.33, basefont = 20, main = "", xlab = "\nTheoretical Quantiles", ylab = "Empirical Quantiles\n")
{
  require(ggplot2)
  
  # compute the sample quantiles and theoretical quantiles
  q.function = eval(parse(text = paste0("q", distribution)))
  d.function = eval(parse(text = paste0("d", distribution)))
  x = na.omit(x)
  ord = order(x)
  n = length(x)
  P = ppoints(length(x))
  df = data.frame(ord.x = x[ord], z = q.function(P, ...))
  
  # compute the quantile line
  Q.x = quantile(df$ord.x, c(probs[1], probs[2]))
  Q.z = q.function(c(probs[1], probs[2]), ...)
  b = diff(Q.x) / diff(Q.z)
  coef = c(Q.x[1] - (b * Q.z[1]), b)
  
  # compute the confidence interval band
  zz = qnorm(1 - (1 - conf) / 2)
  SE = (coef[2] / d.function(df$z, ...)) * sqrt(P * (1 - P) / n)
  fit.value = coef[1] + (coef[2] * df$z)
  df$upper = fit.value + (zz * SE)
  df$lower = fit.value - (zz * SE)
  
  # plot the qqplot
  p = ggplot(df, aes(x = z, y = ord.x)) + 
    geom_point(color = "blue", alpha = alpha) +
    geom_abline(intercept = coef[1], slope = coef[2], size = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
    coord_cartesian(ylim = c(min(df$ord.x), max(df$ord.x))) + 
    labs(x = xlab, y = ylab) +
    theme_bw(base_size = basefont) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # conditional additions
  if(main != "")(p = p + ggtitle(main))
  
  return(p)
}

# ---- plots 6 residual plots -------------------------------------------------------

residplots6 = function(actual, fit, binwidth = NULL, from = NULL, to = NULL, by = NULL, histlabel.y = -10, n = NULL, basefont = 20)
{
  require(ggplot2)
  
  residual = actual - fit 
  DF = data.frame("actual" = actual, "fit" = fit, "residual" = residual)
  
  rvfPlot = ggplot(DF, aes(x = fit, y = residual)) + 
    geom_point(na.rm = TRUE) +
    stat_smooth(method = "loess", se = FALSE, na.rm = TRUE, color = "blue") +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    xlab("Fitted values") +
    ylab("Residuals") +
    ggtitle("Residual vs Fitted Plot") + 
    theme_bw(base_size = basefont) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ggqq = function(x, distribution = "norm", ..., conf = 0.95, probs = c(0.25, 0.75), note = TRUE, alpha = 0.33, main = "", xlab = "\nTheoretical Quantiles", ylab = "Empirical Quantiles\n")
  {
    # compute the sample quantiles and theoretical quantiles
    q.function = eval(parse(text = paste0("q", distribution)))
    d.function = eval(parse(text = paste0("d", distribution)))
    x = na.omit(x)
    ord = order(x)
    n = length(x)
    P = ppoints(length(x))
    DF = data.frame(ord.x = x[ord], z = q.function(P, ...))
    
    # compute the quantile line
    Q.x = quantile(DF$ord.x, c(probs[1], probs[2]))
    Q.z = q.function(c(probs[1], probs[2]), ...)
    b = diff(Q.x) / diff(Q.z)
    coef = c(Q.x[1] - (b * Q.z[1]), b)
    
    # compute the confidence interval band
    zz = qnorm(1 - (1 - conf) / 2)
    SE = (coef[2] / d.function(DF$z, ...)) * sqrt(P * (1 - P) / n)
    fit.value = coef[1] + (coef[2] * DF$z)
    DF$upper = fit.value + (zz * SE)
    DF$lower = fit.value - (zz * SE)
    
    # plot the qqplot
    p = ggplot(DF, aes(x = z, y = ord.x)) + 
      geom_point(color = "black", alpha = alpha) +
      geom_abline(intercept = coef[1], slope = coef[2], size = 1, color = "blue") +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15) +
      coord_cartesian(ylim = c(min(DF$ord.x), max(DF$ord.x))) + 
      labs(x = xlab, y = ylab) +
      theme_bw(base_size = basefont) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    # conditional additions
    if(main != "")(p = p + ggtitle(main))
    
    return(p)
  }
  
  qqPlot = ggqq(residual, 
                alpha = 1,				  
                main = "Normal Q-Q Plot", 
                xlab = "Theoretical Quantiles", 
                ylab = "Residuals")
  
  rvtPlot = ggplot(data.frame("x" = 1:length(DF$residual), "y" = DF$residual), aes(x = x, y = y)) + 
    geom_line(na.rm = TRUE) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    xlab("Obs. Number") +
    ylab("Residuals") +
    ggtitle("Residual Time Series") + 
    theme_bw(base_size = basefont) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  variogramDF = function(x)
  {
    n = length(x) - 2
    
    num = sapply(1:n, function(k)
      sapply(1:(length(x) - k), function(i)
        x[i + k] - x[i]))
    
    num = sapply(1:length(num), function(j)
      var(num[[j]]))
    
    den = var(sapply(1:(length(x) - 1), function(i)
      x[i + 1] - x[i]))
    
    val = num / den
    
    DF = data.frame("Lag" = 1:n, "Variogram" = val)
    
    return(DF)
  }
  
  DFv = variogramDF(x = DF$residual)
  
  varioPlot = ggplot(DFv, aes(x = Lag, y = Variogram)) + 
    geom_point() +
    geom_line(color = "blue") +
    xlab("Lag") +
    ylab("Variogram") +
    ggtitle("Variogram of Residuals") + 
    theme_bw(base_size = basefont) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  test = t.test(DF$residual)
  
  CI = data.frame("x" = test$estimate, 
                  "LCB" = test$conf.int[1], 
                  "UCB" = test$conf.int[2], 
                  row.names = 1)
  
  histPlot = ggplot(DF, aes(x = residual)) +
    geom_histogram(color = "white", fill = "black", binwidth = binwidth) +
    geom_segment(data = CI, aes(x = LCB, xend = LCB, y = 0, yend = Inf), color = "blue") +
    geom_segment(data = CI, aes(x = UCB, xend = UCB, y = 0, yend = Inf), color = "blue") +
    annotate("text", x = CI$x, y = histlabel.y, 
             label = "T-Test C.I.", size = 5, 
             color = "blue", fontface = 2) + 
    ggtitle("Residual Histogram") +
    labs(x = "Residuals", y = "Frequency") +
    theme_bw(base_size = basefont) +
    theme(legend.key.size = unit(.25, "in"),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  if(class(from) != "NULL" & class(to) != "NULL" & class(by) != "NULL") (histPlot = histPlot + scale_x_continuous(breaks = seq(from = from, to = to, by = by)))
  
  ggacf = function(x, n = NULL, conf.level = 0.95, main = "ACF Plot", xlab = "Lag", ylab = "Autocorrelation", basefont = 20) 
  {
    if(class(n) == "NULL")
    {
      n = length(x) - 2
    }
    
    ciline = qnorm((1 - conf.level) / 2) / sqrt(length(x))
    bacf = acf(x, lag.max = n, plot = FALSE)
    bacfdf = with(bacf, data.frame(lag, acf))
    bacfdf = bacfdf[-1,]
    
    p = ggplot(bacfdf, aes(x = lag, y = acf)) + 
      geom_bar(stat = "identity", position = "dodge", fill = "black") +
      geom_hline(yintercept = -ciline, color = "blue", size = 1) +
      geom_hline(yintercept = ciline, color = "blue", size = 1) +
      geom_hline(yintercept = 0, color = "red", size = 1) +
      labs(x = xlab, y = ylab) +
      ggtitle(main) +
      theme_bw(base_size = basefont) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    return(p)
  }
  
  acfPlot = ggacf(x = DF$residual, main = "ACF Plot of Residuals", basefont = basefont, n = n)
  
  return(list("rvfPlot" = rvfPlot, 
              "qqPlot" = qqPlot, 
              "rvtPlot" = rvtPlot, 
              "varioPlot" = varioPlot, 
              "histPlot" = histPlot, 
              "acfPlot" = acfPlot))
}

# ---- plots 4 residual plots -------------------------------------------------------

residplots4 = function(actual, fit, binwidth = NULL, from = NULL, to = NULL, by = NULL, histlabel.y = -10, basefont = 20)
{
  require(ggplot2)
  
  residual = actual - fit 
  DF = data.frame("actual" = actual, "fit" = fit, "residual" = residual)
  
  rvfPlot = ggplot(DF, aes(x = fit, y = residual)) + 
    geom_point(na.rm = TRUE) +
    stat_smooth(method = "loess", se = FALSE, na.rm = TRUE, color = "blue") +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    xlab("Fitted values") +
    ylab("Residuals") +
    ggtitle("Residual vs Fitted Plot") + 
    theme_bw(base_size = basefont) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ggqq = function(x, distribution = "norm", ..., conf = 0.95, probs = c(0.25, 0.75), note = TRUE, alpha = 0.33, main = "", xlab = "\nTheoretical Quantiles", ylab = "Empirical Quantiles\n")
  {
    # compute the sample quantiles and theoretical quantiles
    q.function = eval(parse(text = paste0("q", distribution)))
    d.function = eval(parse(text = paste0("d", distribution)))
    x = na.omit(x)
    ord = order(x)
    n = length(x)
    P = ppoints(length(x))
    DF = data.frame(ord.x = x[ord], z = q.function(P, ...))
    
    # compute the quantile line
    Q.x = quantile(DF$ord.x, c(probs[1], probs[2]))
    Q.z = q.function(c(probs[1], probs[2]), ...)
    b = diff(Q.x) / diff(Q.z)
    coef = c(Q.x[1] - (b * Q.z[1]), b)
    
    # compute the confidence interval band
    zz = qnorm(1 - (1 - conf) / 2)
    SE = (coef[2] / d.function(DF$z, ...)) * sqrt(P * (1 - P) / n)
    fit.value = coef[1] + (coef[2] * DF$z)
    DF$upper = fit.value + (zz * SE)
    DF$lower = fit.value - (zz * SE)
    
    # plot the qqplot
    p = ggplot(DF, aes(x = z, y = ord.x)) + 
      geom_point(color = "black", alpha = alpha) +
      geom_abline(intercept = coef[1], slope = coef[2], size = 1, color = "blue") +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15) +
      coord_cartesian(ylim = c(min(DF$ord.x), max(DF$ord.x))) + 
      labs(x = xlab, y = ylab) +
      theme_bw(base_size = basefont) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    # conditional additions
    if(main != "")(p = p + ggtitle(main))
    
    return(p)
  }
  
  qqPlot = ggqq(residual, 
                alpha = 1,				  
                main = "Normal Q-Q Plot", 
                xlab = "Theoretical Quantiles", 
                ylab = "Residuals")
  
  rvtPlot = ggplot(data.frame("x" = 1:length(DF$residual), "y" = DF$residual), aes(x = x, y = y)) + 
    geom_line(na.rm = TRUE) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    xlab("Obs. Number") +
    ylab("Residuals") +
    ggtitle("Residual Time Series") + 
    theme_bw(base_size = basefont) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  test = t.test(DF$residual)
  
  CI = data.frame("x" = test$estimate, 
                  "LCB" = test$conf.int[1], 
                  "UCB" = test$conf.int[2], 
                  row.names = 1)
  
  histPlot = ggplot(DF, aes(x = residual)) +
    geom_histogram(color = "white", fill = "black", binwidth = binwidth) +
    geom_segment(data = CI, aes(x = LCB, xend = LCB, y = 0, yend = Inf), color = "blue") +
    geom_segment(data = CI, aes(x = UCB, xend = UCB, y = 0, yend = Inf), color = "blue") +
    annotate("text", x = CI$x, y = histlabel.y, 
             label = "T-Test C.I.", size = 5, 
             color = "blue", fontface = 2) + 
    ggtitle("Residual Histogram") +
    labs(x = "Residuals", y = "Frequency") +
    theme_bw(base_size = basefont) +
    theme(legend.key.size = unit(.25, "in"),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  if(class(from) != "NULL" & class(to) != "NULL" & class(by) != "NULL") (histPlot = histPlot + scale_x_continuous(breaks = seq(from = from, to = to, by = by)))
  
  return(list("rvfPlot" = rvfPlot, 
              "qqPlot" = qqPlot, 
              "rvtPlot" = rvtPlot,  
              "histPlot" = histPlot))
}

# ---- builds a square confusion matrix ---------------------------------------------

confusion = function(ytrue, ypred)
{
  require(gtools)
  
  # make predicted and actual vectors into factors, if they aren't already
  if(class(ytrue) != "factor") ytrue = factor(ytrue)
  if(class(ypred) != "factor") ypred = factor(ypred)
  
  # combine their levels into one unique set of levels
  common.levels = mixedsort(unique(c(levels(ytrue), levels(ypred))))
  
  # give each vector the same levels
  ytrue = factor(ytrue, levels = common.levels)
  ypred = factor(ypred, levels = common.levels)
  
  # build the confusion matrix
  output = table("Actual" = ytrue, "Predicted" = ypred)
  
  # return the confusion matrix
  return(output)
}

}

# -----------------------------------------------------------------------------------
# ---- Prepare the Data -------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# ---- clean up the data ----
  
# set the work directory
setwd(mywd)

# import the data
dat = data.table(read.csv("LCforNick.csv", stringsAsFactors = FALSE))
states = data.table(read.csv("USA-State-Data.csv", stringsAsFactors = FALSE))

# check out the data types of each column
types(dat)

# check out the data
dat

# remove the commas from potential_fv, ï..id, annual_inc, revol_bal, total_rev_hi_lim, avg_cur_bal
# relabel ï..id as id
dat[, potential_fv := as.numeric(gsub(",", "", potential_fv))]
dat[, id := as.numeric(gsub(",", "", ï..id))]
dat[, ï..id := NULL]
dat[, annual_inc := as.numeric(gsub(",", "", annual_inc))]
dat[, revol_bal := as.numeric(gsub(",", "", revol_bal))]
dat[, total_rev_hi_lim := as.numeric(gsub(",", "", total_rev_hi_lim))]
dat[, avg_cur_bal := as.numeric(gsub(",", "", avg_cur_bal))]

# order dat by id to maintain a consistent row order
dat = dat[order(id)]

# remove true_ann_return as this is an alternative response variable that we are not interested in
dat[, true_ann_return := NULL]

# identify the response variables of interest
responses = c("pmt_status", "total_return")

# remove any NA's in dat
dat = na.omit(dat)

# make pmt_status into a factor
dat[, pmt_status := factor(pmt_status, levels = c(0, 1))]

# ---- create a multi-class variable for total_return ----

# plot a histogram of total_return to create a categorical version of it
histogram(dat$total_return)

# adjust bin numbers as you like
num.bins = 7
histogram(dat$total_return, nint = num.bins)

# create a table indicating the frequency of observations in each bin
bin.table = table(cut(dat$total_return, breaks = num.bins))

# convert frequencies into percentages
bin.table = bin.table / sum(bin.table)

# look at the bins
bin.table

# lets create a categorical version of total_return
dat[, total_return_class := as.character(cut(dat$total_return, breaks = num.bins))]

# lets combine (-1,-0.749] (-0.749,-0.498] (-0.498,-0.247] and (-0.247,0.00412]
dat[, total_return_class := ifelse(total_return_class %in% c("(-1,-0.749]", "(-0.749,-0.498]", "(-0.498,-0.247]", "(-0.247,0.00412]"), 
                                   "(-1,0.00412]", 
                                   total_return_class)]

# lets split up (0.00412,0.255]
dat[, total_return_class := ifelse(total_return > 0.00412 & total_return <= 0.10, "(0.00412,0.10]",
                                   ifelse(total_return > 0.10 & total_return <= 0.15, "(0.10,0.15]",
                                          ifelse(total_return > 0.15 & total_return <= 0.20, "(0.15,0.20]",
                                                 ifelse(total_return > 0.20 & total_return <= 0.506, "(0.20,0.506]",
                                                        ifelse(total_return_class == "(0.255,0.506]", "(0.20,0.506]",
                                                               total_return_class)))))]

# remove all observation in (0.506,0.759] as this is less than 1% of the data
dat = dat[total_return_class != "(0.506,0.759]"]

# create a table indicating the frequency of observations in each bin
bin.table = table(dat$total_return_class)

# convert frequencies into percentages
bin.table = bin.table / sum(bin.table)

# look at the bins
bin.table * 100

# set up the proper order of the classes
class.order = c("(-1,0.00412]", "(0.00412,0.10]", "(0.10,0.15]", "(0.15,0.20]", "(0.20,0.506]")

# convert total_return_class to a factor
dat[, total_return_class := factor(total_return_class, levels = class.order)]

# determine the best way to normalize total_return
# best.norm = bestNormalize(dat$total_return)
best.norm = yeojohnson(dat$total_return)

# add total_return_norm to dat
dat[, total_return_norm := best.norm$x.t]

# update the response variables of interest
responses = c(responses, "total_return_class", "total_return_norm")

# ---- add state data ----

# get the binary state variables from dat
binary = data.table(dat[,c("id", names(dat)[which(grepl("state_", names(dat)))]), with = FALSE])

# convert binary to long format
binary = melt(binary, id.vars = "id", variable.name = "Label")

# only keep rows where value == 1
binary = binary[value == 1]

# update Label to be the capitalized acronyms of each state
binary[, Label := toupper(gsub("state_", "", Label))]

# compute the count of each label in binary
label.frequency = table(binary$Label)

# make label.frequency into a data table
label.frequency = data.table(Label = names(label.frequency),
                             State_Frequency = as.numeric(label.frequency / sum(label.frequency)))

# update states to not include the State name
states[, State := NULL]

# join states onto binary
setkey(binary, Label)
setkey(states, Label)
binary = states[binary]

# join label.frequency onto binary
setkey(binary, Label)
setkey(label.frequency, Label)
binary = label.frequency[binary]

# remove Label and value from binary
binary[, c("Label", "value") := NULL]

# join binary onto dat
setkey(binary, id)
setkey(dat, id)
dat = binary[dat]

# update dat to not contain any of the binary state indicators
dat = dat[, !names(dat)[which(grepl("state_", names(dat)))], with = FALSE]

# order dat by id
dat = dat[order(id)]

# ---- encode binary variables based on their class frequency ----

# get the binary emp variables from dat
binary = data.table(dat[,c("id", names(dat)[which(grepl("emp_", names(dat)))]), with = FALSE])

# convert binary to long format
binary = melt(binary, id.vars = "id", variable.name = "Label")

# only keep rows where value == 1
binary = binary[value == 1]

# check that each id has an assigned binary value
nrow(binary) == nrow(dat)

# compute the count of each label in binary
label.frequency = table(binary$Label)

# make label.frequency into a data table
label.frequency = data.table(Label = names(label.frequency),
                             emp_Frequency = as.numeric(label.frequency / sum(label.frequency)))

# join label.frequency onto binary
setkey(binary, Label)
setkey(label.frequency, Label)
binary = label.frequency[binary]

# remove Label and value from binary
binary[, c("Label", "value") := NULL]

# join binary onto dat
setkey(binary, id)
setkey(dat, id)
dat = binary[dat]

# order dat by id
dat = dat[order(id)]

# create a mid_local_inc variable
dat[, mid_local_inc := as.integer(1 - high_local_inc - low_local_inc)]

# get the binary local_inc variables from dat
binary = data.table(dat[,c("id", names(dat)[which(grepl("local_inc", names(dat)))]), with = FALSE])

# convert binary to long format
binary = melt(binary, id.vars = "id", variable.name = "Label")

# only keep rows where value == 1
binary = binary[value == 1]

# check that each id has an assigned binary value
nrow(binary) == nrow(dat)

# compute the count of each label in binary
label.frequency = table(binary$Label)

# make label.frequency into a data table
label.frequency = data.table(Label = names(label.frequency),
                             local_inc_Frequency = as.numeric(label.frequency / sum(label.frequency)))

# join label.frequency onto binary
setkey(binary, Label)
setkey(label.frequency, Label)
binary = label.frequency[binary]

# remove Label and value from binary
binary[, c("Label", "value") := NULL]

# join binary onto dat
setkey(binary, id)
setkey(dat, id)
dat = binary[dat]

# order dat by id
dat = dat[order(id)]

# get the binary grade_ variables from dat
binary = data.table(dat[,c("id", names(dat)[which(grepl("grade_", names(dat)))]), with = FALSE])

# remove the subgrade_ variables from binary
binary = data.table(binary[,names(binary)[-which(grepl("subgrade_", names(binary)))], with = FALSE])

# create a grade_na variable
binary[, grade_na := as.integer(1 - grade_a - grade_b - grade_c - grade_d - grade_e - grade_f)]

# convert binary to long format
binary = melt(binary, id.vars = "id", variable.name = "Label")

# only keep rows where value == 1
binary = binary[value == 1]

# check that each id has an assigned binary value
nrow(binary) == nrow(dat)

# compute the count of each label in binary
label.frequency = table(binary$Label)

# make label.frequency into a data table
label.frequency = data.table(Label = names(label.frequency),
                             grade_Frequency = as.numeric(label.frequency / sum(label.frequency)))

# join label.frequency onto binary
setkey(binary, Label)
setkey(label.frequency, Label)
binary = label.frequency[binary]

# remove Label and value from binary
binary[, c("Label", "value") := NULL]

# join binary onto dat
setkey(binary, id)
setkey(dat, id)
dat = binary[dat]

# order dat by id
dat = dat[order(id)]

# get the binary subgrade_a variables from dat
binary = data.table(dat[,c("id", names(dat)[which(grepl("subgrade_a", names(dat)))]), with = FALSE])

# create a subgrade_a0 variable
binary[, subgrade_a0 := as.integer(1 - subgrade_a1 - subgrade_a2 - subgrade_a3 - subgrade_a4 - subgrade_a5)]

# convert binary to long format
binary = melt(binary, id.vars = "id", variable.name = "Label")

# only keep rows where value == 1
binary = binary[value == 1]

# check that each id has an assigned binary value
nrow(binary) == nrow(dat)

# compute the count of each label in binary
label.frequency = table(binary$Label)

# make label.frequency into a data table
label.frequency = data.table(Label = names(label.frequency),
                             subgrade_a_Frequency = as.numeric(label.frequency / sum(label.frequency)))

# join label.frequency onto binary
setkey(binary, Label)
setkey(label.frequency, Label)
binary = label.frequency[binary]

# remove Label and value from binary
binary[, c("Label", "value") := NULL]

# join binary onto dat
setkey(binary, id)
setkey(dat, id)
dat = binary[dat]

# order dat by id
dat = dat[order(id)]

# get the binary subgrade_b variables from dat
binary = data.table(dat[,c("id", names(dat)[which(grepl("subgrade_b", names(dat)))]), with = FALSE])

# create a subgrade_b0 variable
binary[, subgrade_b0 := as.integer(1 - subgrade_b1 - subgrade_b2 - subgrade_b3 - subgrade_b4 - subgrade_b5)]

# convert binary to long format
binary = melt(binary, id.vars = "id", variable.name = "Label")

# only keep rows where value == 1
binary = binary[value == 1]

# check that each id has an assigned binary value
nrow(binary) == nrow(dat)

# compute the count of each label in binary
label.frequency = table(binary$Label)

# make label.frequency into a data table
label.frequency = data.table(Label = names(label.frequency),
                             subgrade_b_Frequency = as.numeric(label.frequency / sum(label.frequency)))

# join label.frequency onto binary
setkey(binary, Label)
setkey(label.frequency, Label)
binary = label.frequency[binary]

# remove Label and value from binary
binary[, c("Label", "value") := NULL]

# join binary onto dat
setkey(binary, id)
setkey(dat, id)
dat = binary[dat]

# order dat by id
dat = dat[order(id)]

# get the binary subgrade_c variables from dat
binary = data.table(dat[,c("id", names(dat)[which(grepl("subgrade_c", names(dat)))]), with = FALSE])

# create a subgrade_c0 variable
binary[, subgrade_c0 := as.integer(1 - subgrade_c1 - subgrade_c2 - subgrade_c3 - subgrade_c4 - subgrade_c5)]

# convert binary to long format
binary = melt(binary, id.vars = "id", variable.name = "Label")

# only keep rows where value == 1
binary = binary[value == 1]

# check that each id has an assigned binary value
nrow(binary) == nrow(dat)

# compute the count of each label in binary
label.frequency = table(binary$Label)

# make label.frequency into a data table
label.frequency = data.table(Label = names(label.frequency),
                             subgrade_c_Frequency = as.numeric(label.frequency / sum(label.frequency)))

# join label.frequency onto binary
setkey(binary, Label)
setkey(label.frequency, Label)
binary = label.frequency[binary]

# remove Label and value from binary
binary[, c("Label", "value") := NULL]

# join binary onto dat
setkey(binary, id)
setkey(dat, id)
dat = binary[dat]

# order dat by id
dat = dat[order(id)]

# get the binary subgrade_d variables from dat
binary = data.table(dat[,c("id", names(dat)[which(grepl("subgrade_d", names(dat)))]), with = FALSE])

# create a subgrade_d0 variable
binary[, subgrade_d0 := as.integer(1 - subgrade_d1 - subgrade_d2 - subgrade_d3 - subgrade_d4 - subgrade_d5)]

# convert binary to long format
binary = melt(binary, id.vars = "id", variable.name = "Label")

# only keep rows where value == 1
binary = binary[value == 1]

# check that each id has an assigned binary value
nrow(binary) == nrow(dat)

# compute the count of each label in binary
label.frequency = table(binary$Label)

# make label.frequency into a data table
label.frequency = data.table(Label = names(label.frequency),
                             subgrade_d_Frequency = as.numeric(label.frequency / sum(label.frequency)))

# join label.frequency onto binary
setkey(binary, Label)
setkey(label.frequency, Label)
binary = label.frequency[binary]

# remove Label and value from binary
binary[, c("Label", "value") := NULL]

# join binary onto dat
setkey(binary, id)
setkey(dat, id)
dat = binary[dat]

# order dat by id
dat = dat[order(id)]

# get the binary subgrade_e variables from dat
binary = data.table(dat[,c("id", names(dat)[which(grepl("subgrade_e", names(dat)))]), with = FALSE])

# create a subgrade_e0 variable
binary[, subgrade_e0 := as.integer(1 - subgrade_e1 - subgrade_e2 - subgrade_e3 - subgrade_e4 - subgrade_e5)]

# convert binary to long format
binary = melt(binary, id.vars = "id", variable.name = "Label")

# only keep rows where value == 1
binary = binary[value == 1]

# check that each id has an assigned binary value
nrow(binary) == nrow(dat)

# compute the count of each label in binary
label.frequency = table(binary$Label)

# make label.frequency into a data table
label.frequency = data.table(Label = names(label.frequency),
                             subgrade_e_Frequency = as.numeric(label.frequency / sum(label.frequency)))

# join label.frequency onto binary
setkey(binary, Label)
setkey(label.frequency, Label)
binary = label.frequency[binary]

# remove Label and value from binary
binary[, c("Label", "value") := NULL]

# join binary onto dat
setkey(binary, id)
setkey(dat, id)
dat = binary[dat]

# order dat by id
dat = dat[order(id)]

# get the binary subgrade_f variables from dat
binary = data.table(dat[,c("id", names(dat)[which(grepl("subgrade_f", names(dat)))]), with = FALSE])

# create a subgrade_f0 variable
binary[, subgrade_f0 := as.integer(1 - subgrade_f1 - subgrade_f2 - subgrade_f3 - subgrade_f4 - subgrade_f5)]

# convert binary to long format
binary = melt(binary, id.vars = "id", variable.name = "Label")

# only keep rows where value == 1
binary = binary[value == 1]

# check that each id has an assigned binary value
nrow(binary) == nrow(dat)

# compute the count of each label in binary
label.frequency = table(binary$Label)

# make label.frequency into a data table
label.frequency = data.table(Label = names(label.frequency),
                             subgrade_f_Frequency = as.numeric(label.frequency / sum(label.frequency)))

# join label.frequency onto binary
setkey(binary, Label)
setkey(label.frequency, Label)
binary = label.frequency[binary]

# remove Label and value from binary
binary[, c("Label", "value") := NULL]

# join binary onto dat
setkey(binary, id)
setkey(dat, id)
dat = binary[dat]

# order dat by id
dat = dat[order(id)]

# get the binary subgrade_g variables from dat
binary = data.table(dat[,c("id", names(dat)[which(grepl("subgrade_g", names(dat)))]), with = FALSE])

# create a subgrade_g0 variable
binary[, subgrade_g0 := as.integer(1 - subgrade_g1 - subgrade_g2 - subgrade_g3 - subgrade_g4 - subgrade_g5)]

# convert binary to long format
binary = melt(binary, id.vars = "id", variable.name = "Label")

# only keep rows where value == 1
binary = binary[value == 1]

# check that each id has an assigned binary value
nrow(binary) == nrow(dat)

# compute the count of each label in binary
label.frequency = table(binary$Label)

# make label.frequency into a data table
label.frequency = data.table(Label = names(label.frequency),
                             subgrade_g_Frequency = as.numeric(label.frequency / sum(label.frequency)))

# join label.frequency onto binary
setkey(binary, Label)
setkey(label.frequency, Label)
binary = label.frequency[binary]

# remove Label and value from binary
binary[, c("Label", "value") := NULL]

# join binary onto dat
setkey(binary, id)
setkey(dat, id)
dat = binary[dat]

# order dat by id
dat = dat[order(id)]

# get the binary housing variables from dat
binary = data.table(dat[,c("id", "rent", "own", "mortgage", "other_home_ownership")])

# convert binary to long format
binary = melt(binary, id.vars = "id", variable.name = "Label")

# only keep rows where value == 1
binary = binary[value == 1]

# check that each id has an assigned binary value
nrow(binary) == nrow(dat)

# compute the count of each label in binary
label.frequency = table(binary$Label)

# make label.frequency into a data table
label.frequency = data.table(Label = names(label.frequency),
                             housing_Frequency = as.numeric(label.frequency / sum(label.frequency)))

# join label.frequency onto binary
setkey(binary, Label)
setkey(label.frequency, Label)
binary = label.frequency[binary]

# remove Label and value from binary
binary[, c("Label", "value") := NULL]

# join binary onto dat
setkey(binary, id)
setkey(dat, id)
dat = binary[dat]

# order dat by id
dat = dat[order(id)]

# get the binary mo variables from dat
binary = data.table(dat[,c("id", "X36mo", "X60mo")])

# convert binary to long format
binary = melt(binary, id.vars = "id", variable.name = "Label")

# only keep rows where value == 1
binary = binary[value == 1]

# check that each id has an assigned binary value
nrow(binary) == nrow(dat)

# compute the count of each label in binary
label.frequency = table(binary$Label)

# make label.frequency into a data table
label.frequency = data.table(Label = names(label.frequency),
                             mo_Frequency = as.numeric(label.frequency / sum(label.frequency)))

# join label.frequency onto binary
setkey(binary, Label)
setkey(label.frequency, Label)
binary = label.frequency[binary]

# remove Label and value from binary
binary[, c("Label", "value") := NULL]

# join binary onto dat
setkey(binary, id)
setkey(dat, id)
dat = binary[dat]

# order dat by id
dat = dat[order(id)]

# get the binary p_ variables from dat
binary = data.table(dat[,c("id", names(dat)[which(grepl("p_", names(dat)))]), with = FALSE])

# remove the emp_ variables from binary
binary = data.table(binary[,names(binary)[-which(grepl("emp_", names(binary)))], with = FALSE])

# convert binary to long format
binary = melt(binary, id.vars = "id", variable.name = "Label")

# only keep rows where value == 1
binary = binary[value == 1]

# check that each id has an assigned binary value
nrow(binary) == nrow(dat)

# compute the count of each label in binary
label.frequency = table(binary$Label)

# make label.frequency into a data table
label.frequency = data.table(Label = names(label.frequency),
                             p_Frequency = as.numeric(label.frequency / sum(label.frequency)))

# join label.frequency onto binary
setkey(binary, Label)
setkey(label.frequency, Label)
binary = label.frequency[binary]

# remove Label and value from binary
binary[, c("Label", "value") := NULL]

# join binary onto dat
setkey(binary, id)
setkey(dat, id)
dat = binary[dat]

# order dat by id
dat = dat[order(id)]

# ---- finalize prep work ----

# create 3 lists for storing model results
predict.pmt_status = list()
predict.total_return = list()
predict.total_return_class = list()

# choose how many folds to use for splitting the data:
# training data set gets ((k.fold - 2) / k.fold) * 100 percent of the data
# validation and testing data sets each get (1 / k.fold) * 100 percent of the data
k.folds = 5

# should we use the normalized response variable?
use.norm = FALSE

# remove some objects
rm(bin.table, class.order, num.bins, binary, states, label.frequency)

# free up RAM
gc()

}

# -----------------------------------------------------------------------------------
# ---- Rank the Indicators with Random Forests --------------------------------------
# -----------------------------------------------------------------------------------

{

# should we rank and remove indicators in the data?
remove.indicators = TRUE

if(remove.indicators)
{
  # ---- remove highly correlated indicators ------------------------------------------
  
  # get the indicators
  indicators = data.table(dat[, !c(responses, "id"), with = FALSE])
  
  # compute correlations
  cors = cor(indicators)
  
  # replace any NA's with 1's
  cors[is.na(cors)] = 1
  
  # find out which indicators are highly correlated
  find.indicators = findCorrelation(cors, cutoff = 0.9, names = TRUE, exact = TRUE)
  find.indicators
  
  # create an original copy of dat
  dat.copy = data.table(dat)
  
  # remove highly correlated indicators
  if(length(find.indicators) > 0) dat = dat[, !find.indicators, with = FALSE]
  
  # pick which way to rank indicators
  use.pmt_status = TRUE
  use.total_return = FALSE
  use.total_return_class = FALSE
  
  # initialize the h2o instance
  h2o.init(nthreads = workers, max_mem_size = "8g")
  
  # remove any objects in the h2o instance
  h2o.removeAll()
  
  # remove the progress bar when model building
  h2o.no_progress()
  
  # ---- rank according to pmt_status -------------------------------------------------
  
  if(use.pmt_status)
  {
    # identify predictors (x) and response (y)
    y = "pmt_status"
    x = names(dat[, !c("id", responses), with = FALSE])
    
    # make x and y into an h2o object
    YX.h2o = as.h2o(dat[, c(y, x), with = FALSE])
    
    # build the fold assignment
    set.seed(42)
    k = 4
    folds = createFolds(y = unname(unlist(dat[, y, with = FALSE])), k = k)
    
    # split up dat into train, valid, and test
    train.rows = unname(unlist(lapply(1:(k - 1), function(f) folds[[f]])))
    train = data.table(dat[train.rows])
    
    valid.rows = unname(unlist(folds[[k]]))
    valid = data.table(dat[valid.rows])
    
    # split up YX.h2o into train, valid, and test
    train.YX.h2o = as.h2o(train[, c(y, x), with = FALSE])
    valid.YX.h2o = as.h2o(valid[, c(y, x), with = FALSE])
    
    # compute the max class weight
    max.class.weight = table(as.vector(YX.h2o[,1]))
    max.class.weight = as.numeric(max(max(max.class.weight) / max.class.weight))
    
    # remove YX.h2o
    rm(YX.h2o)
    
    # set up hyperparameters of interest
    rf.hyper.params = list(ntrees = 250,
                           min_rows = 11,
                           max_depth = 30,
                           # min_rows = c(1, 11, 25),
                           # max_depth = c(20, 40, 60),
                           stopping_metric = "mean_per_class_error")
    
    # lets use a random grid search and specify a time limit and/or model limit
    minutes = 5
    rf.search.criteria = list(strategy = "RandomDiscrete", 
                              max_runtime_secs = minutes * 60, 
                              # max_models = 100, 
                              seed = 42)
    
    # lets run a grid search for a good model, without drop out ratios
    h2o.rm("rf.random.grid")
    rf.random.grid = h2o.grid(algorithm = "randomForest",
                              grid_id = "rf.random.grid",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              # stopping_rounds = 5,
                              # histogram_type = "RoundRobin",
                              nfolds = 4,
                              fold_assignment = "Stratified",
                              seed = 21,
                              balance_classes = TRUE,
                              max_after_balance_size = max.class.weight,
                              hyper_params = rf.hyper.params,
                              search_criteria = rf.search.criteria)
    
    # free up RAM
    gc()
    
    # rank each model in the random grids
    rf.grid = h2o.getGrid("rf.random.grid", sort_by = "auc", decreasing = TRUE)
    rf.grid
    
    # pick the top model from all grid searches
    imp.rf = h2o.getModel(rf.grid@model_ids[[1]])
    
    # extract variable importance
    imp = data.table(imp.rf@model$variable_importances)
    
  } else if(use.total_return)
    
  # ---- rank according to total_return -----------------------------------------------
  
  {
    # identify the response (y)
    if(use.norm)
    {
      y = "total_return_norm"
      
    } else
    {
      y = "total_return"
    }
    
    # identify predictors (x)
    x = names(dat[, !c("id", responses), with = FALSE])
    
    # make x and y into an h2o object
    YX.h2o = as.h2o(dat[, c(y, x), with = FALSE])
    
    # build the fold assignment
    set.seed(42)
    k = 4
    folds = createFolds(y = unname(unlist(dat[, y, with = FALSE])), k = k)
    
    # split up dat into train, valid, and test
    train.rows = unname(unlist(lapply(1:(k - 1), function(f) folds[[f]])))
    train = data.table(dat[train.rows])
    
    valid.rows = unname(unlist(folds[[k]]))
    valid = data.table(dat[valid.rows])
    
    # split up YX.h2o into train, valid, and test
    train.YX.h2o = as.h2o(train[, c(y, x), with = FALSE])
    valid.YX.h2o = as.h2o(valid[, c(y, x), with = FALSE])
    
    # remove YX.h2o
    rm(YX.h2o)
    
    # set up hyperparameters of interest
    rf.hyper.params = list(ntrees = 250,
                           # min_rows = c(1, 11, 25),
                           # max_depth = c(20, 40, 60),
                           min_rows = 11,
                           max_depth = 30,
                           stopping_metric = "rmse")
    
    # lets use a random grid search and specify a time limit and/or model limit
    minutes = 5
    rf.search.criteria = list(strategy = "RandomDiscrete", 
                              max_runtime_secs = minutes * 60, 
                              # max_models = 100, 
                              seed = 42)
    
    # lets run a grid search for a good model, without drop out ratios
    h2o.rm("rf.random.grid")
    rf.random.grid = h2o.grid(algorithm = "randomForest",
                              grid_id = "rf.random.grid",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              # stopping_rounds = 5,
                              # histogram_type = "RoundRobin",
                              nfolds = 4,
                              seed = 21,
                              hyper_params = rf.hyper.params,
                              search_criteria = rf.search.criteria)
    
    # free up RAM
    gc()
    
    # rank each model in the random grids
    rf.grid = h2o.getGrid("rf.random.grid", sort_by = "rmse", decreasing = FALSE)
    rf.grid
    
    # pick the top model from all grid searches
    imp.rf = h2o.getModel(rf.grid@model_ids[[1]])
    
    # extract variable importance
    imp = data.table(imp.rf@model$variable_importances)
    
  } else 
    
  # ---- rank according to total_return_class -----------------------------------------
  
  {
    # identify predictors (x) and response (y)
    y = "total_return_class"
    x = names(dat[, !c("id", responses), with = FALSE])
    
    # make x and y into an h2o object
    YX.h2o = as.h2o(dat[, c(y, x), with = FALSE])
    
    # build the fold assignment
    set.seed(42)
    k = 4
    folds = createFolds(y = unname(unlist(dat[, y, with = FALSE])), k = k)
    
    # split up dat into train, valid, and test
    train.rows = unname(unlist(lapply(1:(k - 1), function(f) folds[[f]])))
    train = data.table(dat[train.rows])
    
    valid.rows = unname(unlist(folds[[k]]))
    valid = data.table(dat[valid.rows])
    
    # split up YX.h2o into train, valid, and test
    train.YX.h2o = as.h2o(train[, c(y, x), with = FALSE])
    valid.YX.h2o = as.h2o(valid[, c(y, x), with = FALSE])
    
    # compute the max class weight
    max.class.weight = table(as.vector(YX.h2o[,1]))
    max.class.weight = as.numeric(max(max(max.class.weight) / max.class.weight))
    
    # remove YX.h2o
    rm(YX.h2o)
    
    # set up hyperparameters of interest
    rf.hyper.params = list(ntrees = 250,
                           # min_rows = c(1, 11, 25),
                           # max_depth = c(20, 40, 60),
                           min_rows = 11,
                           max_depth = 30,
                           stopping_metric = "mean_per_class_error")
    
    # lets use a random grid search and specify a time limit and/or model limit
    minutes = 5
    rf.search.criteria = list(strategy = "RandomDiscrete", 
                              max_runtime_secs = minutes * 60, 
                              # max_models = 100, 
                              seed = 42)
    
    # lets run a grid search for a good model, without drop out ratios
    h2o.rm("rf.random.grid")
    rf.random.grid = h2o.grid(algorithm = "randomForest",
                              grid_id = "rf.random.grid",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              # stopping_rounds = 5,
                              # histogram_type = "RoundRobin",
                              nfolds = 4,
                              fold_assignment = "Stratified",
                              seed = 21,
                              balance_classes = TRUE,
                              max_after_balance_size = max.class.weight,
                              hyper_params = rf.hyper.params,
                              search_criteria = rf.search.criteria)
    
    # free up RAM
    gc()
    
    # rank each model in the random grids
    rf.grid = h2o.getGrid("rf.random.grid", sort_by = "mean_per_class_error", decreasing = FALSE)
    rf.grid
    
    # pick the top model from all grid searches
    imp.rf = h2o.getModel(rf.grid@model_ids[[1]])
    
    # extract variable importance
    imp = data.table(imp.rf@model$variable_importances)
  }
  
  # ---- select indicators ------------------------------------------------------------
  
  # make variable into a factor
  imp[, variable := factor(variable, levels = unique(variable))]
  
  # pick a cut off value
  # this value should show where importance drops the most (look for the center of the "knee")
  cutoff = 0.08
  
  # plot a barplot of variable importance
  imp.plot1 = ggplot(imp, aes(x = variable, y = scaled_importance, fill = scaled_importance, color = scaled_importance)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = cutoff, color = "blue", linetype = "dashed", size = 1.1) + 
    ggtitle("GINI Importance\nRandom Forests") +
    labs(x = "Variable", y = "Scaled Importance") +
    scale_y_continuous(labels = percent) +
    scale_fill_gradient(low = "yellow", high = "red") +
    scale_color_gradient(low = "yellow", high = "red") +
    theme_dark(25) +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          panel.grid.major.x = element_blank())
  
  # plot a density plot of variable importance
  imp.plot2 = ggplot(imp, aes(x = scaled_importance)) +
    geom_density(fill = "cornflowerblue", alpha = 2/3) +
    geom_vline(xintercept = cutoff, color = "red", linetype = "dashed", size = 1.1) + 
    ggtitle("GINI Importance\nRandom Forests") +
    labs(x = "Scaled Importance", y = "Density") +
    scale_x_continuous(labels = percent) +
    coord_flip() + 
    theme_bw(25) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # show the importance plots
  grid.arrange(imp.plot1, imp.plot2, nrow = 1)
  
  # check out imp
  imp
  
  # find which indicators meet the cutoff value
  keep.indicators = as.character(unname(unlist(imp[scaled_importance >= cutoff, .(variable)])))
  
  # update dat to only include the response variables and the indicator variables of interest
  dat = dat[, c("id", responses, keep.indicators), with = FALSE]
  
  # remove some objects
  rm(cors, cutoff, find.indicators, k,
     imp.rf, keep.indicators, indicators, 
     rf.grid, rf.hyper.params, rf.random.grid, rf.search.criteria)
  
  # free up RAM
  gc()
  
} else
{
  # initialize the h2o instance
  h2o.init(nthreads = workers, max_mem_size = "8g")
  
  # remove any objects in the h2o instance
  h2o.removeAll()
  
  # remove the progress bar when model building
  h2o.no_progress()
}

}

# -----------------------------------------------------------------------------------
# ---- Remove Anomolies with Autoencoders -------------------------------------------
# -----------------------------------------------------------------------------------

{

# should we remove the anomolies in the data?
remove.data = TRUE

if(remove.data)
{
  # identify predictors (x)
  x = names(dat[, !c("id", responses), with = FALSE])
  
  # make x into an h2o object
  X.h2o = as.h2o(dat[, x, with = FALSE])
  
  # set up hyperparameters of interest
  nnet.hyper.params = list(hidden = list(c(round((2/3) * length(x), 0), round((2/3)^2 * length(x), 0), round((2/3) * length(x), 0)),
                                         c(round((3/4) * length(x), 0), round((3/4)^2 * length(x), 0), round((3/4) * length(x), 0))),
                           epochs = c(30),
                           activation = "Tanh",
                           l1 = 0,
                           l2 = 0,
                           rho = c(0.9, 0.95, 0.99, 0.999),
                           epsilon = 1e-8,
                           adaptive_rate = TRUE)
  
  # lets use a random grid search and specify a time limit and/or model limit
  minutes = 5
  nnet.search.criteria = list(strategy = "RandomDiscrete", 
                              max_runtime_secs = minutes * 60, 
                              # max_models = 100, 
                              seed = 42)
  
  # run a random grid search for a good model
  h2o.rm("nnet.random.grid")
  nnet.random.grid = h2o.grid(algorithm = "deeplearning",
                              grid_id = "nnet.random.grid",
                              autoencoder = TRUE,
                              x = x,
                              training_frame = X.h2o,
                              seed = 21,
                              hyper_params = nnet.hyper.params,
                              search_criteria = nnet.search.criteria)
  
  # rank each model in the random grid
  nnet.grid = h2o.getGrid("nnet.random.grid", sort_by = "rmse", decreasing = FALSE)
  
  # check out model ranking
  nnet.grid
  
  # get the best model from our grid search
  nnet.mod = h2o.getModel(nnet.grid@model_ids[[1]])
  
  # get the summary table of the grid search
  DT.nnet.grid = data.table(nnet.grid@summary_table)
  
  # compute the reconstruction error for each observation in X
  anomaly.error = unlist(as.data.table(h2o.anomaly(nnet.mod, X.h2o)))
  
  # plot the reconstruction eror to find out a cutoff value
  plot(sort(anomaly.error))
  
  # pick a cutoff value: (1 - cut.prob) of the rows will be removed
  cut.prob = 0.95
  anomaly.sorted = sort(anomaly.error)
  anomaly.cutoff = quantile(x = anomaly.sorted, probs = cut.prob)
  anomaly.cutoff = anomaly.sorted[which.min(abs(anomaly.sorted - anomaly.cutoff))[1]]
  lines(x = rep(anomaly.cutoff, length(anomaly.error)), col = "blue")
  
  # update dat to not contain any observations larger than the reconstruction error cuttoff value
  anomaly.dat = data.table(dat[anomaly.error > anomaly.cutoff])
  dat = dat[anomaly.error <= anomaly.cutoff]
  
  # remove some objects
  rm(anomaly.cutoff, DT.nnet.grid, minutes, nnet.grid, nnet.hyper.params, cut.prob,
     nnet.mod, nnet.random.grid, nnet.search.criteria, X.h2o, anomaly.sorted)
  
  # clean up the data in the h2o cluster
  h2o.removeAll()
  
  # shut down the h2o cluster to relieve RAM
  h2o.shutdown(prompt = FALSE)
  
  # free up RAM
  gc()
  
} else
{
  # clean up the data in the h2o cluster
  h2o.removeAll()
  
  # shut down the h2o cluster to free up RAM
  h2o.shutdown(prompt = FALSE)
  
  # free up RAM
  gc()
}

}

# -----------------------------------------------------------------------------------
# ---- Cluster the Data with Kmeans -------------------------------------------------
# -----------------------------------------------------------------------------------

{

# should we cluster the data into heterogeneous groups?
cluster.data = FALSE

if(cluster.data)
{
  # identify predictors (x)
  x = names(dat[, !c("id", responses), with = FALSE])
  
  # split up dat based on x
  X = data.table(dat[, x, with = FALSE])
  
  # normalize the data set for kmeans clustering
  X = data.table(scale(X))
  
  # build a kmeans model for 2, 3, and 4 groups
  groups = c(2, 3, 4)
  km.mods = lapply(groups, function(k)
  {
    # build the model
    set.seed(42)
    output = kmeans(x = X, 
                    centers = k, 
                    nstart = 10, 
                    iter.max = 1000,
                    algorithm = "Lloyd")
    
    # free memory
    gc()
    
    return(output)
  })
  
  # extract the cluster solutions from each model
  X.clusters = foreach(i = 1:length(groups), .combine = "cbind") %do%
  {
    return(as.numeric(km.mods[[i]]$cluster))
  }
  
  # give each column better names
  colnames(X.clusters) = paste0("Group", groups)
  
  # add each of the cluster solutions to our data set
  X.clusters = data.table(cbind(X.clusters, X))
  
  # convert the cluster solution columns into long format with 2 columns: Model and Cluster
  X.clusters = melt(X.clusters, 
                    measure.vars = paste0("Group", groups), 
                    variable.name = "Model", 
                    value.name = "Cluster")
  
  # compute the distance matrix used by kmeans
  # km.dmat = dist(X)
  
  # lets create cluster plots for each model
  km.plots = lapply(unique(X.clusters$Model), function(i)
  {
    # extract model i from pca.clusters
    DT = data.table(X.clusters[Model == i])
    
    # remove Model from DT
    DT[, Model := NULL]
    
    # build cluster plots
    output = plot.clusters(dat = DT, 
                           cluster.column.name = "Cluster",
                           # distance.matrix = km.dmat,
                           DC.title = paste("Discriminant Coordinate Cluster Plot\nK-Means ", i), 
                           pairs.title = paste("Cluster Pairs Plot\nK-Means ", i), 
                           silhouette.title = paste("Silhouette Width\nK-Means ", i))
    
    return(output)
  })
  
  # look at the silhouette plots of each model
  windows()
  km.plots[[1]]$plot.2D
  
  windows()
  km.plots[[2]]$plot.2D
  
  windows()
  km.plots[[3]]$plot.2D
  
  windows()
  grid.draw(km.plots[[1]]$plot.3D)
  
  windows()
  grid.draw(km.plots[[2]]$plot.3D)
  
  windows()
  grid.draw(km.plots[[3]]$plot.3D)
  
  # update dat to have the Clusters from Group3
  dat[, Cluster := unname(unlist(X.clusters[Model == "Group3", .(Cluster)]))]
  
  # remove some objects
  rm(groups, km.mods, km.plots, X, X.clusters, x, i)
  
  # free up RAM
  gc()
}

}

# -----------------------------------------------------------------------------------
# ---- Build Regression Models ------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# initialize the h2o instance
h2o.init(nthreads = workers, max_mem_size = "8g")

# remove any objects in the h2o instance
h2o.removeAll()

# remove the progress bar when model building
h2o.no_progress()

# ---- build prediction models for pmt_status ---------------------------------------

# determine how many data sets we have to build models for
if(cluster.data)
{
  # get the sub groups of dat
  groups = sort(unique(dat$Cluster))
  
  # also consider using all of dat
  groups = c("All", groups)
  
} else
{
  groups = "All"
}

# determine regression predictive performance for each group
predict.pmt_status.reg = foreach(g = groups) %do%
{
  # ---- split the data ----
  
  # determine which data set to use
  if(cluster.data)
  {
    if(g == "All")
    {
      # update dat to not contain Cluster
      dat.g = data.table(dat[, !"Cluster"])
      dat.g[, pmt_status := factor(pmt_status, levels = 0:1)]
      
    } else
    {
      # update dat to only contain data in Cluster g
      dat.g = data.table(dat[Cluster == g, !"Cluster"])
      dat.g[, pmt_status := factor(pmt_status, levels = 0:1)]
    }
  } else
  {
    # get dat
    dat.g = data.table(dat)
    dat.g[, pmt_status := factor(pmt_status, levels = 0:1)]
  }
  
  # identify predictors (x)
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # consider two-way interactions between the first n indicators in dat.g
  # recall that the indicators are ordered by importance from random forests
  n = floor(length(x) / 3)
  two.way = data.table(dat.g[, x[1:n], with = FALSE])
  
  # should we add two way interactions in the model?
  add.interactions = FALSE
  
  # build interactions if desired
  if(add.interactions)
  {
    # create the two way interactions and add them to dat.g
    two.way = data.table(model.matrix(~.^2, data = two.way)[,-(1:(n + 1))])
    dat.g = cbind(dat.g, two.way)
  }
  
  # identify predictors (x) and response (y)
  y = "pmt_status"
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # build the fold assignment
  set.seed(42)
  k.folds = max(c(k.folds, 3))
  folds = createFolds(y = unname(unlist(dat.g[, y, with = FALSE])), k = k.folds)
  
  # split up dat.g into train, valid, and test
  train.rows = unname(unlist(lapply(1:(k.folds - 2), function(f) folds[[f]])))
  train = data.table(dat.g[train.rows])
  
  valid.rows = unname(unlist(folds[[k.folds - 1]]))
  valid = data.table(dat.g[valid.rows])
  
  test.rows = unname(unlist(folds[[k.folds]]))
  test = data.table(dat.g[test.rows])
  
  # split up YX.h2o into train, valid, and test
  train.YX.h2o = as.h2o(train[, c(y, x), with = FALSE])
  valid.YX.h2o = as.h2o(valid[, c(y, x), with = FALSE])
  test.YX.h2o = as.h2o(test[, c(y, x), with = FALSE])
  
  # ---- grid search for models ----
  
  # set up hyperparameters of interest
  glm.hyper.params = list(lambda = c(1, 0.5, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0),
                          alpha = c(0, 0.5, 1))
  
  # lets use a random grid search and specify a time limit and/or model limit
  minutes = 5
  glm.search.criteria = list(strategy = "RandomDiscrete", 
                             max_runtime_secs = minutes * 60, 
                             # max_models = 100, 
                             seed = 42)
  
  # lets run a grid search for a good model with intercept = FALSE and standardize = FALSE
  h2o.rm("glm.random.gridA")
  glm.random.gridA = h2o.grid(algorithm = "glm",
                              grid_id = "glm.random.gridA",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              early_stopping = TRUE,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              family = "binomial",
                              intercept = FALSE,
                              standardize = FALSE,
                              seed = 21,
                              solver = "COORDINATE_DESCENT",
                              hyper_params = glm.hyper.params,
                              search_criteria = glm.search.criteria)
  
  # lets run a grid search for a good model with intercept = TRUE and standardize = FALSE
  h2o.rm("glm.random.gridB")
  glm.random.gridB = h2o.grid(algorithm = "glm",
                              grid_id = "glm.random.gridB",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              early_stopping = TRUE,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              family = "binomial",
                              intercept = TRUE,
                              standardize = FALSE,
                              seed = 21,
                              solver = "COORDINATE_DESCENT",
                              hyper_params = glm.hyper.params,
                              search_criteria = glm.search.criteria)
  
  # lets run a grid search for a good model with intercept = TRUE and standardize = TRUE
  h2o.rm("glm.random.gridC")
  glm.random.gridC = h2o.grid(algorithm = "glm",
                              grid_id = "glm.random.gridC",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              early_stopping = TRUE,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              family = "binomial",
                              intercept = TRUE,
                              standardize = TRUE,
                              seed = 21,
                              solver = "COORDINATE_DESCENT",
                              hyper_params = glm.hyper.params,
                              search_criteria = glm.search.criteria)
  
  # lets run a grid search for a good model with intercept = FALSE and standardize = TRUE
  h2o.rm("glm.random.gridD")
  glm.random.gridD = h2o.grid(algorithm = "glm",
                              grid_id = "glm.random.gridD",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              early_stopping = TRUE,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              family = "binomial",
                              intercept = FALSE,
                              standardize = TRUE,
                              seed = 21,
                              solver = "COORDINATE_DESCENT",
                              hyper_params = glm.hyper.params,
                              search_criteria = glm.search.criteria)
  
  # free up RAM
  gc()
  
  # rank each model in the random grids
  glm.gridA = h2o.getGrid("glm.random.gridA", sort_by = "auc", decreasing = TRUE)
  glm.gridB = h2o.getGrid("glm.random.gridB", sort_by = "auc", decreasing = TRUE)
  glm.gridC = h2o.getGrid("glm.random.gridC", sort_by = "auc", decreasing = TRUE)
  glm.gridD = h2o.getGrid("glm.random.gridD", sort_by = "auc", decreasing = TRUE)
  
  # combine all the grid tables into one grid table that considers the options for intercept and standardize
  glm.grid = rbind(cbind(data.table(glm.gridA@summary_table), intercept = "FALSE", standardize = "FALSE", search = 1),
                   cbind(data.table(glm.gridB@summary_table), intercept = "TRUE", standardize = "FALSE", search = 2),
                   cbind(data.table(glm.gridC@summary_table), intercept = "TRUE", standardize = "TRUE", search = 3),
                   cbind(data.table(glm.gridD@summary_table), intercept = "FALSE", standardize = "TRUE", search = 4))
  
  # combine all the grid models
  glm.grid.models = list(glm.gridA, glm.gridB, glm.gridC, glm.gridD)
  
  # order grid by auc
  glm.grid = glm.grid[order(as.numeric(auc), decreasing = TRUE)]
  
  # get the summary table of the grid search
  DT.glm.grid = data.table(glm.grid)
  
  # set up the data types for each column in DT.grid for plotting purposes
  DT.glm.grid = DT.glm.grid[, .(alpha = factor(as.numeric(gsub("[", "", gsub("]", "", alpha, fixed = TRUE), fixed = TRUE)), levels = c(1, 0.5, 0)),
                                lambda = factor(as.numeric(gsub("[", "", gsub("]", "", lambda, fixed = TRUE), fixed = TRUE)), levels = c(1, 0.5, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0)),
                                intercept = as.factor(intercept),
                                standardize = as.factor(standardize),
                                model = removePunctuation(gsub("[A-z]+", "", model_ids)),
                                auc = as.numeric(auc),
                                search = as.numeric(search))]
  
  # plot auc v. standardize, lambda, and alpha to see which structure is most robust
  plot.glm.grid = ggplot(DT.glm.grid, aes(x = lambda, y = auc, color = standardize, fill = standardize)) + 
    # geom_boxplot() + 
    geom_jitter(size = 3, alpha = 2/3) + 
    # scale_y_continuous(labels = dollar) + 
    ggtitle("Cross Validation Error") + 
    labs(x = "Strength of Regularization", y = "AUC", color = "Standardize", fill = "Standardize") + 
    facet_wrap(~paste("L1/L2 Distribution:", alpha), nrow = 1) +
    theme_bw(base_size = 25) +
    theme(legend.position = "top", 
          legend.key.size = unit(.25, "in"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
  
  plot.glm.grid
  
  # pick the top model from all grid searches
  # glm.mod = h2o.getModel(glm.grid.models[[DT.glm.grid$search[1]]]@model_ids[[1]])
  
  # find the top k models from all grid searches
  k = 10
  glm.mod.list = lapply(1:k, function(i)
  {
    # get grid search (s), the model id (m), and the position of the model (p)
    s = DT.glm.grid$search[i]
    m = DT.glm.grid$model[i]
    p = which(removePunctuation(gsub("[A-z]+", "", unlist(glm.grid.models[[s]]@model_ids))) == m)
    
    # get the model of interest
    output = glm.grid.models[[s]]@model_ids[[p]]
    
    return(output)
  })
  
  # stack the top k models from all grid searches
  h2o.rm("glm.mod")
  glm.mod = h2o.stackedEnsemble(x = x,
                                y = y,
                                training_frame = train.YX.h2o,
                                validation_frame = valid.YX.h2o,
                                model_id = "glm.mod",
                                base_models = glm.mod.list)
  
  # plot variable importance
  # h2o.std_coef_plot(glm.mod)
  
  # ---- measure model quality ----
  
  # make predictions on each data set
  ynew.train = as.data.frame(predict(glm.mod, newdata = train.YX.h2o))$p1
  ynew.valid = as.data.frame(predict(glm.mod, newdata = valid.YX.h2o))$p1
  ynew.test = as.data.frame(predict(glm.mod, newdata = test.YX.h2o))$p1
  
  # plot the probability predictions on the testing data set
  # densityPlot(ynew.test)
  
  # get the true values from each data set
  ytrue.train = as.data.frame(train.YX.h2o)[,1]
  ytrue.valid = as.data.frame(valid.YX.h2o)[,1]
  ytrue.test = as.data.frame(test.YX.h2o)[,1]
  
  # compute prediction metrics on each data set
  glm.metrics.train = h2o.make_metrics(predicted = as.h2o(ynew.train), 
                                       actuals = as.h2o(factor(ytrue.train, levels = sort(unique(ytrue.train)))))
  
  glm.metrics.valid = h2o.make_metrics(predicted = as.h2o(ynew.valid), 
                                       actuals = as.h2o(factor(ytrue.valid, levels = sort(unique(ytrue.valid)))))
  
  glm.metrics.test = h2o.make_metrics(predicted = as.h2o(ynew.test), 
                                      actuals = as.h2o(factor(ytrue.test, levels = sort(unique(ytrue.test)))))
  
  # free up RAM
  gc()
  
  # check out the threshold metrics
  glm.metrics.train@metrics$max_criteria_and_metric_scores
  glm.metrics.valid@metrics$max_criteria_and_metric_scores
  glm.metrics.test@metrics$max_criteria_and_metric_scores
  
  # pick a threshold for converting probabilities into binary
  threshold = mean(c(tail(glm.metrics.train@metrics$max_criteria_and_metric_scores$threshold, 1),
                     tail(glm.metrics.valid@metrics$max_criteria_and_metric_scores$threshold, 1),
                     tail(glm.metrics.test@metrics$max_criteria_and_metric_scores$threshold, 1)))
  
  # compute a confusion matrix for each data set
  glm.confusion.train = confusion(ytrue = ytrue.train, ypred = as.numeric(ynew.train >= threshold))
  glm.confusion.valid = confusion(ytrue = ytrue.valid, ypred = as.numeric(ynew.valid >= threshold))
  glm.confusion.test = confusion(ytrue = ytrue.test, ypred = as.numeric(ynew.test >= threshold))
  
  # compute Sensitivity (accuracy in predicting 1's) for each data set
  glm.sensitivity.train = glm.confusion.train[2,2] / (glm.confusion.train[2,2] + glm.confusion.train[2,1])
  glm.sensitivity.valid = glm.confusion.valid[2,2] / (glm.confusion.valid[2,2] + glm.confusion.valid[2,1])
  glm.sensitivity.test = glm.confusion.test[2,2] / (glm.confusion.test[2,2] + glm.confusion.test[2,1])
  
  # compute Specificity (accuracy in predicting 0's) for each data set
  glm.specificity.train = glm.confusion.train[1,1] / (glm.confusion.train[1,1] + glm.confusion.train[1,2])
  glm.specificity.valid = glm.confusion.valid[1,1] / (glm.confusion.valid[1,1] + glm.confusion.valid[1,2])
  glm.specificity.test = glm.confusion.test[1,1] / (glm.confusion.test[1,1] + glm.confusion.test[1,2])
  
  # compute Odds Ratio (the odds of success over failure) for each data set
  glm.odds.ratio.train = (glm.confusion.train[1,1] * glm.confusion.train[2,2]) / (glm.confusion.train[2,1] * glm.confusion.train[1,2])
  glm.odds.ratio.valid = (glm.confusion.valid[1,1] * glm.confusion.valid[2,2]) / (glm.confusion.valid[2,1] * glm.confusion.valid[1,2])
  glm.odds.ratio.test = (glm.confusion.test[1,1] * glm.confusion.test[2,2]) / (glm.confusion.test[2,1] * glm.confusion.test[1,2])
  
  # compute Accuracy for each data set
  glm.accuracy.train = (glm.confusion.train[1,1] + glm.confusion.train[2,2]) / (glm.confusion.train[1,1] + glm.confusion.train[1,2] + glm.confusion.train[2,1] + glm.confusion.train[2,2])
  glm.accuracy.valid = (glm.confusion.valid[1,1] + glm.confusion.valid[2,2]) / (glm.confusion.valid[1,1] + glm.confusion.valid[1,2] + glm.confusion.valid[2,1] + glm.confusion.valid[2,2])
  glm.accuracy.test = (glm.confusion.test[1,1] + glm.confusion.test[2,2]) / (glm.confusion.test[1,1] + glm.confusion.test[1,2] + glm.confusion.test[2,1] + glm.confusion.test[2,2])
  
  # compute AUC for each data set
  glm.auc.train = h2o.auc(glm.metrics.train)
  glm.auc.valid = h2o.auc(glm.metrics.valid)
  glm.auc.test = h2o.auc(glm.metrics.test)
  
  # compute Log Loss for each data set
  glm.logloss.train = h2o.logloss(glm.metrics.train)
  glm.logloss.valid = h2o.logloss(glm.metrics.valid)
  glm.logloss.test = h2o.logloss(glm.metrics.test)
  
  # ---- compute kappa ----
  
  # get the total number of observations for each data set
  n.train = sum(glm.confusion.train)
  n.valid = sum(glm.confusion.valid)
  n.test = sum(glm.confusion.test)
  
  # get the vector of correct predictions for each data set 
  dia.train = diag(glm.confusion.train)
  dia.valid = diag(glm.confusion.valid)
  dia.test = diag(glm.confusion.test)
  
  # get the vector of the number of observations per class for each data set
  rsum.train = rowSums(glm.confusion.train)
  rsum.valid = rowSums(glm.confusion.valid)
  rsum.test = rowSums(glm.confusion.test)
  
  # get the vector of the number of predictions per class for each data set
  csum.train = colSums(glm.confusion.train)
  csum.valid = colSums(glm.confusion.valid)
  csum.test = colSums(glm.confusion.test)
  
  # get the proportion of observations per class for each data set
  p.train = rsum.train / n.train
  p.valid = rsum.valid / n.valid
  p.test = rsum.test / n.test
  
  # get the proportion of predcitions per class for each data set
  q.train = csum.train / n.train
  q.valid = csum.valid / n.valid
  q.test = csum.test / n.test
  
  # compute accuracy for each data set
  acc.train = sum(dia.train) / n.train
  acc.valid = sum(dia.valid) / n.valid
  acc.test = sum(dia.test) / n.test
  
  # compute expected accuracy for each data set
  exp.acc.train = sum(p.train * q.train)
  exp.acc.valid = sum(p.valid * q.valid)
  exp.acc.test = sum(p.test * q.test)
  
  # compute kappa for each data set
  kap.train = (acc.train - exp.acc.train) / (1 - exp.acc.train)
  kap.valid = (acc.valid - exp.acc.valid) / (1 - exp.acc.valid)
  kap.test = (acc.test - exp.acc.test) / (1 - exp.acc.test)
  
  # ---- finalize output ----
  
  # build a final metrics table for glm.mod
  glm.mod.table = data.table(Model = "Regression",
                             Group = g,
                             Threshold = threshold,
                             Metric = c("AUC", "Accuracy", "Log_Loss", "Kappa", "Odds_Ratio", "Specificity_0", "Sensitivity_1"),
                             Train = c(glm.auc.train, glm.accuracy.train, glm.logloss.train, kap.train, glm.odds.ratio.train, glm.specificity.train, glm.sensitivity.train),
                             Valid = c(glm.auc.valid, glm.accuracy.valid, glm.logloss.valid, kap.valid, glm.odds.ratio.valid, glm.specificity.valid, glm.sensitivity.valid),
                             Test = c(glm.auc.test, glm.accuracy.test, glm.logloss.test, kap.test, glm.odds.ratio.test, glm.specificity.test, glm.sensitivity.test))
  
  # build a final grid search table
  glm.grid.search = data.table(cbind(Model = rep("Regression", nrow(DT.glm.grid)), 
                                     Group = rep(g, nrow(DT.glm.grid)), 
                                     DT.glm.grid))
  
  # build a list of confusion matrices
  glm.confusion = list("Train" = glm.confusion.train, 
                       "Valid" = glm.confusion.valid, 
                       "Test" = glm.confusion.test)
  
  # build a list of final tables
  glm.list = list(glm.mod.table, glm.grid.search, glm.confusion)
  names(glm.list) = paste0(c("Regression_Metrics_Group_", "Regression_Search_Group_", "Regression_Error_Group_"), g)
  
  # remove some objects
  rm(dat.g, DT.glm.grid, folds, glm.accuracy.test, glm.accuracy.train, glm.accuracy.valid, n,
     glm.auc.test, glm.auc.train, glm.auc.valid, glm.confusion.test, glm.confusion.train, glm.confusion.valid,
     glm.grid, glm.grid.models, glm.gridA, glm.gridB, glm.gridC, glm.gridD, glm.hyper.params, 
     glm.logloss.test, glm.logloss.train, glm.logloss.valid, glm.metrics.test, glm.metrics.train, glm.metrics.valid,
     glm.mod, glm.odds.ratio.test, glm.odds.ratio.train, glm.odds.ratio.valid, glm.random.gridA,
     glm.random.gridB, glm.random.gridC, glm.random.gridD, glm.search.criteria, glm.sensitivity.test,
     glm.sensitivity.train, glm.sensitivity.valid, glm.specificity.test, glm.specificity.train,
     glm.specificity.valid, plot.glm.grid, train, test, valid, k,
     test.YX.h2o, threshold, train.YX.h2o, two.way, valid.YX.h2o,
     x, y, ynew.test, ynew.train, ynew.valid, ytrue.test, ytrue.train, ytrue.valid,
     n.train, dia.train, rsum.train, csum.train, p.train, q.train, acc.train, exp.acc.train, kap.train,
     n.valid, dia.valid, rsum.valid, csum.valid, p.valid, q.valid, acc.valid, exp.acc.valid, kap.valid,
     n.test, dia.test, rsum.test, csum.test, p.test, q.test, acc.test, exp.acc.test, kap.test)
  
  # free up RAM
  gc()
  
  # return the final metrics table
  return(glm.list)
}

# name the results
names(predict.pmt_status.reg) = "h2o.glm"

# combine predict.pmt_status into one table
predict.pmt_status = append(predict.pmt_status, predict.pmt_status.reg)

# remove some objects
rm(glm.mod.table, glm.grid.search, glm.list, g, predict.pmt_status.reg)

# free up RAM
gc()

# ---- build prediction models for total_return ---------------------------------------

# determine how many data sets we have to build models for
if(cluster.data)
{
  # get the sub groups of dat
  groups = sort(unique(dat$Cluster))
  
  # also consider using all of dat
  groups = c("All", groups)
  
} else
{
  groups = "All"
}

# determine regression predictive performance for each group
predict.total_return.reg = foreach(g = groups) %do%
{
  # ---- split the data ----
  
  # determine which data set to use
  if(cluster.data)
  {
    if(g == "All")
    {
      # update dat to not contain Cluster
      dat.g = data.table(dat[, !"Cluster"])
      
    } else
    {
      # update dat to only contain data in Cluster g
      dat.g = data.table(dat[Cluster == g, !"Cluster"])
    }
  } else
  {
    # get dat
    dat.g = data.table(dat)
  }
  
  # identify predictors (x)
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # consider two-way interactions between the first n indicators in dat.g
  # recall that the indicators are ordered by importance from random forests
  n = floor(length(x) / 3)
  two.way = data.table(dat.g[, x[1:n], with = FALSE])
  
  # should we add two way interactions in the model?
  add.interactions = FALSE
  
  # build interactions if desired
  if(add.interactions)
  {
    # create the two way interactions and add them to dat.g
    two.way = data.table(model.matrix(~.^2, data = two.way)[,-(1:(n + 1))])
    dat.g = cbind(dat.g, two.way)
  }
  
  # identify the response (y)
  if(use.norm)
  {
    y = "total_return_norm"
    
  } else
  {
    y = "total_return"
  }
  
  # identify predictors (x)
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # build the fold assignment
  set.seed(42)
  k.folds = max(c(k.folds, 3))
  folds = createFolds(y = unname(unlist(dat.g[, y, with = FALSE])), k = k.folds)
  
  # split up dat.g into train, valid, and test
  train.rows = unname(unlist(lapply(1:(k.folds - 2), function(f) folds[[f]])))
  train = data.table(dat.g[train.rows])
  
  valid.rows = unname(unlist(folds[[k.folds - 1]]))
  valid = data.table(dat.g[valid.rows])
  
  test.rows = unname(unlist(folds[[k.folds]]))
  test = data.table(dat.g[test.rows])
  
  # split up YX.h2o into train, valid, and test
  train.YX.h2o = as.h2o(train[, c(y, x), with = FALSE])
  valid.YX.h2o = as.h2o(valid[, c(y, x), with = FALSE])
  test.YX.h2o = as.h2o(test[, c(y, x), with = FALSE])
  
  # ---- grid search for models ----
  
  # set up hyperparameters of interest
  glm.hyper.params = list(lambda = c(1, 0.5, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0),
                          alpha = c(0, 0.5, 1))
  
  # lets use a random grid search and specify a time limit and/or model limit
  minutes = 5
  glm.search.criteria = list(strategy = "RandomDiscrete", 
                             max_runtime_secs = minutes * 60, 
                             # max_models = 100, 
                             seed = 42)
  
  # lets run a grid search for a good model with intercept = FALSE and standardize = FALSE
  h2o.rm("glm.random.gridA")
  glm.random.gridA = h2o.grid(algorithm = "glm",
                              grid_id = "glm.random.gridA",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              early_stopping = TRUE,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              family = "gaussian",
                              intercept = FALSE,
                              standardize = FALSE,
                              seed = 21,
                              solver = "COORDINATE_DESCENT",
                              hyper_params = glm.hyper.params,
                              search_criteria = glm.search.criteria)
  
  # lets run a grid search for a good model with intercept = TRUE and standardize = FALSE
  h2o.rm("glm.random.gridB")
  glm.random.gridB = h2o.grid(algorithm = "glm",
                              grid_id = "glm.random.gridB",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              early_stopping = TRUE,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              family = "gaussian",
                              intercept = TRUE,
                              standardize = FALSE,
                              seed = 21,
                              solver = "COORDINATE_DESCENT",
                              hyper_params = glm.hyper.params,
                              search_criteria = glm.search.criteria)
  
  # lets run a grid search for a good model with intercept = TRUE and standardize = TRUE
  h2o.rm("glm.random.gridC")
  glm.random.gridC = h2o.grid(algorithm = "glm",
                              grid_id = "glm.random.gridC",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              early_stopping = TRUE,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              family = "gaussian",
                              intercept = TRUE,
                              standardize = TRUE,
                              seed = 21,
                              solver = "COORDINATE_DESCENT",
                              hyper_params = glm.hyper.params,
                              search_criteria = glm.search.criteria)
  
  # lets run a grid search for a good model with intercept = FALSE and standardize = TRUE
  h2o.rm("glm.random.gridD")
  glm.random.gridD = h2o.grid(algorithm = "glm",
                              grid_id = "glm.random.gridD",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              early_stopping = TRUE,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              family = "gaussian",
                              intercept = FALSE,
                              standardize = TRUE,
                              seed = 21,
                              solver = "COORDINATE_DESCENT",
                              hyper_params = glm.hyper.params,
                              search_criteria = glm.search.criteria)
  
  # free up RAM
  gc()
  
  # rank each model in the random grids
  glm.gridA = h2o.getGrid("glm.random.gridA", sort_by = "rmse", decreasing = FALSE)
  glm.gridB = h2o.getGrid("glm.random.gridB", sort_by = "rmse", decreasing = FALSE)
  glm.gridC = h2o.getGrid("glm.random.gridC", sort_by = "rmse", decreasing = FALSE)
  glm.gridD = h2o.getGrid("glm.random.gridD", sort_by = "rmse", decreasing = FALSE)
  
  # combine all the grid tables into one grid table that considers the options for intercept and standardize
  glm.grid = rbind(cbind(data.table(glm.gridA@summary_table), intercept = "FALSE", standardize = "FALSE", search = 1),
                   cbind(data.table(glm.gridB@summary_table), intercept = "TRUE", standardize = "FALSE", search = 2),
                   cbind(data.table(glm.gridC@summary_table), intercept = "TRUE", standardize = "TRUE", search = 3),
                   cbind(data.table(glm.gridD@summary_table), intercept = "FALSE", standardize = "TRUE", search = 4))
  
  # combine all the grid models
  glm.grid.models = list(glm.gridA, glm.gridB, glm.gridC, glm.gridD)
  
  # order grid by rmse
  glm.grid = glm.grid[order(as.numeric(rmse), decreasing = FALSE)]
  
  # get the summary table of the grid search
  DT.glm.grid = data.table(glm.grid)
  
  # set up the data types for each column in DT.grid for plotting purposes
  DT.glm.grid = DT.glm.grid[, .(alpha = factor(as.numeric(gsub("[", "", gsub("]", "", alpha, fixed = TRUE), fixed = TRUE)), levels = c(1, 0.5, 0)),
                                lambda = factor(as.numeric(gsub("[", "", gsub("]", "", lambda, fixed = TRUE), fixed = TRUE)), levels = c(1, 0.5, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0)),
                                intercept = as.factor(intercept),
                                standardize = as.factor(standardize),
                                model = removePunctuation(gsub("[A-z]+", "", model_ids)),
                                rmse = as.numeric(rmse),
                                search = as.numeric(search))]
  
  # plot rmse v. standardize, lambda, and alpha to see which structure is most robust
  plot.glm.grid = ggplot(DT.glm.grid, aes(x = lambda, y = rmse, color = standardize, fill = standardize)) + 
    # geom_boxplot() + 
    geom_jitter(size = 3, alpha = 2/3) + 
    # scale_y_continuous(labels = dollar) + 
    ggtitle("Cross Validation Error") + 
    labs(x = "Strength of Regularization", y = "RMSE", color = "Standardize", fill = "Standardize") + 
    facet_wrap(~paste("L1/L2 Distribution:", alpha), nrow = 1) +
    theme_bw(base_size = 25) +
    theme(legend.position = "top", 
          legend.key.size = unit(.25, "in"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
  
  plot.glm.grid
  
  # pick the top model from all grid searches
  # glm.mod = h2o.getModel(glm.grid.models[[DT.glm.grid$search[1]]]@model_ids[[1]])
  
  # find the top k models from all grid searches
  k = 10
  glm.mod.list = lapply(1:k, function(i)
  {
    # get grid search (s), the model id (m), and the position of the model (p)
    s = DT.glm.grid$search[i]
    m = DT.glm.grid$model[i]
    p = which(removePunctuation(gsub("[A-z]+", "", unlist(glm.grid.models[[s]]@model_ids))) == m)
    
    # get the model of interest
    output = glm.grid.models[[s]]@model_ids[[p]]
    
    return(output)
  })
  
  # stack the top k models from all grid searches
  h2o.rm("glm.mod")
  glm.mod = h2o.stackedEnsemble(x = x,
                                y = y,
                                training_frame = train.YX.h2o,
                                validation_frame = valid.YX.h2o,
                                model_id = "glm.mod",
                                base_models = glm.mod.list)
  
  # plot variable importance
  # names(h2o.std_coef_plot(glm.mod))
  
  # ---- measure model quality ----
  
  # make predictions on each data set
  ynew.train = as.data.frame(predict(glm.mod, newdata = train.YX.h2o))$predict
  ynew.valid = as.data.frame(predict(glm.mod, newdata = valid.YX.h2o))$predict
  ynew.test = as.data.frame(predict(glm.mod, newdata = test.YX.h2o))$predict
  
  # get the true values from each data set
  ytrue.train = as.data.frame(train.YX.h2o)[,1]
  ytrue.valid = as.data.frame(valid.YX.h2o)[,1]
  ytrue.test = as.data.frame(test.YX.h2o)[,1]
  
  # transform the predicted and actual values back to original units if needed
  if(use.norm)
  {
    # map predicted values to the original units of the response variable
    ynew.train = predict(best.norm, newdata = ynew.train, inverse = TRUE)
    ynew.valid = predict(best.norm, newdata = ynew.valid, inverse = TRUE)
    ynew.test = predict(best.norm, newdata = ynew.test, inverse = TRUE)
    
    # map actual values to the original units of the response variable
    ytrue.train = predict(best.norm, newdata = ytrue.train, inverse = TRUE)
    ytrue.valid = predict(best.norm, newdata = ytrue.valid, inverse = TRUE)
    ytrue.test = predict(best.norm, newdata = ytrue.test, inverse = TRUE)
  }
  
  # compute prediction metrics on each data set
  glm.metrics.train = h2o.make_metrics(predicted = as.h2o(ynew.train), 
                                       actuals = as.h2o(ytrue.train))
  
  glm.metrics.valid = h2o.make_metrics(predicted = as.h2o(ynew.valid), 
                                       actuals = as.h2o(ytrue.valid))
  
  glm.metrics.test = h2o.make_metrics(predicted = as.h2o(ynew.test), 
                                      actuals = as.h2o(ytrue.test))
  
  # free up RAM
  gc()
  
  # check out the prediction metrics
  glm.metrics.train
  glm.metrics.valid
  glm.metrics.test
  
  # fit the normal distribution to the residuals
  glm.resid.train = fitdist(data = ytrue.train - ynew.train, distr = "norm")
  glm.resid.valid = fitdist(data = ytrue.valid - ynew.valid, distr = "norm") 
  glm.resid.test = fitdist(data = ytrue.test - ynew.test, distr = "norm") 
  
  # ---- finalize output ----
  
  # build a final metrics table for glm.mod
  glm.mod.table = data.table(Model = "Regression",
                             Group = g,
                             Metric = c("R2", "RMSE"),
                             Train = c(glm.metrics.train@metrics$r2, glm.metrics.train@metrics$RMSE),
                             Valid = c(glm.metrics.valid@metrics$r2, glm.metrics.valid@metrics$RMSE),
                             Test = c(glm.metrics.test@metrics$r2, glm.metrics.test@metrics$RMSE))
  
  # build a final grid search table
  glm.grid.search = data.table(cbind(Model = rep("Regression", nrow(DT.glm.grid)), 
                                     Group = rep(g, nrow(DT.glm.grid)), 
                                     DT.glm.grid))
  
  # build a list of residual plots
  glm.resid = list("Train" = glm.resid.train, 
                   "Valid" = glm.resid.valid, 
                   "Test" = glm.resid.test)
  
  # build a list of final tables
  glm.list = list(glm.mod.table, glm.grid.search, glm.resid)
  names(glm.list) = paste0(c("Regression_Metrics_Group_", "Regression_Search_Group_", "Regression_Error_Group_"), g)
  
  # remove some objects
  rm(dat.g, DT.glm.grid, folds, n,
     glm.grid, glm.grid.models, glm.gridA, glm.gridB, glm.gridC, glm.gridD, glm.hyper.params, 
     glm.metrics.test, glm.metrics.train, glm.metrics.valid,
     glm.mod, glm.random.gridA, glm.resid.test, glm.resid.train, glm.resid.valid,
     glm.random.gridB, glm.random.gridC, glm.random.gridD, glm.search.criteria,
     plot.glm.grid, train, test, valid, 
     test.YX.h2o, train.YX.h2o, two.way, valid.YX.h2o,
     x, y, ynew.test, ynew.train, ynew.valid, ytrue.test, ytrue.train, ytrue.valid)
  
  # free up RAM
  gc()
  
  # return the final metrics table
  return(glm.list)
}

# name the results
names(predict.total_return.reg) = "h2o.glm"

# combine predict.total_return into one table
predict.total_return = append(predict.total_return, predict.total_return.reg)

# remove some objects
rm(glm.mod.table, glm.grid.search, glm.list, g, predict.total_return.reg)

# free up RAM
gc()

# ---- build prediction models for total_return_class ---------------------------------------

# determine how many data sets we have to build models for
if(cluster.data)
{
  # get the sub groups of dat
  groups = sort(unique(dat$Cluster))
  
  # also consider using all of dat
  groups = c("All", groups)
  
} else
{
  groups = "All"
}

# determine regression predictive performance for each group
predict.total_return_class.reg = foreach(g = groups) %do%
{
  # ---- split the data ----
  
  # determine which data set to use
  if(cluster.data)
  {
    if(g == "All")
    {
      # update dat to not contain Cluster
      dat.g = data.table(dat[, !"Cluster"])
      
    } else
    {
      # update dat to only contain data in Cluster g
      dat.g = data.table(dat[Cluster == g, !"Cluster"])
    }
  } else
  {
    # get dat
    dat.g = data.table(dat)
  }
  
  # identify predictors (x)
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # consider two-way interactions between the first n indicators in dat.g
  # recall that the indicators are ordered by importance from random forests
  n = floor(length(x) / 3)
  two.way = data.table(dat.g[, x[1:n], with = FALSE])
  
  # should we add two way interactions in the model?
  add.interactions = FALSE
  
  # build interactions if desired
  if(add.interactions)
  {
    # create the two way interactions and add them to dat.g
    two.way = data.table(model.matrix(~.^2, data = two.way)[,-(1:(n + 1))])
    dat.g = cbind(dat.g, two.way)
  }
  
  # identify predictors (x) and response (y)
  y = "total_return_class"
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # build the fold assignment
  set.seed(42)
  k.folds = max(c(k.folds, 3))
  folds = createFolds(y = unname(unlist(dat.g[, y, with = FALSE])), k = k.folds)
  
  # split up dat.g into train, valid, and test
  train.rows = unname(unlist(lapply(1:(k.folds - 2), function(f) folds[[f]])))
  train = data.table(dat.g[train.rows])
  
  valid.rows = unname(unlist(folds[[k.folds - 1]]))
  valid = data.table(dat.g[valid.rows])
  
  test.rows = unname(unlist(folds[[k.folds]]))
  test = data.table(dat.g[test.rows])
  
  # split up YX.h2o into train, valid, and test
  train.YX.h2o = as.h2o(train[, c(y, x), with = FALSE])
  valid.YX.h2o = as.h2o(valid[, c(y, x), with = FALSE])
  test.YX.h2o = as.h2o(test[, c(y, x), with = FALSE])
  
  # ---- grid search for models ----
  
  # set up hyperparameters of interest
  glm.hyper.params = list(lambda = c(1, 0.5, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0),
                          alpha = c(0, 0.5, 1))
  
  # lets use a random grid search and specify a time limit and/or model limit
  minutes = 5
  glm.search.criteria = list(strategy = "RandomDiscrete", 
                             max_runtime_secs = minutes * 60, 
                             # max_models = 100, 
                             seed = 42)
  
  # lets run a grid search for a good model with intercept = FALSE and standardize = FALSE
  h2o.rm("glm.random.gridA")
  glm.random.gridA = h2o.grid(algorithm = "glm",
                              grid_id = "glm.random.gridA",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              early_stopping = TRUE,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              family = "multinomial",
                              intercept = FALSE,
                              standardize = FALSE,
                              seed = 21,
                              solver = "COORDINATE_DESCENT",
                              hyper_params = glm.hyper.params,
                              search_criteria = glm.search.criteria)
  
  # lets run a grid search for a good model with intercept = TRUE and standardize = FALSE
  h2o.rm("glm.random.gridB")
  glm.random.gridB = h2o.grid(algorithm = "glm",
                              grid_id = "glm.random.gridB",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              early_stopping = TRUE,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              family = "multinomial",
                              intercept = TRUE,
                              standardize = FALSE,
                              seed = 21,
                              solver = "COORDINATE_DESCENT",
                              hyper_params = glm.hyper.params,
                              search_criteria = glm.search.criteria)
  
  # lets run a grid search for a good model with intercept = TRUE and standardize = TRUE
  h2o.rm("glm.random.gridC")
  glm.random.gridC = h2o.grid(algorithm = "glm",
                              grid_id = "glm.random.gridC",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              early_stopping = TRUE,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              family = "multinomial",
                              intercept = TRUE,
                              standardize = TRUE,
                              seed = 21,
                              solver = "COORDINATE_DESCENT",
                              hyper_params = glm.hyper.params,
                              search_criteria = glm.search.criteria)
  
  # lets run a grid search for a good model with intercept = FALSE and standardize = TRUE
  h2o.rm("glm.random.gridD")
  glm.random.gridD = h2o.grid(algorithm = "glm",
                              grid_id = "glm.random.gridD",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              early_stopping = TRUE,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              family = "multinomial",
                              intercept = FALSE,
                              standardize = TRUE,
                              seed = 21,
                              solver = "COORDINATE_DESCENT",
                              hyper_params = glm.hyper.params,
                              search_criteria = glm.search.criteria)
  
  # free up RAM
  gc()
  
  # rank each model in the random grids
  glm.gridA = h2o.getGrid("glm.random.gridA", sort_by = "mean_per_class_error", decreasing = FALSE)
  glm.gridB = h2o.getGrid("glm.random.gridB", sort_by = "mean_per_class_error", decreasing = FALSE)
  glm.gridC = h2o.getGrid("glm.random.gridC", sort_by = "mean_per_class_error", decreasing = FALSE)
  glm.gridD = h2o.getGrid("glm.random.gridD", sort_by = "mean_per_class_error", decreasing = FALSE)
  
  # combine all the grid tables into one grid table that considers the options for intercept and standardize
  glm.grid = rbind(cbind(data.table(glm.gridA@summary_table), intercept = "FALSE", standardize = "FALSE", search = 1),
                   cbind(data.table(glm.gridB@summary_table), intercept = "TRUE", standardize = "FALSE", search = 2),
                   cbind(data.table(glm.gridC@summary_table), intercept = "TRUE", standardize = "TRUE", search = 3),
                   cbind(data.table(glm.gridD@summary_table), intercept = "FALSE", standardize = "TRUE", search = 4))
  
  # combine all the grid models
  glm.grid.models = list(glm.gridA, glm.gridB, glm.gridC, glm.gridD)
  
  # order grid by mean_per_class_error
  glm.grid = glm.grid[order(as.numeric(mean_per_class_error), decreasing = FALSE)]
  
  # get the summary table of the grid search
  DT.glm.grid = data.table(glm.grid)
  
  # set up the data types for each column in DT.grid for plotting purposes
  DT.glm.grid = DT.glm.grid[, .(alpha = factor(as.numeric(gsub("[", "", gsub("]", "", alpha, fixed = TRUE), fixed = TRUE)), levels = c(1, 0.5, 0)),
                                lambda = factor(as.numeric(gsub("[", "", gsub("]", "", lambda, fixed = TRUE), fixed = TRUE)), levels = c(1, 0.5, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0)),
                                intercept = as.factor(intercept),
                                standardize = as.factor(standardize),
                                model = removePunctuation(gsub("[A-z]+", "", model_ids)),
                                mean_per_class_error = as.numeric(mean_per_class_error),
                                search = as.numeric(search))]
  
  # plot mean_per_class_error v. standardize, lambda, and alpha to see which structure is most robust
  plot.glm.grid = ggplot(DT.glm.grid, aes(x = lambda, y = mean_per_class_error, color = standardize, fill = standardize)) + 
    # geom_boxplot() + 
    geom_jitter(size = 3, alpha = 2/3) + 
    # scale_y_continuous(labels = dollar) + 
    ggtitle("Cross Validation Error") + 
    labs(x = "Strength of Regularization", y = "Log Loss", color = "Standardize", fill = "Standardize") + 
    facet_wrap(~paste("L1/L2 Distribution:", alpha), nrow = 1) +
    theme_bw(base_size = 25) +
    theme(legend.position = "top", 
          legend.key.size = unit(.25, "in"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
  
  plot.glm.grid
  
  # pick the top model from all grid searches
  # glm.mod = h2o.getModel(glm.grid.models[[DT.glm.grid$search[1]]]@model_ids[[1]])
  
  # find the top k models from all grid searches
  k = 10
  glm.mod.list = lapply(1:k, function(i)
  {
    # get grid search (s), the model id (m), and the position of the model (p)
    s = DT.glm.grid$search[i]
    m = DT.glm.grid$model[i]
    p = which(removePunctuation(gsub("[A-z]+", "", unlist(glm.grid.models[[s]]@model_ids))) == m)
    
    # get the model of interest
    output = glm.grid.models[[s]]@model_ids[[p]]
    
    return(output)
  })
  
  # stack the top k models from all grid searches
  h2o.rm("glm.mod")
  glm.mod = h2o.stackedEnsemble(x = x,
                                y = y,
                                training_frame = train.YX.h2o,
                                validation_frame = valid.YX.h2o,
                                model_id = "glm.mod",
                                base_models = glm.mod.list)
  
  # ---- measure model quality ----
  
  # make predictions on each data set
  ynew.train = as.matrix(predict(glm.mod, newdata = train.YX.h2o)[,-1])
  ynew.valid = as.matrix(predict(glm.mod, newdata = valid.YX.h2o)[,-1])
  ynew.test = as.matrix(predict(glm.mod, newdata = test.YX.h2o)[,-1])
  
  # get the true values from each data set
  ytrue.train = as.data.frame(train.YX.h2o[,1])
  ytrue.valid = as.data.frame(valid.YX.h2o[,1])
  ytrue.test = as.data.frame(test.YX.h2o[,1])
  
  # ---- compute multi log loss ----
  
  # build a matrix indicating the true class values for each data set
  ytrue.train.mat = model.matrix(~., data = ytrue.train)[,-1]
  ytrue.train.mat = cbind(1 - rowSums(ytrue.train.mat), ytrue.train.mat)
  
  ytrue.valid.mat = model.matrix(~., data = ytrue.valid)[,-1]
  ytrue.valid.mat = cbind(1 - rowSums(ytrue.valid.mat), ytrue.valid.mat)
  
  ytrue.test.mat = model.matrix(~., data = ytrue.test)[,-1]
  ytrue.test.mat = cbind(1 - rowSums(ytrue.test.mat), ytrue.test.mat)
  
  # compute the multi-class logarithmic loss for each data set
  mll.train = MultiLogLoss(y_pred = ynew.train, y_true = ytrue.train.mat)
  mll.valid = MultiLogLoss(y_pred = ynew.valid, y_true = ytrue.valid.mat)
  mll.test = MultiLogLoss(y_pred = ynew.test, y_true = ytrue.test.mat)
  
  # free up RAM
  gc()
  
  # ---- compute kappa ----
  
  # get the predicted classes and actual classes for each data set
  ynew.train.code = apply(ynew.train, 1, which.max) - 1
  ytrue.train.code = factor(as.numeric(ytrue.train[,1]) - 1, levels = 0:(ncol(ytrue.train.mat) - 1))
  
  ynew.valid.code = apply(ynew.valid, 1, which.max) - 1
  ytrue.valid.code = factor(as.numeric(ytrue.valid[,1]) - 1, levels = 0:(ncol(ytrue.valid.mat) - 1))
  
  ynew.test.code = apply(ynew.test, 1, which.max) - 1
  ytrue.test.code = factor(as.numeric(ytrue.test[,1]) - 1, levels = 0:(ncol(ytrue.test.mat) - 1))
  
  # build a square confusion matrix for each data set
  conf.train = confusion(ytrue = ytrue.train.code, ypred = ynew.train.code)
  conf.valid = confusion(ytrue = ytrue.valid.code, ypred = ynew.valid.code)
  conf.test = confusion(ytrue = ytrue.test.code, ypred = ynew.test.code)
  
  # get the total number of observations for each data set
  n.train = sum(conf.train)
  n.valid = sum(conf.valid)
  n.test = sum(conf.test)
  
  # get the vector of correct predictions for each data set 
  dia.train = diag(conf.train)
  dia.valid = diag(conf.valid)
  dia.test = diag(conf.test)
  
  # get the vector of the number of observations per class for each data set
  rsum.train = rowSums(conf.train)
  rsum.valid = rowSums(conf.valid)
  rsum.test = rowSums(conf.test)
  
  # get the vector of the number of predictions per class for each data set
  csum.train = colSums(conf.train)
  csum.valid = colSums(conf.valid)
  csum.test = colSums(conf.test)
  
  # get the proportion of observations per class for each data set
  p.train = rsum.train / n.train
  p.valid = rsum.valid / n.valid
  p.test = rsum.test / n.test
  
  # get the proportion of predcitions per class for each data set
  q.train = csum.train / n.train
  q.valid = csum.valid / n.valid
  q.test = csum.test / n.test
  
  # compute accuracy for each data set
  acc.train = sum(dia.train) / n.train
  acc.valid = sum(dia.valid) / n.valid
  acc.test = sum(dia.test) / n.test
  
  # compute expected accuracy for each data set
  exp.acc.train = sum(p.train * q.train)
  exp.acc.valid = sum(p.valid * q.valid)
  exp.acc.test = sum(p.test * q.test)
  
  # compute kappa for each data set
  kap.train = (acc.train - exp.acc.train) / (1 - exp.acc.train)
  kap.valid = (acc.valid - exp.acc.valid) / (1 - exp.acc.valid)
  kap.test = (acc.test - exp.acc.test) / (1 - exp.acc.test)
  
  # ---- compute one-vs-all metrics ----
  
  # compute a binary confusion matrix for each class, for each data set
  one.v.all.train = lapply(1:nrow(conf.train), function(i)
  {
    # get the four entries of a binary confusion matrix
    v = c(conf.train[i,i], 
          rsum.train[i] - conf.train[i,i], 
          csum.train[i] - conf.train[i,i], 
          n.train - rsum.train[i] - csum.train[i] + conf.train[i,i]);
    
    # build the confusion matrix
    return(matrix(v, nrow = 2, byrow = TRUE))
  })
  
  one.v.all.valid = lapply(1:nrow(conf.valid), function(i)
  {
    # get the four entries of a binary confusion matrix
    v = c(conf.valid[i,i], 
          rsum.valid[i] - conf.valid[i,i], 
          csum.valid[i] - conf.valid[i,i], 
          n.valid - rsum.valid[i] - csum.valid[i] + conf.valid[i,i]);
    
    # build the confusion matrix
    return(matrix(v, nrow = 2, byrow = TRUE))
  })
  
  one.v.all.test = lapply(1:nrow(conf.test), function(i)
  {
    # get the four entries of a binary confusion matrix
    v = c(conf.test[i,i], 
          rsum.test[i] - conf.test[i,i], 
          csum.test[i] - conf.test[i,i], 
          n.test - rsum.test[i] - csum.test[i] + conf.test[i,i]);
    
    # build the confusion matrix
    return(matrix(v, nrow = 2, byrow = TRUE))
  })
  
  # sum up all of the matrices for each data set
  one.v.all.train = Reduce('+', one.v.all.train)
  one.v.all.valid = Reduce('+', one.v.all.valid)
  one.v.all.test = Reduce('+', one.v.all.test)
  
  # compute the micro average accuracy for each data set
  micro.acc.train = sum(diag(one.v.all.train)) / sum(one.v.all.train)
  micro.acc.valid = sum(diag(one.v.all.valid)) / sum(one.v.all.valid)
  micro.acc.test = sum(diag(one.v.all.test)) / sum(one.v.all.test)
  
  # get the macro accuracy for each data set
  macro.acc.train = acc.train
  macro.acc.valid = acc.valid
  macro.acc.test = acc.test
  
  # ---- finalize output ----
  
  # build a final metrics table for glm.mod
  glm.mod.table = data.table(Model = "Regression",
                             Group = g,
                             Metric = c("Log_Loss", "Kappa", "Macro_Accuracy", "Micro_Accuracy"),
                             Train = c(mll.train, kap.train, macro.acc.train, micro.acc.train),
                             Valid = c(mll.valid, kap.valid, macro.acc.valid, micro.acc.valid),
                             Test = c(mll.test, kap.test, macro.acc.test, micro.acc.test))
  
  # build a final grid search table
  glm.grid.search = data.table(cbind(Model = rep("Regression", nrow(DT.glm.grid)), 
                                     Group = rep(g, nrow(DT.glm.grid)), 
                                     DT.glm.grid))
  
  # update the row and column names of the confusion matrices
  colnames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
  rownames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
  
  colnames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
  rownames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
  
  colnames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
  rownames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
  
  # build a list of confusion matrices
  glm.confusion = list("Train" = conf.train, 
                       "Valid" = conf.valid, 
                       "Test" = conf.test)
  
  # build a list of final tables
  glm.list = list(glm.mod.table, glm.grid.search, glm.confusion)
  names(glm.list) = paste0(c("Regression_Metrics_Group_", "Regression_Search_Group_", "Regression_Error_Group_"), g)
  
  # remove some objects
  rm(dat.g, DT.glm.grid, folds, n,
     glm.grid, glm.grid.models, glm.gridA, glm.gridB, glm.gridC, glm.gridD, glm.hyper.params, 
     dia.train, ynew.train.code, ytrue.train.code, ytrue.train.mat, 
     dia.valid, ynew.valid.code, ytrue.valid.code, ytrue.valid.mat, 
     dia.test, ynew.test.code, ytrue.test.code, ytrue.test.mat, 
     p.train, q.train, acc.train, exp.acc.train,
     p.valid, q.valid, acc.valid, exp.acc.valid,
     p.test, q.test, acc.test, exp.acc.test,
     conf.train, rsum.train, csum.train, n.train, 
     conf.valid, rsum.valid, csum.valid, n.valid, 
     conf.test, rsum.test, csum.test, n.test, 
     one.v.all.train, one.v.all.valid, one.v.all.test,
     mll.train, kap.train, macro.acc.train, micro.acc.train,
     mll.valid, kap.valid, macro.acc.valid, micro.acc.valid,
     mll.test, kap.test, macro.acc.test, micro.acc.test,
     glm.mod, glm.random.gridA, glm.random.gridB, glm.random.gridC, glm.random.gridD, glm.search.criteria,
     plot.glm.grid, train, test, valid, 
     test.YX.h2o, train.YX.h2o, two.way, valid.YX.h2o,
     x, y, YX.h2o, ynew.test, ynew.train, ynew.valid, ytrue.test, ytrue.train, ytrue.valid)
  
  # free up RAM
  gc()
  
  # return the final tables
  return(glm.list)
}

# name the results
names(predict.total_return_class.reg) = "h2o.glm"

# combine predict.total_return_class into one table
predict.total_return_class = append(predict.total_return_class, predict.total_return_class.reg)

# remove some objects
rm(glm.mod.table, glm.grid.search, glm.list, g, predict.total_return_class.reg)


# clean up the data in the h2o cluster
h2o.removeAll()

# shut down the h2o cluster to free up RAM
h2o.shutdown(prompt = FALSE)

# free up RAM
gc()

}

# -----------------------------------------------------------------------------------
# ---- Build Neural Network Models --------------------------------------------------
# -----------------------------------------------------------------------------------

{

# initialize the h2o instance
h2o.init(nthreads = workers, max_mem_size = "8g")

# remove any objects in the h2o instance
h2o.removeAll()

# remove the progress bar when model building
h2o.no_progress()

# ---- build prediction models for pmt_status ---------------------------------------

# determine how many data sets we have to build models for
if(cluster.data)
{
  # get the sub groups of dat
  groups = sort(unique(dat$Cluster))
  
  # also consider using all of dat
  groups = c("All", groups)
  
} else
{
  groups = "All"
}

# determine Neural Network predictive performance for each group
predict.pmt_status.nnet = foreach(g = groups) %do%
{
  # ---- split the data ----
  
  # determine which data set to use
  if(cluster.data)
  {
    if(g == "All")
    {
      # update dat to not contain Cluster
      dat.g = data.table(dat[, !"Cluster"])
      
    } else
    {
      # update dat to only contain data in Cluster g
      dat.g = data.table(dat[Cluster == g, !"Cluster"])
    }
  } else
  {
    # get dat
    dat.g = data.table(dat)
  }
  
  # identify predictors (x)
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # consider two-way interactions between the first n indicators in dat.g
  # recall that the indicators are ordered by importance from random forests
  n = floor(length(x) / 3)
  two.way = data.table(dat.g[, x[1:n], with = FALSE])
  
  # should we add two way interactions in the model?
  add.interactions = FALSE
  
  # build interactions if desired
  if(add.interactions)
  {
    # create the two way interactions and add them to dat.g
    two.way = data.table(model.matrix(~.^2, data = two.way)[,-(1:(n + 1))])
    dat.g = cbind(dat.g, two.way)
  }
  
  # identify predictors (x) and response (y)
  y = "pmt_status"
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # build the fold assignment
  set.seed(42)
  k.folds = max(c(k.folds, 3))
  folds = createFolds(y = unname(unlist(dat.g[, y, with = FALSE])), k = k.folds)
  
  # split up dat.g into train, valid, and test
  train.rows = unname(unlist(lapply(1:(k.folds - 2), function(f) folds[[f]])))
  train = data.table(dat.g[train.rows])
  
  valid.rows = unname(unlist(folds[[k.folds - 1]]))
  valid = data.table(dat.g[valid.rows])
  
  test.rows = unname(unlist(folds[[k.folds]]))
  test = data.table(dat.g[test.rows])
  
  # split up YX.h2o into train, valid, and test
  train.YX.h2o = as.h2o(train[, c(y, x), with = FALSE])
  valid.YX.h2o = as.h2o(valid[, c(y, x), with = FALSE])
  test.YX.h2o = as.h2o(test[, c(y, x), with = FALSE])
  
  # compute the max class weight
  max.class.weight = table(unname(unlist(dat.g[, y, with = FALSE])))
  max.class.weight = as.numeric(max(max(max.class.weight) / max.class.weight))
  
  # compute small, medium, and large hidden layer sizes relative to the input layer
  small = round(length(x) * (1 + 0.5), 0)
  medium = round(length(x) * (1 + 0.5)^2, 0)
  large = round(length(x) * (1 + 0.5)^3, 0)
  
  # ---- grid search for models ----
  
  # set up hyperparameters of interest, with drop out ratios
  nnet.hyper.params = list(adaptive_rate = TRUE,
                           
                           # -- initial tuning --
                           # epochs = c(10, 1e2, 1e3, 1e4),
                           hidden = list(small, medium, large, c(small, small), c(medium, medium), c(large, large), 
                                         c(small, small, small), c(medium, medium, medium), c(large, large, large)),
                           # activation = c("RectifierWithDropout", "TanhWithDropout"),
                           # input_dropout_ratio = c(0, 0.1, 0.2),
                           # l1 = c(0, 1e-3, 1e-5),
                           # l2 = c(0, 1e-3, 1e-5),
                           rho = c(0.9, 0.95, 0.99, 0.999),
                           # epsilon = c(1e-10, 1e-8),
                           
                           # -- re-fine tuning --
                           epochs = 10,
                           # hidden = list(small, medium, large),
                           activation = "RectifierWithDropout",
                           input_dropout_ratio = 0,
                           l1 = c(0, 1e-5),
                           l2 = c(0, 1e-5),
                           # rho = 0.99,
                           epsilon = 1e-8,
                           
                           # -- scoring function -- 
                           stopping_metric = "mean_per_class_error")
  
  # lets use a random grid search and specify a time limit and/or model limit
  minutes = 20
  nnet.search.criteria = list(strategy = "RandomDiscrete", 
                              max_runtime_secs = minutes * 60, 
                              # max_models = 100, 
                              seed = 42)
  
  # lets run a grid search for a good model, without drop out ratios
  h2o.rm("nnet.random.grid")
  nnet.random.grid = h2o.grid(algorithm = "deeplearning",
                              grid_id = "nnet.random.grid",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              mini_batch_size = 10,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              seed = 21,
                              variable_importances = FALSE,
                              balance_classes = TRUE,
                              max_after_balance_size = max.class.weight,
                              hyper_params = nnet.hyper.params,
                              search_criteria = nnet.search.criteria)
  
  # free up RAM
  gc()
  
  # rank each model in the random grids
  nnet.grid = h2o.getGrid("nnet.random.grid", sort_by = "auc", decreasing = TRUE)
  
  # get the summary table of the grid search
  DT.nnet.grid = data.table(nnet.grid@summary_table)
  DT.nnet.grid
  
  # set up the data types for each column in DT.grid for plotting purposes
  DT.nnet.grid = DT.nnet.grid[, .(activation = as.factor(activation),
                                  epochs = as.numeric(epochs),
                                  hidden = as.factor(hidden),
                                  input_dropout_ratio = as.factor(input_dropout_ratio),
                                  rho = as.factor(rho),
                                  epsilon = as.factor(epsilon),
                                  l1 = as.factor(l1),
                                  l2 = as.factor(l2),
                                  model = removePunctuation(gsub("[A-z]+", "", model_ids)),
                                  auc = as.numeric(auc))]
  
  # plot auc v. hidden and activation to see which structure is most robust
  plot.nnet.grid = ggplot(DT.nnet.grid, aes(x = hidden, y = auc, color = activation, fill = activation)) + 
    # geom_boxplot() + 
    geom_jitter(size = 3, alpha = 2/3) + 
    # scale_y_continuous(labels = dollar) + 
    ggtitle("Cross Validation Error") + 
    labs(x = "Hidden Layer Structure", y = "AUC", color = "Activation Function", fill = "Activation Function") + 
    facet_wrap(~paste("Learning Rate:", rho), nrow = 2) +
    theme_bw(base_size = 25) +
    theme(legend.position = "top", 
          legend.key.size = unit(.25, "in"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
  
  plot.nnet.grid
  
  # pick the top model from all grid searches
  # nnet.mod = h2o.getModel(nnet.grid@model_ids[[1]])
  
  # find the top k models from all grid searches
  k = 10
  nnet.mod.list = lapply(1:k, function(i) nnet.grid@model_ids[[i]])
  
  # stack the top k models from all grid searches
  h2o.rm("nnet.mod")
  nnet.mod = h2o.stackedEnsemble(x = x,
                                y = y,
                                training_frame = train.YX.h2o,
                                validation_frame = valid.YX.h2o,
                                model_id = "nnet.mod",
                                base_models = nnet.mod.list)
  
  # plot variable importance
  # h2o.varimp_plot(nnet.mod)
  
  # ---- measure model quality ----
  
  # make predictions on each data set
  ynew.train = as.data.frame(predict(nnet.mod, newdata = train.YX.h2o))$p1
  ynew.valid = as.data.frame(predict(nnet.mod, newdata = valid.YX.h2o))$p1
  ynew.test = as.data.frame(predict(nnet.mod, newdata = test.YX.h2o))$p1
  
  # get the true values from each data set
  ytrue.train = as.data.frame(train.YX.h2o)[,1]
  ytrue.valid = as.data.frame(valid.YX.h2o)[,1]
  ytrue.test = as.data.frame(test.YX.h2o)[,1]
  
  # plot the probability predictions on the testing data set
  # densityPlot(ynew.test)
  
  # compute prediction metrics on each data set
  nnet.metrics.train = h2o.make_metrics(predicted = as.h2o(ynew.train), 
                                       actuals = as.h2o(factor(ytrue.train, levels = sort(unique(ytrue.train)))))
  
  nnet.metrics.valid = h2o.make_metrics(predicted = as.h2o(ynew.valid), 
                                       actuals = as.h2o(factor(ytrue.valid, levels = sort(unique(ytrue.valid)))))
  
  nnet.metrics.test = h2o.make_metrics(predicted = as.h2o(ynew.test), 
                                      actuals = as.h2o(factor(ytrue.test, levels = sort(unique(ytrue.test)))))
  
  # free up RAM
  gc()
  
  # check out the threshold metrics
  nnet.metrics.train@metrics$max_criteria_and_metric_scores
  nnet.metrics.valid@metrics$max_criteria_and_metric_scores
  nnet.metrics.test@metrics$max_criteria_and_metric_scores
  
  # pick a threshold for converting probabilities into binary
  threshold = mean(c(tail(nnet.metrics.train@metrics$max_criteria_and_metric_scores$threshold, 1),
                     tail(nnet.metrics.valid@metrics$max_criteria_and_metric_scores$threshold, 1),
                     tail(nnet.metrics.test@metrics$max_criteria_and_metric_scores$threshold, 1)))
  
  # compute a confusion matrix for each data set
  nnet.confusion.train = confusion(ytrue = ytrue.train, ypred = as.numeric(ynew.train >= threshold))
  nnet.confusion.valid = confusion(ytrue = ytrue.valid, ypred = as.numeric(ynew.valid >= threshold))
  nnet.confusion.test = confusion(ytrue = ytrue.test, ypred = as.numeric(ynew.test >= threshold))
  
  # compute Sensitivity (accuracy in predicting 1's) for each data set
  nnet.sensitivity.train = nnet.confusion.train[2,2] / (nnet.confusion.train[2,2] + nnet.confusion.train[2,1])
  nnet.sensitivity.valid = nnet.confusion.valid[2,2] / (nnet.confusion.valid[2,2] + nnet.confusion.valid[2,1])
  nnet.sensitivity.test = nnet.confusion.test[2,2] / (nnet.confusion.test[2,2] + nnet.confusion.test[2,1])
  
  # compute Specificity (accuracy in predicting 0's) for each data set
  nnet.specificity.train = nnet.confusion.train[1,1] / (nnet.confusion.train[1,1] + nnet.confusion.train[1,2])
  nnet.specificity.valid = nnet.confusion.valid[1,1] / (nnet.confusion.valid[1,1] + nnet.confusion.valid[1,2])
  nnet.specificity.test = nnet.confusion.test[1,1] / (nnet.confusion.test[1,1] + nnet.confusion.test[1,2])
  
  # compute Odds Ratio (the odds of success over failure) for each data set
  nnet.odds.ratio.train = (nnet.confusion.train[1,1] * nnet.confusion.train[2,2]) / (nnet.confusion.train[2,1] * nnet.confusion.train[1,2])
  nnet.odds.ratio.valid = (nnet.confusion.valid[1,1] * nnet.confusion.valid[2,2]) / (nnet.confusion.valid[2,1] * nnet.confusion.valid[1,2])
  nnet.odds.ratio.test = (nnet.confusion.test[1,1] * nnet.confusion.test[2,2]) / (nnet.confusion.test[2,1] * nnet.confusion.test[1,2])
  
  # compute Accuracy for each data set
  nnet.accuracy.train = (nnet.confusion.train[1,1] + nnet.confusion.train[2,2]) / (nnet.confusion.train[1,1] + nnet.confusion.train[1,2] + nnet.confusion.train[2,1] + nnet.confusion.train[2,2])
  nnet.accuracy.valid = (nnet.confusion.valid[1,1] + nnet.confusion.valid[2,2]) / (nnet.confusion.valid[1,1] + nnet.confusion.valid[1,2] + nnet.confusion.valid[2,1] + nnet.confusion.valid[2,2])
  nnet.accuracy.test = (nnet.confusion.test[1,1] + nnet.confusion.test[2,2]) / (nnet.confusion.test[1,1] + nnet.confusion.test[1,2] + nnet.confusion.test[2,1] + nnet.confusion.test[2,2])
  
  # compute AUC for each data set
  nnet.auc.train = h2o.auc(nnet.metrics.train)
  nnet.auc.valid = h2o.auc(nnet.metrics.valid)
  nnet.auc.test = h2o.auc(nnet.metrics.test)
  
  # compute Log Loss for each data set
  nnet.logloss.train = h2o.logloss(nnet.metrics.train)
  nnet.logloss.valid = h2o.logloss(nnet.metrics.valid)
  nnet.logloss.test = h2o.logloss(nnet.metrics.test)
  
  # ---- compute kappa ----
  
  # get the total number of observations for each data set
  n.train = sum(nnet.confusion.train)
  n.valid = sum(nnet.confusion.valid)
  n.test = sum(nnet.confusion.test)
  
  # get the vector of correct predictions for each data set 
  dia.train = diag(nnet.confusion.train)
  dia.valid = diag(nnet.confusion.valid)
  dia.test = diag(nnet.confusion.test)
  
  # get the vector of the number of observations per class for each data set
  rsum.train = rowSums(nnet.confusion.train)
  rsum.valid = rowSums(nnet.confusion.valid)
  rsum.test = rowSums(nnet.confusion.test)
  
  # get the vector of the number of predictions per class for each data set
  csum.train = colSums(nnet.confusion.train)
  csum.valid = colSums(nnet.confusion.valid)
  csum.test = colSums(nnet.confusion.test)
  
  # get the proportion of observations per class for each data set
  p.train = rsum.train / n.train
  p.valid = rsum.valid / n.valid
  p.test = rsum.test / n.test
  
  # get the proportion of predcitions per class for each data set
  q.train = csum.train / n.train
  q.valid = csum.valid / n.valid
  q.test = csum.test / n.test
  
  # compute accuracy for each data set
  acc.train = sum(dia.train) / n.train
  acc.valid = sum(dia.valid) / n.valid
  acc.test = sum(dia.test) / n.test
  
  # compute expected accuracy for each data set
  exp.acc.train = sum(p.train * q.train)
  exp.acc.valid = sum(p.valid * q.valid)
  exp.acc.test = sum(p.test * q.test)
  
  # compute kappa for each data set
  kap.train = (acc.train - exp.acc.train) / (1 - exp.acc.train)
  kap.valid = (acc.valid - exp.acc.valid) / (1 - exp.acc.valid)
  kap.test = (acc.test - exp.acc.test) / (1 - exp.acc.test)
  
  # ---- finalize output ----
  
  # build a final metrics table for nnet.mod
  nnet.mod.table = data.table(Model = "Neural Network",
                             Group = g,
                             Threshold = threshold,
                             Metric = c("AUC", "Accuracy", "Log_Loss", "Kappa", "Odds_Ratio", "Specificity_0", "Sensitivity_1"),
                             Train = c(nnet.auc.train, nnet.accuracy.train, nnet.logloss.train, kap.train, nnet.odds.ratio.train, nnet.specificity.train, nnet.sensitivity.train),
                             Valid = c(nnet.auc.valid, nnet.accuracy.valid, nnet.logloss.valid, kap.valid, nnet.odds.ratio.valid, nnet.specificity.valid, nnet.sensitivity.valid),
                             Test = c(nnet.auc.test, nnet.accuracy.test, nnet.logloss.test, kap.test, nnet.odds.ratio.test, nnet.specificity.test, nnet.sensitivity.test))
  
  # build a final grid search table
  nnet.grid.search = data.table(cbind(Model = rep("Neural Network", nrow(DT.nnet.grid)), 
                                     Group = rep(g, nrow(DT.nnet.grid)), 
                                     DT.nnet.grid))
  
  # build a list of confusion matrices
  nnet.confusion = list("Train" = nnet.confusion.train, 
                        "Valid" = nnet.confusion.valid, 
                        "Test" = nnet.confusion.test)
  
  # build a list of final tables
  nnet.list = list(nnet.mod.table, nnet.grid.search, nnet.confusion)
  names(nnet.list) = paste0(c("Neural_Network_Metrics_Group_", "Neural_Network_Search_Group_", "Neural_Network_Error_Group_"), g)
  
  # remove some objects
  rm(dat.g, DT.nnet.grid, folds, nnet.accuracy.test, nnet.accuracy.train, nnet.accuracy.valid,
     nnet.auc.test, nnet.auc.train, nnet.auc.valid, nnet.confusion.test, nnet.confusion.train, nnet.confusion.valid,
     nnet.grid, nnet.hyper.params, valid.rows, nnet.confusion, nnet.mod.table, nnet.grid.search,
     nnet.logloss.test, nnet.logloss.train, nnet.logloss.valid, nnet.metrics.test, nnet.metrics.train, nnet.metrics.valid,
     nnet.mod, nnet.odds.ratio.test, nnet.odds.ratio.train, nnet.odds.ratio.valid,
     nnet.search.criteria, nnet.sensitivity.test, nnet.mod.list,
     nnet.sensitivity.train, nnet.sensitivity.valid, nnet.specificity.test, nnet.specificity.train,
     nnet.specificity.valid, max.class.weight, plot.nnet.grid, train, test, valid, 
     small, test.rows, train.rows, k,
     test.YX.h2o, threshold, train.YX.h2o, valid.YX.h2o, nnet.random.grid,
     x, y, ynew.test, ynew.train, ynew.valid, ytrue.test, ytrue.train, ytrue.valid,
     n.train, dia.train, rsum.train, csum.train, p.train, q.train, acc.train, exp.acc.train, kap.train,
     n.valid, dia.valid, rsum.valid, csum.valid, p.valid, q.valid, acc.valid, exp.acc.valid, kap.valid,
     n.test, dia.test, rsum.test, csum.test, p.test, q.test, acc.test, exp.acc.test, kap.test)
  
  # free up RAM
  gc()
  
  # return the final metrics table
  return(nnet.list)
}

# name the results
names(predict.pmt_status.nnet) = "h2o.deeplearning"

# combine predict.pmt_status into one table
predict.pmt_status = append(predict.pmt_status, predict.pmt_status.nnet)

# remove some objects
rm(nnet.list, g, predict.pmt_status.nnet)

# free up RAM
gc()

# ---- build prediction models for total_return ---------------------------------------

# determine how many data sets we have to build models for
if(cluster.data)
{
  # get the sub groups of dat
  groups = sort(unique(dat$Cluster))
  
  # also consider using all of dat
  groups = c("All", groups)
  
} else
{
  groups = "All"
}

# determine Neural Network predictive performance for each group
predict.total_return.nnet = foreach(g = groups) %do%
{
  # ---- split the data ----
  
  # determine which data set to use
  if(cluster.data)
  {
    if(g == "All")
    {
      # update dat to not contain Cluster
      dat.g = data.table(dat[, !"Cluster"])
      
    } else
    {
      # update dat to only contain data in Cluster g
      dat.g = data.table(dat[Cluster == g, !"Cluster"])
    }
  } else
  {
    # get dat
    dat.g = data.table(dat)
  }
  
  # identify predictors (x)
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # consider two-way interactions between the first n indicators in dat.g
  # recall that the indicators are ordered by importance from random forests
  n = floor(length(x) / 3)
  two.way = data.table(dat.g[, x[1:n], with = FALSE])
  
  # should we add two way interactions in the model?
  add.interactions = FALSE
  
  # build interactions if desired
  if(add.interactions)
  {
    # create the two way interactions and add them to dat.g
    two.way = data.table(model.matrix(~.^2, data = two.way)[,-(1:(n + 1))])
    dat.g = cbind(dat.g, two.way)
  }
  
  # identify predictors (x) and response (y)
  y = "total_return"
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # build the fold assignment
  set.seed(42)
  k.folds = max(c(k.folds, 3))
  folds = createFolds(y = unname(unlist(dat.g[, y, with = FALSE])), k = k.folds)
  
  # split up dat.g into train, valid, and test
  train.rows = unname(unlist(lapply(1:(k.folds - 2), function(f) folds[[f]])))
  train = data.table(dat.g[train.rows])
  
  valid.rows = unname(unlist(folds[[k.folds - 1]]))
  valid = data.table(dat.g[valid.rows])
  
  test.rows = unname(unlist(folds[[k.folds]]))
  test = data.table(dat.g[test.rows])
  
  # split up YX.h2o into train, valid, and test
  train.YX.h2o = as.h2o(train[, c(y, x), with = FALSE])
  valid.YX.h2o = as.h2o(valid[, c(y, x), with = FALSE])
  test.YX.h2o = as.h2o(test[, c(y, x), with = FALSE])
  
  # compute small, medium, and large hidden layer sizes relative to the input layer
  small = round(length(x) * (1 + 0.5), 0)
  medium = round(length(x) * (1 + 0.5)^2, 0)
  large = round(length(x) * (1 + 0.5)^3, 0)
  
  # ---- grid search for models ----
  
  # set up hyperparameters of interest, with drop out ratios
  nnet.hyper.params = list(adaptive_rate = TRUE,
                           
                           # -- initial tuning --
                           # epochs = c(10, 1e2, 1e3, 1e4),
                           hidden = list(small, medium, large, c(small, small), c(medium, medium), c(large, large), 
                                         c(small, small, small), c(medium, medium, medium), c(large, large, large)),
                           # activation = c("RectifierWithDropout", "TanhWithDropout"),
                           # input_dropout_ratio = c(0, 0.1, 0.2),
                           # l1 = c(0, 1e-3, 1e-5),
                           # l2 = c(0, 1e-3, 1e-5),
                           rho = c(0.9, 0.95, 0.99, 0.999),
                           # epsilon = c(1e-10, 1e-8),
                           
                           # -- re-fine tuning --
                           epochs = 10,
                           # hidden = list(small, medium, large),
                           activation = "RectifierWithDropout",
                           input_dropout_ratio = 0,
                           l1 = c(0, 1e-5),
                           l2 = c(0, 1e-5),
                           # rho = 0.99,
                           epsilon = 1e-8,
                           
                           # -- scoring function -- 
                           stopping_metric = "RMSE")
  
  # lets use a random grid search and specify a time limit and/or model limit
  minutes = 20
  nnet.search.criteria = list(strategy = "RandomDiscrete", 
                              max_runtime_secs = minutes * 60, 
                              # max_models = 100, 
                              seed = 42)
  
  # lets run a grid search for a good model, without drop out ratios
  h2o.rm("nnet.random.grid")
  nnet.random.grid = h2o.grid(algorithm = "deeplearning",
                              grid_id = "nnet.random.grid",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              mini_batch_size = 10,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              seed = 21,
                              variable_importances = FALSE,
                              hyper_params = nnet.hyper.params,
                              search_criteria = nnet.search.criteria)
  
  # free up RAM
  gc()
  
  # rank each model in the random grids
  nnet.grid = h2o.getGrid("nnet.random.grid", sort_by = "rmse", decreasing = FALSE)
  
  # get the summary table of the grid search
  DT.nnet.grid = data.table(nnet.grid@summary_table)
  DT.nnet.grid
  
  # set up the data types for each column in DT.grid for plotting purposes
  DT.nnet.grid = DT.nnet.grid[, .(activation = as.factor(activation),
                                  epochs = as.numeric(epochs),
                                  hidden = as.factor(hidden),
                                  input_dropout_ratio = as.factor(input_dropout_ratio),
                                  rho = as.factor(rho),
                                  epsilon = as.factor(epsilon),
                                  l1 = as.factor(l1),
                                  l2 = as.factor(l2),
                                  model = removePunctuation(gsub("[A-z]+", "", model_ids)),
                                  rmse = as.numeric(rmse))]
  
  # plot rmse v. hidden and activation to see which structure is most robust
  plot.nnet.grid = ggplot(DT.nnet.grid, aes(x = hidden, y = rmse, color = activation, fill = activation)) + 
    # geom_boxplot() + 
    geom_jitter(size = 3, alpha = 2/3) + 
    # scale_y_continuous(labels = dollar) + 
    ggtitle("Cross Validation Error") + 
    labs(x = "Hidden Layer Structure", y = "RMSE", color = "Activation Function", fill = "Activation Function") + 
    facet_wrap(~paste("Learning Rate:", rho), nrow = 2) +
    theme_bw(base_size = 25) +
    theme(legend.position = "top", 
          legend.key.size = unit(.25, "in"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
  
  plot.nnet.grid
  
  # pick the top model from all grid searches
  # nnet.mod = h2o.getModel(nnet.grid@model_ids[[1]])
  
  # find the top k models from all grid searches
  k = 10
  nnet.mod.list = lapply(1:k, function(i) nnet.grid@model_ids[[i]])
  
  # stack the top k models from all grid searches
  h2o.rm("nnet.mod")
  nnet.mod = h2o.stackedEnsemble(x = x,
                                 y = y,
                                 training_frame = train.YX.h2o,
                                 validation_frame = valid.YX.h2o,
                                 model_id = "nnet.mod",
                                 base_models = nnet.mod.list)
  
  # plot variable importance
  # h2o.varimp_plot(nnet.mod)
  
  # ---- measure model quality ----
  
  # make predictions on each data set
  ynew.train = as.data.frame(predict(nnet.mod, newdata = train.YX.h2o))$predict
  ynew.valid = as.data.frame(predict(nnet.mod, newdata = valid.YX.h2o))$predict
  ynew.test = as.data.frame(predict(nnet.mod, newdata = test.YX.h2o))$predict
  
  # get the true values from each data set
  ytrue.train = as.data.frame(train.YX.h2o)[,1]
  ytrue.valid = as.data.frame(valid.YX.h2o)[,1]
  ytrue.test = as.data.frame(test.YX.h2o)[,1]
  
  # transform the predicted and actual values back to original units if needed
  if(use.norm)
  {
    # map predicted values to the original units of the response variable
    ynew.train = predict(best.norm, newdata = ynew.train, inverse = TRUE)
    ynew.valid = predict(best.norm, newdata = ynew.valid, inverse = TRUE)
    ynew.test = predict(best.norm, newdata = ynew.test, inverse = TRUE)
    
    # map actual values to the original units of the response variable
    ytrue.train = predict(best.norm, newdata = ytrue.train, inverse = TRUE)
    ytrue.valid = predict(best.norm, newdata = ytrue.valid, inverse = TRUE)
    ytrue.test = predict(best.norm, newdata = ytrue.test, inverse = TRUE)
  }
  
  # compute prediction metrics on each data set
  nnet.metrics.train = h2o.make_metrics(predicted = as.h2o(ynew.train), 
                                       actuals = as.h2o(ytrue.train))
  
  nnet.metrics.valid = h2o.make_metrics(predicted = as.h2o(ynew.valid), 
                                       actuals = as.h2o(ytrue.valid))
  
  nnet.metrics.test = h2o.make_metrics(predicted = as.h2o(ynew.test), 
                                      actuals = as.h2o(ytrue.test))
  
  # free up RAM
  gc()
  
  # check out the prediction metrics
  nnet.metrics.train
  nnet.metrics.valid
  nnet.metrics.test
  
  # fit the normal distribution to the residuals
  nnet.resid.train = fitdist(data = ytrue.train - ynew.train, distr = "norm")
  nnet.resid.valid = fitdist(data = ytrue.valid - ynew.valid, distr = "norm") 
  nnet.resid.test = fitdist(data = ytrue.test - ynew.test, distr = "norm") 
  
  # ---- finalize output ----
  
  # build a final metrics table for nnet.mod
  nnet.mod.table = data.table(Model = "Neural Network",
                             Group = g,
                             Metric = c("R2", "RMSE"),
                             Train = c(nnet.metrics.train@metrics$r2, nnet.metrics.train@metrics$RMSE),
                             Valid = c(nnet.metrics.valid@metrics$r2, nnet.metrics.valid@metrics$RMSE),
                             Test = c(nnet.metrics.test@metrics$r2, nnet.metrics.test@metrics$RMSE))
  
  # build a final grid search table
  nnet.grid.search = data.table(cbind(Model = rep("Neural Network", nrow(DT.nnet.grid)), 
                                     Group = rep(g, nrow(DT.nnet.grid)), 
                                     DT.nnet.grid))
  
  # build a list of residual plots
  nnet.resid = list("Train" = nnet.resid.train, 
                   "Valid" = nnet.resid.valid, 
                   "Test" = nnet.resid.test)
  
  # build a list of final tables
  nnet.list = list(nnet.mod.table, nnet.grid.search, nnet.resid)
  names(nnet.list) = paste0(c("Neural_Network_Metrics_Group_", "Neural_Network_Search_Group_", "Neural_Network_Error_Group_"), g)
  
  # remove some objects
  rm(dat.g, DT.nnet.grid, folds, n,
     nnet.grid, nnet.hyper.params, 
     nnet.metrics.test, nnet.metrics.train, nnet.metrics.valid,
     nnet.mod, nnet.resid.test, nnet.resid.train, nnet.resid.valid,
     nnet.search.criteria, small, medium, large, k,
     plot.nnet.grid, train, test, valid, 
     nnet.random.grid, nnet.mod.table, nnet.grid.search, nnet.resid,
     test.YX.h2o, train.YX.h2o, two.way, valid.YX.h2o,
     x, y, ynew.test, ynew.train, ynew.valid, ytrue.test, ytrue.train, ytrue.valid)
  
  # free up RAM
  gc()
  
  # return the final metrics table
  return(nnet.list)
}

# name the results
names(predict.total_return.nnet) = "h2o.deeplearning"

# combine predict.total_return into one table
predict.total_return = append(predict.total_return, predict.total_return.nnet)

# remove some objects
rm(nnet.list, g, predict.total_return.nnet)

# free up RAM
gc()

# ---- build prediction models for total_return_class ---------------------------------------

# determine how many data sets we have to build models for
if(cluster.data)
{
  # get the sub groups of dat
  groups = sort(unique(dat$Cluster))
  
  # also consider using all of dat
  groups = c("All", groups)
  
} else
{
  groups = "All"
}

# determine Neural Network predictive performance for each group
predict.total_return_class.nnet = foreach(g = groups) %do%
{
  # ---- split the data ----
  
  # determine which data set to use
  if(cluster.data)
  {
    if(g == "All")
    {
      # update dat to not contain Cluster
      dat.g = data.table(dat[, !"Cluster"])
      
    } else
    {
      # update dat to only contain data in Cluster g
      dat.g = data.table(dat[Cluster == g, !"Cluster"])
    }
  } else
  {
    # get dat
    dat.g = data.table(dat)
  }
  
  # identify predictors (x)
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # consider two-way interactions between the first n indicators in dat.g
  # recall that the indicators are ordered by importance from random forests
  n = floor(length(x) / 3)
  two.way = data.table(dat.g[, x[1:n], with = FALSE])
  
  # should we add two way interactions in the model?
  add.interactions = FALSE
  
  # build interactions if desired
  if(add.interactions)
  {
    # create the two way interactions and add them to dat.g
    two.way = data.table(model.matrix(~.^2, data = two.way)[,-(1:(n + 1))])
    dat.g = cbind(dat.g, two.way)
  }
  
  # identify predictors (x) and response (y)
  y = "total_return_class"
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # build the fold assignment
  set.seed(42)
  k.folds = max(c(k.folds, 3))
  folds = createFolds(y = unname(unlist(dat.g[, y, with = FALSE])), k = k.folds)
  
  # split up dat.g into train, valid, and test
  train.rows = unname(unlist(lapply(1:(k.folds - 2), function(f) folds[[f]])))
  train = data.table(dat.g[train.rows])
  
  valid.rows = unname(unlist(folds[[k.folds - 1]]))
  valid = data.table(dat.g[valid.rows])
  
  test.rows = unname(unlist(folds[[k.folds]]))
  test = data.table(dat.g[test.rows])
  
  # split up YX.h2o into train, valid, and test
  train.YX.h2o = as.h2o(train[, c(y, x), with = FALSE])
  valid.YX.h2o = as.h2o(valid[, c(y, x), with = FALSE])
  test.YX.h2o = as.h2o(test[, c(y, x), with = FALSE])
  
  # compute the max class weight
  max.class.weight = table(unname(unlist(dat.g[, y, with = FALSE])))
  max.class.weight = as.numeric(max(max(max.class.weight) / max.class.weight))
  
  # compute small, medium, and large hidden layer sizes relative to the input layer
  small = round(length(x) * (1 + 0.5), 0)
  medium = round(length(x) * (1 + 0.5)^2, 0)
  large = round(length(x) * (1 + 0.5)^3, 0)
  
  # ---- grid search for models ----
  
  # set up hyperparameters of interest, with drop out ratios
  nnet.hyper.params = list(adaptive_rate = TRUE,
                           
                           # -- initial tuning --
                           # epochs = c(10, 1e2, 1e3, 1e4),
                           hidden = list(small, medium, large, c(small, small), c(medium, medium), c(large, large), 
                                         c(small, small, small), c(medium, medium, medium), c(large, large, large)),
                           # activation = c("RectifierWithDropout", "TanhWithDropout"),
                           # input_dropout_ratio = c(0, 0.1, 0.2),
                           # l1 = c(0, 1e-3, 1e-5),
                           # l2 = c(0, 1e-3, 1e-5),
                           rho = c(0.9, 0.95, 0.99, 0.999),
                           # epsilon = c(1e-10, 1e-8),
                           
                           # -- re-fine tuning --
                           epochs = 10,
                           # hidden = list(small, medium, large),
                           activation = "RectifierWithDropout",
                           input_dropout_ratio = 0,
                           l1 = c(0, 1e-5),
                           l2 = c(0, 1e-5),
                           # rho = 0.99,
                           epsilon = 1e-8,
                           
                           # -- scoring function -- 
                           stopping_metric = "mean_per_class_error")
  
  # lets use a random grid search and specify a time limit and/or model limit
  minutes = 20
  nnet.search.criteria = list(strategy = "RandomDiscrete", 
                              max_runtime_secs = minutes * 60, 
                              # max_models = 100, 
                              seed = 42)
  
  # lets run a grid search for a good model, without drop out ratios
  h2o.rm("nnet.random.grid")
  nnet.random.grid = h2o.grid(algorithm = "deeplearning",
                              grid_id = "nnet.random.grid",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              mini_batch_size = 10,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              seed = 21,
                              variable_importances = FALSE,
                              balance_classes = TRUE,
                              max_after_balance_size = max.class.weight,
                              hyper_params = nnet.hyper.params,
                              search_criteria = nnet.search.criteria)
  
  # free up RAM
  gc()
  
  # rank each model in the random grids
  nnet.grid = h2o.getGrid("nnet.random.grid", sort_by = "mean_per_class_error", decreasing = FALSE)
  
  # get the summary table of the grid search
  DT.nnet.grid = data.table(nnet.grid@summary_table)
  DT.nnet.grid
  
  # set up the data types for each column in DT.grid for plotting purposes
  DT.nnet.grid = DT.nnet.grid[, .(activation = as.factor(activation),
                                  epochs = as.numeric(epochs),
                                  hidden = as.factor(hidden),
                                  input_dropout_ratio = as.factor(input_dropout_ratio),
                                  rho = as.factor(rho),
                                  epsilon = as.factor(epsilon),
                                  l1 = as.factor(l1),
                                  l2 = as.factor(l2),
                                  model = removePunctuation(gsub("[A-z]+", "", model_ids)),
                                  mean_per_class_error = as.numeric(mean_per_class_error))]
  
  # plot mean_per_class_error v. hidden and activation to see which structure is most robust
  plot.nnet.grid = ggplot(DT.nnet.grid, aes(x = hidden, y = mean_per_class_error, color = activation, fill = activation)) + 
    # geom_boxplot() + 
    geom_jitter(size = 3, alpha = 2/3) + 
    # scale_y_continuous(labels = dollar) + 
    ggtitle("Cross Validation Error") + 
    labs(x = "Hidden Layer Structure", y = "Mean Per Class Error", color = "Activation Function", fill = "Activation Function") + 
    facet_wrap(~paste("Learning Rate:", rho), nrow = 2) +
    theme_bw(base_size = 25) +
    theme(legend.position = "top", 
          legend.key.size = unit(.25, "in"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
  
  plot.nnet.grid
  
  # pick the top model from all grid searches
  # nnet.mod = h2o.getModel(nnet.grid@model_ids[[1]])
  
  # find the top k models from all grid searches
  k = 10
  nnet.mod.list = lapply(1:k, function(i) nnet.grid@model_ids[[i]])
  
  # stack the top k models from all grid searches
  h2o.rm("nnet.mod")
  nnet.mod = h2o.stackedEnsemble(x = x,
                                 y = y,
                                 training_frame = train.YX.h2o,
                                 validation_frame = valid.YX.h2o,
                                 model_id = "nnet.mod",
                                 base_models = nnet.mod.list)
  
  # plot variable importance
  # h2o.varimp_plot(nnet.mod)
  
  # ---- measure model quality ----
  
  # make predictions on each data set
  ynew.train = as.matrix(predict(nnet.mod, newdata = train.YX.h2o)[,-1])
  ynew.valid = as.matrix(predict(nnet.mod, newdata = valid.YX.h2o)[,-1])
  ynew.test = as.matrix(predict(nnet.mod, newdata = test.YX.h2o)[,-1])
  
  # get the true values from each data set
  ytrue.train = as.data.frame(train.YX.h2o[,1])
  ytrue.valid = as.data.frame(valid.YX.h2o[,1])
  ytrue.test = as.data.frame(test.YX.h2o[,1])
  
  # ---- compute multi log loss ----
  
  # build a matrix indicating the true class values for each data set
  ytrue.train.mat = model.matrix(~., data = ytrue.train)[,-1]
  ytrue.train.mat = cbind(1 - rowSums(ytrue.train.mat), ytrue.train.mat)
  
  ytrue.valid.mat = model.matrix(~., data = ytrue.valid)[,-1]
  ytrue.valid.mat = cbind(1 - rowSums(ytrue.valid.mat), ytrue.valid.mat)
  
  ytrue.test.mat = model.matrix(~., data = ytrue.test)[,-1]
  ytrue.test.mat = cbind(1 - rowSums(ytrue.test.mat), ytrue.test.mat)
  
  # compute the multi-class logarithmic loss for each data set
  mll.train = MultiLogLoss(y_pred = ynew.train, y_true = ytrue.train.mat)
  mll.valid = MultiLogLoss(y_pred = ynew.valid, y_true = ytrue.valid.mat)
  mll.test = MultiLogLoss(y_pred = ynew.test, y_true = ytrue.test.mat)
  
  # free up RAM
  gc()
  
  # ---- compute kappa ----
  
  # get the predicted classes and actual classes for each data set
  ynew.train.code = apply(ynew.train, 1, which.max) - 1
  ytrue.train.code = factor(as.numeric(ytrue.train[,1]) - 1, levels = 0:(ncol(ytrue.train.mat) - 1))
  
  ynew.valid.code = apply(ynew.valid, 1, which.max) - 1
  ytrue.valid.code = factor(as.numeric(ytrue.valid[,1]) - 1, levels = 0:(ncol(ytrue.valid.mat) - 1))
  
  ynew.test.code = apply(ynew.test, 1, which.max) - 1
  ytrue.test.code = factor(as.numeric(ytrue.test[,1]) - 1, levels = 0:(ncol(ytrue.test.mat) - 1))
  
  # build a square confusion matrix for each data set
  conf.train = confusion(ytrue = ytrue.train.code, ypred = ynew.train.code)
  conf.valid = confusion(ytrue = ytrue.valid.code, ypred = ynew.valid.code)
  conf.test = confusion(ytrue = ytrue.test.code, ypred = ynew.test.code)
  
  # get the total number of observations for each data set
  n.train = sum(conf.train)
  n.valid = sum(conf.valid)
  n.test = sum(conf.test)
  
  # get the vector of correct predictions for each data set 
  dia.train = diag(conf.train)
  dia.valid = diag(conf.valid)
  dia.test = diag(conf.test)
  
  # get the vector of the number of observations per class for each data set
  rsum.train = rowSums(conf.train)
  rsum.valid = rowSums(conf.valid)
  rsum.test = rowSums(conf.test)
  
  # get the vector of the number of predictions per class for each data set
  csum.train = colSums(conf.train)
  csum.valid = colSums(conf.valid)
  csum.test = colSums(conf.test)
  
  # get the proportion of observations per class for each data set
  p.train = rsum.train / n.train
  p.valid = rsum.valid / n.valid
  p.test = rsum.test / n.test
  
  # get the proportion of predcitions per class for each data set
  q.train = csum.train / n.train
  q.valid = csum.valid / n.valid
  q.test = csum.test / n.test
  
  # compute accuracy for each data set
  acc.train = sum(dia.train) / n.train
  acc.valid = sum(dia.valid) / n.valid
  acc.test = sum(dia.test) / n.test
  
  # compute expected accuracy for each data set
  exp.acc.train = sum(p.train * q.train)
  exp.acc.valid = sum(p.valid * q.valid)
  exp.acc.test = sum(p.test * q.test)
  
  # compute kappa for each data set
  kap.train = (acc.train - exp.acc.train) / (1 - exp.acc.train)
  kap.valid = (acc.valid - exp.acc.valid) / (1 - exp.acc.valid)
  kap.test = (acc.test - exp.acc.test) / (1 - exp.acc.test)
  
  # ---- compute one-vs-all metrics ----
  
  # compute a binary confusion matrix for each class, for each data set
  one.v.all.train = lapply(1:nrow(conf.train), function(i)
  {
    # get the four entries of a binary confusion matrix
    v = c(conf.train[i,i], 
          rsum.train[i] - conf.train[i,i], 
          csum.train[i] - conf.train[i,i], 
          n.train - rsum.train[i] - csum.train[i] + conf.train[i,i]);
    
    # build the confusion matrix
    return(matrix(v, nrow = 2, byrow = TRUE))
  })
  
  one.v.all.valid = lapply(1:nrow(conf.valid), function(i)
  {
    # get the four entries of a binary confusion matrix
    v = c(conf.valid[i,i], 
          rsum.valid[i] - conf.valid[i,i], 
          csum.valid[i] - conf.valid[i,i], 
          n.valid - rsum.valid[i] - csum.valid[i] + conf.valid[i,i]);
    
    # build the confusion matrix
    return(matrix(v, nrow = 2, byrow = TRUE))
  })
  
  one.v.all.test = lapply(1:nrow(conf.test), function(i)
  {
    # get the four entries of a binary confusion matrix
    v = c(conf.test[i,i], 
          rsum.test[i] - conf.test[i,i], 
          csum.test[i] - conf.test[i,i], 
          n.test - rsum.test[i] - csum.test[i] + conf.test[i,i]);
    
    # build the confusion matrix
    return(matrix(v, nrow = 2, byrow = TRUE))
  })
  
  # sum up all of the matrices for each data set
  one.v.all.train = Reduce('+', one.v.all.train)
  one.v.all.valid = Reduce('+', one.v.all.valid)
  one.v.all.test = Reduce('+', one.v.all.test)
  
  # compute the micro average accuracy for each data set
  micro.acc.train = sum(diag(one.v.all.train)) / sum(one.v.all.train)
  micro.acc.valid = sum(diag(one.v.all.valid)) / sum(one.v.all.valid)
  micro.acc.test = sum(diag(one.v.all.test)) / sum(one.v.all.test)
  
  # get the macro accuracy for each data set
  macro.acc.train = acc.train
  macro.acc.valid = acc.valid
  macro.acc.test = acc.test
  
  # ---- finalize output ----
  
  # build a final metrics table for nnet.mod
  nnet.mod.table = data.table(Model = "Neural Network",
                             Group = g,
                             Metric = c("Log_Loss", "Kappa", "Macro_Accuracy", "Micro_Accuracy"),
                             Train = c(mll.train, kap.train, macro.acc.train, micro.acc.train),
                             Valid = c(mll.valid, kap.valid, macro.acc.valid, micro.acc.valid),
                             Test = c(mll.test, kap.test, macro.acc.test, micro.acc.test))
  
  # build a final grid search table
  nnet.grid.search = data.table(cbind(Model = rep("Neural Network", nrow(DT.nnet.grid)), 
                                     Group = rep(g, nrow(DT.nnet.grid)), 
                                     DT.nnet.grid))
  
  # build a list of final tables
  nnet.list = list(nnet.mod.table, nnet.grid.search)
  names(nnet.list) = paste0(c("Neural_Network_Metrics_Group_", "Neural_Network_Search_Group_"), g)
  
  # remove some objects
  rm(dat.g, DT.nnet.grid, folds, n,
     nnet.grid, nnet.hyper.params, 
     dia.train, ynew.train.code, ytrue.train.code, ytrue.train.mat, 
     dia.valid, ynew.valid.code, ytrue.valid.code, ytrue.valid.mat, 
     dia.test, ynew.test.code, ytrue.test.code, ytrue.test.mat, 
     p.train, q.train, acc.train, exp.acc.train,
     p.valid, q.valid, acc.valid, exp.acc.valid,
     p.test, q.test, acc.test, exp.acc.test, max.class.weight,
     conf.train, rsum.train, csum.train, n.train, 
     conf.valid, rsum.valid, csum.valid, n.valid, 
     conf.test, rsum.test, csum.test, n.test, nnet.mod.table, nnet.grid.search,
     one.v.all.train, one.v.all.valid, one.v.all.test,
     mll.train, kap.train, macro.acc.train, micro.acc.train,
     mll.valid, kap.valid, macro.acc.valid, micro.acc.valid,
     mll.test, kap.test, macro.acc.test, micro.acc.test,
     nnet.mod, nnet.search.criteria, nnet.mod.list, nnet.random.grid,
     plot.nnet.grid, train, test, valid, small, medium, large,
     test.YX.h2o, train.YX.h2o, two.way, valid.YX.h2o, minutes, k,
     x, y, ynew.test, ynew.train, ynew.valid, ytrue.test, ytrue.train, ytrue.valid)
  
  # free up RAM
  gc()
  
  # return the final tables
  return(nnet.list)
}

# name the results
names(predict.total_return_class.nnet) = "h2o.deeplearning"

# combine predict.total_return into one table
predict.total_return_class = append(predict.total_return_class, predict.total_return_class.nnet)

# remove some objects
rm(nnet.list, g, predict.total_return_class.nnet)

# clean up the data in the h2o cluster
h2o.removeAll()

# shut down the h2o cluster to free up RAM
h2o.shutdown(prompt = FALSE)

# free up RAM
gc()

}

# -----------------------------------------------------------------------------------
# ---- Build Gradient Boosting Models -----------------------------------------------
# -----------------------------------------------------------------------------------

{

# initialize the h2o instance
h2o.init(nthreads = workers, max_mem_size = "8g")

# remove any objects in the h2o instance
h2o.removeAll()

# remove the progress bar when model building
h2o.no_progress()

# should we use extreme gradient boosting?
use.xgb = h2o.xgboost.available()

# ---- build prediction models for pmt_status ---------------------------------------

# determine how many data sets we have to build models for
if(cluster.data)
{
  # get the sub groups of dat
  groups = sort(unique(dat$Cluster))
  
  # also consider using all of dat
  groups = c("All", groups)
  
} else
{
  groups = "All"
}

# determine Gradient Boosting predictive performance for each group
predict.pmt_status.gb = foreach(g = groups) %do%
{
  # ---- split the data ----
  
  # determine which data set to use
  if(cluster.data)
  {
    if(g == "All")
    {
      # update dat to not contain Cluster
      dat.g = data.table(dat[, !"Cluster"])
      
    } else
    {
      # update dat to only contain data in Cluster g
      dat.g = data.table(dat[Cluster == g, !"Cluster"])
    }
  } else
  {
    # get dat
    dat.g = data.table(dat)
  }
  
  # identify predictors (x)
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # consider two-way interactions between the first n indicators in dat.g
  # recall that the indicators are ordered by importance from random forests
  n = floor(length(x) / 3)
  two.way = data.table(dat.g[, x[1:n], with = FALSE])
  
  # should we add two way interactions in the model?
  add.interactions = FALSE
  
  # build interactions if desired
  if(add.interactions)
  {
    # create the two way interactions and add them to dat.g
    two.way = data.table(model.matrix(~.^2, data = two.way)[,-(1:(n + 1))])
    dat.g = cbind(dat.g, two.way)
  }
  
  # identify predictors (x) and response (y)
  y = "pmt_status"
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # account for imbalance between classes
  if(use.xgb)
  {
    # the classes are imbalanced so lets define the class.weights parameter where a class with more observations is weighted less
    class.weights = table(unname(unlist(dat.g[, y, with = FALSE])))
    class.weights = max(class.weights) / class.weights
    
    # make class.weights into a table so we can join it onto dat
    class.weights = data.table(names(class.weights), 
                               as.numeric(class.weights))
    
    # update the column names of class.weights
    setnames(class.weights, c(y, "class.weights"))
    
    # join class.weights onto dat
    setkeyv(dat.g, y)
    setkeyv(class.weights, y)
    dat.g = data.table(class.weights[dat.g])
    
    # order dat.g by id
    dat.g = dat.g[order(id)]
    
  } else
  {
    # compute the max class weight
    max.class.weight = table(unname(unlist(dat.g[, y, with = FALSE])))
    max.class.weight = as.numeric(max(max(max.class.weight) / max.class.weight))
  }
  
  # build the fold assignment
  set.seed(42)
  k.folds = max(c(k.folds, 3))
  folds = createFolds(y = unname(unlist(dat.g[, y, with = FALSE])), k = k.folds)
  
  # split up dat.g into train, valid, and test
  train.rows = unname(unlist(lapply(1:(k.folds - 2), function(f) folds[[f]])))
  train = data.table(dat.g[train.rows])
  
  valid.rows = unname(unlist(folds[[k.folds - 1]]))
  valid = data.table(dat.g[valid.rows])
  
  test.rows = unname(unlist(folds[[k.folds]]))
  test = data.table(dat.g[test.rows])
  
  # split up YX.h2o into train, valid, and test
  train.YX.h2o = as.h2o(train[, c(y, x), with = FALSE])
  valid.YX.h2o = as.h2o(valid[, c(y, x), with = FALSE])
  test.YX.h2o = as.h2o(test[, c(y, x), with = FALSE])
  
  # ---- grid search for models ----
  
  # set up hyperparameters of interest
  gb.hyper.params = list(
    
    # -- initial tuning --
    ntrees = 250,
    learn_rate = c(0.025, 0.1, 0.3),
    max_depth = c(5, 15, 25), 
    min_rows = c(11, 25, 41),
    # sample_rate = c(0.7, 1),
    # col_sample_rate = c(0.7, 1),
    # min_split_improvement = c(0, 1e-5),
    
    # -- re-fine tuning --
    # ntrees = 100,
    # learn_rate = c(0.1, 0.3),
    # max_depth = 5, 
    # min_rows = 11,
    sample_rate = 0.7,
    col_sample_rate = 0.7,
    min_split_improvement = 1e-5,
    
    # -- scoring function -- 
    stopping_metric = "mean_per_class_error")
  
  # lets use a random grid search and specify a time limit and/or model limit
  minutes = 20
  gb.search.criteria = list(strategy = "RandomDiscrete", 
                            max_runtime_secs = minutes * 60, 
                            # max_models = 100, 
                            seed = 42)
  
  # build models
  if(use.xgb)
  {
    # lets run a grid search for a good model
    h2o.rm("gb.random.grid")
    gb.random.grid = h2o.grid(algorithm = "xgboost",
                              grid_id = "gb.random.grid",
                              distribution = "bernoulli",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              stopping_rounds = 5,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              seed = 21,
                              weights_column = "class.weights",
                              hyper_params = gb.hyper.params,
                              search_criteria = gb.search.criteria)
    
    # free up RAM
    gc()
    
  } else
  {
    # lets run a grid search for a good model, without drop out ratios
    h2o.rm("gb.random.grid")
    gb.random.grid = h2o.grid(algorithm = "gbm",
                              grid_id = "gb.random.grid",
                              distribution = "bernoulli",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              stopping_rounds = 5,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              seed = 21,
                              balance_classes = TRUE,
                              max_after_balance_size = max.class.weight,
                              hyper_params = gb.hyper.params,
                              search_criteria = gb.search.criteria)
    
    # free up RAM
    gc()
  }
  
  # rank each model in the random grids
  gb.grid = h2o.getGrid("gb.random.grid", sort_by = "auc", decreasing = TRUE)
  
  # get the summary table of the grid search
  DT.gb.grid = data.table(gb.grid@summary_table)
  DT.gb.grid
  
  # set up the data types for each column in DT.grid for plotting purposes
  DT.gb.grid = DT.gb.grid[, .(ntrees = as.numeric(ntrees),
                              learn_rate = as.numeric(learn_rate),
                              max_depth = as.factor(max_depth),
                              min_rows = as.factor(min_rows),
                              sample_rate = as.numeric(sample_rate),
                              col_sample_rate = as.numeric(col_sample_rate),
                              min_split_improvement = as.numeric(min_split_improvement),
                              model = removePunctuation(gsub("[A-z]+", "", model_ids)),
                              auc = as.numeric(auc))]
  
  # plot auc v. min_rows and activation to see which structure is most robust
  plot.gb.grid = ggplot(DT.gb.grid, aes(x = min_rows, y = auc, color = ntrees, fill = ntrees)) + 
    # geom_boxplot() + 
    geom_jitter(size = 3, alpha = 2/3) + 
    # scale_y_continuous(labels = dollar) + 
    ggtitle("Cross Validation Error") + 
    labs(x = "Minimum Child Weight", y = "AUC", color = "Number of Trees", fill = "Number of Trees") + 
    facet_wrap(~paste("Max Depth:", max_depth), nrow = 2) +
    theme_bw(base_size = 25) +
    theme(legend.position = "top", 
          legend.key.size = unit(.25, "in"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
  
  plot.gb.grid
  
  # pick the top model from all grid searches
  # gb.mod = h2o.getModel(gb.grid@model_ids[[1]])
  
  # find the top k models from all grid searches
  k = 10
  gb.mod.list = lapply(1:k, function(i) gb.grid@model_ids[[i]])
  
  # stack the top k models from all grid searches
  h2o.rm("gb.mod")
  gb.mod = h2o.stackedEnsemble(x = x,
                                 y = y,
                                 training_frame = train.YX.h2o,
                                 validation_frame = valid.YX.h2o,
                                 model_id = "gb.mod",
                                 base_models = gb.mod.list)
  
  # plot variable importance
  # h2o.varimp_plot(gb.mod)
  
  # ---- measure model quality ----
  
  # make predictions on each data set
  ynew.train = as.data.frame(predict(gb.mod, newdata = train.YX.h2o))$p1
  ynew.valid = as.data.frame(predict(gb.mod, newdata = valid.YX.h2o))$p1
  ynew.test = as.data.frame(predict(gb.mod, newdata = test.YX.h2o))$p1
  
  # get the true values from each data set
  ytrue.train = as.data.frame(train.YX.h2o)[,1]
  ytrue.valid = as.data.frame(valid.YX.h2o)[,1]
  ytrue.test = as.data.frame(test.YX.h2o)[,1]
  
  # plot the probability predictions on the testing data set
  # densityPlot(ynew.test)
  
  # compute prediction metrics on each data set
  gb.metrics.train = h2o.make_metrics(predicted = as.h2o(ynew.train), 
                                        actuals = as.h2o(factor(ytrue.train, levels = sort(unique(ytrue.train)))))
  
  gb.metrics.valid = h2o.make_metrics(predicted = as.h2o(ynew.valid), 
                                        actuals = as.h2o(factor(ytrue.valid, levels = sort(unique(ytrue.valid)))))
  
  gb.metrics.test = h2o.make_metrics(predicted = as.h2o(ynew.test), 
                                       actuals = as.h2o(factor(ytrue.test, levels = sort(unique(ytrue.test)))))
  
  # free up RAM
  gc()
  
  # check out the threshold metrics
  gb.metrics.train@metrics$max_criteria_and_metric_scores
  gb.metrics.valid@metrics$max_criteria_and_metric_scores
  gb.metrics.test@metrics$max_criteria_and_metric_scores
  
  # pick a threshold for converting probabilities into binary
  threshold = mean(c(tail(gb.metrics.train@metrics$max_criteria_and_metric_scores$threshold, 1),
                     tail(gb.metrics.valid@metrics$max_criteria_and_metric_scores$threshold, 1),
                     tail(gb.metrics.test@metrics$max_criteria_and_metric_scores$threshold, 1)))
  
  # compute a confusion matrix for each data set
  gb.confusion.train = confusion(ytrue = ytrue.train, ypred = as.numeric(ynew.train >= threshold))
  gb.confusion.valid = confusion(ytrue = ytrue.valid, ypred = as.numeric(ynew.valid >= threshold))
  gb.confusion.test = confusion(ytrue = ytrue.test, ypred = as.numeric(ynew.test >= threshold))
  
  # compute Sensitivity (accuracy in predicting 1's) for each data set
  gb.sensitivity.train = gb.confusion.train[2,2] / (gb.confusion.train[2,2] + gb.confusion.train[2,1])
  gb.sensitivity.valid = gb.confusion.valid[2,2] / (gb.confusion.valid[2,2] + gb.confusion.valid[2,1])
  gb.sensitivity.test = gb.confusion.test[2,2] / (gb.confusion.test[2,2] + gb.confusion.test[2,1])
  
  # compute Specificity (accuracy in predicting 0's) for each data set
  gb.specificity.train = gb.confusion.train[1,1] / (gb.confusion.train[1,1] + gb.confusion.train[1,2])
  gb.specificity.valid = gb.confusion.valid[1,1] / (gb.confusion.valid[1,1] + gb.confusion.valid[1,2])
  gb.specificity.test = gb.confusion.test[1,1] / (gb.confusion.test[1,1] + gb.confusion.test[1,2])
  
  # compute Odds Ratio (the odds of success over failure) for each data set
  gb.odds.ratio.train = (gb.confusion.train[1,1] * gb.confusion.train[2,2]) / (gb.confusion.train[2,1] * gb.confusion.train[1,2])
  gb.odds.ratio.valid = (gb.confusion.valid[1,1] * gb.confusion.valid[2,2]) / (gb.confusion.valid[2,1] * gb.confusion.valid[1,2])
  gb.odds.ratio.test = (gb.confusion.test[1,1] * gb.confusion.test[2,2]) / (gb.confusion.test[2,1] * gb.confusion.test[1,2])
  
  # compute Accuracy for each data set
  gb.accuracy.train = (gb.confusion.train[1,1] + gb.confusion.train[2,2]) / (gb.confusion.train[1,1] + gb.confusion.train[1,2] + gb.confusion.train[2,1] + gb.confusion.train[2,2])
  gb.accuracy.valid = (gb.confusion.valid[1,1] + gb.confusion.valid[2,2]) / (gb.confusion.valid[1,1] + gb.confusion.valid[1,2] + gb.confusion.valid[2,1] + gb.confusion.valid[2,2])
  gb.accuracy.test = (gb.confusion.test[1,1] + gb.confusion.test[2,2]) / (gb.confusion.test[1,1] + gb.confusion.test[1,2] + gb.confusion.test[2,1] + gb.confusion.test[2,2])
  
  # compute AUC for each data set
  gb.auc.train = h2o.auc(gb.metrics.train)
  gb.auc.valid = h2o.auc(gb.metrics.valid)
  gb.auc.test = h2o.auc(gb.metrics.test)
  
  # compute Log Loss for each data set
  gb.logloss.train = h2o.logloss(gb.metrics.train)
  gb.logloss.valid = h2o.logloss(gb.metrics.valid)
  gb.logloss.test = h2o.logloss(gb.metrics.test)
  
  # ---- compute kappa ----
  
  # get the total number of observations for each data set
  n.train = sum(gb.confusion.train)
  n.valid = sum(gb.confusion.valid)
  n.test = sum(gb.confusion.test)
  
  # get the vector of correct predictions for each data set 
  dia.train = diag(gb.confusion.train)
  dia.valid = diag(gb.confusion.valid)
  dia.test = diag(gb.confusion.test)
  
  # get the vector of the number of observations per class for each data set
  rsum.train = rowSums(gb.confusion.train)
  rsum.valid = rowSums(gb.confusion.valid)
  rsum.test = rowSums(gb.confusion.test)
  
  # get the vector of the number of predictions per class for each data set
  csum.train = colSums(gb.confusion.train)
  csum.valid = colSums(gb.confusion.valid)
  csum.test = colSums(gb.confusion.test)
  
  # get the proportion of observations per class for each data set
  p.train = rsum.train / n.train
  p.valid = rsum.valid / n.valid
  p.test = rsum.test / n.test
  
  # get the proportion of predcitions per class for each data set
  q.train = csum.train / n.train
  q.valid = csum.valid / n.valid
  q.test = csum.test / n.test
  
  # compute accuracy for each data set
  acc.train = sum(dia.train) / n.train
  acc.valid = sum(dia.valid) / n.valid
  acc.test = sum(dia.test) / n.test
  
  # compute expected accuracy for each data set
  exp.acc.train = sum(p.train * q.train)
  exp.acc.valid = sum(p.valid * q.valid)
  exp.acc.test = sum(p.test * q.test)
  
  # compute kappa for each data set
  kap.train = (acc.train - exp.acc.train) / (1 - exp.acc.train)
  kap.valid = (acc.valid - exp.acc.valid) / (1 - exp.acc.valid)
  kap.test = (acc.test - exp.acc.test) / (1 - exp.acc.test)
  
  # ---- finalize output ----
  
  # build a final metrics table for gb.mod
  gb.mod.table = data.table(Model = "Gradient Boosting",
                              Group = g,
                              Threshold = threshold,
                              Metric = c("AUC", "Accuracy", "Log_Loss", "Kappa", "Odds_Ratio", "Specificity_0", "Sensitivity_1"),
                              Train = c(gb.auc.train, gb.accuracy.train, gb.logloss.train, kap.train, gb.odds.ratio.train, gb.specificity.train, gb.sensitivity.train),
                              Valid = c(gb.auc.valid, gb.accuracy.valid, gb.logloss.valid, kap.valid, gb.odds.ratio.valid, gb.specificity.valid, gb.sensitivity.valid),
                              Test = c(gb.auc.test, gb.accuracy.test, gb.logloss.test, kap.test, gb.odds.ratio.test, gb.specificity.test, gb.sensitivity.test))
  
  # build a final grid search table
  gb.grid.search = data.table(cbind(Model = rep("Gradient Boosting", nrow(DT.gb.grid)), 
                                      Group = rep(g, nrow(DT.gb.grid)), 
                                      DT.gb.grid))
  
  # build a list of confusion matrices
  gb.confusion = list("Train" = gb.confusion.train, 
                        "Valid" = gb.confusion.valid, 
                        "Test" = gb.confusion.test)
  
  # build a list of final tables
  gb.list = list(gb.mod.table, gb.grid.search, gb.confusion)
  names(gb.list) = paste0(c("Gradient_Boosting_Metrics_Group_", "Gradient_Boosting_Search_Group_", "Gradient_Boosting_Error_Group_"), g)
  
  # remove some objects
  rm(dat.g, DT.gb.grid, folds, gb.accuracy.test, gb.accuracy.train, gb.accuracy.valid,
     gb.auc.test, gb.auc.train, gb.auc.valid, gb.confusion.test, gb.confusion.train, gb.confusion.valid,
     gb.grid, gb.hyper.params, valid.rows, gb.confusion, gb.mod.table, gb.grid.search,
     gb.logloss.test, gb.logloss.train, gb.logloss.valid, gb.metrics.test, gb.metrics.train, gb.metrics.valid,
     gb.mod, gb.odds.ratio.test, gb.odds.ratio.train, gb.odds.ratio.valid,
     gb.search.criteria, gb.sensitivity.test, gb.mod.list,
     gb.sensitivity.train, gb.sensitivity.valid, gb.specificity.test, gb.specificity.train,
     gb.specificity.valid, max.class.weight, plot.gb.grid, train, test, valid, 
     small, test.rows, train.rows, k,
     test.YX.h2o, threshold, train.YX.h2o, valid.YX.h2o, gb.random.grid,
     x, y, ynew.test, ynew.train, ynew.valid, ytrue.test, ytrue.train, ytrue.valid,
     n.train, dia.train, rsum.train, csum.train, p.train, q.train, acc.train, exp.acc.train, kap.train,
     n.valid, dia.valid, rsum.valid, csum.valid, p.valid, q.valid, acc.valid, exp.acc.valid, kap.valid,
     n.test, dia.test, rsum.test, csum.test, p.test, q.test, acc.test, exp.acc.test, kap.test)
  
  # free up RAM
  gc()
  
  # return the final metrics table
  return(gb.list)
}

# name the results
names(predict.pmt_status.gb) = "h2o.gb"

# combine predict.pmt_status into one table
predict.pmt_status = append(predict.pmt_status, predict.pmt_status.gb)

# remove some objects
rm(gb.list, g, predict.pmt_status.gb)

# free up RAM
gc()

# ---- build prediction models for total_return ---------------------------------------

# determine how many data sets we have to build models for
if(cluster.data)
{
  # get the sub groups of dat
  groups = sort(unique(dat$Cluster))
  
  # also consider using all of dat
  groups = c("All", groups)
  
} else
{
  groups = "All"
}

# determine Gradient Boosting predictive performance for each group
predict.total_return.gb = foreach(g = groups) %do%
{
  # ---- split the data ----
  
  # determine which data set to use
  if(cluster.data)
  {
    if(g == "All")
    {
      # update dat to not contain Cluster
      dat.g = data.table(dat[, !"Cluster"])
      
    } else
    {
      # update dat to only contain data in Cluster g
      dat.g = data.table(dat[Cluster == g, !"Cluster"])
    }
  } else
  {
    # get dat
    dat.g = data.table(dat)
  }
  
  # identify predictors (x)
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # consider two-way interactions between the first n indicators in dat.g
  # recall that the indicators are ordered by importance from random forests
  n = floor(length(x) / 3)
  two.way = data.table(dat.g[, x[1:n], with = FALSE])
  
  # should we add two way interactions in the model?
  add.interactions = FALSE
  
  # build interactions if desired
  if(add.interactions)
  {
    # create the two way interactions and add them to dat.g
    two.way = data.table(model.matrix(~.^2, data = two.way)[,-(1:(n + 1))])
    dat.g = cbind(dat.g, two.way)
  }
  
  # identify predictors (x) and response (y)
  y = "total_return"
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # account for imbalance between classes
  if(use.xgb)
  {
    # the classes are imbalanced so lets define the class.weights parameter where a class with more observations is weighted less
    class.weights = table(unname(unlist(dat.g[, y, with = FALSE])))
    class.weights = max(class.weights) / class.weights
    
    # make class.weights into a table so we can join it onto dat
    class.weights = data.table(names(class.weights), 
                               as.numeric(class.weights))
    
    # update the column names of class.weights
    setnames(class.weights, c(y, "class.weights"))
    
    # join class.weights onto dat
    setkeyv(dat.g, y)
    setkeyv(class.weights, y)
    dat.g = data.table(class.weights[dat.g])
    
    # order dat.g by id
    dat.g = dat.g[order(id)]
    
  } else
  {
    # compute the max class weight
    max.class.weight = table(unname(unlist(dat.g[, y, with = FALSE])))
    max.class.weight = as.numeric(max(max(max.class.weight) / max.class.weight))
  }
  
  # build the fold assignment
  set.seed(42)
  k.folds = max(c(k.folds, 3))
  folds = createFolds(y = unname(unlist(dat.g[, y, with = FALSE])), k = k.folds)
  
  # split up dat.g into train, valid, and test
  train.rows = unname(unlist(lapply(1:(k.folds - 2), function(f) folds[[f]])))
  train = data.table(dat.g[train.rows])
  
  valid.rows = unname(unlist(folds[[k.folds - 1]]))
  valid = data.table(dat.g[valid.rows])
  
  test.rows = unname(unlist(folds[[k.folds]]))
  test = data.table(dat.g[test.rows])
  
  # split up YX.h2o into train, valid, and test
  train.YX.h2o = as.h2o(train[, c(y, x), with = FALSE])
  valid.YX.h2o = as.h2o(valid[, c(y, x), with = FALSE])
  test.YX.h2o = as.h2o(test[, c(y, x), with = FALSE])
  
  # ---- grid search for models ----
  
  # set up hyperparameters of interest
  gb.hyper.params = list(
    
    # -- initial tuning --
    ntrees = 250,
    learn_rate = c(0.025, 0.1, 0.3),
    max_depth = c(5, 15, 25), 
    min_rows = c(11, 25, 41),
    # sample_rate = c(0.7, 1),
    # col_sample_rate = c(0.7, 1),
    # min_split_improvement = c(0, 1e-5),
    
    # -- re-fine tuning --
    # ntrees = 100,
    # learn_rate = c(0.1, 0.3),
    # max_depth = 5, 
    # min_rows = 11,
    sample_rate = 0.7,
    col_sample_rate = 0.7,
    min_split_improvement = 1e-5,
    
    # -- scoring function -- 
    stopping_metric = "RMSE")
  
  # lets use a random grid search and specify a time limit and/or model limit
  minutes = 20
  gb.search.criteria = list(strategy = "RandomDiscrete", 
                            max_runtime_secs = minutes * 60, 
                            # max_models = 100, 
                            seed = 42)
  
  # build models
  if(use.xgb)
  {
    # lets run a grid search for a good model, without drop out ratios
    h2o.rm("gb.random.grid")
    gb.random.grid = h2o.grid(algorithm = "xgboost",
                              grid_id = "gb.random.grid",
                              distribution = "gaussian",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              stopping_rounds = 5,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              seed = 21,
                              hyper_params = gb.hyper.params,
                              search_criteria = gb.search.criteria)
    
    # free up RAM
    gc()
    
  } else
  {
    # lets run a grid search for a good model, without drop out ratios
    h2o.rm("gb.random.grid")
    gb.random.grid = h2o.grid(algorithm = "gbm",
                              grid_id = "gb.random.grid",
                              distribution = "gaussian",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              stopping_rounds = 5,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              seed = 21,
                              hyper_params = gb.hyper.params,
                              search_criteria = gb.search.criteria)
    
    # free up RAM
    gc()
  }
  
  # rank each model in the random grids
  gb.grid = h2o.getGrid("gb.random.grid", sort_by = "rmse", decreasing = FALSE)
  
  # get the summary table of the grid search
  DT.gb.grid = data.table(gb.grid@summary_table)
  DT.gb.grid
  
  # set up the data types for each column in DT.grid for plotting purposes
  DT.gb.grid = DT.gb.grid[, .(ntrees = as.numeric(ntrees),
                              learn_rate = as.numeric(learn_rate),
                              max_depth = as.factor(max_depth),
                              min_rows = as.factor(min_rows),
                              sample_rate = as.numeric(sample_rate),
                              col_sample_rate = as.numeric(col_sample_rate),
                              min_split_improvement = as.numeric(min_split_improvement),
                              model = removePunctuation(gsub("[A-z]+", "", model_ids)),
                              rmse = as.numeric(rmse))]
  
  # plot rmse v. min_rows and activation to see which structure is most robust
  plot.gb.grid = ggplot(DT.gb.grid, aes(x = min_rows, y = rmse, color = ntrees, fill = ntrees)) + 
    # geom_boxplot() + 
    geom_jitter(size = 3, alpha = 2/3) + 
    # scale_y_continuous(labels = dollar) + 
    ggtitle("Cross Validation Error") + 
    labs(x = "Minimum Child Weight", y = "RMSE", color = "Number of Trees", fill = "Number of Trees") + 
    facet_wrap(~paste("Max Depth:", max_depth), nrow = 2) +
    theme_bw(base_size = 25) +
    theme(legend.position = "top", 
          legend.key.size = unit(.25, "in"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
  
  plot.gb.grid
  
  # pick the top model from all grid searches
  # gb.mod = h2o.getModel(gb.grid@model_ids[[1]])
  
  # find the top k models from all grid searches
  k = 10
  gb.mod.list = lapply(1:k, function(i) gb.grid@model_ids[[i]])
  
  # stack the top k models from all grid searches
  h2o.rm("gb.mod")
  gb.mod = h2o.stackedEnsemble(x = x,
                                 y = y,
                                 training_frame = train.YX.h2o,
                                 validation_frame = valid.YX.h2o,
                                 model_id = "gb.mod",
                                 base_models = gb.mod.list)
  
  # plot variable importance
  # h2o.varimp_plot(gb.mod)
  
  # ---- measure model quality ----
  
  # make predictions on each data set
  ynew.train = as.data.frame(predict(gb.mod, newdata = train.YX.h2o))$predict
  ynew.valid = as.data.frame(predict(gb.mod, newdata = valid.YX.h2o))$predict
  ynew.test = as.data.frame(predict(gb.mod, newdata = test.YX.h2o))$predict
  
  # get the true values from each data set
  ytrue.train = as.data.frame(train.YX.h2o)[,1]
  ytrue.valid = as.data.frame(valid.YX.h2o)[,1]
  ytrue.test = as.data.frame(test.YX.h2o)[,1]
  
  # transform the predicted and actual values back to original units if needed
  if(use.norm)
  {
    # map predicted values to the original units of the response variable
    ynew.train = predict(best.norm, newdata = ynew.train, inverse = TRUE)
    ynew.valid = predict(best.norm, newdata = ynew.valid, inverse = TRUE)
    ynew.test = predict(best.norm, newdata = ynew.test, inverse = TRUE)
    
    # map actual values to the original units of the response variable
    ytrue.train = predict(best.norm, newdata = ytrue.train, inverse = TRUE)
    ytrue.valid = predict(best.norm, newdata = ytrue.valid, inverse = TRUE)
    ytrue.test = predict(best.norm, newdata = ytrue.test, inverse = TRUE)
  }
  
  # compute prediction metrics on each data set
  gb.metrics.train = h2o.make_metrics(predicted = as.h2o(ynew.train), 
                                        actuals = as.h2o(ytrue.train))
  
  gb.metrics.valid = h2o.make_metrics(predicted = as.h2o(ynew.valid), 
                                        actuals = as.h2o(ytrue.valid))
  
  gb.metrics.test = h2o.make_metrics(predicted = as.h2o(ynew.test), 
                                       actuals = as.h2o(ytrue.test))
  
  # free up RAM
  gc()
  
  # check out the prediction metrics
  gb.metrics.train
  gb.metrics.valid
  gb.metrics.test
  
  # fit the normal distribution to the residuals
  gb.resid.train = fitdist(data = ytrue.train - ynew.train, distr = "norm")
  gb.resid.valid = fitdist(data = ytrue.valid - ynew.valid, distr = "norm") 
  gb.resid.test = fitdist(data = ytrue.test - ynew.test, distr = "norm") 
  
  # ---- finalize output ----
  
  # build a final metrics table for gb.mod
  gb.mod.table = data.table(Model = "Gradient Boosting",
                              Group = g,
                              Metric = c("R2", "RMSE"),
                              Train = c(gb.metrics.train@metrics$r2, gb.metrics.train@metrics$RMSE),
                              Valid = c(gb.metrics.valid@metrics$r2, gb.metrics.valid@metrics$RMSE),
                              Test = c(gb.metrics.test@metrics$r2, gb.metrics.test@metrics$RMSE))
  
  # build a final grid search table
  gb.grid.search = data.table(cbind(Model = rep("Gradient Boosting", nrow(DT.gb.grid)), 
                                      Group = rep(g, nrow(DT.gb.grid)), 
                                      DT.gb.grid))
  
  # build a list of residual plots
  gb.resid = list("Train" = gb.resid.train, 
                    "Valid" = gb.resid.valid, 
                    "Test" = gb.resid.test)
  
  # build a list of final tables
  gb.list = list(gb.mod.table, gb.grid.search, gb.resid)
  names(gb.list) = paste0(c("Gradient_Boosting_Metrics_Group_", "Gradient_Boosting_Search_Group_", "Gradient_Boosting_Error_Group_"), g)
  
  # remove some objects
  rm(dat.g, DT.gb.grid, folds, n,
     gb.grid, gb.hyper.params, 
     gb.metrics.test, gb.metrics.train, gb.metrics.valid,
     gb.mod, gb.resid.test, gb.resid.train, gb.resid.valid,
     gb.search.criteria, small, medium, large, k,
     plot.gb.grid, train, test, valid, 
     gb.random.grid, gb.mod.table, gb.grid.search, gb.resid,
     test.YX.h2o, train.YX.h2o, two.way, valid.YX.h2o,
     x, y, ynew.test, ynew.train, ynew.valid, ytrue.test, ytrue.train, ytrue.valid)
  
  # free up RAM
  gc()
  
  # return the final metrics table
  return(gb.list)
}

# name the results
names(predict.total_return.gb) = "h2o.gb"

# combine predict.total_return into one table
predict.total_return = append(predict.total_return, predict.total_return.gb)

# remove some objects
rm(gb.list, g, predict.total_return.gb)

# free up RAM
gc()

# ---- build prediction models for total_return_class ---------------------------------------

# determine how many data sets we have to build models for
if(cluster.data)
{
  # get the sub groups of dat
  groups = sort(unique(dat$Cluster))
  
  # also consider using all of dat
  groups = c("All", groups)
  
} else
{
  groups = "All"
}

# determine Gradient Boosting predictive performance for each group
predict.total_return_class.gb = foreach(g = groups) %do%
{
  # ---- split the data ----
  
  # determine which data set to use
  if(cluster.data)
  {
    if(g == "All")
    {
      # update dat to not contain Cluster
      dat.g = data.table(dat[, !"Cluster"])
      
    } else
    {
      # update dat to only contain data in Cluster g
      dat.g = data.table(dat[Cluster == g, !"Cluster"])
    }
  } else
  {
    # get dat
    dat.g = data.table(dat)
  }
  
  # identify predictors (x)
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # consider two-way interactions between the first n indicators in dat.g
  # recall that the indicators are ordered by importance from random forests
  n = floor(length(x) / 3)
  two.way = data.table(dat.g[, x[1:n], with = FALSE])
  
  # should we add two way interactions in the model?
  add.interactions = FALSE
  
  # build interactions if desired
  if(add.interactions)
  {
    # create the two way interactions and add them to dat.g
    two.way = data.table(model.matrix(~.^2, data = two.way)[,-(1:(n + 1))])
    dat.g = cbind(dat.g, two.way)
  }
  
  # identify predictors (x) and response (y)
  y = "total_return_class"
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # account for imbalance between classes
  if(use.xgb)
  {
    # the classes are imbalanced so lets define the class.weights parameter where a class with more observations is weighted less
    class.weights = table(unname(unlist(dat.g[, y, with = FALSE])))
    class.weights = max(class.weights) / class.weights
    
    # make class.weights into a table so we can join it onto dat
    class.weights = data.table(names(class.weights), 
                               as.numeric(class.weights))
    
    # update the column names of class.weights
    setnames(class.weights, c(y, "class.weights"))
    
    # join class.weights onto dat
    setkeyv(dat.g, y)
    setkeyv(class.weights, y)
    dat.g = data.table(class.weights[dat.g])
    
    # order dat.g by id
    dat.g = dat.g[order(id)]
    
  } else
  {
    # compute the max class weight
    max.class.weight = table(unname(unlist(dat.g[, y, with = FALSE])))
    max.class.weight = as.numeric(max(max(max.class.weight) / max.class.weight))
  }
  
  # build the fold assignment
  set.seed(42)
  k.folds = max(c(k.folds, 3))
  folds = createFolds(y = unname(unlist(dat.g[, y, with = FALSE])), k = k.folds)
  
  # split up dat.g into train, valid, and test
  train.rows = unname(unlist(lapply(1:(k.folds - 2), function(f) folds[[f]])))
  train = data.table(dat.g[train.rows])
  
  valid.rows = unname(unlist(folds[[k.folds - 1]]))
  valid = data.table(dat.g[valid.rows])
  
  test.rows = unname(unlist(folds[[k.folds]]))
  test = data.table(dat.g[test.rows])
  
  # split up YX.h2o into train, valid, and test
  train.YX.h2o = as.h2o(train[, c(y, x), with = FALSE])
  valid.YX.h2o = as.h2o(valid[, c(y, x), with = FALSE])
  test.YX.h2o = as.h2o(test[, c(y, x), with = FALSE])
  
  # ---- grid search for models ----
  
  # set up hyperparameters of interest
  gb.hyper.params = list(
    
    # -- initial tuning --
    ntrees = 250,
    learn_rate = c(0.025, 0.1, 0.3),
    max_depth = c(5, 15, 25), 
    min_rows = c(11, 25, 41),
    # sample_rate = c(0.7, 1),
    # col_sample_rate = c(0.7, 1),
    # min_split_improvement = c(0, 1e-5),
    
    # -- re-fine tuning --
    # ntrees = 100,
    # learn_rate = c(0.1, 0.3),
    # max_depth = 5, 
    # min_rows = 11,
    sample_rate = 0.7,
    col_sample_rate = 0.7,
    min_split_improvement = 1e-5,
    
    # -- scoring function -- 
    stopping_metric = "mean_per_class_error")
  
  # lets use a random grid search and specify a time limit and/or model limit
  minutes = 15
  gb.search.criteria = list(strategy = "RandomDiscrete", 
                            max_runtime_secs = minutes * 60, 
                            # max_models = 100, 
                            seed = 42)
  
  # build models
  if(use.xgb)
  {
    # lets run a grid search for a good model, without drop out ratios
    h2o.rm("gb.random.grid")
    gb.random.grid = h2o.grid(algorithm = "xgboost",
                              grid_id = "gb.random.grid",
                              distribution = "multinomial",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              stopping_rounds = 5,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              seed = 21,
                              weights_column = "class.weights",
                              hyper_params = gb.hyper.params,
                              search_criteria = gb.search.criteria)
    
    # free up RAM
    gc()
    
  } else
  {
    # lets run a grid search for a good model, without drop out ratios
    h2o.rm("gb.random.grid")
    gb.random.grid = h2o.grid(algorithm = "gbm",
                              grid_id = "gb.random.grid",
                              distribution = "multinomial",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                              stopping_rounds = 5,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              seed = 21,
                              balance_classes = TRUE,
                              max_after_balance_size = max.class.weight,
                              hyper_params = gb.hyper.params,
                              search_criteria = gb.search.criteria)
    
    # free up RAM
    gc()
  }
  
  # rank each model in the random grids
  gb.grid = h2o.getGrid("gb.random.grid", sort_by = "mean_per_class_error", decreasing = FALSE)
  
  # get the summary table of the grid search
  DT.gb.grid = data.table(gb.grid@summary_table)
  DT.gb.grid
  
  # set up the data types for each column in DT.grid for plotting purposes
  DT.gb.grid = DT.gb.grid[, .(ntrees = as.numeric(ntrees),
                              learn_rate = as.numeric(learn_rate),
                              max_depth = as.factor(max_depth),
                              min_rows = as.factor(min_rows),
                              sample_rate = as.numeric(sample_rate),
                              col_sample_rate = as.numeric(col_sample_rate),
                              min_split_improvement = as.numeric(min_split_improvement),
                              model = removePunctuation(gsub("[A-z]+", "", model_ids)),
                              mean_per_class_error = as.numeric(mean_per_class_error))]
  
  # plot mean_per_class_error v. min_rows and activation to see which structure is most robust
  plot.gb.grid = ggplot(DT.gb.grid, aes(x = min_rows, y = mean_per_class_error, color = ntrees, fill = ntrees)) + 
    # geom_boxplot() + 
    geom_jitter(size = 3, alpha = 2/3) + 
    # scale_y_continuous(labels = dollar) + 
    ggtitle("Cross Validation Error") + 
    labs(x = "Minimum Child Weight", y = "Mean Per Class Error", color = "Number of Trees", fill = "Number of Trees") + 
    facet_wrap(~paste("Max Depth:", max_depth), nrow = 2) +
    theme_bw(base_size = 25) +
    theme(legend.position = "top", 
          legend.key.size = unit(.25, "in"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
  
  plot.gb.grid
  
  # pick the top model from all grid searches
  # gb.mod = h2o.getModel(gb.grid@model_ids[[1]])
  
  # find the top k models from all grid searches
  k = 10
  gb.mod.list = lapply(1:k, function(i) gb.grid@model_ids[[i]])
  
  # stack the top k models from all grid searches
  h2o.rm("gb.mod")
  gb.mod = h2o.stackedEnsemble(x = x,
                                 y = y,
                                 training_frame = train.YX.h2o,
                                 validation_frame = valid.YX.h2o,
                                 model_id = "gb.mod",
                                 base_models = gb.mod.list)
  
  # plot variable importance
  # h2o.varimp_plot(gb.mod)
  
  # ---- measure model quality ----
  
  # make predictions on each data set
  ynew.train = as.matrix(predict(gb.mod, newdata = train.YX.h2o)[,-1])
  ynew.valid = as.matrix(predict(gb.mod, newdata = valid.YX.h2o)[,-1])
  ynew.test = as.matrix(predict(gb.mod, newdata = test.YX.h2o)[,-1])
  
  # get the true values from each data set
  ytrue.train = as.data.frame(train.YX.h2o[,1])
  ytrue.valid = as.data.frame(valid.YX.h2o[,1])
  ytrue.test = as.data.frame(test.YX.h2o[,1])
  
  # ---- compute multi log loss ----
  
  # build a matrix indicating the true class values for each data set
  ytrue.train.mat = model.matrix(~., data = ytrue.train)[,-1]
  ytrue.train.mat = cbind(1 - rowSums(ytrue.train.mat), ytrue.train.mat)
  
  ytrue.valid.mat = model.matrix(~., data = ytrue.valid)[,-1]
  ytrue.valid.mat = cbind(1 - rowSums(ytrue.valid.mat), ytrue.valid.mat)
  
  ytrue.test.mat = model.matrix(~., data = ytrue.test)[,-1]
  ytrue.test.mat = cbind(1 - rowSums(ytrue.test.mat), ytrue.test.mat)
  
  # compute the multi-class logarithmic loss for each data set
  mll.train = MultiLogLoss(y_pred = ynew.train, y_true = ytrue.train.mat)
  mll.valid = MultiLogLoss(y_pred = ynew.valid, y_true = ytrue.valid.mat)
  mll.test = MultiLogLoss(y_pred = ynew.test, y_true = ytrue.test.mat)
  
  # free up RAM
  gc()
  
  # ---- compute kappa ----
  
  # get the predicted classes and actual classes for each data set
  ynew.train.code = apply(ynew.train, 1, which.max) - 1
  ytrue.train.code = factor(as.numeric(ytrue.train[,1]) - 1, levels = 0:(ncol(ytrue.train.mat) - 1))
  
  ynew.valid.code = apply(ynew.valid, 1, which.max) - 1
  ytrue.valid.code = factor(as.numeric(ytrue.valid[,1]) - 1, levels = 0:(ncol(ytrue.valid.mat) - 1))
  
  ynew.test.code = apply(ynew.test, 1, which.max) - 1
  ytrue.test.code = factor(as.numeric(ytrue.test[,1]) - 1, levels = 0:(ncol(ytrue.test.mat) - 1))
  
  # build a square confusion matrix for each data set
  conf.train = confusion(ytrue = ytrue.train.code, ypred = ynew.train.code)
  conf.valid = confusion(ytrue = ytrue.valid.code, ypred = ynew.valid.code)
  conf.test = confusion(ytrue = ytrue.test.code, ypred = ynew.test.code)
  
  # get the total number of observations for each data set
  n.train = sum(conf.train)
  n.valid = sum(conf.valid)
  n.test = sum(conf.test)
  
  # get the vector of correct predictions for each data set 
  dia.train = diag(conf.train)
  dia.valid = diag(conf.valid)
  dia.test = diag(conf.test)
  
  # get the vector of the number of observations per class for each data set
  rsum.train = rowSums(conf.train)
  rsum.valid = rowSums(conf.valid)
  rsum.test = rowSums(conf.test)
  
  # get the vector of the number of predictions per class for each data set
  csum.train = colSums(conf.train)
  csum.valid = colSums(conf.valid)
  csum.test = colSums(conf.test)
  
  # get the proportion of observations per class for each data set
  p.train = rsum.train / n.train
  p.valid = rsum.valid / n.valid
  p.test = rsum.test / n.test
  
  # get the proportion of predcitions per class for each data set
  q.train = csum.train / n.train
  q.valid = csum.valid / n.valid
  q.test = csum.test / n.test
  
  # compute accuracy for each data set
  acc.train = sum(dia.train) / n.train
  acc.valid = sum(dia.valid) / n.valid
  acc.test = sum(dia.test) / n.test
  
  # compute expected accuracy for each data set
  exp.acc.train = sum(p.train * q.train)
  exp.acc.valid = sum(p.valid * q.valid)
  exp.acc.test = sum(p.test * q.test)
  
  # compute kappa for each data set
  kap.train = (acc.train - exp.acc.train) / (1 - exp.acc.train)
  kap.valid = (acc.valid - exp.acc.valid) / (1 - exp.acc.valid)
  kap.test = (acc.test - exp.acc.test) / (1 - exp.acc.test)
  
  # ---- compute one-vs-all metrics ----
  
  # compute a binary confusion matrix for each class, for each data set
  one.v.all.train = lapply(1:nrow(conf.train), function(i)
  {
    # get the four entries of a binary confusion matrix
    v = c(conf.train[i,i], 
          rsum.train[i] - conf.train[i,i], 
          csum.train[i] - conf.train[i,i], 
          n.train - rsum.train[i] - csum.train[i] + conf.train[i,i]);
    
    # build the confusion matrix
    return(matrix(v, nrow = 2, byrow = TRUE))
  })
  
  one.v.all.valid = lapply(1:nrow(conf.valid), function(i)
  {
    # get the four entries of a binary confusion matrix
    v = c(conf.valid[i,i], 
          rsum.valid[i] - conf.valid[i,i], 
          csum.valid[i] - conf.valid[i,i], 
          n.valid - rsum.valid[i] - csum.valid[i] + conf.valid[i,i]);
    
    # build the confusion matrix
    return(matrix(v, nrow = 2, byrow = TRUE))
  })
  
  one.v.all.test = lapply(1:nrow(conf.test), function(i)
  {
    # get the four entries of a binary confusion matrix
    v = c(conf.test[i,i], 
          rsum.test[i] - conf.test[i,i], 
          csum.test[i] - conf.test[i,i], 
          n.test - rsum.test[i] - csum.test[i] + conf.test[i,i]);
    
    # build the confusion matrix
    return(matrix(v, nrow = 2, byrow = TRUE))
  })
  
  # sum up all of the matrices for each data set
  one.v.all.train = Reduce('+', one.v.all.train)
  one.v.all.valid = Reduce('+', one.v.all.valid)
  one.v.all.test = Reduce('+', one.v.all.test)
  
  # compute the micro average accuracy for each data set
  micro.acc.train = sum(diag(one.v.all.train)) / sum(one.v.all.train)
  micro.acc.valid = sum(diag(one.v.all.valid)) / sum(one.v.all.valid)
  micro.acc.test = sum(diag(one.v.all.test)) / sum(one.v.all.test)
  
  # get the macro accuracy for each data set
  macro.acc.train = acc.train
  macro.acc.valid = acc.valid
  macro.acc.test = acc.test
  
  # ---- finalize output ----
  
  # build a final metrics table for gb.mod
  gb.mod.table = data.table(Model = "Gradient Boosting",
                              Group = g,
                              Metric = c("Log_Loss", "Kappa", "Macro_Accuracy", "Micro_Accuracy"),
                              Train = c(mll.train, kap.train, macro.acc.train, micro.acc.train),
                              Valid = c(mll.valid, kap.valid, macro.acc.valid, micro.acc.valid),
                              Test = c(mll.test, kap.test, macro.acc.test, micro.acc.test))
  
  # build a final grid search table
  gb.grid.search = data.table(cbind(Model = rep("Gradient Boosting", nrow(DT.gb.grid)), 
                                      Group = rep(g, nrow(DT.gb.grid)), 
                                      DT.gb.grid))
  
  # build a list of final tables
  gb.list = list(gb.mod.table, gb.grid.search)
  names(gb.list) = paste0(c("Gradient_Boosting_Metrics_Group_", "Gradient_Boosting_Search_Group_"), g)
  
  # remove some objects
  rm(dat.g, DT.gb.grid, folds, n,
     gb.grid, gb.hyper.params, 
     dia.train, ynew.train.code, ytrue.train.code, ytrue.train.mat, 
     dia.valid, ynew.valid.code, ytrue.valid.code, ytrue.valid.mat, 
     dia.test, ynew.test.code, ytrue.test.code, ytrue.test.mat, 
     p.train, q.train, acc.train, exp.acc.train,
     p.valid, q.valid, acc.valid, exp.acc.valid,
     p.test, q.test, acc.test, exp.acc.test, max.class.weight,
     conf.train, rsum.train, csum.train, n.train, 
     conf.valid, rsum.valid, csum.valid, n.valid, 
     conf.test, rsum.test, csum.test, n.test, gb.mod.table, gb.grid.search,
     one.v.all.train, one.v.all.valid, one.v.all.test,
     mll.train, kap.train, macro.acc.train, micro.acc.train,
     mll.valid, kap.valid, macro.acc.valid, micro.acc.valid,
     mll.test, kap.test, macro.acc.test, micro.acc.test,
     gb.mod, gb.search.criteria, gb.mod.list, gb.random.grid,
     plot.gb.grid, train, test, valid, small, medium, large,
     test.YX.h2o, train.YX.h2o, two.way, valid.YX.h2o, minutes, k,
     x, y, ynew.test, ynew.train, ynew.valid, ytrue.test, ytrue.train, ytrue.valid)
  
  # free up RAM
  gc()
  
  # return the final tables
  return(gb.list)
}

# name the results
names(predict.total_return_class.gb) = "h2o.gb"

# combine predict.total_return into one table
predict.total_return_class = append(predict.total_return_class, predict.total_return_class.gb)

# remove some objects
rm(gb.list, g, predict.total_return_class.gb)

# clean up the data in the h2o cluster
h2o.removeAll()

# shut down the h2o cluster to free up RAM
h2o.shutdown(prompt = FALSE)

# free up RAM
gc()

}

# -----------------------------------------------------------------------------------
# ---- Build Random Forest Models ---------------------------------------------------
# -----------------------------------------------------------------------------------

{

# initialize the h2o instance
h2o.init(nthreads = workers, max_mem_size = "8g")

# remove any objects in the h2o instance
h2o.removeAll()

# remove the progress bar when model building
h2o.no_progress()

# ---- build prediction models for pmt_status ---------------------------------------

# determine how many data sets we have to build models for
if(cluster.data)
{
  # get the sub groups of dat
  groups = sort(unique(dat$Cluster))
  
  # also consider using all of dat
  groups = c("All", groups)
  
} else
{
  groups = "All"
}

# determine Random Forest predictive performance for each group
predict.pmt_status.rf = foreach(g = groups) %do%
{
  # ---- split the data ----
  
  # determine which data set to use
  if(cluster.data)
  {
    if(g == "All")
    {
      # update dat to not contain Cluster
      dat.g = data.table(dat[, !"Cluster"])
      
    } else
    {
      # update dat to only contain data in Cluster g
      dat.g = data.table(dat[Cluster == g, !"Cluster"])
    }
  } else
  {
    # get dat
    dat.g = data.table(dat)
  }
  
  # identify predictors (x)
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # consider two-way interactions between the first n indicators in dat.g
  # recall that the indicators are ordered by importance from random forests
  n = floor(length(x) / 3)
  two.way = data.table(dat.g[, x[1:n], with = FALSE])
  
  # should we add two way interactions in the model?
  add.interactions = FALSE
  
  # build interactions if desired
  if(add.interactions)
  {
    # create the two way interactions and add them to dat.g
    two.way = data.table(model.matrix(~.^2, data = two.way)[,-(1:(n + 1))])
    dat.g = cbind(dat.g, two.way)
  }
  
  # identify predictors (x) and response (y)
  y = "pmt_status"
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # build the fold assignment
  set.seed(42)
  k.folds = max(c(k.folds, 3))
  folds = createFolds(y = unname(unlist(dat.g[, y, with = FALSE])), k = k.folds)
  
  # split up dat.g into train, valid, and test
  train.rows = unname(unlist(lapply(1:(k.folds - 2), function(f) folds[[f]])))
  train = data.table(dat.g[train.rows])
  
  valid.rows = unname(unlist(folds[[k.folds - 1]]))
  valid = data.table(dat.g[valid.rows])
  
  test.rows = unname(unlist(folds[[k.folds]]))
  test = data.table(dat.g[test.rows])
  
  # split up YX.h2o into train, valid, and test
  train.YX.h2o = as.h2o(train[, c(y, x), with = FALSE])
  valid.YX.h2o = as.h2o(valid[, c(y, x), with = FALSE])
  test.YX.h2o = as.h2o(test[, c(y, x), with = FALSE])
  
  # compute the max class weight
  max.class.weight = table(unname(unlist(dat.g[, y, with = FALSE])))
  max.class.weight = as.numeric(max(max(max.class.weight) / max.class.weight))
  
  # ---- grid search for models ----
  
  # set up hyperparameters of interest
  rf.hyper.params = list(
    
    # -- initial tuning --
    ntrees = 250,
    min_rows = c(1, 11, 25, 41),
    max_depth = c(20, 40, 60),
    # min_split_improvement = c(0, 1e-5),
    
    # -- re-fine tuning --
    # ntrees = 250,
    # min_rows = 11,
    # max_depth = 30,
    min_split_improvement = 1e-5,
    
    # -- scoring function -- 
    stopping_metric = "mean_per_class_error")
  
  # lets use a random grid search and specify a time limit and/or model limit
  minutes = 20
  rf.search.criteria = list(strategy = "RandomDiscrete", 
                              max_runtime_secs = minutes * 60, 
                              # max_models = 100, 
                              seed = 42)
  
  # lets run a grid search for a good model, without drop out ratios
  h2o.rm("rf.random.grid")
  rf.random.grid = h2o.grid(algorithm = "randomForest",
                              grid_id = "rf.random.grid",
                              y = y,
                              x = x,
                              training_frame = train.YX.h2o,
                              validation_frame = valid.YX.h2o,
                            stopping_rounds = 5,
                              nfolds = 5,
                              keep_cross_validation_predictions = TRUE,
                              fold_assignment = "Modulo",
                              seed = 21,
                              balance_classes = TRUE,
                              max_after_balance_size = max.class.weight,
                              hyper_params = rf.hyper.params,
                              search_criteria = rf.search.criteria)
  
  # free up RAM
  gc()
  
  # rank each model in the random grids
  rf.grid = h2o.getGrid("rf.random.grid", sort_by = "auc", decreasing = TRUE)
  
  # get the summary table of the grid search
  DT.rf.grid = data.table(rf.grid@summary_table)
  DT.rf.grid
  
  # set up the data types for each column in DT.grid for plotting purposes
  DT.rf.grid = DT.rf.grid[, .(ntrees = as.numeric(ntrees),
                              min_rows = as.factor(min_rows),
                              max_depth = as.factor(max_depth),
                              min_split_improvement = as.factor(min_split_improvement),
                                  model = removePunctuation(gsub("[A-z]+", "", model_ids)),
                                  auc = as.numeric(auc))]
  
  # plot auc v. max_depth and ntrees to see which structure is most robust
  plot.rf.grid = ggplot(DT.rf.grid, aes(x = max_depth, y = auc, color = ntrees, fill = ntrees)) + 
    # geom_boxplot() + 
    geom_jitter(size = 3, alpha = 2/3) + 
    # scale_y_continuous(labels = dollar) + 
    ggtitle("Cross Validation Error") + 
    labs(x = "Max Depth", y = "AUC", color = "Trees", fill = "Trees") + 
    facet_wrap(~paste("Min Rows:", min_rows), nrow = 2) +
    theme_bw(base_size = 25) +
    theme(legend.position = "top", 
          legend.key.size = unit(.25, "in"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
  
  plot.rf.grid
  
  # pick the top model from all grid searches
  # rf.mod = h2o.getModel(rf.grid@model_ids[[1]])
  
  # find the top k models from all grid searches
  k = 5
  rf.mod.list = lapply(1:k, function(i) rf.grid@model_ids[[i]])
  
  # stack the top k models from all grid searches
  h2o.rm("rf.mod")
  rf.mod = h2o.stackedEnsemble(x = x,
                                 y = y,
                                 training_frame = train.YX.h2o,
                                 validation_frame = valid.YX.h2o,
                                 model_id = "rf.mod",
                                 base_models = rf.mod.list)
  
  # plot variable importance
  # h2o.varimp_plot(rf.mod)
  
  # ---- measure model quality ----
  
  # make predictions on each data set
  ynew.train = as.data.frame(predict(rf.mod, newdata = train.YX.h2o))$p1
  ynew.valid = as.data.frame(predict(rf.mod, newdata = valid.YX.h2o))$p1
  ynew.test = as.data.frame(predict(rf.mod, newdata = test.YX.h2o))$p1
  
  # get the true values from each data set
  ytrue.train = as.data.frame(train.YX.h2o)[,1]
  ytrue.valid = as.data.frame(valid.YX.h2o)[,1]
  ytrue.test = as.data.frame(test.YX.h2o)[,1]
  
  # plot the probability predictions on the testing data set
  # densityPlot(ynew.test)
  
  # compute prediction metrics on each data set
  rf.metrics.train = h2o.make_metrics(predicted = as.h2o(ynew.train), 
                                        actuals = as.h2o(factor(ytrue.train, levels = sort(unique(ytrue.train)))))
  
  rf.metrics.valid = h2o.make_metrics(predicted = as.h2o(ynew.valid), 
                                        actuals = as.h2o(factor(ytrue.valid, levels = sort(unique(ytrue.valid)))))
  
  rf.metrics.test = h2o.make_metrics(predicted = as.h2o(ynew.test), 
                                       actuals = as.h2o(factor(ytrue.test, levels = sort(unique(ytrue.test)))))
  
  # free up RAM
  gc()
  
  # check out the threshold metrics
  rf.metrics.train@metrics$max_criteria_and_metric_scores
  rf.metrics.valid@metrics$max_criteria_and_metric_scores
  rf.metrics.test@metrics$max_criteria_and_metric_scores
  
  # pick a threshold for converting probabilities into binary
  threshold = mean(c(tail(rf.metrics.train@metrics$max_criteria_and_metric_scores$threshold, 1),
                     tail(rf.metrics.valid@metrics$max_criteria_and_metric_scores$threshold, 1),
                     tail(rf.metrics.test@metrics$max_criteria_and_metric_scores$threshold, 1)))
  
  # compute a confusion matrix for each data set
  rf.confusion.train = confusion(ytrue = ytrue.train, ypred = as.numeric(ynew.train >= threshold))
  rf.confusion.valid = confusion(ytrue = ytrue.valid, ypred = as.numeric(ynew.valid >= threshold))
  rf.confusion.test = confusion(ytrue = ytrue.test, ypred = as.numeric(ynew.test >= threshold))
  
  # compute Sensitivity (accuracy in predicting 1's) for each data set
  rf.sensitivity.train = rf.confusion.train[2,2] / (rf.confusion.train[2,2] + rf.confusion.train[2,1])
  rf.sensitivity.valid = rf.confusion.valid[2,2] / (rf.confusion.valid[2,2] + rf.confusion.valid[2,1])
  rf.sensitivity.test = rf.confusion.test[2,2] / (rf.confusion.test[2,2] + rf.confusion.test[2,1])
  
  # compute Specificity (accuracy in predicting 0's) for each data set
  rf.specificity.train = rf.confusion.train[1,1] / (rf.confusion.train[1,1] + rf.confusion.train[1,2])
  rf.specificity.valid = rf.confusion.valid[1,1] / (rf.confusion.valid[1,1] + rf.confusion.valid[1,2])
  rf.specificity.test = rf.confusion.test[1,1] / (rf.confusion.test[1,1] + rf.confusion.test[1,2])
  
  # compute Odds Ratio (the odds of success over failure) for each data set
  rf.odds.ratio.train = (rf.confusion.train[1,1] * rf.confusion.train[2,2]) / (rf.confusion.train[2,1] * rf.confusion.train[1,2])
  rf.odds.ratio.valid = (rf.confusion.valid[1,1] * rf.confusion.valid[2,2]) / (rf.confusion.valid[2,1] * rf.confusion.valid[1,2])
  rf.odds.ratio.test = (rf.confusion.test[1,1] * rf.confusion.test[2,2]) / (rf.confusion.test[2,1] * rf.confusion.test[1,2])
  
  # compute Accuracy for each data set
  rf.accuracy.train = (rf.confusion.train[1,1] + rf.confusion.train[2,2]) / (rf.confusion.train[1,1] + rf.confusion.train[1,2] + rf.confusion.train[2,1] + rf.confusion.train[2,2])
  rf.accuracy.valid = (rf.confusion.valid[1,1] + rf.confusion.valid[2,2]) / (rf.confusion.valid[1,1] + rf.confusion.valid[1,2] + rf.confusion.valid[2,1] + rf.confusion.valid[2,2])
  rf.accuracy.test = (rf.confusion.test[1,1] + rf.confusion.test[2,2]) / (rf.confusion.test[1,1] + rf.confusion.test[1,2] + rf.confusion.test[2,1] + rf.confusion.test[2,2])
  
  # compute AUC for each data set
  rf.auc.train = h2o.auc(rf.metrics.train)
  rf.auc.valid = h2o.auc(rf.metrics.valid)
  rf.auc.test = h2o.auc(rf.metrics.test)
  
  # compute Log Loss for each data set
  rf.logloss.train = h2o.logloss(rf.metrics.train)
  rf.logloss.valid = h2o.logloss(rf.metrics.valid)
  rf.logloss.test = h2o.logloss(rf.metrics.test)
  
  # ---- compute kappa ----
  
  # get the total number of observations for each data set
  n.train = sum(rf.confusion.train)
  n.valid = sum(rf.confusion.valid)
  n.test = sum(rf.confusion.test)
  
  # get the vector of correct predictions for each data set 
  dia.train = diag(rf.confusion.train)
  dia.valid = diag(rf.confusion.valid)
  dia.test = diag(rf.confusion.test)
  
  # get the vector of the number of observations per class for each data set
  rsum.train = rowSums(rf.confusion.train)
  rsum.valid = rowSums(rf.confusion.valid)
  rsum.test = rowSums(rf.confusion.test)
  
  # get the vector of the number of predictions per class for each data set
  csum.train = colSums(rf.confusion.train)
  csum.valid = colSums(rf.confusion.valid)
  csum.test = colSums(rf.confusion.test)
  
  # get the proportion of observations per class for each data set
  p.train = rsum.train / n.train
  p.valid = rsum.valid / n.valid
  p.test = rsum.test / n.test
  
  # get the proportion of predcitions per class for each data set
  q.train = csum.train / n.train
  q.valid = csum.valid / n.valid
  q.test = csum.test / n.test
  
  # compute accuracy for each data set
  acc.train = sum(dia.train) / n.train
  acc.valid = sum(dia.valid) / n.valid
  acc.test = sum(dia.test) / n.test
  
  # compute expected accuracy for each data set
  exp.acc.train = sum(p.train * q.train)
  exp.acc.valid = sum(p.valid * q.valid)
  exp.acc.test = sum(p.test * q.test)
  
  # compute kappa for each data set
  kap.train = (acc.train - exp.acc.train) / (1 - exp.acc.train)
  kap.valid = (acc.valid - exp.acc.valid) / (1 - exp.acc.valid)
  kap.test = (acc.test - exp.acc.test) / (1 - exp.acc.test)
  
  # ---- finalize output ----
  
  # build a final metrics table for rf.mod
  rf.mod.table = data.table(Model = "Random Forest",
                              Group = g,
                              Threshold = threshold,
                              Metric = c("AUC", "Accuracy", "Log_Loss", "Kappa", "Odds_Ratio", "Specificity_0", "Sensitivity_1"),
                              Train = c(rf.auc.train, rf.accuracy.train, rf.logloss.train, kap.train, rf.odds.ratio.train, rf.specificity.train, rf.sensitivity.train),
                              Valid = c(rf.auc.valid, rf.accuracy.valid, rf.logloss.valid, kap.valid, rf.odds.ratio.valid, rf.specificity.valid, rf.sensitivity.valid),
                              Test = c(rf.auc.test, rf.accuracy.test, rf.logloss.test, kap.test, rf.odds.ratio.test, rf.specificity.test, rf.sensitivity.test))
  
  # build a final grid search table
  rf.grid.search = data.table(cbind(Model = rep("Random Forest", nrow(DT.rf.grid)), 
                                      Group = rep(g, nrow(DT.rf.grid)), 
                                      DT.rf.grid))
  
  # build a list of confusion matrices
  rf.confusion = list("Train" = rf.confusion.train, 
                        "Valid" = rf.confusion.valid, 
                        "Test" = rf.confusion.test)
  
  # build a list of final tables
  rf.list = list(rf.mod.table, rf.grid.search, rf.confusion)
  names(rf.list) = paste0(c("Random_Forest_Metrics_Group_", "Random_Forest_Search_Group_", "Random_Forest_Error_Group_"), g)
  
  # remove some objects
  rm(dat.g, DT.rf.grid, folds, rf.accuracy.test, rf.accuracy.train, rf.accuracy.valid,
     rf.auc.test, rf.auc.train, rf.auc.valid, rf.confusion.test, rf.confusion.train, rf.confusion.valid,
     rf.grid, rf.hyper.params, valid.rows, rf.confusion, rf.mod.table, rf.grid.search,
     rf.logloss.test, rf.logloss.train, rf.logloss.valid, rf.metrics.test, rf.metrics.train, rf.metrics.valid,
     rf.mod, rf.odds.ratio.test, rf.odds.ratio.train, rf.odds.ratio.valid,
     rf.search.criteria, rf.sensitivity.test, rf.mod.list,
     rf.sensitivity.train, rf.sensitivity.valid, rf.specificity.test, rf.specificity.train,
     rf.specificity.valid, max.class.weight, plot.rf.grid, train, test, valid, 
     small, test.rows, train.rows, k,
     test.YX.h2o, threshold, train.YX.h2o, valid.YX.h2o, rf.random.grid,
     x, y, ynew.test, ynew.train, ynew.valid, ytrue.test, ytrue.train, ytrue.valid,
     n.train, dia.train, rsum.train, csum.train, p.train, q.train, acc.train, exp.acc.train, kap.train,
     n.valid, dia.valid, rsum.valid, csum.valid, p.valid, q.valid, acc.valid, exp.acc.valid, kap.valid,
     n.test, dia.test, rsum.test, csum.test, p.test, q.test, acc.test, exp.acc.test, kap.test)
  
  # free up RAM
  gc()
  
  # return the final metrics table
  return(rf.list)
}

# name the results
names(predict.pmt_status.rf) = "h2o.randomForest"

# combine predict.pmt_status into one table
predict.pmt_status = append(predict.pmt_status, predict.pmt_status.rf)

# remove some objects
rm(rf.list, g, predict.pmt_status.rf)

# free up RAM
gc()

# ---- build prediction models for total_return ---------------------------------------

# determine how many data sets we have to build models for
if(cluster.data)
{
  # get the sub groups of dat
  groups = sort(unique(dat$Cluster))
  
  # also consider using all of dat
  groups = c("All", groups)
  
} else
{
  groups = "All"
}

# determine Random Forest predictive performance for each group
predict.total_return.rf = foreach(g = groups) %do%
{
  # ---- split the data ----
  
  # determine which data set to use
  if(cluster.data)
  {
    if(g == "All")
    {
      # update dat to not contain Cluster
      dat.g = data.table(dat[, !"Cluster"])
      
    } else
    {
      # update dat to only contain data in Cluster g
      dat.g = data.table(dat[Cluster == g, !"Cluster"])
    }
  } else
  {
    # get dat
    dat.g = data.table(dat)
  }
  
  # identify predictors (x)
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # consider two-way interactions between the first n indicators in dat.g
  # recall that the indicators are ordered by importance from random forests
  n = floor(length(x) / 3)
  two.way = data.table(dat.g[, x[1:n], with = FALSE])
  
  # should we add two way interactions in the model?
  add.interactions = FALSE
  
  # build interactions if desired
  if(add.interactions)
  {
    # create the two way interactions and add them to dat.g
    two.way = data.table(model.matrix(~.^2, data = two.way)[,-(1:(n + 1))])
    dat.g = cbind(dat.g, two.way)
  }
  
  # identify predictors (x) and response (y)
  y = "total_return"
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # build the fold assignment
  set.seed(42)
  k.folds = max(c(k.folds, 3))
  folds = createFolds(y = unname(unlist(dat.g[, y, with = FALSE])), k = k.folds)
  
  # split up dat.g into train, valid, and test
  train.rows = unname(unlist(lapply(1:(k.folds - 2), function(f) folds[[f]])))
  train = data.table(dat.g[train.rows])
  
  valid.rows = unname(unlist(folds[[k.folds - 1]]))
  valid = data.table(dat.g[valid.rows])
  
  test.rows = unname(unlist(folds[[k.folds]]))
  test = data.table(dat.g[test.rows])
  
  # split up YX.h2o into train, valid, and test
  train.YX.h2o = as.h2o(train[, c(y, x), with = FALSE])
  valid.YX.h2o = as.h2o(valid[, c(y, x), with = FALSE])
  test.YX.h2o = as.h2o(test[, c(y, x), with = FALSE])
  
  # ---- grid search for models ----
  
  # set up hyperparameters of interest, with drop out ratios
  rf.hyper.params = list(
    
    # -- initial tuning --
    ntrees = 250,
    min_rows = c(1, 11, 25, 41),
    max_depth = c(20, 40, 60),
    # min_split_improvement = c(0, 1e-5),
    
    # -- re-fine tuning --
    # ntrees = 250,
    # min_rows = 11,
    # max_depth = 30,
    min_split_improvement = 1e-5,
    
    # -- scoring function -- 
    stopping_metric = "RMSE")
  
  # lets use a random grid search and specify a time limit and/or model limit
  minutes = 20
  rf.search.criteria = list(strategy = "RandomDiscrete", 
                              max_runtime_secs = minutes * 60, 
                              # max_models = 100, 
                              seed = 42)
  
  # lets run a grid search for a good model, without drop out ratios
  h2o.rm("rf.random.grid")
  rf.random.grid = h2o.grid(algorithm = "randomForest",
                            grid_id = "rf.random.grid",
                            y = y,
                            x = x,
                            training_frame = train.YX.h2o,
                            validation_frame = valid.YX.h2o,
                            stopping_rounds = 5,
                            nfolds = 5,
                            keep_cross_validation_predictions = TRUE,
                            fold_assignment = "Modulo",
                            seed = 21,
                            hyper_params = rf.hyper.params,
                            search_criteria = rf.search.criteria)
  
  # free up RAM
  gc()
  
  # rank each model in the random grids
  rf.grid = h2o.getGrid("rf.random.grid", sort_by = "rmse", decreasing = FALSE)
  
  # get the summary table of the grid search
  DT.rf.grid = data.table(rf.grid@summary_table)
  DT.rf.grid
  
  # set up the data types for each column in DT.grid for plotting purposes
  DT.rf.grid = DT.rf.grid[, .(ntrees = as.numeric(ntrees),
                              min_rows = as.factor(min_rows),
                              max_depth = as.factor(max_depth),
                              min_split_improvement = as.factor(min_split_improvement),
                              model = removePunctuation(gsub("[A-z]+", "", model_ids)),
                              rmse = as.numeric(rmse))]
  
  # plot rmse v. max_depth and ntrees to see which structure is most robust
  plot.rf.grid = ggplot(DT.rf.grid, aes(x = max_depth, y = rmse, color = ntrees, fill = ntrees)) + 
    # geom_boxplot() + 
    geom_jitter(size = 3, alpha = 2/3) + 
    # scale_y_continuous(labels = dollar) + 
    ggtitle("Cross Validation Error") + 
    labs(x = "Max Depth", y = "RMSE", color = "Trees", fill = "Trees") + 
    facet_wrap(~paste("Min Rows:", min_rows), nrow = 2) +
    theme_bw(base_size = 25) +
    theme(legend.position = "top", 
          legend.key.size = unit(.25, "in"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
  
  plot.rf.grid
  
  # pick the top model from all grid searches
  # rf.mod = h2o.getModel(rf.grid@model_ids[[1]])
  
  # find the top k models from all grid searches
  k = 5
  rf.mod.list = lapply(1:k, function(i) rf.grid@model_ids[[i]])
  
  # stack the top k models from all grid searches
  h2o.rm("rf.mod")
  rf.mod = h2o.stackedEnsemble(x = x,
                                 y = y,
                                 training_frame = train.YX.h2o,
                                 validation_frame = valid.YX.h2o,
                                 model_id = "rf.mod",
                                 base_models = rf.mod.list)
  
  # plot variable importance
  # h2o.varimp_plot(rf.mod)
  
  # ---- measure model quality ----
  
  # make predictions on each data set
  ynew.train = as.data.frame(predict(rf.mod, newdata = train.YX.h2o))$predict
  ynew.valid = as.data.frame(predict(rf.mod, newdata = valid.YX.h2o))$predict
  ynew.test = as.data.frame(predict(rf.mod, newdata = test.YX.h2o))$predict
  
  # get the true values from each data set
  ytrue.train = as.data.frame(train.YX.h2o)[,1]
  ytrue.valid = as.data.frame(valid.YX.h2o)[,1]
  ytrue.test = as.data.frame(test.YX.h2o)[,1]
  
  # transform the predicted and actual values back to original units if needed
  if(use.norm)
  {
    # map predicted values to the original units of the response variable
    ynew.train = predict(best.norm, newdata = ynew.train, inverse = TRUE)
    ynew.valid = predict(best.norm, newdata = ynew.valid, inverse = TRUE)
    ynew.test = predict(best.norm, newdata = ynew.test, inverse = TRUE)
    
    # map actual values to the original units of the response variable
    ytrue.train = predict(best.norm, newdata = ytrue.train, inverse = TRUE)
    ytrue.valid = predict(best.norm, newdata = ytrue.valid, inverse = TRUE)
    ytrue.test = predict(best.norm, newdata = ytrue.test, inverse = TRUE)
  }
  
  # compute prediction metrics on each data set
  rf.metrics.train = h2o.make_metrics(predicted = as.h2o(ynew.train), 
                                        actuals = as.h2o(ytrue.train))
  
  rf.metrics.valid = h2o.make_metrics(predicted = as.h2o(ynew.valid), 
                                        actuals = as.h2o(ytrue.valid))
  
  rf.metrics.test = h2o.make_metrics(predicted = as.h2o(ynew.test), 
                                       actuals = as.h2o(ytrue.test))
  
  # free up RAM
  gc()
  
  # check out the prediction metrics
  rf.metrics.train
  rf.metrics.valid
  rf.metrics.test
  
  # fit the normal distribution to the residuals
  rf.resid.train = fitdist(data = ytrue.train - ynew.train, distr = "norm")
  rf.resid.valid = fitdist(data = ytrue.valid - ynew.valid, distr = "norm") 
  rf.resid.test = fitdist(data = ytrue.test - ynew.test, distr = "norm") 
  
  # ---- finalize output ----
  
  # build a final metrics table for rf.mod
  rf.mod.table = data.table(Model = "Random Forest",
                              Group = g,
                              Metric = c("R2", "RMSE"),
                              Train = c(rf.metrics.train@metrics$r2, rf.metrics.train@metrics$RMSE),
                              Valid = c(rf.metrics.valid@metrics$r2, rf.metrics.valid@metrics$RMSE),
                              Test = c(rf.metrics.test@metrics$r2, rf.metrics.test@metrics$RMSE))
  
  # build a final grid search table
  rf.grid.search = data.table(cbind(Model = rep("Random Forest", nrow(DT.rf.grid)), 
                                      Group = rep(g, nrow(DT.rf.grid)), 
                                      DT.rf.grid))
  
  # build a list of residual plots
  rf.resid = list("Train" = rf.resid.train, 
                    "Valid" = rf.resid.valid, 
                    "Test" = rf.resid.test)
  
  # build a list of final tables
  rf.list = list(rf.mod.table, rf.grid.search, rf.resid)
  names(rf.list) = paste0(c("Random_Forest_Metrics_Group_", "Random_Forest_Search_Group_", "Random_Forest_Error_Group_"), g)
  
  # remove some objects
  rm(dat.g, DT.rf.grid, folds, n,
     rf.grid, rf.hyper.params, 
     rf.metrics.test, rf.metrics.train, rf.metrics.valid,
     rf.mod, rf.resid.test, rf.resid.train, rf.resid.valid,
     rf.search.criteria, small, medium, large, k,
     plot.rf.grid, train, test, valid, 
     rf.random.grid, rf.mod.table, rf.grid.search, rf.resid,
     test.YX.h2o, train.YX.h2o, two.way, valid.YX.h2o,
     x, y, ynew.test, ynew.train, ynew.valid, ytrue.test, ytrue.train, ytrue.valid)
  
  # free up RAM
  gc()
  
  # return the final metrics table
  return(rf.list)
}

# name the results
names(predict.total_return.rf) = "h2o.randomForest"

# combine predict.total_return into one table
predict.total_return = append(predict.total_return, predict.total_return.rf)

# remove some objects
rm(rf.list, g, predict.total_return.rf)

# free up RAM
gc()

# ---- build prediction models for total_return_class ---------------------------------------

# determine how many data sets we have to build models for
if(cluster.data)
{
  # get the sub groups of dat
  groups = sort(unique(dat$Cluster))
  
  # also consider using all of dat
  groups = c("All", groups)
  
} else
{
  groups = "All"
}

# determine Random Forest predictive performance for each group
predict.total_return_class.rf = foreach(g = groups) %do%
{
  # ---- split the data ----
  
  # determine which data set to use
  if(cluster.data)
  {
    if(g == "All")
    {
      # update dat to not contain Cluster
      dat.g = data.table(dat[, !"Cluster"])
      
    } else
    {
      # update dat to only contain data in Cluster g
      dat.g = data.table(dat[Cluster == g, !"Cluster"])
    }
  } else
  {
    # get dat
    dat.g = data.table(dat)
  }
  
  # identify predictors (x)
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # consider two-way interactions between the first n indicators in dat.g
  # recall that the indicators are ordered by importance from random forests
  n = floor(length(x) / 3)
  two.way = data.table(dat.g[, x[1:n], with = FALSE])
  
  # should we add two way interactions in the model?
  add.interactions = FALSE
  
  # build interactions if desired
  if(add.interactions)
  {
    # create the two way interactions and add them to dat.g
    two.way = data.table(model.matrix(~.^2, data = two.way)[,-(1:(n + 1))])
    dat.g = cbind(dat.g, two.way)
  }
  
  # identify predictors (x) and response (y)
  y = "total_return_class"
  x = names(dat.g[, !c("id", responses), with = FALSE])
  
  # build the fold assignment
  set.seed(42)
  k.folds = max(c(k.folds, 3))
  folds = createFolds(y = unname(unlist(dat.g[, y, with = FALSE])), k = k.folds)
  
  # split up dat.g into train, valid, and test
  train.rows = unname(unlist(lapply(1:(k.folds - 2), function(f) folds[[f]])))
  train = data.table(dat.g[train.rows])
  
  valid.rows = unname(unlist(folds[[k.folds - 1]]))
  valid = data.table(dat.g[valid.rows])
  
  test.rows = unname(unlist(folds[[k.folds]]))
  test = data.table(dat.g[test.rows])
  
  # split up YX.h2o into train, valid, and test
  train.YX.h2o = as.h2o(train[, c(y, x), with = FALSE])
  valid.YX.h2o = as.h2o(valid[, c(y, x), with = FALSE])
  test.YX.h2o = as.h2o(test[, c(y, x), with = FALSE])
  
  # compute the max class weight
  max.class.weight = table(unname(unlist(dat.g[, y, with = FALSE])))
  max.class.weight = as.numeric(max(max(max.class.weight) / max.class.weight))
  
  # ---- grid search for models ----
  
  # set up hyperparameters of interest, with drop out ratios
  rf.hyper.params = list(
    
    # -- initial tuning --
    ntrees = 250,
    min_rows = c(1, 11, 25, 41),
    max_depth = c(20, 40, 60),
    # min_split_improvement = c(0, 1e-5),
    
    # -- re-fine tuning --
    # ntrees = 250,
    # min_rows = 11,
    # max_depth = 30,
    min_split_improvement = 1e-5,
    
    # -- scoring function -- 
    stopping_metric = "mean_per_class_error")
  
  # lets use a random grid search and specify a time limit and/or model limit
  minutes = 20
  rf.search.criteria = list(strategy = "RandomDiscrete", 
                              max_runtime_secs = minutes * 60, 
                              # max_models = 100, 
                              seed = 42)
  
  # lets run a grid search for a good model, without drop out ratios
  h2o.rm("rf.random.grid")
  rf.random.grid = h2o.grid(algorithm = "randomForest",
                            grid_id = "rf.random.grid",
                            y = y,
                            x = x,
                            training_frame = train.YX.h2o,
                            validation_frame = valid.YX.h2o,
                            stopping_rounds = 5,
                            nfolds = 5,
                            keep_cross_validation_predictions = TRUE,
                            fold_assignment = "Modulo",
                            seed = 21,
                            balance_classes = TRUE,
                            max_after_balance_size = max.class.weight,
                            hyper_params = rf.hyper.params,
                            search_criteria = rf.search.criteria)
  
  # free up RAM
  gc()
  
  # rank each model in the random grids
  rf.grid = h2o.getGrid("rf.random.grid", sort_by = "mean_per_class_error", decreasing = FALSE)
  
  # get the summary table of the grid search
  DT.rf.grid = data.table(rf.grid@summary_table)
  DT.rf.grid
  
  # set up the data types for each column in DT.grid for plotting purposes
  DT.rf.grid = DT.rf.grid[, .(ntrees = as.numeric(ntrees),
                              min_rows = as.factor(min_rows),
                              max_depth = as.factor(max_depth),
                              min_split_improvement = as.factor(min_split_improvement),
                              model = removePunctuation(gsub("[A-z]+", "", model_ids)),
                              mean_per_class_error = as.numeric(mean_per_class_error))]
  
  # plot mean_per_class_error v. max_depth and ntrees to see which structure is most robust
  plot.rf.grid = ggplot(DT.rf.grid, aes(x = max_depth, y = mean_per_class_error, color = ntrees, fill = ntrees)) + 
    # geom_boxplot() + 
    geom_jitter(size = 3, alpha = 2/3) + 
    # scale_y_continuous(labels = dollar) + 
    ggtitle("Cross Validation Error") + 
    labs(x = "Max Depth", y = "Mean Per Class Error", color = "Trees", fill = "Trees") + 
    facet_wrap(~paste("Min Rows:", min_rows), nrow = 2) +
    theme_bw(base_size = 25) +
    theme(legend.position = "top", 
          legend.key.size = unit(.25, "in"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
  
  plot.rf.grid
  
  # pick the top model from all grid searches
  # rf.mod = h2o.getModel(rf.grid@model_ids[[1]])
  
  # find the top k models from all grid searches
  k = 5
  rf.mod.list = lapply(1:k, function(i) rf.grid@model_ids[[i]])
  
  # stack the top k models from all grid searches
  h2o.rm("rf.mod")
  rf.mod = h2o.stackedEnsemble(x = x,
                                 y = y,
                                 training_frame = train.YX.h2o,
                                 validation_frame = valid.YX.h2o,
                                 model_id = "rf.mod",
                                 base_models = rf.mod.list)
  
  # plot variable importance
  # h2o.varimp_plot(rf.mod)
  
  # ---- measure model quality ----
  
  # make predictions on each data set
  ynew.train = as.matrix(predict(rf.mod, newdata = train.YX.h2o)[,-1])
  ynew.valid = as.matrix(predict(rf.mod, newdata = valid.YX.h2o)[,-1])
  ynew.test = as.matrix(predict(rf.mod, newdata = test.YX.h2o)[,-1])
  
  # get the true values from each data set
  ytrue.train = as.data.frame(train.YX.h2o[,1])
  ytrue.valid = as.data.frame(valid.YX.h2o[,1])
  ytrue.test = as.data.frame(test.YX.h2o[,1])
  
  # ---- compute multi log loss ----
  
  # build a matrix indicating the true class values for each data set
  ytrue.train.mat = model.matrix(~., data = ytrue.train)[,-1]
  ytrue.train.mat = cbind(1 - rowSums(ytrue.train.mat), ytrue.train.mat)
  
  ytrue.valid.mat = model.matrix(~., data = ytrue.valid)[,-1]
  ytrue.valid.mat = cbind(1 - rowSums(ytrue.valid.mat), ytrue.valid.mat)
  
  ytrue.test.mat = model.matrix(~., data = ytrue.test)[,-1]
  ytrue.test.mat = cbind(1 - rowSums(ytrue.test.mat), ytrue.test.mat)
  
  # compute the multi-class logarithmic loss for each data set
  mll.train = MultiLogLoss(y_pred = ynew.train, y_true = ytrue.train.mat)
  mll.valid = MultiLogLoss(y_pred = ynew.valid, y_true = ytrue.valid.mat)
  mll.test = MultiLogLoss(y_pred = ynew.test, y_true = ytrue.test.mat)
  
  # free up RAM
  gc()
  
  # ---- compute kappa ----
  
  # get the predicted classes and actual classes for each data set
  ynew.train.code = apply(ynew.train, 1, which.max) - 1
  ytrue.train.code = factor(as.numeric(ytrue.train[,1]) - 1, levels = 0:(ncol(ytrue.train.mat) - 1))
  
  ynew.valid.code = apply(ynew.valid, 1, which.max) - 1
  ytrue.valid.code = factor(as.numeric(ytrue.valid[,1]) - 1, levels = 0:(ncol(ytrue.valid.mat) - 1))
  
  ynew.test.code = apply(ynew.test, 1, which.max) - 1
  ytrue.test.code = factor(as.numeric(ytrue.test[,1]) - 1, levels = 0:(ncol(ytrue.test.mat) - 1))
  
  # build a square confusion matrix for each data set
  conf.train = confusion(ytrue = ytrue.train.code, ypred = ynew.train.code)
  conf.valid = confusion(ytrue = ytrue.valid.code, ypred = ynew.valid.code)
  conf.test = confusion(ytrue = ytrue.test.code, ypred = ynew.test.code)
  
  # get the total number of observations for each data set
  n.train = sum(conf.train)
  n.valid = sum(conf.valid)
  n.test = sum(conf.test)
  
  # get the vector of correct predictions for each data set 
  dia.train = diag(conf.train)
  dia.valid = diag(conf.valid)
  dia.test = diag(conf.test)
  
  # get the vector of the number of observations per class for each data set
  rsum.train = rowSums(conf.train)
  rsum.valid = rowSums(conf.valid)
  rsum.test = rowSums(conf.test)
  
  # get the vector of the number of predictions per class for each data set
  csum.train = colSums(conf.train)
  csum.valid = colSums(conf.valid)
  csum.test = colSums(conf.test)
  
  # get the proportion of observations per class for each data set
  p.train = rsum.train / n.train
  p.valid = rsum.valid / n.valid
  p.test = rsum.test / n.test
  
  # get the proportion of predcitions per class for each data set
  q.train = csum.train / n.train
  q.valid = csum.valid / n.valid
  q.test = csum.test / n.test
  
  # compute accuracy for each data set
  acc.train = sum(dia.train) / n.train
  acc.valid = sum(dia.valid) / n.valid
  acc.test = sum(dia.test) / n.test
  
  # compute expected accuracy for each data set
  exp.acc.train = sum(p.train * q.train)
  exp.acc.valid = sum(p.valid * q.valid)
  exp.acc.test = sum(p.test * q.test)
  
  # compute kappa for each data set
  kap.train = (acc.train - exp.acc.train) / (1 - exp.acc.train)
  kap.valid = (acc.valid - exp.acc.valid) / (1 - exp.acc.valid)
  kap.test = (acc.test - exp.acc.test) / (1 - exp.acc.test)
  
  # ---- compute one-vs-all metrics ----
  
  # compute a binary confusion matrix for each class, for each data set
  one.v.all.train = lapply(1:nrow(conf.train), function(i)
  {
    # get the four entries of a binary confusion matrix
    v = c(conf.train[i,i], 
          rsum.train[i] - conf.train[i,i], 
          csum.train[i] - conf.train[i,i], 
          n.train - rsum.train[i] - csum.train[i] + conf.train[i,i]);
    
    # build the confusion matrix
    return(matrix(v, nrow = 2, byrow = TRUE))
  })
  
  one.v.all.valid = lapply(1:nrow(conf.valid), function(i)
  {
    # get the four entries of a binary confusion matrix
    v = c(conf.valid[i,i], 
          rsum.valid[i] - conf.valid[i,i], 
          csum.valid[i] - conf.valid[i,i], 
          n.valid - rsum.valid[i] - csum.valid[i] + conf.valid[i,i]);
    
    # build the confusion matrix
    return(matrix(v, nrow = 2, byrow = TRUE))
  })
  
  one.v.all.test = lapply(1:nrow(conf.test), function(i)
  {
    # get the four entries of a binary confusion matrix
    v = c(conf.test[i,i], 
          rsum.test[i] - conf.test[i,i], 
          csum.test[i] - conf.test[i,i], 
          n.test - rsum.test[i] - csum.test[i] + conf.test[i,i]);
    
    # build the confusion matrix
    return(matrix(v, nrow = 2, byrow = TRUE))
  })
  
  # sum up all of the matrices for each data set
  one.v.all.train = Reduce('+', one.v.all.train)
  one.v.all.valid = Reduce('+', one.v.all.valid)
  one.v.all.test = Reduce('+', one.v.all.test)
  
  # compute the micro average accuracy for each data set
  micro.acc.train = sum(diag(one.v.all.train)) / sum(one.v.all.train)
  micro.acc.valid = sum(diag(one.v.all.valid)) / sum(one.v.all.valid)
  micro.acc.test = sum(diag(one.v.all.test)) / sum(one.v.all.test)
  
  # get the macro accuracy for each data set
  macro.acc.train = acc.train
  macro.acc.valid = acc.valid
  macro.acc.test = acc.test
  
  # ---- finalize output ----
  
  # build a final metrics table for rf.mod
  rf.mod.table = data.table(Model = "Random Forest",
                              Group = g,
                              Metric = c("Log_Loss", "Kappa", "Macro_Accuracy", "Micro_Accuracy"),
                              Train = c(mll.train, kap.train, macro.acc.train, micro.acc.train),
                              Valid = c(mll.valid, kap.valid, macro.acc.valid, micro.acc.valid),
                              Test = c(mll.test, kap.test, macro.acc.test, micro.acc.test))
  
  # build a final grid search table
  rf.grid.search = data.table(cbind(Model = rep("Random Forest", nrow(DT.rf.grid)), 
                                      Group = rep(g, nrow(DT.rf.grid)), 
                                      DT.rf.grid))
  
  # build a list of final tables
  rf.list = list(rf.mod.table, rf.grid.search)
  names(rf.list) = paste0(c("Random_Forest_Metrics_Group_", "Random_Forest_Search_Group_"), g)
  
  # remove some objects
  rm(dat.g, DT.rf.grid, folds, n,
     rf.grid, rf.hyper.params, 
     dia.train, ynew.train.code, ytrue.train.code, ytrue.train.mat, 
     dia.valid, ynew.valid.code, ytrue.valid.code, ytrue.valid.mat, 
     dia.test, ynew.test.code, ytrue.test.code, ytrue.test.mat, 
     p.train, q.train, acc.train, exp.acc.train,
     p.valid, q.valid, acc.valid, exp.acc.valid,
     p.test, q.test, acc.test, exp.acc.test, max.class.weight,
     conf.train, rsum.train, csum.train, n.train, 
     conf.valid, rsum.valid, csum.valid, n.valid, 
     conf.test, rsum.test, csum.test, n.test, rf.mod.table, rf.grid.search,
     one.v.all.train, one.v.all.valid, one.v.all.test,
     mll.train, kap.train, macro.acc.train, micro.acc.train,
     mll.valid, kap.valid, macro.acc.valid, micro.acc.valid,
     mll.test, kap.test, macro.acc.test, micro.acc.test,
     rf.mod, rf.search.criteria, rf.mod.list, rf.random.grid,
     plot.rf.grid, train, test, valid, small, medium, large,
     test.YX.h2o, train.YX.h2o, two.way, valid.YX.h2o, minutes, k,
     x, y, ynew.test, ynew.train, ynew.valid, ytrue.test, ytrue.train, ytrue.valid)
  
  # free up RAM
  gc()
  
  # return the final tables
  return(rf.list)
}

# name the results
names(predict.total_return_class.rf) = "h2o.randomForest"

# combine predict.total_return into one table
predict.total_return_class = append(predict.total_return_class, predict.total_return_class.rf)

# remove some objects
rm(rf.list, g, predict.total_return_class.rf)

# clean up the data in the h2o cluster
h2o.removeAll()

# shut down the h2o cluster to free up RAM
h2o.shutdown(prompt = FALSE)

# free up RAM
gc()

}



























