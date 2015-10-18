renderWave = function(expr, env=parent.frame(), quoted=FALSE) {
    # Convert the expression + environment into a function
    func <- exprToFunction(expr, env, quoted)
    
    function() {
      val <- func()
      inputw(wave = wav, f = wav@samp.rate)
    }
  }


# output$timeSeries1 <- renderTimeSeries({
#   ts(matrix(rnorm(300), 100, 3), start=c(1961, 1), frequency=12)
# })
