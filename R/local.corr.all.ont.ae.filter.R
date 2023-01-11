local.corr.all.ont.ae.filter = function (s, t.p, ae.index.v) 
{
  well.list <- unique(s[[t.p]]$cw)
  temp = list()
  r.vector = list()
  r.vector[[length(well.list)]] = 0
  names(r.vector) = well.list
  for (cur.well in well.list) {
    index.ch = which(s[[t.p]]$cw == cur.well)
    # restricting the list of channels to active electrodes
    index.ch = index.ch[ae.index.v[index.ch]]
    # if none of the electrodes in cur.well are active, then the correlation for that well is undefined
    if(length(index.ch) %in% c(0,1)) {
      # if 0 active elctrodes in well - 
      # In create_burst_ont_Data, we only select [well.names] (wells with active electrodes)
      # So wells with 0 active electrodes will be dropped regardless
      # If there is only 1 AE, there can be no correlation.
      # So r will be set to 0
      r.vector[[cur.well]] = 0
      next
    }
    start.ms = round(s[[t.p]]$rec.time[1] * 10^2, digits = 0)
    stop.ms = round(s[[t.p]]$rec.time[2] * 10^2, digits = 0)
    total.ms = stop.ms - start.ms
    dig.s = list()
    dig.s[[length(index.ch)]] = rep(0, total.ms)
    for (cur.index in c(1:length(index.ch))) {
      cur.ch = index.ch[cur.index] # since index.ch is restricted to AEs, the current channel (cur.ch) will only include active electrodes
      dig.s[[cur.index]] = rep(0, total.ms)
      #s[[t.p]]$spikes is a list of every electrode that has at least one spike
      # s[[t.p]]$spikes[A1_11] is a list of all of the start times of spikes on A1_11
      # ordered and grouped by well
      # so current.s.ms is a list of the all the start times of every spike in the currrent electrode ("channel") 
      current.s.ms = round(s[[t.p]]$spikes[[cur.ch]] * 
                             10^2, digits = 0) - start.ms
      for (j in current.s.ms) {
        dig.s[[cur.index]][j] = 1
      }
    }
    c.vector = c()
    ch.names = c()
    for (cur.ch in c(1:length(index.ch))) {
      temp.corr = c()
      # correlation calculated here using cor()
      for (other.ch in setdiff(c(1:length(index.ch)), cur.ch)) {
        temp.corr = c(temp.corr, cor(dig.s[[cur.ch]], 
                                     dig.s[[other.ch]]))
      }
      c.vector = c(c.vector, mean(temp.corr, na.rm = T))
    }
    names(c.vector) = s[[t.p]]$channels[index.ch]
    r.vector[[cur.well]] = c.vector
  }
  r.vector
}
