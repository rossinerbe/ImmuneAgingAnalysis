#function to get higher order cell types from mixture result object

higherOrderImmuneTypes = function(mixtureData){
  mixtureData$Bcells = apply(mixtureData, 1, function(x) {sum(as.numeric(x[1:3]))})/apply(mixtureData, 1, function(x) {sum(as.numeric(x[1:22]))})
  mixtureData$Tcells = apply(mixtureData, 1, function(x) {sum(as.numeric(x[4:10]))})/apply(mixtureData, 1, function(x) {sum(as.numeric(x[1:22]))})
  mixtureData$NKcells = apply(mixtureData, 1, function(x) {sum(as.numeric(x[11:12]))})/apply(mixtureData, 1, function(x) {sum(as.numeric(x[1:22]))})
  mixtureData$MiscMyeloid = apply(mixtureData, 1, function(x) {sum(as.numeric(x[c(13,19:22)]))})/apply(mixtureData, 1, function(x) {sum(as.numeric(x[1:22]))})
  mixtureData$Macrophages = apply(mixtureData, 1, function(x) {sum(as.numeric(x[14:16]))})/apply(mixtureData, 1, function(x) {sum(as.numeric(x[1:22]))})
  mixtureData$Dendriticcells = apply(mixtureData, 1, function(x) {sum(as.numeric(x[17:18]))})/apply(mixtureData, 1, function(x) {sum(as.numeric(x[1:22]))})
  return(mixtureData)
}