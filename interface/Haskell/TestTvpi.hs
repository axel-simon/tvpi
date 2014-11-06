import Tvpi

main = do
  tvpi1 <- new 5
  tvpi2 <- copy tvpi1
  var0 <- createVariable tvpi1 5
  var1 <- createUnboundedVariable tvpi1
  var2 <- createUnboundedVariable tvpi1
  addInequalities tvpi1 var0 var1
    [Inequality 1 2 3, Inequality 5 6 2]
  dump tvpi1 ["var"++show n | n <- [0..2]]