
context("Results consistent with original version of the code")

drug.case = 'FLUN'
drug.control = 'FLU'

load('data/flu1.RData')
load('data/flu2.RData')
load('data/group.RData')

test_that("Example data type 1 and method = aeks", {
  all.equal(enrich(df = flu1, dd.group = group, drug.case = drug.case, 
                   drug.control = drug.control, method = 'aeks', n_iter=1000, 
                   seed = 2277)$Final_result, 
            load(file = 'tests/data/ks1.RData'))
})

test_that("Example data type 2 and method = aeks", {
  all.equal(enrich(df = flu2, dd.group = group, drug.case = drug.case, 
                   drug.control = drug.control, method = 'aeks', n_iter=1000, 
                   seed = 2277)$Final_result, 
            load(file = 'tests/data/ks2.RData'))
})

test_that("Example data type 1 and method = aefisher", {
  all.equal(enrich(df = flu1, dd.group = group, drug.case = drug.case, 
                   drug.control = drug.control, method = 'aeks', n_iter=1000, 
                   seed = 2277)$Final_result, 
            load(file = 'tests/data/fr1.RData'))
})

test_that("Example data type 2 and method = aefisher", {
  all.equal(enrich(df = flu2, dd.group = group, drug.case = drug.case, 
                   drug.control = drug.control, method = 'aeks', n_iter=1000, 
                   seed = 2277)$Final_result, 
            load(file = 'tests/data/fr2.RData'))
})
