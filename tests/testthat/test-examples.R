
context("Results consistent with original version of the code")

drug.case = 'FLUN'
drug.control = 'FLU'

tmp <- tempfile()

test_that("Example data type 1 and method = aeks", {
  expect_known_output(
    enrich(df = flu1, dd.group = group, drug.case = drug.case, 
           drug.control = drug.control, method = 'aeks', n_iter=1000, 
           seed = 2277)$Final_result, tmp)
})

test_that("Example data type 2 and method = aeks", {
  expect_known_output(
    enrich(df = flu2, dd.group = group, drug.case = drug.case, 
           drug.control = drug.control, method = 'aeks', n_iter=1000, 
           seed = 2277)$Final_result, tmp)
})

test_that("Example data type 1 and method = aefisher", {
  expect_known_output(
    enrich(df = flu1, dd.group = group, drug.case = drug.case, 
           drug.control = drug.control, method = 'aeks', n_iter=1000, 
           seed = 2277)$Final_result, tmp)
})

test_that("Example data type 2 and method = aefisher", {
  expect_known_output(
    enrich(df = flu2, dd.group = group, drug.case = drug.case, 
           drug.control = drug.control, method = 'aeks', n_iter=1000, 
           seed = 2277)$Final_result, tmp)
})
